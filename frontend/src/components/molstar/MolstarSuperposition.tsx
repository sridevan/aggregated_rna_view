import { useEffect, useRef, useState, type ReactNode } from 'react'
import { PDBeMolstarPlugin } from 'pdbe-molstar/lib/viewer'
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder'
import { StateTransforms } from 'molstar/lib/mol-plugin-state/transforms'
import { createStructureRepresentationParams } from 'molstar/lib/mol-plugin-state/helpers/structure-representation-params'
import { PluginCommands } from 'molstar/lib/mol-plugin/commands'
import { setSubtreeVisibility } from 'molstar/lib/mol-plugin/behavior/static/state'
import { Script } from 'molstar/lib/mol-script/script'
import { StructureSelection, StructureElement, StructureProperties, type Structure } from 'molstar/lib/mol-model/structure'
import { Mat4, Vec3, Vec4 } from 'molstar/lib/mol-math/linear-algebra'
import { Color } from 'molstar/lib/mol-util/color'
import { OrderedSet } from 'molstar/lib/mol-data/int'
import type { OrientationTransform } from '../../types/backend'
import { isApplicable, toMat4 } from '../../lib/orientation'
import 'pdbe-molstar/build/pdbe-molstar-light.css'
import { chainBytes, surroundingsBytes, residueLabel, sceneFetchCount, type SceneApi, type SceneAtom, type SceneEntry, type SceneHooks } from '../../lib/molstarScene'

export type { SceneApi, SceneAtom, SceneEntry, SceneHooks }

/*
 * Mol* scene of one or two RNA chains superposed in the family frame (Structures tab, Interactions compare). Same engine, camera
 * convention and look as the Family Overview viewer (PdbeMolstarViewer): one plugin per mount, white background, fixed view direction
 * in the family frame, uniform cartoons. Each chain is fetched once as BinaryCIF from the PDBe ModelServer (label_asym_id, model 1;
 * cached per chain for the page session), parsed in Mol*, and placed by a TransformStructureConformation node holding its
 * orientation.json transform (x' = R x + t, the family frame), optionally preceded by an extra rigid transform (e.g. a local
 * base-anchored fit). Changing a transform, the emphasis or the focused residues never refetches or re-parses a chain.
 */
const WHITE = { r: 255, g: 255, b: 255 }
/** Show ligands: as in the Family Overview viewer (PdbeMolstarViewer), small-molecule ligand instances with any atom within this distance of the chain */
const LIGAND_RADIUS = 4.0
const LIGAND_EXCLUDED_CLASSES = ['additive', 'metal_ion']
const VIEW_UP = Vec3.create(0, 1, 0), VIEW_DIR = Vec3.create(0, 0, -1)   // fixed camera orientation in the family frame (as the Family Overview)

const hexColor = (h: string) => Color(parseInt(h.replace('#', ''), 16))
const orientationMat = (t: OrientationTransform | null) => (isApplicable(t) ? toMat4(t!) : null) ?? Mat4.identity()

interface Nodes { entry: SceneEntry; data: string; structure: string; transform: string; component: string; repr: string; focus: string | null }
/** ligands near one chain: a separate ModelServer surroundings structure, placed with the chain's matrix */
interface LigNodes { data: string; transform: string; shown: string[] }

export function MolstarSuperposition({ entries, emphasis = null, extra = {}, focus = null, onApi, onStatus, onHooks, ligandClasses, children, testId, className = '', style, dataAttrs = {} }: {
  entries: SceneEntry[]
  /** seq_name drawn thicker and opaque while the others are faded (only with more than one structure) */
  emphasis?: string | null
  /** extra rigid transform (column-major 4x4) applied in the family frame after the orientation, per seq_name */
  extra?: Record<string, Mat4 | null>
  /** residues (label_seq_id) drawn as sticks in their structure's colour; the camera then focuses on them */
  focus?: Record<string, number[]> | null
  onApi?: (api: SceneApi | null) => void
  onStatus?: (s: 'loading' | 'ready' | 'error') => void
  /** test hooks after every scene change (loaded chains, fetch count, transforms, residue probe) */
  onHooks?: (h: SceneHooks) => void
  /** family ligand classification (CCD -> class) for Show ligands: ions and crystallisation additives are left out */
  ligandClasses?: () => Promise<Record<string, string>>
  /** overlays drawn over the canvas (labels) */
  children?: ReactNode
  testId: string; className?: string; style?: React.CSSProperties
  dataAttrs?: Record<string, string>
}) {
  const host = useRef<HTMLDivElement>(null), outer = useRef<HTMLDivElement>(null)
  const plugin = useRef<PDBeMolstarPlugin | null>(null)
  const queue = useRef<Promise<void>>(Promise.resolve())
  const nodes = useRef(new Map<string, Nodes>())
  const ligNodes = useRef(new Map<string, LigNodes>())
  const [showLigands, setShowLigands] = useState(false)
  const showLigRef = useRef(showLigands); showLigRef.current = showLigands
  const [ligInfo, setLigInfo] = useState<Record<string, string[]> | null>(null)   // shown ligand instances per seq_name (after the first Show)
  const [ligBusy, setLigBusy] = useState(false)
  const classesRef = useRef(ligandClasses); classesRef.current = ligandClasses
  const classes = useRef<Promise<Record<string, string>> | null>(null)
  const [gen, setGen] = useState(0)
  const [status, setStatusState] = useState<'loading' | 'ready' | 'error'>('loading')
  const [loadedKey, setLoadedKey] = useState<string | null>(null)   // entries the scene was last built for: until it matches, the scene is loading
  const [err, setErr] = useState<string | null>(null)
  const [isFull, setIsFull] = useState(false)
  const [hover, setHover] = useState<{ text: string; x: number; y: number; flip: boolean } | null>(null)
  const pointer = useRef({ x: 0, y: 0, flip: false })
  const statusRef = useRef(onStatus); statusRef.current = onStatus
  const apiRef = useRef(onApi); apiRef.current = onApi
  const hooksRef = useRef(onHooks); hooksRef.current = onHooks
  const setStatus = (s: 'loading' | 'ready' | 'error') => { setStatusState(s); statusRef.current?.(s) }

  const key = entries.map((e) => `${e.seqName}:${e.color}`).join(',')
  const extraKey = JSON.stringify(entries.map((e) => (extra[e.seqName] ? Array.from(extra[e.seqName]!) : null)))
  const focusKey = JSON.stringify(focus)
  const emphasised = emphasis && entries.length > 1 && entries.some((e) => e.seqName === emphasis) ? emphasis : null

  // ---------------------------------------------------------------- plugin (one per mount)
  useEffect(() => {
    const el = host.current
    if (!el) return
    let alive = true
    const p = new PDBeMolstarPlugin()
    // a fresh mount point per plugin: React's StrictMode remount must not call createRoot() again on the element of the torn-down plugin
    const mount = document.createElement('div'); mount.style.position = 'absolute'; mount.style.inset = '0'
    el.appendChild(mount)
    queue.current = queue.current.then(async () => {
      if (!alive) return
      await p.render(mount, { hideControls: true, sequencePanel: false, bgColor: WHITE, pdbeLink: false, expanded: false, loadingOverlay: false, subscribeEvents: false,
        hideCanvasControls: ['expand', 'controlToggle', 'controlInfo', 'animation', 'trajectory', 'selection'] })
      p.plugin.canvas3d?.setProps({ camera: { manualReset: true } })
      await p.canvas.setBgColor(WHITE)
      if (!alive) return
      plugin.current = p
      // residue hover label from Mol*'s hover loci (author chain ids in full, also multi-letter ones)
      p.plugin.behaviors.interaction.hover.subscribe(({ current }) => {
        const l = current.loci
        if (!StructureElement.Loci.is(l) || StructureElement.Loci.isEmpty(l)) { setHover(null); return }
        const loc = StructureElement.Loci.getFirstLocation(l)
        if (!loc) { setHover(null); return }
        const pdb = loc.unit.model.entryId
        setHover({ text: residueLabel(pdb, { chain: StructureProperties.chain.auth_asym_id(loc), resn: StructureProperties.atom.auth_comp_id(loc), resi: String(StructureProperties.residue.auth_seq_id(loc)), icode: StructureProperties.residue.pdbx_PDB_ins_code(loc) || '' }), ...pointer.current })
      })
      setGen((g) => g + 1)
    }).catch((e: unknown) => { if (alive) { setErr(String(e)); setStatus('error') } })
    const onFull = () => { setIsFull(document.fullscreenElement === outer.current); setTimeout(() => { try { p.plugin.canvas3d?.handleResize(); void fit() } catch { /* torn down */ } }, 80) }
    document.addEventListener('fullscreenchange', onFull)
    // Mol* only follows window resizes: follow the container too (panels that grow with a drill-down, stacked layouts)
    const ro = new ResizeObserver(() => { try { p.plugin?.canvas3d?.handleResize() } catch { /* torn down */ } })
    ro.observe(el)
    return () => {
      alive = false; document.removeEventListener('fullscreenchange', onFull); ro.disconnect()
      plugin.current = null; nodes.current.clear(); ligNodes.current.clear(); apiRef.current?.(null)
      queue.current = queue.current.catch(() => undefined).then(() => (p.plugin ? p.clear() : undefined)).catch(() => undefined).then(() => { try { p.plugin?.dispose() } catch { /* never rendered */ } mount.remove() })
    }
  }, [])   // eslint-disable-line react-hooks/exhaustive-deps

  // ---------------------------------------------------------------- helpers
  const structureData = (ref: string): Structure | undefined => plugin.current?.plugin.state.data.cells.get(ref)?.obj?.data as Structure | undefined
  const chainExpr = (asym: string) => MS.struct.generator.atomGroups({ 'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), asym]) })
  /** cartoon of a chain: emphasised = thicker; faded (another structure is emphasised) or residues focused = semi-transparent */
  const cartoonParams = (color: string, mode: 'normal' | 'emphasis' | 'faded') => createStructureRepresentationParams(plugin.current!.plugin, undefined, {
    type: 'cartoon', typeParams: { sizeFactor: mode === 'emphasis' ? 0.32 : 0.22, alpha: focusRef.current ? 0.3 : mode === 'faded' ? 0.35 : 1 }, color: 'uniform', colorParams: { value: hexColor(color) } })
  const modeOf = (sn: string) => (!emphasisRef.current ? 'normal' : sn === emphasisRef.current ? 'emphasis' : 'faded') as 'normal' | 'emphasis' | 'faded'
  const emphasisRef = useRef(emphasised); emphasisRef.current = emphasised
  const matrixFor = (e: SceneEntry) => { const o = orientationMat(e.transform), x = extraRef.current[e.seqName]; return x ? Mat4.mul(Mat4(), x, o) : o }
  const extraRef = useRef(extra); extraRef.current = extra

  /** loci of the given residues (or the whole chain) in each structure's displayed (transformed) coordinates */
  const lociOf = (sel?: Record<string, number[]> | null) => [...nodes.current.values()].flatMap((n) => {
    const data = structureData(n.transform); if (!data) return []
    const seqs = sel ? sel[n.entry.seqName] : undefined
    if (sel && (!seqs || !seqs.length)) return []
    const expr = seqs ? MS.struct.generator.atomGroups({ 'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), n.entry.labelAsymId]),
      'residue-test': MS.core.set.has([MS.set(...seqs), MS.struct.atomProperty.macromolecular.label_seq_id()]) }) : chainExpr(n.entry.labelAsymId)
    const loci = StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(expr, data))
    return StructureElement.Loci.isEmpty(loci) ? [] : [loci]
  })
  /** fixed view direction in the family frame (as the Family Overview); the camera fits the given loci with a small margin */
  const fit = async (sel?: Record<string, number[]> | null) => {
    const p = plugin.current, ctx = p?.plugin
    if (!ctx?.canvas3d) return
    const all = lociOf(sel)
    if (!all.length) return
    const pts: Vec3[] = []
    for (const l of all) for (const e of l.elements) { const v = Vec3(); OrderedSet.forEach(e.indices, (i) => { e.unit.conformation.position(e.unit.elements[i], v); pts.push(Vec3.clone(v)) }) }
    const c = Vec3(); for (const v of pts) Vec3.add(c, c, v); Vec3.scale(c, c, 1 / pts.length)
    let r = 0; for (const v of pts) r = Math.max(r, Vec3.distance(c, v))
    const cam = ctx.canvas3d.camera, vp = cam.viewport, aspect = vp.width > 0 && vp.height > 0 ? vp.width / vp.height : 1
    const half = Math.min(cam.state.fov / 2, Math.atan(Math.tan(cam.state.fov / 2) * aspect))
    const radius = Math.max(r * (sel ? 1.6 : 1.08), sel ? 9 : 5)
    const distance = radius / Math.sin(half)
    const dir = Vec3.normalize(Vec3(), VIEW_DIR)
    await PluginCommands.Camera.SetSnapshot(ctx, { snapshot: { target: c, position: Vec3.sub(Vec3(), c, Vec3.scale(Vec3(), dir, distance)), up: Vec3.clone(VIEW_UP), radius, radiusMax: radius }, durationMs: sel ? 250 : 0 })
  }
  const exposeHooks = () => {
    hooksRef.current?.({
      loaded: [...nodes.current.keys()], fetches: sceneFetchCount(),
      transformed: [...nodes.current.values()].map((n) => isApplicable(n.entry.transform)),
      // screen position (page coordinates) and label of the i-th residue's O3' atom (Mol*'s nucleic-acid cartoon trace atom, so the point lies on the tube)
      probe: (seqName: string, i: number) => {
        const n = nodes.current.get(seqName), ctx = plugin.current?.plugin, data = n && structureData(n.transform)
        if (!n || !ctx?.canvas3d || !data) return null
        const loci = StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(MS.struct.generator.atomGroups({
          'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), n.entry.labelAsymId]), 'atom-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_atom_id(), "O3'"]) }), data))
        let k = 0, found: { pos: Vec3; label: string } | null = null
        StructureElement.Loci.forEachLocation(loci, (loc) => {
          if (found || k++ !== i) return
          const pos = Vec3.create(StructureProperties.atom.x(loc), StructureProperties.atom.y(loc), StructureProperties.atom.z(loc))
          found = { pos, label: residueLabel(n.entry.pdbId, { chain: StructureProperties.chain.auth_asym_id(loc), resn: StructureProperties.atom.auth_comp_id(loc), resi: String(StructureProperties.residue.auth_seq_id(loc)), icode: StructureProperties.residue.pdbx_PDB_ins_code(loc) || '' }) }
        })
        if (!found) return null
        const f = found as { pos: Vec3; label: string }
        const out = ctx.canvas3d.camera.project(Vec4(), f.pos)
        const canvas = host.current?.querySelector('canvas'), rect = canvas?.getBoundingClientRect()
        if (!canvas || !rect) return null
        const sx = rect.width / ctx.canvas3d.camera.viewport.width, sy = rect.height / ctx.canvas3d.camera.viewport.height
        return { x: rect.left + out[0] * sx, y: rect.top + (ctx.canvas3d.camera.viewport.height - out[1]) * sy, label: f.label }
      } })
  }

  // ---------------------------------------------------------------- scene: add / remove chains
  useEffect(() => {
    const p = plugin.current
    if (!p) return
    let alive = true
    setStatus('loading'); setErr(null); setHover(null)
    queue.current = queue.current.catch(() => undefined).then(async () => {
      const ctx = p.plugin
      const want = new Map(entries.map((e) => [e.seqName, e]))
      for (const [sn, n] of [...nodes.current]) {
        if (!want.has(sn) || want.get(sn)!.color !== n.entry.color) {
          await PluginCommands.State.RemoveObject(ctx, { state: ctx.state.data, ref: n.data }); nodes.current.delete(sn)
          const l = ligNodes.current.get(sn); if (l) { await PluginCommands.State.RemoveObject(ctx, { state: ctx.state.data, ref: l.data }); ligNodes.current.delete(sn) }
        }
      }
      const bytes = await Promise.all(entries.map((e) => (nodes.current.has(e.seqName) ? null : chainBytes(e))))
      if (!alive) return
      for (const [k, e] of entries.entries()) {
        if (nodes.current.has(e.seqName)) continue
        const b = ctx.build().toRoot()
        const data = b.apply(StateTransforms.Data.RawData, { data: bytes[k]!, label: e.label })
        const structure = data.apply(StateTransforms.Data.ParseCif).apply(StateTransforms.Model.TrajectoryFromMmCif).apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 }).apply(StateTransforms.Model.StructureFromModel)
        const transform = structure.apply(StateTransforms.Model.TransformStructureConformation, { transform: { name: 'matrix' as const, params: { data: matrixFor(e), transpose: false } } })
        const component = transform.apply(StateTransforms.Model.StructureSelectionFromExpression, { expression: chainExpr(e.labelAsymId), label: e.label })
        const repr = component.apply(StateTransforms.Representation.StructureRepresentation3D, cartoonParams(e.color, modeOf(e.seqName)))
        await b.commit()
        nodes.current.set(e.seqName, { entry: e, data: data.ref, structure: structure.ref, transform: transform.ref, component: component.ref, repr: repr.ref, focus: null })
      }
      if (!alive) return
      await applyFocus()
      await applyLigands()
      await fit(focusRef.current)
      exposeHooks()
      apiRef.current?.({ residueAtoms })
      setLoadedKey(key); setStatus('ready')
    }).catch((e: unknown) => { if (alive) { setErr(String((e as Error)?.message ?? e)); setStatus('error') } })
    return () => { alive = false }
  }, [key, gen])   // eslint-disable-line react-hooks/exhaustive-deps

  /** atoms of one residue in the family frame (orientation only): read from the untransformed model, then mapped with the orientation */
  const residueAtoms = (seqName: string, labelSeqId: number): SceneAtom[] => {
    const n = nodes.current.get(seqName), data = n && structureData(n.structure)
    if (!n || !data) return []
    const loci = StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(MS.struct.generator.atomGroups({
      'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), n.entry.labelAsymId]), 'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_seq_id(), labelSeqId]) }), data))
    const o = orientationMat(n.entry.transform), out: SceneAtom[] = []
    StructureElement.Loci.forEachLocation(loci, (loc) => {
      if (StructureProperties.atom.type_symbol(loc) === 'H') return
      const v = Vec3.transformMat4(Vec3(), Vec3.create(StructureProperties.atom.x(loc), StructureProperties.atom.y(loc), StructureProperties.atom.z(loc)), o)
      out.push({ name: StructureProperties.atom.label_atom_id(loc), x: v[0], y: v[1], z: v[2], seq: labelSeqId, resn: StructureProperties.atom.label_comp_id(loc) })
    })
    return out
  }

  // ---------------------------------------------------------------- transforms (extra rigid fit): update the matrices only
  useEffect(() => {
    const p = plugin.current
    if (!p || status !== 'ready') return
    queue.current = queue.current.catch(() => undefined).then(async () => {
      const b = p.plugin.state.data.build()
      for (const n of nodes.current.values()) {
        b.to(n.transform).update({ transform: { name: 'matrix' as const, params: { data: matrixFor(n.entry), transpose: false } } })
        const l = ligNodes.current.get(n.entry.seqName); if (l) b.to(l.transform).update({ transform: { name: 'matrix' as const, params: { data: matrixFor(n.entry), transpose: false } } })
      }
      await b.commit()
      await fit(focusRef.current)
      exposeHooks()
    })
  }, [extraKey])   // eslint-disable-line react-hooks/exhaustive-deps

  // ---------------------------------------------------------------- emphasis: representation parameters only
  useEffect(() => {
    const p = plugin.current
    if (!p || status !== 'ready') return
    queue.current = queue.current.catch(() => undefined).then(async () => {
      const b = p.plugin.state.data.build()
      for (const n of nodes.current.values()) b.to(n.repr).update(cartoonParams(n.entry.color, modeOf(n.entry.seqName)))
      await b.commit()
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), [`${testId}Emphasis`]: emphasisRef.current }
    })
  }, [emphasised, status])   // eslint-disable-line react-hooks/exhaustive-deps

  // ---------------------------------------------------------------- focused residues: sticks in each structure's colour
  const focusRef = useRef(focus); focusRef.current = focus
  const applyFocus = async () => {
    const p = plugin.current
    if (!p) return
    const ctx = p.plugin
    for (const n of nodes.current.values()) {
      if (n.focus) { await PluginCommands.State.RemoveObject(ctx, { state: ctx.state.data, ref: n.focus }); n.focus = null }
      const seqs = focusRef.current?.[n.entry.seqName]
      setSubtreeVisibility(ctx.state.data, n.component, false)   // visible (third argument = hidden)
      if (!seqs?.length) continue
      const b = ctx.build()
      const sel = b.to(n.transform).apply(StateTransforms.Model.StructureSelectionFromExpression, { expression: MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), n.entry.labelAsymId]),
        'residue-test': MS.core.set.has([MS.set(...seqs), MS.struct.atomProperty.macromolecular.label_seq_id()]) }), label: `${n.entry.label} focus` })
      sel.apply(StateTransforms.Representation.StructureRepresentation3D, createStructureRepresentationParams(ctx, undefined, { type: 'ball-and-stick', typeParams: { sizeFactor: 0.3 }, color: 'uniform', colorParams: { value: hexColor(n.entry.color) } }))
      await b.commit()
      n.focus = sel.ref
    }
    const b = ctx.state.data.build()
    for (const n of nodes.current.values()) b.to(n.repr).update(cartoonParams(n.entry.color, modeOf(n.entry.seqName)))
    await b.commit()
  }
  useEffect(() => {
    const p = plugin.current
    if (!p || status !== 'ready') return
    queue.current = queue.current.catch(() => undefined).then(async () => { await applyFocus(); await fit(focusRef.current); exposeHooks() })
  }, [focusKey])   // eslint-disable-line react-hooks/exhaustive-deps

  // Mol*'s own canvas buttons and highlight box are hidden (class selectors above): this viewer has its own controls and residue label
  // the entries change one render before the scene effect sets 'loading': never report 'ready' for entries that are not loaded yet
  const shownStatus = status === 'ready' && loadedKey !== key ? 'loading' : status
  // ---------------------------------------------------------------- Show ligands
  /** small-molecule residues (not ions, oligosaccharides or lipids, not the excluded classes) with any atom within LIGAND_RADIUS of the chain */
  const nearChainLigands = (asym: string, exclude: string[]) => MS.struct.filter.within({
    0: MS.struct.generator.atomGroups({
      'entity-test': MS.core.logic.and([MS.core.rel.eq([MS.ammp('entityType'), 'non-polymer']), MS.core.logic.not([MS.core.str.match([MS.re('(ion|oligosaccharide|lipid)', 'i'), MS.ammp('entitySubtype')])])]),
      ...(exclude.length ? { 'residue-test': MS.core.logic.not([MS.core.set.has([MS.set(...exclude), MS.ammp('label_comp_id')])]) } : {}),
      'group-by': MS.ammp('residueKey') }),
    target: chainExpr(asym), 'max-radius': LIGAND_RADIUS })
  /** fetch and build the ligands of every displayed chain once (only after Show ligands was ticked), then show or hide them */
  const applyLigands = async () => {
    const p = plugin.current
    if (!p) return
    const ctx = p.plugin, show = showLigRef.current
    if (show) {
      const missing = [...nodes.current.values()].filter((n) => !ligNodes.current.has(n.entry.seqName))
      if (missing.length) {
        setLigBusy(true)
        try {
          classes.current ??= (classesRef.current?.() ?? Promise.resolve({})).catch(() => ({}))
          const cls = await classes.current
          const exclude = Object.entries(cls).filter(([, c]) => LIGAND_EXCLUDED_CLASSES.includes(c)).map(([k]) => k)
          const bytes = await Promise.all(missing.map((n) => surroundingsBytes(n.entry, LIGAND_RADIUS)))
          for (const [k, n] of missing.entries()) {
            if (plugin.current !== p || !nodes.current.has(n.entry.seqName)) continue
            const b = ctx.build().toRoot()
            const data = b.apply(StateTransforms.Data.RawData, { data: bytes[k], label: `${n.entry.label} ligands` })
            const transform = data.apply(StateTransforms.Data.ParseCif).apply(StateTransforms.Model.TrajectoryFromMmCif).apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 })
              .apply(StateTransforms.Model.StructureFromModel).apply(StateTransforms.Model.TransformStructureConformation, { transform: { name: 'matrix' as const, params: { data: matrixFor(n.entry), transpose: false } } })
            const expression = nearChainLigands(n.entry.labelAsymId, exclude)
            transform.apply(StateTransforms.Model.StructureSelectionFromExpression, { expression, label: `${n.entry.label} ligands` })
              // sticks; carbons in the structure's colour, other elements CPK (as the Family Overview viewer)
              .apply(StateTransforms.Representation.StructureRepresentation3D, createStructureRepresentationParams(ctx, undefined, {
                type: 'ball-and-stick', typeParams: { sizeFactor: 0.22, sizeAspectRatio: 1 }, color: 'element-symbol',
                colorParams: { carbonColor: { name: 'uniform', params: { value: hexColor(n.entry.color), saturation: 0, lightness: 0 } }, lightness: 0 } }))
            await b.commit()
            const st = structureData(transform.ref), shown = new Set<string>()
            if (st) StructureElement.Loci.forEachLocation(StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(expression, st)), (l) => {
              shown.add(`${StructureProperties.atom.label_comp_id(l)} ${StructureProperties.chain.auth_asym_id(l)}${StructureProperties.residue.auth_seq_id(l)}${StructureProperties.residue.pdbx_PDB_ins_code(l) || ''}`) })
            ligNodes.current.set(n.entry.seqName, { data: data.ref, transform: transform.ref, shown: [...shown].sort() })
          }
        } finally { setLigBusy(false) }
      }
      setLigInfo(Object.fromEntries([...nodes.current.keys()].map((sn) => [sn, ligNodes.current.get(sn)?.shown ?? []])))
    }
    for (const l of ligNodes.current.values()) setSubtreeVisibility(ctx.state.data, l.data, !show)
    // test hook: Show ligands state, chains whose ligands are built, and the ligand stick style
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), [`${testId}Ligands`]: { show, built: [...ligNodes.current.keys()],
      visible: [...ligNodes.current.keys()].filter(() => show), shown: Object.fromEntries([...ligNodes.current].map(([k, l]) => [k, l.shown])) } }
  }
  useEffect(() => {
    const p = plugin.current
    if (!p || status !== 'ready') return
    queue.current = queue.current.catch(() => undefined).then(() => applyLigands()).catch((e: unknown) => { setErr(String((e as Error)?.message ?? e)) })
  }, [showLigands])   // eslint-disable-line react-hooks/exhaustive-deps

  const toggleFull = () => { const el = outer.current; if (!el) return; if (document.fullscreenElement) void document.exitFullscreen?.(); else void el.requestFullscreen?.() }
  const recenter = () => { queue.current = queue.current.catch(() => undefined).then(() => fit(focusRef.current)) }
  const shownLig = ligInfo ? entries.map((e) => [e.seqName, ligInfo[e.seqName]] as const).filter(([, v]) => v) : []
  const ligCodes = [...new Set(shownLig.flatMap(([, v]) => v!.map((x) => x.split(' ')[0])))]
  const noneNearby = showLigands && !ligBusy && ligInfo !== null && shownLig.length > 0 && ligCodes.length === 0
  const ligTitle = `Show ligand instances within ${LIGAND_RADIUS} Å of displayed RNAs as sticks. Excludes water, ions and crystallisation additives.`
    + (showLigands && ligCodes.length ? `\nShown: ${ligCodes.join(', ')}` : '')
  return (
    <div ref={outer} className={`relative flex flex-col overflow-hidden rounded border border-stone-200 bg-white [&_.msp-highlight-toast-wrapper]:hidden [&_.msp-viewport-controls]:hidden ${className}`} style={style} data-testid={testId} data-status={shownStatus}
      data-structures={entries.map((e) => e.seqName).join(',')} data-colors={entries.map((e) => e.color).join(',')} data-emphasis={emphasised ?? ''} data-engine="molstar" {...dataAttrs}>
      {/* controls in a toolbar above the canvas (as the Family Overview viewer), never over the structures */}
      <div className="flex shrink-0 flex-wrap items-center gap-1.5 border-b border-stone-200 bg-stone-50 px-2 py-1.5 text-xs" data-testid="viewer-toolbar">
        <label className="flex items-center gap-1 rounded border border-stone-300 bg-white px-2 py-0.5 text-stone-700" title={ligTitle} data-testid="viewer-show-ligands-control"
          data-shown={shownLig.map(([sn, v]) => `${sn}:${v!.join('|')}`).join(';')}>
          <input type="checkbox" checked={showLigands} disabled={!entries.length || shownStatus !== 'ready'} onChange={(e) => setShowLigands(e.target.checked)} data-testid="viewer-show-ligands" />
          Show ligands{ligBusy && <span className="text-[10px] text-stone-500"> (loading…)</span>}{noneNearby && <span className="text-[10px] text-stone-500" data-testid="viewer-ligands-none"> (none nearby)</span>}
        </label>
        {shownStatus === 'loading' && <span className="px-1 py-0.5 text-stone-600">Loading…</span>}
        <span className="ml-auto flex gap-1">
          <button type="button" className="rounded border border-stone-300 bg-white px-2 py-0.5 text-stone-700 hover:bg-stone-50" onClick={recenter} data-testid="viewer-centre">Centre view</button>
          <button type="button" className="rounded border border-stone-300 bg-white px-2 py-0.5 text-stone-700 hover:bg-stone-50" onClick={toggleFull} data-testid="viewer-fullscreen">{isFull ? 'Exit fullscreen' : 'Fullscreen'}</button>
        </span>
      </div>
      <div className="relative min-h-0 flex-1" data-testid={`${testId}-canvas-area`}>
        <div ref={host} className="absolute inset-0 bg-white" onMouseMove={(e) => { const r = e.currentTarget.getBoundingClientRect(); pointer.current = { x: e.clientX - r.left, y: e.clientY - r.top, flip: e.clientX - r.left > r.width - 200 } }} onMouseLeave={() => setHover(null)} />
        {shownStatus === 'loading' && <div className="pointer-events-none absolute inset-0 z-[5] bg-white" />}
        {hover && <div role="tooltip" className="pointer-events-none absolute z-20 whitespace-nowrap rounded border border-stone-300 bg-white/95 px-1.5 py-0.5 font-mono text-[11px] text-stone-800 shadow" style={{ left: hover.flip ? undefined : hover.x + 12, right: hover.flip ? `calc(100% - ${hover.x - 12}px)` : undefined, top: hover.y + 12 }} data-testid="viewer-residue-hover">{hover.text}</div>}
        {children}
        {err && <div className="absolute inset-x-2 bottom-2 z-10 rounded bg-red-50 px-2 py-1 text-xs text-red-800" data-testid="viewer-error">{err}</div>}
        {!entries.length && <div className="absolute inset-0 z-10 flex items-center justify-center text-sm text-stone-500">No structure displayed.</div>}
      </div>
    </div>
  )
}
