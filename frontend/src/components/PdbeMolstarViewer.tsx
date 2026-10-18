import { useEffect, useRef, useState } from 'react'
import { PDBeMolstarPlugin } from 'pdbe-molstar/lib/viewer'
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder'
import { StateTransforms } from 'molstar/lib/mol-plugin-state/transforms'
import { createStructureRepresentationParams } from 'molstar/lib/mol-plugin-state/helpers/structure-representation-params'
import { PluginCommands } from 'molstar/lib/mol-plugin/commands'
import { setSubtreeVisibility } from 'molstar/lib/mol-plugin/behavior/static/state'
import { StructureFocusRepresentation } from 'molstar/lib/mol-plugin/behavior/dynamic/selection/structure-focus-representation'
import { Script } from 'molstar/lib/mol-script/script'
import { StructureSelection, StructureElement, StructureProperties } from 'molstar/lib/mol-model/structure'
import { Vec3 } from 'molstar/lib/mol-math/linear-algebra'
import { Color } from 'molstar/lib/mol-util/color'
import { OrderedSet } from 'molstar/lib/mol-data/int'
import type { Structure } from 'molstar/lib/mol-model/structure'
import type { OrientationTransform } from '../types/backend'
import { isApplicable, toMat4 } from '../lib/orientation'
import 'pdbe-molstar/build/pdbe-molstar-light.css'

export interface MolstarSelectionItem { label_seq_id: number }
export type ViewMode = 'region' | 'complete'
interface Props {
  pdbId: string
  labelAsymId: string
  authAsymId: string
  region?: { start: number; end: number }
  selection: MolstarSelectionItem[]
  /** centre the camera on the selection whenever it changes (single-residue / group exploration); the view direction is kept */
  focusOnSelection?: boolean
  /** precomputed rigid-body transform into the family reference frame (orientation.json); null/insufficient = shown as deposited */
  orientation?: OrientationTransform | null
  onResidueClick?: (labelSeqId: number) => void
  /** changes on every explicit deselection: clear highlight and focus, restore the default whole-structure (region) camera */
  resetKey?: number
  /** one fixed colour (hex) for the polymer / Rfam-region representation of every structure; undefined = PDBe Mol* default colouring */
  structureColor?: number
  /** "Clear selection" in the viewer controls: clears the shared selection (the parent then restores the default 2D and 3D views) */
  onClearSelection?: () => void
  canClear?: boolean
  className?: string
  initialMode?: ViewMode
}
interface MolstarClickEvent extends Event { eventData?: { seq_id?: number; auth_seq_id?: number; auth_asym_id?: string; entity_id?: string; ins_code?: string } }
/** selection (nucleotide or whole loop): violet #7C3AED, distinct from the slate structure colour; identical for every instance */
const SELECT_HEX = 0x7c3aed
const SELECT_COLOR = { r: 0x7c, g: 0x3a, b: 0xed }
const WHITE = { r: 255, g: 255, b: 255 }
const REGION_LABEL = 'Rfam-mapped region'
const ORIENT_TAG = 'rna-explorer-orientation'
/** simplified display: the focus representation shows only the selected nucleotides (ball-and-stick); no 5 A surroundings, no contact lines */
const FOCUS_COMPONENTS = ['target'] as const
/** water is hidden (visibility only: the component and the coordinates stay in the structure) */
const isWater = (c: { key?: string; cell: { obj?: { label?: string } } }) => c.key === 'water' || c.cell.obj?.label === 'Water'
const VIEW_UP = Vec3.create(0, 1, 0), VIEW_DIR = Vec3.create(0, 0, -1)   // fixed camera orientation in the reference frame

/**
 * Reusable PDBe Mol* wrapper. One plugin per mount; structures are switched with `visual.update(..., fullLoad=false)` (load the new
 * entry, then remove the old one) so the canvas and its white background are never re-initialised. A precomputed rigid-body transform
 * (TransformStructureConformation, inserted below the structure node) moves each structure into the family reference frame; the
 * camera is then set from a fixed view direction, so all instances appear in the same orientation. Selection colours residues; with
 * `focusOnSelection` the camera is also centred on them (direction kept, context margin), and clearing never moves the camera.
 */
export function PdbeMolstarViewer({ pdbId, labelAsymId, authAsymId, region, selection, focusOnSelection = false, orientation, onResidueClick, resetKey = 0, structureColor, onClearSelection, canClear = false, className = '', initialMode = 'region' }: Props) {
  const host = useRef<HTMLDivElement>(null)
  const outer = useRef<HTMLDivElement>(null)
  const plugin = useRef<PDBeMolstarPlugin | null>(null)
  const queue = useRef<Promise<void>>(Promise.resolve())
  const regionRef = useRef<string | null>(null)
  const [status, setStatus] = useState<'loading' | 'ready' | 'error'>('loading')
  const [loadedPdb, setLoadedPdb] = useState<string | null>(null)
  const [pluginGen, setPluginGen] = useState(0)        // bumped when a plugin instance exists: drives the load effect deterministically
  const [retry, setRetry] = useState(0)
  const [errorMsg, setErrorMsg] = useState<string | null>(null)
  const [mode, setMode] = useState<ViewMode>(region ? initialMode : 'complete')
  const clickRef = useRef(onResidueClick); clickRef.current = onResidueClick
  const authRef = useRef(authAsymId); authRef.current = authAsymId
  const selRef = useRef<{ label_asym_id: string; label_seq_id: number; color: typeof SELECT_COLOR }[]>([])
  const applied = isApplicable(orientation) ? orientation : null

  const structureOf = (p: PDBeMolstarPlugin) => {
    const all = p.plugin.managers.structure.hierarchy.current.structures
    return [...all].sort((a, b) => b.components.length - a.components.length)[0]
  }
  const visibleStructureData = (p: PDBeMolstarPlugin): Structure | undefined => {
    const s = structureOf(p)
    return (s?.transform?.cell.obj?.data as Structure | undefined) ?? (s?.cell.obj?.data as Structure | undefined)
  }

  /** Insert (or update) the conformation transform below the structure node, mirroring Mol*'s own superposition helper. */
  const applyOrientation = async (p: PDBeMolstarPlugin, t: OrientationTransform | null) => {
    const ctx = p.plugin
    const s = structureOf(p)
    if (!s) return
    const existing = ctx.state.data.selectQ((q) => q.byRef(s.cell.transform.ref).subtree().withTransformer(StateTransforms.Model.TransformStructureConformation))[0]
    const m = t ? toMat4(t) : null
    if (!m) { if (existing) await PluginCommands.State.RemoveObject(ctx, { state: ctx.state.data, ref: existing.transform.ref }); return }
    const params = { transform: { name: 'matrix' as const, params: { data: m, transpose: false } } }
    const b = existing ? ctx.state.data.build().to(existing).update(params) : ctx.state.data.build().to(s.cell.transform.ref).insert(StateTransforms.Model.TransformStructureConformation, params, { tags: ORIENT_TAG })
    await ctx.runTask(ctx.state.data.updateTree(b))
  }

  const applyMode = async (p: PDBeMolstarPlugin, m: ViewMode) => {
    const ctx = p.plugin
    for (let i = 0; i < 60 && ctx.managers.structure.hierarchy.current.structures.length === 0; i++) await new Promise((r) => setTimeout(r, 50))
    const s = structureOf(p)
    if (!s) return
    const comps = () => ctx.managers.structure.hierarchy.current.structures.flatMap((x) => x.components)
    if (regionRef.current && !ctx.state.data.cells.has(regionRef.current)) regionRef.current = null   // region of a structure that has been removed
    if (m === 'region' && region) {
      if (!regionRef.current) {
        const expr = MS.struct.generator.atomGroups({
          'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), labelAsymId]),
          'residue-test': MS.core.rel.inRange([MS.struct.atomProperty.macromolecular.label_seq_id(), region.start, region.end]),
        })
        const parentRef = s.transform?.cell.transform.ref ?? s.cell.transform.ref   // below the orientation transform so the region is transformed too
        // the structure (or its transformed copy) may not have resolved yet right after a switch: wait for it rather than skipping the region
        let structureData = visibleStructureData(p)
        for (let i = 0; i < 60 && !structureData; i++) { await new Promise((r) => setTimeout(r, 50)); structureData = visibleStructureData(p) }
        if (structureData) {
          const update = ctx.build()
          const sel = update.to(parentRef).apply(StateTransforms.Model.StructureSelectionFromExpression, { expression: expr, label: REGION_LABEL })
          // cartoon only (Mol* draws the bases as rings/ladders): the selected nucleotide then stands out in its atom-level focus representation
          sel.apply(StateTransforms.Representation.StructureRepresentation3D, createStructureRepresentationParams(ctx, structureData, { type: 'cartoon', ...uniformColor() }))
          await update.commit()
          regionRef.current = sel.ref
        } else {
          console.warn('[PdbeMolstarViewer] region component could not be built; showing the complete structure')
        }
      }
      // never leave an empty scene: without a region component the default components stay visible (Mol* drops camera updates when nothing is visible)
      for (const c of comps()) setSubtreeVisibility(ctx.state.data, c.cell.transform.ref, regionRef.current ? c.cell.transform.ref !== regionRef.current : isWater(c))
    } else {
      for (const c of comps()) setSubtreeVisibility(ctx.state.data, c.cell.transform.ref, c.cell.transform.ref === regionRef.current || isWater(c))
    }
    // the same fixed colour on the polymer component (Complete structure) for every instance; ligands, ions and other components keep their colouring
    if (structureColor != null) {
      const polymers = comps().filter((c) => c.key === 'polymer' || c.cell.obj?.label === 'Polymer')
      if (polymers.length) await ctx.managers.structure.component.updateRepresentationsTheme(polymers, { color: 'uniform', colorParams: { value: Color(structureColor) } })
    }
    const regionCell = regionRef.current ? ctx.state.data.cells.get(regionRef.current) : undefined
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), molstarColours: representationColours(p) }
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), molstarComponents: [...comps().map((c) => ({ label: c.cell.obj?.label, hidden: c.cell.state.isHidden })), { label: REGION_LABEL, hidden: regionCell?.state.isHidden ?? null, present: !!regionCell }] }
    await resetView(p, m)
  }

  /**
   * Deterministic, tight camera: the visible atoms (region or whole chain) are projected onto the fixed view axes (dir, up) of the
   * reference frame; the camera is centred on the projected bounding box and placed at the distance at which that box fills the
   * canvas (given the camera fov and the actual canvas aspect ratio) with a small margin. No arbitrary zoom factor; hidden
   * components never enter the calculation; the view direction is identical for every instance.
   */
  const resetView = async (p: PDBeMolstarPlugin, m: ViewMode) => {
    const ctx = p.plugin
    const data = visibleStructureData(p)
    if (!data || !ctx.canvas3d) return
    const expr = m === 'region' && region
      ? MS.struct.generator.atomGroups({ 'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), labelAsymId]), 'residue-test': MS.core.rel.inRange([MS.struct.atomProperty.macromolecular.label_seq_id(), region.start, region.end]) })
      : MS.struct.generator.all()
    const loci = StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(expr, data))
    const sphere = StructureElement.Loci.getBoundary(loci).sphere
    // view axes
    const dir = Vec3.normalize(Vec3(), VIEW_DIR)
    const right = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), dir, VIEW_UP))
    const up = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), right, dir))
    // projected coordinates of every visible atom in the view frame (x: right, y: up, z: towards the camera = -dir)
    const xs: number[] = [], ys: number[] = [], zs: number[] = []
    const v = Vec3(), d = Vec3()
    for (const e of loci.elements) {
      const unit = e.unit
      OrderedSet.forEach(e.indices, (i) => {
        unit.conformation.position(unit.elements[i], v); Vec3.sub(d, v, sphere.center)
        xs.push(Vec3.dot(d, right)); ys.push(Vec3.dot(d, up)); zs.push(-Vec3.dot(d, dir))
      })
    }
    const n = xs.length
    if (!n) return
    const cam = ctx.canvas3d.camera
    const fov = cam.state.fov
    const vp = cam.viewport
    const aspect = vp.width > 0 && vp.height > 0 ? vp.width / vp.height : 1
    const tanHalf = Math.tan(fov / 2)
    const MARGIN = 1.06                                   // 6 % breathing space so sticks at the silhouette are never clipped
    const OVERLAY_PX = 40                                 // top controls / bottom badge overlay the canvas: keep the molecule out of those bands
    const usableH = Math.max(vp.height - 2 * OVERLAY_PX, vp.height * 0.5)
    const tanHalfY = tanHalf * usableH / vp.height        // effective vertical half-angle of the usable band
    const xmin = Math.min(...xs), xmax = Math.max(...xs), ymin = Math.min(...ys), ymax = Math.max(...ys)
    const cx = (xmin + xmax) / 2, cy = (ymin + ymax) / 2
    // exact perspective fit: an atom at (x, y) displaced z towards the camera is inside the frustum when |y| <= (d - z)·tan(fov/2)
    // and |x| <= (d - z)·tan(fov/2)·aspect, so the camera distance is the maximum over atoms of the required d
    let distance = 0
    for (let i = 0; i < n; i++) {
      const need = Math.max(Math.abs(ys[i] - cy) * MARGIN / tanHalfY, Math.abs(xs[i] - cx) * MARGIN / (tanHalf * aspect)) + zs[i]
      if (need > distance) distance = need
    }
    distance = Math.max(distance, 5)
    const target = Vec3.scaleAndAdd(Vec3(), Vec3.scaleAndAdd(Vec3(), Vec3.clone(sphere.center), right, cx), up, cy)
    const position = Vec3.sub(Vec3(), target, Vec3.scale(Vec3(), dir, distance))
    // Camera.Reset re-derives the distance from the bounding radius (Camera.targetDistance), discarding the position; SetSnapshot applies ours as is
    const snapshot = { target, position, up: Vec3.clone(VIEW_UP), radius: Math.max(sphere.radius, 5), radiusMax: Math.max(sphere.radius, 5) }
    const halfW = (xmax - xmin) / 2, halfH = (ymax - ymin) / 2
    // camera resets requested by the loader are resolved on the next frames; let them pass, apply ours, then verify once
    const frames = (n: number) => new Promise<void>((res) => { const step = (k: number) => (k <= 0 ? res() : requestAnimationFrame(() => step(k - 1))); step(n) })
    await frames(2)
    await PluginCommands.Camera.SetSnapshot(ctx, { snapshot, durationMs: 0 })
    await frames(2)
    const st = ctx.canvas3d.camera.state
    if (Math.abs(Vec3.distance(st.position, st.target) - distance) > 0.01 * distance) await PluginCommands.Camera.SetSnapshot(ctx, { snapshot, durationMs: 0 })
    const occupancy = { width: halfW / (distance * tanHalf * aspect), height: halfH / (distance * tanHalf) }   // projected extent at the box centre depth / visible extent
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), lastView: { center: [...sphere.center], radius: sphere.radius, up: [...VIEW_UP], dir: [...VIEW_DIR], atoms: n, distance, aspect, fov, occupancy } }
  }

  useEffect(() => {
    const el = host.current, wrap = outer.current
    if (!el || !wrap) return
    let alive = true
    let p: PDBeMolstarPlugin | null = null
    const handler = (e: Event) => {
      const d = (e as MolstarClickEvent).eventData
      if (!d || d.seq_id === undefined || d.auth_asym_id !== authRef.current) return
      clickRef.current?.(d.seq_id)
    }
    wrap.addEventListener('PDB.molstar.click', handler)
    const timer = window.setTimeout(() => {
      if (!alive) return
      p = new PDBeMolstarPlugin()
      plugin.current = p
      const inst = p
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), molstar: inst, getCamera: () => { const c = inst.plugin.canvas3d?.camera.state; return c ? { position: [...c.position], target: [...c.target], up: [...c.up], radius: c.radius } : null } }
      setPluginGen((g) => g + 1)
    }, 0)
    return () => {
      alive = false
      window.clearTimeout(timer)
      wrap.removeEventListener('PDB.molstar.click', handler)
      const inst = p
      plugin.current = null
      queue.current = queue.current.catch(() => undefined).then(() => inst?.clear()).catch(() => undefined).then(() => { el.innerHTML = '' })
    }
  }, [])

  // load / switch structure without re-initialising the canvas, then orientation + mode
  useEffect(() => {
    const p = plugin.current, el = host.current
    if (!p || !el) return
    let alive = true
    setStatus('loading'); setErrorMsg(null)
    // watchdog: never leave an indefinite blank viewer
    const watchdog = window.setTimeout(() => { if (alive) { setErrorMsg(`Loading ${pdbId.toUpperCase()} timed out after 60 s`); setStatus('error') } }, 60_000)
    queue.current = queue.current.catch(() => undefined).then(async () => {
      if (!alive || plugin.current !== p) return
      regionRef.current = null   // reset inside the queue: a still-running previous load must not leave its region ref behind
      if (loadedPdb === null && !el.querySelector('canvas')) {
        await p.render(el, { moleculeId: pdbId, subscribeEvents: true, hideControls: true, sequencePanel: false, bgColor: WHITE, pdbeLink: false, expanded: false, loadingOverlay: false,
          hideCanvasControls: ['expand', 'controlToggle', 'controlInfo', 'animation', 'trajectory'] })
        await p.canvas.setBgColor(WHITE)
        p.plugin.canvas3d?.setProps({ camera: { manualReset: true } })   // no automatic camera resets on state changes: the view is set explicitly below
        await p.canvas.setBgColor(WHITE)                                   // setProps may re-derive renderer defaults: re-assert the white background
        await setFocusComponents(p)
      } else {
        const ctx = p.plugin
        const old = ctx.managers.structure.hierarchy.current.trajectories.map((t) => t.cell.transform.ref)
        p.plugin.managers.structure.focus.clear()                     // no focus representation of the previous structure survives the switch
        await p.visual.update({ moleculeId: pdbId }, false)          // fullLoad=false: no plugin.clear(), canvas settings untouched
        const oldRefs = ctx.managers.structure.hierarchy.current.trajectories.filter((t) => old.includes(t.cell.transform.ref))
        if (oldRefs.length) await ctx.managers.structure.hierarchy.remove(oldRefs, false)
        await p.canvas.setBgColor(WHITE)
      }
      await applyOrientation(p, applied)
      await applyMode(p, region ? mode : 'complete')
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), structuresLoaded: p.plugin.managers.structure.hierarchy.current.structures.length }
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), orientationApplied: applied ? { pdbId, rmsd: applied.rmsd, n_pairs: applied.n_pairs, fit_quality: applied.fit_quality } : null, bgColor: p.plugin.canvas3d?.props.renderer.backgroundColor }
      if (alive && plugin.current === p) { window.clearTimeout(watchdog); setLoadedPdb(pdbId); setStatus('ready') }
    }).catch((err: unknown) => { if (alive) { console.error('[PdbeMolstarViewer] load', err); window.clearTimeout(watchdog); setErrorMsg(`Could not load ${pdbId.toUpperCase()}: ${err instanceof Error ? err.message : String(err)}`); setStatus('error') } })
    return () => { alive = false; window.clearTimeout(watchdog) }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pdbId, labelAsymId, region?.start, region?.end, applied, pluginGen, retry])

  useEffect(() => {
    const p = plugin.current
    if (!p || status !== 'ready' || loadedPdb !== pdbId) return
    queue.current = queue.current.catch(() => undefined).then(() => applyMode(p, mode)).catch((err: unknown) => console.error('[PdbeMolstarViewer] mode', err))
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [mode])

  const key = selection.map((s) => s.label_seq_id).sort((a, b) => a - b).join(',')
  useEffect(() => {
    const p = plugin.current
    if (!p || status !== 'ready' || loadedPdb !== pdbId) return
    let alive = true
    queue.current = queue.current.catch(() => undefined).then(async () => {
      if (!alive || plugin.current !== p) return
      await p.visual.clearSelection()
      const data = key ? key.split(',').map((n) => ({ label_asym_id: labelAsymId, label_seq_id: +n, color: SELECT_COLOR })) : []
      selRef.current = data
      // atom-level close-up: Mol*'s structure focus restricted to the selected nucleotides (ball-and-stick, no surroundings or contact lines; see
      // FOCUS_COMPONENTS) replaces any previous focus, so repeated selections never accumulate representations; cleared when nothing is selected
      const loci = data.length ? lociFor(p, data) : null
      if (loci) p.plugin.managers.structure.focus.setFromLoci(loci); else p.plugin.managers.structure.focus.clear()
      if (data.length) await p.visual.select({ data })
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), lastMolstarSelection: { pdbId, labelAsymId, labelSeqIds: data.map((d) => d.label_seq_id), color: SELECT_HEX }, molstarMode: mode, focusRepresentation: focusedSeqIds(p) }
      if (loci && focusOnSelection) await focusOn(p, data)   // clearing never touches the camera
    }).catch((err: unknown) => console.error('[PdbeMolstarViewer] selection', err))
    return () => { alive = false }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [key, labelAsymId, pdbId, status, loadedPdb])

  // explicit deselection (same nucleotide/loop clicked again, Clear selection): no highlight, no focus, default camera of the current mode
  const lastReset = useRef(resetKey)
  useEffect(() => {
    if (resetKey === lastReset.current) return
    lastReset.current = resetKey
    const p = plugin.current
    if (!p || status !== 'ready' || loadedPdb !== pdbId) return
    queue.current = queue.current.catch(() => undefined).then(async () => {
      if (plugin.current !== p) return
      p.plugin.managers.structure.focus.clear()
      await p.visual.clearSelection()
      selRef.current = []
      await resetView(p, region ? mode : 'complete')
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), focusRepresentation: focusedSeqIds(p), molstarReset: Number(window.__rnaExplorer?.molstarReset ?? 0) + 1 }
    }).catch((err: unknown) => console.error('[PdbeMolstarViewer] reset', err))
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [resetKey])

  const uniformColor = () => (structureColor != null ? { color: 'uniform' as const, colorParams: { value: Color(structureColor) } } : {})
  /** colour theme of every representation in the scene (test hook): parent component label, theme name, uniform value, visibility */
  const representationColours = (p: PDBeMolstarPlugin) => p.plugin.managers.structure.hierarchy.current.structures.flatMap((st) => st.components.flatMap((c) => c.representations.map((r) => {
    const ct = (r.cell.transform.params as { colorTheme?: { name: string; params?: { value?: number } } } | undefined)?.colorTheme
    return { component: c.cell.obj?.label ?? '', theme: ct?.name ?? null, value: ct?.params?.value ?? null, hidden: !!c.cell.state.isHidden }
  })))

  /** Restrict Mol*'s focus representation to the target (selected nucleotides, heavy atoms): no surroundings, no non-covalent contact lines. */
  const setFocusComponents = async (p: PDBeMolstarPlugin) => {
    const ctx = p.plugin
    const cell = ctx.state.behaviors.selectQ((q) => q.ofTransformer(StructureFocusRepresentation))[0]
    if (cell) await ctx.state.updateBehavior(StructureFocusRepresentation, (old) => {
      old.components = [...FOCUS_COMPONENTS]; old.ignoreHydrogens = true            // heavy atoms only: less clutter on large loops
      // contrasting target: uniform dark selection colour, slightly thicker ball-and-stick than Mol*'s default
      old.targetParams = createStructureRepresentationParams(ctx, undefined, { type: 'ball-and-stick', typeParams: { sizeFactor: 0.32 }, color: 'uniform', colorParams: { value: Color(SELECT_HEX) } })
    })
    const now = ctx.state.behaviors.selectQ((q) => q.ofTransformer(StructureFocusRepresentation))[0]
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), focusComponents: (now?.transform.params as { components?: string[] } | undefined)?.components ?? null, focusIgnoreHydrogens: (now?.transform.params as { ignoreHydrogens?: boolean } | undefined)?.ignoreHydrogens ?? null }
  }

  /** Loci of the given residues in the displayed (transformed) structure; null when none of them is present (no coordinates). */
  const lociFor = (p: PDBeMolstarPlugin, items: { label_asym_id: string; label_seq_id: number }[]) => {
    const data = visibleStructureData(p)
    if (!data || !items.length) return null
    const expr = MS.struct.generator.atomGroups({
      'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), items[0].label_asym_id]),
      'residue-test': MS.core.set.has([MS.set(...items.map((i) => i.label_seq_id)), MS.struct.atomProperty.macromolecular.label_seq_id()]),
    })
    const loci = StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(expr, data))
    return StructureElement.Loci.isEmpty(loci) ? null : loci
  }
  /** label_seq_ids of the residues in Mol*'s current focus representation (test hook). */
  const focusedSeqIds = (p: PDBeMolstarPlugin): number[] => {
    const cur = p.plugin.managers.structure.focus.current
    if (!cur) return []
    const ids = new Set<number>()
    StructureElement.Loci.forEachLocation(cur.loci, (l) => { ids.add(StructureProperties.residue.label_seq_id(l)) })
    return [...ids].sort((a, b) => a - b)
  }
  /** Close-up of the selected residues, keeping the current view direction: the nucleotide fills the view with its immediate neighbours as context. */
  const focusOn = async (p: PDBeMolstarPlugin, items: { label_asym_id: string; label_seq_id: number }[]) => {
    const ctx = p.plugin
    const loci = lociFor(p, items)
    if (!loci || !ctx.canvas3d) return
    const sphere = StructureElement.Loci.getBoundary(loci).sphere
    // focusLoci keeps the camera direction/up and moves target + distance (smooth transition); the margin keeps the 5 A surroundings in view
    ctx.managers.camera.focusLoci(loci, { extraRadius: 7, minRadius: 12, durationMs: 250 })
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), lastFocus: items.map((d) => d.label_seq_id), lastFocusCenter: [...sphere.center], lastFocusRadius: sphere.radius }
  }
  const retryLoad = () => {
    const p = plugin.current
    if (!p) return
    queue.current = queue.current.catch(() => undefined).then(() => p.clear()).catch(() => undefined)
    setLoadedPdb(null); setRetry((r) => r + 1)
  }
  const resetOrientation = () => { const p = plugin.current; if (p) queue.current = queue.current.catch(() => undefined).then(() => resetView(p, region ? mode : 'complete')) }

  return (
    <div ref={outer} className={`relative overflow-hidden rounded border border-stone-200 bg-white ${className}`} data-testid="molstar-viewer" data-pdb={status === 'ready' ? loadedPdb ?? pdbId : pdbId} data-status={status} data-mode={region ? mode : 'complete'} data-oriented={applied ? 'true' : 'false'}>
      <div ref={host} className="absolute inset-0 bg-white" />
      {status !== 'ready' && <div className="pointer-events-none absolute inset-0 z-[5] bg-white" data-testid="molstar-loading-overlay" />}
      <div className="absolute left-2 top-2 z-10 flex flex-wrap items-center gap-1 text-xs">
        {region && (
          <div className="flex overflow-hidden rounded border border-stone-300 bg-white/90" role="group" aria-label="View mode">
            <button data-testid="mode-region" onClick={() => setMode('region')} aria-pressed={mode === 'region'} className={`px-2 py-0.5 ${mode === 'region' ? 'bg-blue-600 text-white' : 'text-stone-700'}`}>Rfam region</button>
            <button data-testid="mode-complete" onClick={() => setMode('complete')} aria-pressed={mode === 'complete'} className={`px-2 py-0.5 ${mode === 'complete' ? 'bg-blue-600 text-white' : 'text-stone-700'}`}>Complete structure</button>
          </div>
        )}
        <button data-testid="reset-orientation" onClick={resetOrientation} className="rounded border border-stone-300 bg-white/90 px-2 py-0.5 text-stone-700">Reset orientation / Fit</button>
        {onClearSelection && <button data-testid="molstar-clear-selection" onClick={onClearSelection} disabled={!canClear} className="rounded border border-stone-300 bg-white/90 px-2 py-0.5 text-stone-700 disabled:opacity-40">Clear selection</button>}
        {status === 'loading' && <span className="rounded bg-white/90 px-2 py-0.5 text-stone-600">Loading {pdbId.toUpperCase()}…</span>}
        {status === 'error' && <span className="flex items-center gap-2 rounded bg-red-50 px-2 py-0.5 text-red-700" data-testid="molstar-error">{errorMsg ?? `Could not load ${pdbId.toUpperCase()}`}<button data-testid="molstar-retry" onClick={retryLoad} className="rounded border border-red-300 bg-white px-1.5 py-0.5 text-red-700">Retry</button></span>}
      </div>
    </div>
  )
}
