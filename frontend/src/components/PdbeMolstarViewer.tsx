import { useEffect, useRef, useState } from 'react'
import { PDBeMolstarPlugin } from 'pdbe-molstar/lib/viewer'
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder'
import { StateTransforms } from 'molstar/lib/mol-plugin-state/transforms'
import { createStructureRepresentationParams } from 'molstar/lib/mol-plugin-state/helpers/structure-representation-params'
import { PluginCommands } from 'molstar/lib/mol-plugin/commands'
import { setSubtreeVisibility } from 'molstar/lib/mol-plugin/behavior/static/state'
import { Script } from 'molstar/lib/mol-script/script'
import { StructureSelection, StructureElement } from 'molstar/lib/mol-model/structure'
import { Vec3 } from 'molstar/lib/mol-math/linear-algebra'
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
  /** precomputed rigid-body transform into the family reference frame (orientation.json); null/insufficient = shown as deposited */
  orientation?: OrientationTransform | null
  onResidueClick?: (labelSeqId: number) => void
  className?: string
  initialMode?: ViewMode
}
interface MolstarClickEvent extends Event { eventData?: { seq_id?: number; auth_seq_id?: number; auth_asym_id?: string; entity_id?: string; ins_code?: string } }
const SELECT_COLOR = { r: 98, g: 80, b: 214 }
const WHITE = { r: 255, g: 255, b: 255 }
const REGION_LABEL = 'Rfam-mapped region'
const ORIENT_TAG = 'rna-explorer-orientation'
const VIEW_UP = Vec3.create(0, 1, 0), VIEW_DIR = Vec3.create(0, 0, -1)   // fixed camera orientation in the reference frame

/**
 * Reusable PDBe Mol* wrapper. One plugin per mount; structures are switched with `visual.update(..., fullLoad=false)` (load the new
 * entry, then remove the old one) so the canvas and its white background are never re-initialised. A precomputed rigid-body transform
 * (TransformStructureConformation, inserted below the structure node) moves each structure into the family reference frame; the
 * camera is then set from a fixed view direction, so all instances appear in the same orientation. Selection never moves the camera.
 */
export function PdbeMolstarViewer({ pdbId, labelAsymId, authAsymId, region, selection, orientation, onResidueClick, className = '', initialMode = 'region' }: Props) {
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
    if (m === 'region' && region) {
      if (!regionRef.current) {
        const expr = MS.struct.generator.atomGroups({
          'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), labelAsymId]),
          'residue-test': MS.core.rel.inRange([MS.struct.atomProperty.macromolecular.label_seq_id(), region.start, region.end]),
        })
        const parentRef = s.transform?.cell.transform.ref ?? s.cell.transform.ref   // below the orientation transform so the region is transformed too
        const structureData = visibleStructureData(p)
        if (structureData) {
          const update = ctx.build()
          const sel = update.to(parentRef).apply(StateTransforms.Model.StructureSelectionFromExpression, { expression: expr, label: REGION_LABEL })
          sel.apply(StateTransforms.Representation.StructureRepresentation3D, createStructureRepresentationParams(ctx, structureData, { type: 'cartoon' }))
          sel.apply(StateTransforms.Representation.StructureRepresentation3D, createStructureRepresentationParams(ctx, structureData, { type: 'ball-and-stick', typeParams: { sizeFactor: 0.2 } }))
          await update.commit()
          regionRef.current = sel.ref
        }
      }
      for (const c of comps()) setSubtreeVisibility(ctx.state.data, c.cell.transform.ref, c.cell.transform.ref !== regionRef.current)
    } else {
      for (const c of comps()) setSubtreeVisibility(ctx.state.data, c.cell.transform.ref, c.cell.transform.ref === regionRef.current)
    }
    const regionCell = regionRef.current ? ctx.state.data.cells.get(regionRef.current) : undefined
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
    regionRef.current = null
    // watchdog: never leave an indefinite blank viewer
    const watchdog = window.setTimeout(() => { if (alive) { setErrorMsg(`Loading ${pdbId.toUpperCase()} timed out after 60 s`); setStatus('error') } }, 60_000)
    queue.current = queue.current.catch(() => undefined).then(async () => {
      if (!alive || plugin.current !== p) return
      if (loadedPdb === null && !el.querySelector('canvas')) {
        await p.render(el, { moleculeId: pdbId, subscribeEvents: true, hideControls: true, sequencePanel: false, bgColor: WHITE, pdbeLink: false, expanded: false, loadingOverlay: false,
          hideCanvasControls: ['expand', 'controlToggle', 'controlInfo', 'animation', 'trajectory'] })
        await p.canvas.setBgColor(WHITE)
        p.plugin.canvas3d?.setProps({ camera: { manualReset: true } })   // no automatic camera resets on state changes: the view is set explicitly below
        await p.canvas.setBgColor(WHITE)                                   // setProps may re-derive renderer defaults: re-assert the white background
      } else {
        const ctx = p.plugin
        const old = ctx.managers.structure.hierarchy.current.trajectories.map((t) => t.cell.transform.ref)
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
      if (data.length) await p.visual.select({ data })
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), lastMolstarSelection: { pdbId, labelAsymId, labelSeqIds: data.map((d) => d.label_seq_id) }, molstarMode: mode }
    }).catch((err: unknown) => console.error('[PdbeMolstarViewer] selection', err))
    return () => { alive = false }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [key, labelAsymId, pdbId, status, loadedPdb])

  const focusSelection = () => {
    const p = plugin.current
    if (!p || !selRef.current.length) return
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), lastFocus: selRef.current.map((d) => d.label_seq_id) }
    queue.current = queue.current.catch(() => undefined).then(() => p.visual.focus(selRef.current))
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
        <button data-testid="focus-selection" onClick={focusSelection} disabled={!selection.length} className="rounded border border-stone-300 bg-white/90 px-2 py-0.5 text-stone-700 disabled:opacity-40">Focus selection</button>
        {status === 'loading' && <span className="rounded bg-white/90 px-2 py-0.5 text-stone-600">Loading {pdbId.toUpperCase()}…</span>}
        {status === 'error' && <span className="flex items-center gap-2 rounded bg-red-50 px-2 py-0.5 text-red-700" data-testid="molstar-error">{errorMsg ?? `Could not load ${pdbId.toUpperCase()}`}<button data-testid="molstar-retry" onClick={retryLoad} className="rounded border border-red-300 bg-white px-1.5 py-0.5 text-red-700">Retry</button></span>}
      </div>
      <div className="absolute bottom-2 left-2 z-10 rounded bg-white/90 px-2 py-0.5 text-[11px] text-stone-600" data-testid="orientation-badge">
        {orientation == null ? 'Orientation: as deposited (no transform available)'
          : applied ? `Reference frame: superposed on the family reference, RMSD ${applied.rmsd?.toFixed(2)} Å over ${applied.n_pairs} C4′ pairs (${applied.fit_quality})${applied.is_identity ? ' · reference structure' : ''}`
          : `Orientation: as deposited (fit ${orientation.fit_quality}: ${orientation.n_pairs} shared positions with coordinates)`}
      </div>
    </div>
  )
}
