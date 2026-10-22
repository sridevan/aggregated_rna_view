import { useEffect, useRef, useState } from 'react'
import { PDBeMolstarPlugin } from 'pdbe-molstar/lib/viewer'
import { StateTransforms } from 'molstar/lib/mol-plugin-state/transforms'
import { createStructureRepresentationParams } from 'molstar/lib/mol-plugin-state/helpers/structure-representation-params'
import { StructureElement, StructureProperties } from 'molstar/lib/mol-model/structure'
import { PluginCommands } from 'molstar/lib/mol-plugin/commands'
import { Vec3 } from 'molstar/lib/mol-math/linear-algebra'
import { InfoTip } from '../InfoTip'
import 'pdbe-molstar/build/pdbe-molstar-light.css'
import { checkIdealSdf, IDEAL_INFO, idealSdfUrls, principalAxes } from '../../lib/idealLigand'

/*
 * Ideal 3D structure of one chemical component (CCD), on demand in a modal. Coordinates: the CCD ideal conformer as SDF from RCSB
 * (files.rcsb.org/ligands/download/{CCD}_ideal.sdf, the documented download of the RCSB ligand pages); PDBe's identical
 * pdbechem_v2/{CCD}_ideal.sdf only when the RCSB request fails. Both allow cross-origin requests (checked 2026-10-22 for all 72 gallery
 * CCDs). Model (experimental) coordinates are never substituted. The CCD records ideal coordinates as all zeros when they are missing
 * (e.g. MIY, TAC): reported as unavailable, not rendered. One Mol* plugin per open modal, disposed on close; a new CCD clears the
 * previous structure first and a request token drops stale loads.
 */
/** Fetch the ideal SDF (RCSB, then PDBe on a network or HTTP failure). Returns the text and its source, or null when both fail. */
async function fetchIdealSdf(ccd: string, signal: AbortSignal): Promise<{ text: string; url: string } | null> {
  for (const url of idealSdfUrls(ccd)) {
    try { const r = await fetch(url, { signal }); if (r.ok) return { text: await r.text(), url } } catch (e) { if (signal.aborted) throw e }
  }
  return null
}

type W = { __rnaExplorer?: Record<string, unknown> }
/** centred and fitted with a margin, looking along the smallest principal axis with the largest one horizontal (a planar molecule is seen face-on) */
function fitView(ctx: PDBeMolstarPlugin['plugin'], xyz: number[][] | null, durationMs: number) {
  const cam = ctx.canvas3d?.camera
  if (!xyz || !cam) return
  const { centre, axes } = principalAxes(xyz)
  const r = Math.max(1.5, ...xyz.map((p) => Math.hypot(p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]))) + 0.8   // + atom radius margin
  const vp = cam.viewport, aspect = vp.width > 0 && vp.height > 0 ? vp.width / vp.height : 1
  const half = Math.min(cam.state.fov / 2, Math.atan(Math.tan(cam.state.fov / 2) * aspect))
  const dir = Vec3.create(axes[2][0], axes[2][1], axes[2][2]), target = Vec3.create(centre[0], centre[1], centre[2])
  const radius = r, distance = radius / Math.sin(half)
  void PluginCommands.Camera.SetSnapshot(ctx, { snapshot: { target, position: Vec3.sub(Vec3(), target, Vec3.scale(Vec3(), dir, distance)), up: Vec3.create(axes[1][0], axes[1][1], axes[1][2]), radius, radiusMax: radius }, durationMs })
}
const WHITE = { r: 255, g: 255, b: 255 }

export function LigandIdealViewer({ ccd, name, pageUrl, onClose }: { ccd: string; name: string; pageUrl: string; onClose: () => void }) {
  const host = useRef<HTMLDivElement>(null), closeBtn = useRef<HTMLButtonElement>(null), panel = useRef<HTMLDivElement>(null)
  const plugin = useRef<PDBeMolstarPlugin | null>(null), coords = useRef<number[][] | null>(null)
  const queue = useRef<Promise<unknown>>(Promise.resolve())
  const [ready, setReady] = useState(false)
  const [status, setStatus] = useState<{ state: 'loading' } | { state: 'ready'; nAtoms: number; source: string } | { state: 'error'; message: string }>({ state: 'loading' })
  const [hover, setHover] = useState<string | null>(null)

  // modal behaviour: focus inside on open, Escape closes, focus returns to the opener, page scroll locked while open (set up once per open)
  const closeRef = useRef(onClose); closeRef.current = onClose
  useEffect(() => {
    const opener = document.activeElement as HTMLElement | null
    closeBtn.current?.focus()
    const onKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') { e.stopPropagation(); closeRef.current() }
      if (e.key === 'Tab' && panel.current) {   // keep keyboard focus inside the dialog
        const f = [...panel.current.querySelectorAll<HTMLElement>('a[href], button:not([disabled]), [tabindex]:not([tabindex="-1"])')]
        if (f.length && e.shiftKey && document.activeElement === f[0]) { f[f.length - 1].focus(); e.preventDefault() }
        else if (f.length && !e.shiftKey && document.activeElement === f[f.length - 1]) { f[0].focus(); e.preventDefault() }
      }
    }
    document.addEventListener('keydown', onKey)
    const overflow = document.body.style.overflow; document.body.style.overflow = 'hidden'
    return () => { document.removeEventListener('keydown', onKey); document.body.style.overflow = overflow; opener?.focus?.() }
  }, [])

  // one plugin per open modal (created on demand, disposed on close)
  useEffect(() => {
    const el = host.current; if (!el) return
    let alive = true
    const p = new PDBeMolstarPlugin()
    const mount = document.createElement('div'); mount.style.position = 'absolute'; mount.style.inset = '0'; el.appendChild(mount)
    queue.current = queue.current.then(async () => {
      if (!alive) return
      await p.render(mount, { hideControls: true, sequencePanel: false, bgColor: WHITE, pdbeLink: false, expanded: false, loadingOverlay: false, subscribeEvents: false,
        hideCanvasControls: ['expand', 'controlToggle', 'controlInfo', 'animation', 'trajectory', 'selection'] })
      await p.canvas.setBgColor(WHITE)
      if (!alive) return
      plugin.current = p
      // atom element and number on hover
      p.plugin.behaviors.interaction.hover.subscribe(({ current }) => {
        const l = current.loci
        const loc = StructureElement.Loci.is(l) && !StructureElement.Loci.isEmpty(l) ? StructureElement.Loci.getFirstLocation(l) : undefined
        // SDF has no atom names: element symbol and the atom's number in the file
        setHover(loc ? `${StructureProperties.atom.type_symbol(loc)} · atom ${StructureProperties.atom.id(loc)}` : null)
      })
      setReady(true)
    }).catch((e: unknown) => { if (alive) setStatus({ state: 'error', message: `3D viewer could not start: ${String((e as Error)?.message ?? e)}` }) })
    const ro = new ResizeObserver(() => { try { p.plugin?.canvas3d?.handleResize() } catch { /* torn down */ } }); ro.observe(el)
    return () => {
      alive = false; ro.disconnect(); plugin.current = null
      queue.current = queue.current.catch(() => undefined).then(() => (p.plugin ? p.clear() : undefined)).catch(() => undefined)
        .then(() => { try { p.plugin?.dispose() } catch { /* never rendered */ } mount.remove() })
      const w = window as unknown as W; if (w.__rnaExplorer) w.__rnaExplorer.ligandViewer = null
    }
  }, [])

  // load the selected CCD: previous structure cleared first; a newer request or closing aborts this one
  useEffect(() => {
    if (!ready) return
    const ac = new AbortController()
    setStatus({ state: 'loading' }); setHover(null)
    const expose = (v: Record<string, unknown>) => { const w = window as unknown as W; w.__rnaExplorer = { ...(w.__rnaExplorer ?? {}), ligandViewer: { ccd, ...v } } }
    expose({ status: 'loading' })
    queue.current = queue.current.then(async () => {
      const p = plugin.current; if (!p || ac.signal.aborted) return
      await p.clear(); coords.current = null
      const got = await fetchIdealSdf(ccd, ac.signal)
      if (ac.signal.aborted) return
      if (!got) { setStatus({ state: 'error', message: `Ideal coordinates for ${ccd} could not be retrieved from RCSB or PDBe. Check the connection and try again.` }); expose({ status: 'error', reason: 'network' }); return }
      const chk = checkIdealSdf(got.text)
      if (!chk.ok) {
        setStatus({ state: 'error', message: chk.reason === 'no_ideal_coordinates'
          ? `The Chemical Component Dictionary has no ideal coordinates for ${ccd}, so no ideal 3D structure is shown.`
          : `The ideal coordinate file for ${ccd} could not be read.` })
        expose({ status: 'error', reason: chk.reason }); return
      }
      const ctx = p.plugin
      const b = ctx.build().toRoot()
      b.apply(StateTransforms.Data.RawData, { data: got.text, label: `${ccd} ideal` })
        .apply(StateTransforms.Model.TrajectoryFromSDF).apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 }).apply(StateTransforms.Model.StructureFromModel)
        .apply(StateTransforms.Representation.StructureRepresentation3D, createStructureRepresentationParams(ctx, undefined, { type: 'ball-and-stick', color: 'element-symbol', colorParams: { carbonColor: { name: 'element-symbol', params: {} } } }))
      await b.commit()
      if (ac.signal.aborted) { await p.clear(); return }
      coords.current = chk.xyz
      ctx.canvas3d?.setProps({ camera: { manualReset: true } })
      fitView(ctx, chk.xyz, 0)
      const source = got.url.includes('rcsb.org') ? 'RCSB' : 'PDBe'
      setStatus({ state: 'ready', nAtoms: chk.nAtoms, source }); expose({ status: 'ready', nAtoms: chk.nAtoms, source: got.url })
    }).catch((e: unknown) => { if (!ac.signal.aborted) { setStatus({ state: 'error', message: `The ideal 3D structure of ${ccd} could not be displayed: ${String((e as Error)?.message ?? e)}` }); expose({ status: 'error', reason: 'render' }) } })
    return () => ac.abort()
  }, [ccd, ready])

  const reset = () => { const p = plugin.current; if (p) fitView(p.plugin, coords.current, 250) }
  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-stone-900/40 p-4" onMouseDown={(e) => { if (e.target === e.currentTarget) onClose() }} data-testid="ligand-viewer-backdrop">
      <div ref={panel} role="dialog" aria-modal="true" aria-labelledby="ligand-viewer-title" className="flex max-h-full w-full max-w-3xl flex-col rounded border border-stone-300 bg-white shadow-xl" data-testid="ligand-viewer" data-ccd={ccd} data-status={status.state}>
        <header className="flex items-start gap-2 border-b border-stone-200 px-4 py-2.5">
          <div className="min-w-0 flex-1">
            <h2 id="ligand-viewer-title" className="line-clamp-2 break-words text-lg font-semibold leading-snug text-stone-900" title={`${ccd} — ${name}`} data-testid="ligand-viewer-title"><span className="font-mono">{ccd}</span> — {name}</h2>
            <p className="flex items-center gap-1 text-xs italic text-stone-600" data-testid="ligand-viewer-subtitle">Ideal 3D coordinates from the PDB Chemical Component Dictionary<InfoTip label="About ideal coordinates" text={IDEAL_INFO} /></p>
          </div>
          <button ref={closeBtn} type="button" onClick={onClose} aria-label="Close the ideal 3D structure" className="rounded border border-stone-300 bg-white px-2 py-0.5 text-xs text-stone-700 hover:bg-stone-50 focus:outline-none focus:ring-2 focus:ring-blue-500" data-testid="ligand-viewer-close">Close</button>
        </header>
        <div className="relative h-[min(60vh,480px)] w-full bg-white" data-testid="ligand-viewer-canvas">
          <div ref={host} className="absolute inset-0" />
          {status.state === 'loading' && <div className="absolute inset-0 flex items-center justify-center bg-white/70 text-sm text-stone-600" data-testid="ligand-viewer-loading">Loading ideal coordinates…</div>}
          {status.state === 'error' && <div className="absolute inset-0 flex items-center justify-center bg-white p-6 text-center text-sm text-stone-700" role="alert" data-testid="ligand-viewer-error">{status.message}</div>}
          {hover && status.state === 'ready' && <div className="pointer-events-none absolute left-2 top-2 rounded border border-stone-200 bg-white/90 px-1.5 py-0.5 text-[11px] text-stone-800" data-testid="ligand-viewer-hover">{hover}</div>}
        </div>
        <footer className="flex flex-wrap items-center gap-2 border-t border-stone-200 px-4 py-2 text-xs text-stone-600">
          <button type="button" onClick={reset} disabled={status.state !== 'ready'} className="rounded border border-stone-300 bg-white px-2 py-0.5 text-stone-700 hover:bg-stone-50 disabled:opacity-40" data-testid="ligand-viewer-reset">Reset view</button>
          <span className="text-[11px] text-stone-500">Drag to rotate · scroll to zoom · right-drag to move{status.state === 'ready' ? ` · ${status.nAtoms} atoms (${status.source})` : ''}</span>
          <a href={pageUrl} target="_blank" rel="noreferrer" className="ml-auto text-blue-700 hover:underline" data-testid="ligand-viewer-link">PDBe chemical component page ↗</a>
        </footer>
      </div>
    </div>
  )
}
