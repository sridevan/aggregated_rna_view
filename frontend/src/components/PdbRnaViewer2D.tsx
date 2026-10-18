import { useEffect, useRef, useState } from 'react'

interface Props {
  pdbId: string
  entityId: number
  /** author chain id (the PDBe RNA viewer and its topology data are keyed by author chain) */
  authAsymId: string
  /** label_seq_ids to show as selected (the viewer's API is label_seq_id based) */
  selection: number[]
  onResidueClick?: (labelSeqId: number) => void
  /** changes on every explicit deselection: clear the highlight and restore the default fitted view (zoom and pan reset) */
  resetKey?: number
  className?: string
}
interface RnaViewerEvent extends Event { eventData?: { pdbId: string; entityId: string; label_seq_ids?: number[]; label_seq_id?: number } }

/** Make the whole topology visible: drop the viewer's initial zoom transform, set the viewBox to the drawing's bounding box (with padding)
 * and let the SVG scale to the container (the viewer sets fixed pixel dimensions and a d3 zoom transform sized for its own container). */
function fitSvg(el: HTMLElement): void {
  const svg = el.querySelector('svg.rnaTopoSvg') as SVGSVGElement | null
  if (!svg) return
  for (const g of svg.querySelectorAll<SVGGElement>('g[transform]')) g.removeAttribute('transform')
  let b: DOMRect
  try { b = svg.getBBox() } catch { return }
  if (!b.width || !b.height) return
  const pad = Math.max(b.width, b.height) * 0.08
  svg.setAttribute('viewBox', `${b.x - pad} ${b.y - pad} ${b.width + 2 * pad} ${b.height + 2 * pad}`)
  svg.setAttribute('preserveAspectRatio', 'xMidYMid meet')
  svg.style.width = '100%'; svg.style.height = '100%'
  svg.removeAttribute('width'); svg.removeAttribute('height')
  // keep d3-zoom's stored transform consistent with the removed group transforms (identity), so the next wheel/drag starts from the fitted view
  const z = (svg as unknown as { __zoom?: { constructor: new (k: number, x: number, y: number) => unknown } }).__zoom
  if (z) (svg as unknown as { __zoom: unknown }).__zoom = new z.constructor(1, 0, 0)
}

/** Chains longer than this get the selection brought into focus in 2D; shorter ones keep the whole fitted topology (already readable). */
export const FOCUS_MIN_NT = 200
/** Minimum focus window, in median nucleotide spacings, so a single nucleotide is shown with its helix and surroundings, not alone. */
const FOCUS_MIN_SPAN = 26
const FOCUS_MS = 350

export interface Box { x: number; y: number; w: number; h: number }
/** Window that frames `sel` (bounding box of the selected nucleotides) with context: at least `minSide` on each side, plus 25% padding
 *  around larger selections, centred on the selection and never larger than the whole drawing `full`. */
export function focusBox(sel: Box, minSide: number, full: Box): Box {
  const w = Math.min(full.w, Math.max(minSide, sel.w * 1.25)), h = Math.min(full.h, Math.max(minSide, sel.h * 1.25))
  return { x: sel.x + sel.w / 2 - w / 2, y: sel.y + sel.h / 2 - h / 2, w, h }
}
export const contains = (outer: Box, inner: Box) => inner.x >= outer.x && inner.y >= outer.y && inner.x + inner.w <= outer.x + outer.w && inner.y + inner.h <= outer.y + outer.h
const parseVb = (svg: SVGSVGElement): Box | null => { const v = (svg.getAttribute('viewBox') ?? '').split(/\s+/).map(Number); return v.length === 4 && v.every(Number.isFinite) ? { x: v[0], y: v[1], w: v[2], h: v[3] } : null }
const setVb = (svg: SVGSVGElement, b: Box) => svg.setAttribute('viewBox', `${b.x} ${b.y} ${b.w} ${b.h}`)

/** Drawing-space region visible now: the viewBox, undone by any d3-zoom transform the user applied with wheel/drag (on the outermost group). */
function visibleRegion(svg: SVGSVGElement): Box | null {
  const vb = parseVb(svg); if (!vb) return null
  const g = svg.querySelector<SVGGElement>('g[transform]')
  const m = g?.transform.baseVal.consolidate()?.matrix
  if (!m || !m.a) return vb
  return { x: (vb.x - m.e) / m.a, y: (vb.y - m.f) / m.d, w: vb.w / m.a, h: vb.h / m.d }
}

/** Bounding box (drawing space) of the given label_seq_ids, and the median spacing between consecutive nucleotides (the drawing's length scale). */
function nucleotideGeometry(svg: SVGSVGElement, pdbId: string, ids: number[]): { sel: Box | null; spacing: number } {
  const centre = (id: number): Box | null => {
    const els = svg.getElementsByClassName(`rnaview_${pdbId}_${id}`)
    let b: Box | null = null
    for (const el of els) {
      let r: DOMRect; try { r = (el as SVGGraphicsElement).getBBox() } catch { continue }
      if (!r.width && !r.height) continue
      b = b ? { x: Math.min(b.x, r.x), y: Math.min(b.y, r.y), w: Math.max(b.x + b.w, r.x + r.width) - Math.min(b.x, r.x), h: Math.max(b.y + b.h, r.y + r.height) - Math.min(b.y, r.y) } : { x: r.x, y: r.y, w: r.width, h: r.height }
    }
    return b
  }
  let sel: Box | null = null
  for (const id of ids) {
    const b = centre(id); if (!b) continue
    sel = sel ? { x: Math.min(sel.x, b.x), y: Math.min(sel.y, b.y), w: Math.max(sel.x + sel.w, b.x + b.w) - Math.min(sel.x, b.x), h: Math.max(sel.y + sel.h, b.y + b.h) - Math.min(sel.y, b.y) } : b
  }
  // spacing from up to 40 consecutive pairs around the selection
  const d: number[] = []
  const start = Math.max(1, (ids[0] ?? 1) - 20)
  let prev = centre(start)
  for (let i = start + 1; i < start + 41; i++) {
    const c = centre(i)
    if (c && prev) d.push(Math.hypot(c.x + c.w / 2 - prev.x - prev.w / 2, c.y + c.h / 2 - prev.y - prev.h / 2))
    prev = c
  }
  d.sort((a, b) => a - b)
  return { sel, spacing: d.length ? d[Math.floor(d.length / 2)] : 0 }
}

/** Wrapper around the vendored PDBe RNA Viewer. One instance at a time; React state is authoritative (subscribeEvents is off). */
export function PdbRnaViewer2D({ pdbId, entityId, authAsymId, selection, onResidueClick, resetKey = 0, className = '' }: Props) {
  const host = useRef<HTMLDivElement>(null)
  const inst = useRef<PdbRnaViewerPluginInstance | null>(null)
  const ready = useRef<PdbRnaViewerPluginInstance | null>(null)   // the instance whose render() has completed
  const boxRef = useRef<HTMLDivElement | null>(null)              // container of the current render (for re-fitting after a selection)
  const [status, setStatus] = useState<'loading' | 'ready' | 'error' | 'unavailable'>('loading')
  const nNt = useRef(0)                                           // nucleotides drawn in the current topology
  const fullBox = useRef<Box | null>(null)                        // the fitted (whole topology) viewBox
  const anim = useRef(0)
  const [focused, setFocused] = useState(false)
  const clickRef = useRef(onResidueClick)
  clickRef.current = onResidueClick

  useEffect(() => {
    const el = host.current
    if (!el) return
    if (!window.PdbRnaViewerPlugin) { setStatus('unavailable'); return }
    let alive = true
    setStatus('loading')
    el.innerHTML = ''
    // each render gets its own child container: a render that completes after this effect was cleaned up only touches its detached box
    const box = document.createElement('div'); box.className = 'absolute inset-0'; el.appendChild(box); boxRef.current = box
    const v = new window.PdbRnaViewerPlugin()
    inst.current = v
    ready.current = null
    // (fetches below are cheap and idempotent; StrictMode double-mount just repeats them)
    const handler = (e: Event) => {
      const d = (e as RnaViewerEvent).eventData
      const id = d?.label_seq_ids?.[0] ?? d?.label_seq_id
      if (id !== undefined) clickRef.current?.(+id)
    }
    el.addEventListener('PDB.RNA.viewer.click', handler)
    // Fetch the two PDBe static inputs ourselves: a missing topology or FR3D file means "no 2D available" (the plugin throws on missing FR3D data).
    const base = 'https://www.ebi.ac.uk/pdbe/static/entry'
    Promise.all([fetch(`${base}/${pdbId}_${entityId}_${authAsymId.toUpperCase()}.json`), fetch(`${base}/${pdbId}_${authAsymId.toUpperCase()}_basepair.json`)])
      .then(async ([a, b]) => {
        if (!a.ok || !b.ok) { if (alive) setStatus('unavailable'); return }
        const [apiData, FR3DData] = await Promise.all([a.json(), b.json()])
        await v.render(box, { pdbId, entityId: String(entityId), chainId: authAsymId, subscribeEvents: false, apiData, FR3DData })
        if (!alive) return
        fitSvg(box)
        const svgEl = box.querySelector<SVGSVGElement>('svg.rnaTopoSvg')
        fullBox.current = svgEl ? parseVb(svgEl) : null
        nNt.current = new Set([...box.querySelectorAll('.rnaviewEle')].flatMap((e) => [...e.classList].filter((c) => c.startsWith(`rnaview_${pdbId}_`)))).size
        setFocused(false)
        if (box.querySelector('svg')) { ready.current = v; setStatus('ready') } else setStatus('unavailable')
      })
      .catch((err: unknown) => { console.error('[PdbRnaViewer2D]', err); if (alive) setStatus('error') })
    return () => {
      alive = false
      el.removeEventListener('PDB.RNA.viewer.click', handler)
      if (ready.current === v) { try { v.clearSelection(undefined, true) } catch { /* ignore */ } }
      inst.current = null
      ready.current = null
      box.remove()
    }
  }, [pdbId, entityId, authAsymId])

  const key = [...selection].sort((a, b) => a - b).join(',')
  useEffect(() => {
    const v = inst.current
    if (!v || status !== 'ready' || ready.current !== v) return   // only drive an instance whose render() completed
    try {
      const U = (window as unknown as { UiActionsService?: { selectNucleotide: (...a: unknown[]) => void; unSelectNucleotide?: (...a: unknown[]) => void } }).UiActionsService
      const o = (v as unknown as { options?: { pdbId: string; entityId: string | number } }).options
      // clear the previous colours WITHOUT the viewer's zoom reset: its clearSelection() is a 'click' deselect, which also resets zoom and pan
      if (U?.unSelectNucleotide && o) U.unSelectNucleotide(o.pdbId, o.entityId, undefined, 'select', false, undefined, undefined, true)
      else v.clearSelection(undefined, true)           // doNotReport: no echo event
      const ids = key ? key.split(',').map(Number) : []
      // colour all selected nucleotides at once WITHOUT the viewer's zoom-to-nucleotide (event type other than 'click'), so the current
      // 2D zoom and pan are preserved (whole loops, discontinuous segments included); falls back to the range API if unavailable
      if (ids.length && U?.selectNucleotide && o) U.selectNucleotide(o.pdbId, o.entityId, ids, 'select', false, undefined, undefined, undefined, true, true)
      else { for (const id of ids) v.selectResidueRange(id, id, undefined, true); if (ids.length && boxRef.current) fitSvg(boxRef.current) }
      if (typeof window !== 'undefined') window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), rna2dSelected: ids }
      if (ids.length) focusOn(ids)
    } catch (err) { console.error('[PdbRnaViewer2D] selection', err) }
  }, [key, status])

  /** Long chains: bring the selection (with its helix and surroundings) into focus with a short eased pan/zoom. Nothing moves when the selection
   *  is already inside a zoomed-in view (e.g. a click in the 2D itself); the whole fitted view always zooms in. Wheel/drag keep working afterwards. */
  function focusOn(ids: number[]) {
    const svg = boxRef.current?.querySelector<SVGSVGElement>('svg.rnaTopoSvg'), full = fullBox.current
    if (!svg || !full || nNt.current <= FOCUS_MIN_NT) return
    const { sel, spacing } = nucleotideGeometry(svg, pdbId, ids)
    if (!sel || !spacing) return
    const from = visibleRegion(svg); if (!from) return
    const zoomedIn = from.w < full.w * 0.98 || from.h < full.h * 0.98
    if (zoomedIn && contains(from, sel)) return
    const to = focusBox(sel, spacing * FOCUS_MIN_SPAN, full)
    // fold the user's wheel/drag transform into the viewBox (no visual jump), reset d3-zoom to identity, then animate the viewBox
    for (const g of svg.querySelectorAll<SVGGElement>('g[transform]')) g.removeAttribute('transform')
    const z = (svg as unknown as { __zoom?: { constructor: new (k: number, x: number, y: number) => unknown } }).__zoom
    if (z) (svg as unknown as { __zoom: unknown }).__zoom = new z.constructor(1, 0, 0)
    setVb(svg, from)
    cancelAnimationFrame(anim.current)
    const t0 = performance.now()
    const step = (t: number) => {
      const u = Math.min(1, (t - t0) / FOCUS_MS), e = u < 0.5 ? 2 * u * u : 1 - (-2 * u + 2) ** 2 / 2
      setVb(svg, { x: from.x + (to.x - from.x) * e, y: from.y + (to.y - from.y) * e, w: from.w + (to.w - from.w) * e, h: from.h + (to.h - from.h) * e })
      if (u < 1) anim.current = requestAnimationFrame(step)
    }
    anim.current = requestAnimationFrame(step)
    setFocused(true)
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), rna2dFocus: { ids, target: to, full } }
  }
  useEffect(() => () => cancelAnimationFrame(anim.current), [])

  // explicit deselection: default view (whole topology fitted, zoom/pan reset); the selection effect above has already cleared the colours
  const lastReset = useRef(resetKey)
  useEffect(() => {
    if (resetKey === lastReset.current) return
    lastReset.current = resetKey
    const v = inst.current
    if (!v || status !== 'ready' || ready.current !== v || !boxRef.current) return
    try { v.clearSelection(undefined, true) } catch { /* ignore */ }   // 'click' deselect: the viewer's own zoom reset (to identity)
    cancelAnimationFrame(anim.current); setFocused(false)
    fitSvg(boxRef.current)                                             // whole topology fitted, d3 zoom state = identity
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), rna2dReset: Number(window.__rnaExplorer?.rna2dReset ?? 0) + 1 }
  }, [resetKey, status])

  return (
    <div className={`relative rounded border border-stone-200 bg-white ${className}`} data-testid="rna2d-viewer" data-pdb={pdbId} data-status={status} data-focused={String(focused)}>
      <div ref={host} className="absolute inset-0" />
      {status === 'loading' && <div className="pointer-events-none absolute left-2 top-2 rounded bg-white/80 px-2 py-0.5 text-xs text-stone-600">Loading 2D for {pdbId.toUpperCase()} {authAsymId}…</div>}
      {status === 'unavailable' && <div className="absolute inset-0 flex items-center justify-center p-4 text-center text-sm text-stone-500">No PDBe 2D topology available for {pdbId.toUpperCase()} chain {authAsymId}.</div>}
      {status === 'error' && <div className="absolute left-2 top-2 rounded bg-red-50 px-2 py-0.5 text-xs text-red-700">2D viewer failed to load.</div>}
    </div>
  )
}
