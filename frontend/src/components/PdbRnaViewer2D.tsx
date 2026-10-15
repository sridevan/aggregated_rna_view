import { useEffect, useRef, useState } from 'react'

interface Props {
  pdbId: string
  entityId: number
  /** author chain id (the PDBe RNA viewer and its topology data are keyed by author chain) */
  authAsymId: string
  /** label_seq_ids to show as selected (the viewer's API is label_seq_id based) */
  selection: number[]
  onResidueClick?: (labelSeqId: number) => void
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
}

/** Wrapper around the vendored PDBe RNA Viewer. One instance at a time; React state is authoritative (subscribeEvents is off). */
export function PdbRnaViewer2D({ pdbId, entityId, authAsymId, selection, onResidueClick, className = '' }: Props) {
  const host = useRef<HTMLDivElement>(null)
  const inst = useRef<PdbRnaViewerPluginInstance | null>(null)
  const ready = useRef<PdbRnaViewerPluginInstance | null>(null)   // the instance whose render() has completed
  const [status, setStatus] = useState<'loading' | 'ready' | 'error' | 'unavailable'>('loading')
  const clickRef = useRef(onResidueClick)
  clickRef.current = onResidueClick

  useEffect(() => {
    const el = host.current
    if (!el) return
    if (!window.PdbRnaViewerPlugin) { setStatus('unavailable'); return }
    let alive = true
    setStatus('loading')
    el.innerHTML = ''
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
        await v.render(el, { pdbId, entityId: String(entityId), chainId: authAsymId, subscribeEvents: false, apiData, FR3DData })
        if (!alive) return
        fitSvg(el)
        if (el.querySelector('svg')) { ready.current = v; setStatus('ready') } else setStatus('unavailable')
      })
      .catch((err: unknown) => { console.error('[PdbRnaViewer2D]', err); if (alive) setStatus('error') })
    return () => {
      alive = false
      el.removeEventListener('PDB.RNA.viewer.click', handler)
      if (ready.current === v) { try { v.clearSelection(undefined, true) } catch { /* ignore */ } }
      inst.current = null
      ready.current = null
      el.innerHTML = ''
    }
  }, [pdbId, entityId, authAsymId])

  const key = [...selection].sort((a, b) => a - b).join(',')
  useEffect(() => {
    const v = inst.current
    if (!v || status !== 'ready' || ready.current !== v) return   // only drive an instance whose render() completed
    try {
      v.clearSelection(undefined, true)           // doNotReport: no echo event
      for (const id of key ? key.split(',').map(Number) : []) v.selectResidueRange(id, id, undefined, true)   // range form accumulates (multipleSelection)
    } catch (err) { console.error('[PdbRnaViewer2D] selection', err) }
  }, [key, status])

  return (
    <div className={`relative rounded border border-stone-200 bg-white ${className}`} data-testid="rna2d-viewer" data-pdb={pdbId} data-status={status}>
      <div ref={host} className="absolute inset-0" />
      {status === 'loading' && <div className="pointer-events-none absolute left-2 top-2 rounded bg-white/80 px-2 py-0.5 text-xs text-stone-600">Loading 2D for {pdbId.toUpperCase()} {authAsymId}…</div>}
      {status === 'unavailable' && <div className="absolute inset-0 flex items-center justify-center p-4 text-center text-sm text-stone-500">No PDBe 2D topology available for {pdbId.toUpperCase()} chain {authAsymId}.</div>}
      {status === 'error' && <div className="absolute left-2 top-2 rounded bg-red-50 px-2 py-0.5 text-xs text-red-700">2D viewer failed to load.</div>}
    </div>
  )
}
