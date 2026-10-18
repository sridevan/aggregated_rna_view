import { useEffect, useMemo, useRef, useState } from 'react'
import type { InstanceRow } from '../../types/backend'
import { DIAGONAL, fmt, instanceLabel, viridis, type DissimilarityMatrix } from '../../lib/structureGroups'
import { displayOrganism } from '../../lib/format'
import { HoverCard } from '../InfoTip'

/*
 * Structural Dissimilarity mode: the full geometric discrepancy matrix in the stored TPPL seriation order (rows and columns share it).
 * Same viridis scale (0 to the family maximum) and grey diagonal as the clustering heatmap. Interactions: hover any cell for the pair and
 * its value (a diagonal cell: that RNA instance's details); click a diagonal cell to show that instance alone in 3D; click a cell BELOW the diagonal to superpose the pair in the 3D view
 * (upper-triangle cells are hover-only). Hover and clicks are handled once on the SVG (no per-cell handlers), and the cells are memoised,
 * so hovering does not rebuild the matrix.
 */
export const LABEL_LIMIT = 50          // Y-axis instance labels only up to this many instances
const LABEL_W = 64

const cellAt = (svg: SVGSVGElement, e: React.MouseEvent, labelW: number, cell: number, n: number) => {
  const b = svg.getBoundingClientRect()
  const c = Math.floor((e.clientX - b.left - labelW) / cell), r = Math.floor((e.clientY - b.top) / cell)
  return r >= 0 && r < n && c >= 0 && c < n ? { r, c, anchor: new DOMRect(b.left + labelW + c * cell, b.top + r * cell, cell, cell) } : null
}

/** PDB ID and chain, title, method, resolution and organism of one RNA instance (absent fields omitted); used on hover and in the pinned card */
function InstanceDetails({ seqName, row, prefix }: { seqName: string; row: InstanceRow | undefined; prefix: string }) {
  const org = row ? displayOrganism(row) : null
  return <>
    <p className="font-semibold" data-testid={`${prefix}-id`}>{instanceLabel(seqName)}</p>
    {row?.title && <p data-testid={`${prefix}-title`}>{row.title}</p>}
    {row?.experimental_method && <p data-testid={`${prefix}-method`}><span className="text-stone-500">Method:</span> {row.experimental_method}</p>}
    {row?.resolution != null && <p data-testid={`${prefix}-resolution`}><span className="text-stone-500">Resolution:</span> {row.resolution} Å</p>}
    {org?.name && <p data-testid={`${prefix}-organism`}><span className="text-stone-500">Organism:</span> <span className={org.isSpecies ? 'italic' : ''}>{org.name}</span></p>}
  </>
}

export function DissimilarityHeatmap({ m, nShared, rowOf, shown, onSelectPair, onSelectOne }: {
  m: DissimilarityMatrix; nShared: number; rowOf: (seqName: string) => InstanceRow | undefined
  /** the structures in the 3D view: one (checked on the diagonal) or a pair (its cell outlined, both checked on the diagonal) */
  shown: string[]; onSelectPair: (a: string, b: string) => void
  /** a diagonal click: show that instance alone */
  onSelectOne: (seqName: string) => void
}) {
  const boxRef = useRef<HTMLDivElement>(null), svgRef = useRef<SVGSVGElement>(null)
  const [W, setW] = useState(560)
  useEffect(() => {
    const el = boxRef.current
    if (!el || typeof ResizeObserver === 'undefined') return
    const fit = () => { const w = el.getBoundingClientRect().width; if (w > 0) setW(Math.round(w)) }
    fit(); const ro = new ResizeObserver(fit); ro.observe(el); return () => ro.disconnect()
  }, [])
  const n = m.ids.length
  const labels = n <= LABEL_LIMIT
  const labelW = labels ? LABEL_W : 0
  const cell = Math.max(3, Math.min(18, Math.floor((W - labelW - 2) / n)))
  const [hover, setHover] = useState<{ r: number; c: number; anchor: DOMRect } | null>(null)
  const pr = shown.length === 2 ? shown.map((sn) => m.ids.indexOf(sn)) : null
  const checks = shown.map((sn) => m.ids.indexOf(sn)).filter((k) => k >= 0)
  const sel = pr && pr[0] >= 0 && pr[1] >= 0 && pr[0] !== pr[1] ? { r: Math.max(pr[0], pr[1]), c: Math.min(pr[0], pr[1]) } : null

  const cells = useMemo(() => m.values.map((row, r) => row.map((v, c) => (
    <rect key={`${r}-${c}`} x={labelW + c * cell} y={r * cell} width={cell} height={cell} fill={r === c ? DIAGONAL : v == null ? 'url(#dm-missing)' : viridis(m.max > 0 ? v / m.max : 0)}
      data-testid={`dm-cell-${r}-${c}`} data-value={r === c || v == null ? '' : v.toFixed(6)} className={r === c || r > c ? 'cursor-pointer' : ''} />
  ))), [m, cell, labelW])

  const onMove = (e: React.MouseEvent) => { const h = svgRef.current && cellAt(svgRef.current, e, labelW, cell, n); setHover((x) => (h ? (x && x.r === h.r && x.c === h.c ? x : h) : null)) }
  const onClick = (e: React.MouseEvent) => {
    const h = svgRef.current && cellAt(svgRef.current, e, labelW, cell, n)
    if (!h) return
    if (h.r === h.c) onSelectOne(m.ids[h.r])                          // the diagonal: that instance alone
    else if (h.r > h.c) onSelectPair(m.ids[h.r], m.ids[h.c])           // below the diagonal: superpose the pair
  }
  const hv = hover ? m.values[hover.r][hover.c] : null
  return (
    <div ref={boxRef} className="w-full" data-testid="dissimilarity-heatmap" data-n={n} data-labels={String(labels)}>
      <div className="overflow-x-auto">
        <svg ref={svgRef} width={labelW + n * cell} height={n * cell} role="img" aria-label={`Geometric discrepancy matrix of ${n} RNA instances in seriation order`} onMouseMove={onMove} onMouseLeave={() => setHover(null)} onClick={onClick} data-testid="dm-svg" data-cell={cell}>
          <defs>
            <pattern id="dm-missing" width="4" height="4" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width="4" height="4" fill="#f5f5f4" /><line x1="0" y1="0" x2="0" y2="4" stroke="#a8a29e" strokeWidth="1.5" /></pattern>
          </defs>
          {labels && m.ids.map((sn, r) => {
            const on = checks.includes(r)
            return <text key={sn} x={labelW - 4} y={(r + 0.5) * cell} textAnchor="end" dominantBaseline="central" className={`text-[10px] ${on ? 'fill-blue-700 font-bold' : 'fill-stone-700'}`} data-testid={`dm-row-${sn}`}>{instanceLabel(sn)}</text>
          })}
          {cells}
          {sel && (
            <g pointerEvents="none" data-testid="dm-selection" data-row={sel.r} data-col={sel.c}>
              <rect x={labelW} y={sel.r * cell} width={n * cell} height={cell} fill="none" stroke="#1d4ed8" strokeOpacity={0.6} strokeWidth={1} />
              <rect x={labelW + sel.c * cell} y={0} width={cell} height={n * cell} fill="none" stroke="#1d4ed8" strokeOpacity={0.6} strokeWidth={1} />
              <rect x={labelW + sel.r * cell} y={sel.c * cell} width={cell} height={cell} fill="none" stroke="#1d4ed8" strokeWidth={1} strokeDasharray="2 1" />
              <rect x={labelW + sel.c * cell - 1} y={sel.r * cell - 1} width={cell + 2} height={cell + 2} fill="none" stroke="#111827" strokeWidth={2} data-testid="dm-selected-cell" />
            </g>
          )}
          {/* the displayed structures, checked on the diagonal (the representative alone by default, or both instances of a pair) */}
          <g pointerEvents="none">
            {checks.map((k) => (
              <g key={k} data-testid={`dm-diag-check-${k}`}>
                <rect x={labelW + k * cell} y={k * cell} width={cell} height={cell} fill="#1d4ed8" />
                {cell >= 7 && <path d={`M${labelW + k * cell + cell * 0.22},${k * cell + cell * 0.52} L${labelW + k * cell + cell * 0.42},${k * cell + cell * 0.72} L${labelW + k * cell + cell * 0.8},${k * cell + cell * 0.3}`} fill="none" stroke="#fff" strokeWidth={Math.max(1.2, cell / 7)} strokeLinecap="round" strokeLinejoin="round" />}
              </g>
            ))}
          </g>
        </svg>
      </div>
      <div className="mt-2 flex flex-wrap items-center gap-x-2 gap-y-1 text-[11px] text-stone-600" data-testid="dissimilarity-legend">
        <span>0</span>
        <span className="inline-block h-2.5 w-40 rounded-[2px]" style={{ background: `linear-gradient(to right, ${[0, 0.25, 0.5, 0.75, 1].map((t) => viridis(t)).join(',')})` }} />
        <span>{fmt(m.max)}</span>
        <span className="ml-1">More similar ← <b>geometric discrepancy</b> → More different</span>
        <span className="ml-1 inline-flex items-center gap-1"><span className="inline-block h-2.5 w-2.5 rounded-[2px]" style={{ background: DIAGONAL }} />instance (hover for details, click to show)</span>
        {m.nMissing > 0 && <span className="inline-flex items-center gap-1" data-testid="dm-legend-missing"><svg width="10" height="10"><rect width="10" height="10" fill="url(#dm-missing)" /></svg>not measured ({m.nMissing} pairs)</span>}
      </div>
      <HoverCard anchor={hover ? hover.anchor : null} testId="dm-tooltip">
        {hover && (hover.r === hover.c
          ? <div className="space-y-0.5"><InstanceDetails seqName={m.ids[hover.r]} row={rowOf(m.ids[hover.r])} prefix="dm-tooltip" /></div>
          : <>
              <p className="font-semibold"><span data-testid="dm-tooltip-a">{instanceLabel(m.ids[hover.r])}</span> × <span data-testid="dm-tooltip-b">{instanceLabel(m.ids[hover.c])}</span></p>
              <p data-testid="dm-tooltip-value">{hv == null ? 'Geometric discrepancy: not measured' : `Geometric discrepancy ${hv.toFixed(3)}`}</p>
              <p className="text-stone-500">{nShared} shared Rfam positions compared</p>
              {hover.r > hover.c && <p className="text-stone-500">Click to superpose this pair in 3D</p>}
            </>)}
      </HoverCard>
    </div>
  )
}
