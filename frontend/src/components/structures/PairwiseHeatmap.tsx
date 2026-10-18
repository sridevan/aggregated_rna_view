import { useEffect, useRef, useState } from 'react'
import type { Clustering, TreeNode } from '../../lib/structureGroups'
import { DIAGONAL, fmt, instanceLabel, viridis } from '../../lib/structureGroups'
import { Dendrogram, type DendroGroup } from './Dendrogram'

/*
 * Instance-level pairwise view (<= detail_max instances): dendrogram in the left gutter, row labels, the matrix in dendrogram order, column
 * labels on top, a viridis colour bar. Adapted from antibody_antigen_interface/app-react/src/states/DissimilarityHeatmap.jsx (viridis, grey
 * diagonal, labelled bar) without its drag-to-drill and measure switch. The colour scale is fixed per metric from 0 to the family-wide
 * maximum dissimilarity, so a group's drill-down is coloured on the same scale as any other view of that metric.
 */
const LABEL_W = 64, TOP = 58, GUTTER = 150

export function PairwiseHeatmap({ c, root, groups, selectedId, onSelectGroup, selectedInstance, onSelectInstance, nPositions, compareLeaves = [] }: {
  c: Clustering; root: TreeNode; groups: DendroGroup[]; selectedId: number | null; onSelectGroup: (id: number) => void
  selectedInstance: string | null; onSelectInstance: (seqName: string) => void; nPositions: number; compareLeaves?: { leaf: number; color: string; label: string }[]
}) {
  const boxRef = useRef<HTMLDivElement>(null)
  const [W, setW] = useState(560)
  useEffect(() => {
    const el = boxRef.current
    if (!el || typeof ResizeObserver === 'undefined') return
    const fit = () => { const w = el.getBoundingClientRect().width; if (w > 0) setW(Math.round(w)) }
    fit(); const ro = new ResizeObserver(fit); ro.observe(el); return () => ro.disconnect()
  }, [])
  const order = root.leaves
  const n = order.length
  const cell = Math.max(8, Math.min(40, Math.floor((W - GUTTER - LABEL_W - 4) / n)))
  const maxD = c.metric.distance_range.max
  const sim = c.metric.kind === 'similarity'
  const q = sim ? '1 − TM-score' : 'geometric discrepancy'
  const selIdx = selectedInstance ? c.ids.indexOf(selectedInstance) : -1
  const sizeW = GUTTER + LABEL_W + n * cell
  const heightPx = TOP + n * cell
  const describe = (i: number, j: number) => {
    const a = instanceLabel(c.ids[i]), b = instanceLabel(c.ids[j])
    if (i === j) return `${a}: itself`
    return sim ? `${a} × ${b}: TM-score ${c.value[i][j].toFixed(3)} (1 − TM ${c.dist[i][j].toFixed(3)}; ${nPositions} positions)` : `${a} × ${b}: geometric discrepancy ${c.dist[i][j].toFixed(3)} (${nPositions} positions)`
  }
  return (
    <div ref={boxRef} className="w-full" data-testid="pairwise-heatmap">
      <div className="overflow-x-auto">
        <div className="relative" style={{ width: sizeW, height: heightPx }}>
          <div className="absolute left-0 top-0"><Dendrogram root={root} groups={groups} selectedId={selectedId} onSelectGroup={onSelectGroup} cell={cell} top={TOP} width={GUTTER}
            maxHeight={root.height || c.maxHeight} cutHeight={groups.length > 1 ? c.cutHeight : null} floor={c.floor} quantity={sim ? '1 − TM' : 'discrepancy'} compareLeaves={compareLeaves} markLeaf={selIdx >= 0 ? selIdx : null} /></div>
          <svg className="absolute top-0" style={{ left: GUTTER }} width={LABEL_W + n * cell} height={heightPx} role="img" aria-label={`Pairwise ${q} matrix of ${n} RNA instances`}>
            {order.map((i, r) => {
              const isRep = groups.some((g) => c.ids[i] === c.groups.find((x) => x.group_id === g.id)?.representative)
              return (
                <g key={`r${i}`}>
                  <text x={LABEL_W - 4} y={TOP + (r + 0.5) * cell} textAnchor="end" dominantBaseline="central" className={`cursor-pointer text-[10px] ${i === selIdx ? 'fill-blue-700 font-bold' : 'fill-stone-700'}`}
                    onClick={() => onSelectInstance(c.ids[i])} data-testid={`hm-row-${c.ids[i]}`}>{instanceLabel(c.ids[i])}{isRep ? ' ◆' : ''}<title>{`${c.ids[i]}${isRep ? ' (group representative, medoid)' : ''}. Click to select this instance.`}</title></text>
                  <text transform={`translate(${LABEL_W + (r + 0.5) * cell},${TOP - 4}) rotate(-60)`} className={`text-[9px] ${i === selIdx ? 'fill-blue-700 font-bold' : 'fill-stone-600'}`}>{instanceLabel(c.ids[i])}</text>
                  {order.map((j, col) => (
                    <rect key={j} x={LABEL_W + col * cell} y={TOP + r * cell} width={cell} height={cell} fill={i === j ? DIAGONAL : viridis(maxD > 0 ? c.dist[i][j] / maxD : 0)} stroke="#fff" strokeWidth={cell >= 14 ? 0.5 : 0}
                      data-testid={`hm-cell-${r}-${col}`} data-value={i === j ? '' : (sim ? c.value[i][j] : c.dist[i][j]).toFixed(4)}><title>{describe(i, j)}</title></rect>
                  ))}
                </g>
              )
            })}
            {selIdx >= 0 && order.includes(selIdx) && (() => { const r = order.indexOf(selIdx); return <><rect x={LABEL_W} y={TOP + r * cell} width={n * cell} height={cell} fill="none" stroke="#1d4ed8" strokeWidth={1.5} /><rect x={LABEL_W + r * cell} y={TOP} width={cell} height={n * cell} fill="none" stroke="#1d4ed8" strokeWidth={1.5} /></> })()}
          </svg>
        </div>
      </div>
      <div className="mt-2 flex items-center gap-2 text-[11px] text-stone-600" data-testid="heatmap-legend">
        <span>0</span>
        <span className="inline-block h-2.5 w-40 rounded-[2px]" style={{ background: `linear-gradient(to right, ${[0, 0.25, 0.5, 0.75, 1].map((t) => viridis(t)).join(',')})` }} />
        <span>{fmt(maxD)}</span>
        <span className="ml-1">More similar ← <b>{q}</b> → More different</span>
        <span className="ml-1 inline-block h-2.5 w-2.5 rounded-[2px]" style={{ background: DIAGONAL }} /> <span>self</span>
      </div>
    </div>
  )
}
