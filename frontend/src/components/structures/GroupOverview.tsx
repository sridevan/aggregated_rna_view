import { useEffect, useRef, useState } from 'react'
import type { Clustering, TreeNode } from '../../lib/structureGroups'
import { fmt, instanceLabel, pathBetween } from '../../lib/structureGroups'
import './structures.css'

/*
 * Collapsed structural-group dendrogram for a clustering too large for the instance-level view (> detail_max instances). Ported from
 * antibody_antigen_interface/app-react/src/groups/GroupOverview.jsx: the tree drawn down to the groups, one wedge per group, clamped row
 * heights with a proportional bar, the dashed cut line, every row selects its group. Added for RNA: the duplicate-floor line, the
 * representative id on each row, and a mark on the row holding the selected RNA instance.
 */
const MIN_ROW = 34, MAX_ROW = 132, GAP = 6, MIN_ROWS_TOTAL = 320, TOP = 40, BOTTOM = 24, STUB = 10, STRIP_W = 14
const INK = '#15191f', ABOVE = '#9aa2ad'

export function GroupOverview({ c, selectedId, onSelectGroup, quantity, selectedInstance, compare = null }: {
  c: Clustering; selectedId: number | null; onSelectGroup: (id: number) => void; quantity: string; selectedInstance: string | null
  /** representative comparison: the two group ids; their route through the stored tree is drawn in ink (as in the reference demo) */
  compare?: [number, number] | null
}) {
  const boxRef = useRef<HTMLDivElement>(null)
  const [W, setW] = useState(540)
  useEffect(() => {
    const el = boxRef.current
    if (!el || typeof ResizeObserver === 'undefined') return
    const fit = () => { const w = el.getBoundingClientRect().width; if (w > 0) setW(Math.round(w)) }
    fit(); const ro = new ResizeObserver(fit); ro.observe(el); return () => ro.disconnect()
  }, [])
  const { root, groups, cutHeight, maxHeight, floor } = c
  const rows = [...groups].sort((a, b) => a.from - b.from)
  const largest = Math.max(...groups.map((g) => g.n_instances)), total = root.leaves.length
  const natural = (n: number) => MIN_ROW + (MAX_ROW - MIN_ROW) * (largest > 1 ? (n - 1) / (largest - 1) : 0)
  const sum = rows.reduce((t, g) => t + natural(g.n_instances), 0) + GAP * (rows.length - 1)
  const stretch = Math.max(1, MIN_ROWS_TOTAL / sum)
  let y = TOP
  const at = new Map<number, { y: number; h: number; mid: number }>()
  for (const g of rows) { const h = Math.round(natural(g.n_instances) * stretch); at.set(g.group_id, { y, h, mid: y + h / 2 }); y += h + GAP }
  const bottom = y - GAP, height = bottom + BOTTOM
  const xBase = Math.round(W * 0.38), xZero = xBase - STUB, xRoot = 10
  const x = (h: number) => (maxHeight > 0 ? xZero - (h / maxHeight) * (xZero - xRoot) : xZero)
  const xStrip = xBase + 6, xText = xStrip + STRIP_W + 10, barMax = Math.max(40, W - xText - 12)
  const topOf = new Map(groups.map((g) => [g.node, g]))
  const lines: React.ReactNode[] = [], wedges: React.ReactNode[] = [], routeLines: React.ReactNode[] = []
  const sides = (compare ?? []).map((id) => groups.find((g) => g.group_id === id)).filter((g): g is NonNullable<typeof g> => !!g)
  const route = sides.length === 2 ? pathBetween(root, sides[0].node, sides[1].node) : null
  const compared = new Set(sides.map((g) => g.group_id))
  const yOf = (node: TreeNode): number => {
    const g = topOf.get(node.id)
    if (g) {
      const r = at.get(g.group_id)!, on = g.group_id === selectedId
      if (route?.has(node.id)) routeLines.push(<line key={`rg${g.group_id}`} x1={node.left ? x(node.height) : xZero} x2={xBase} y1={r.mid} y2={r.mid} />)
      if (node.left) wedges.push(<polygon key={`w${g.group_id}`} className={'sg-ov-wedge' + (on ? ' on' : '')} points={`${x(node.height)},${r.mid} ${xBase},${r.y + 3} ${xBase},${r.y + r.h - 3}`} />)
      else lines.push(<line key={`s${g.group_id}`} x1={xZero} x2={xBase} y1={r.mid} y2={r.mid} stroke={on ? INK : ABOVE} strokeWidth={on ? 1.8 : 1.2} />)
      return r.mid
    }
    const yl = yOf(node.left!), yr = yOf(node.right!), xn = x(node.height)
    if (route?.has(node.id)) {
      for (const [ch, yc] of [[node.left!, yl], [node.right!, yr]] as const) {
        if (!route.has(ch.id)) continue
        routeLines.push(<line key={`rh${ch.id}`} x1={xn} x2={x(ch.height)} y1={yc} y2={yc} />)
        routeLines.push(<line key={`rv${ch.id}`} x1={xn} x2={xn} y1={yc} y2={(yl + yr) / 2} />)
      }
    }
    lines.push(<line key={`v${node.id}`} x1={xn} x2={xn} y1={yl} y2={yr} stroke={ABOVE} strokeWidth={1.2} />)
    for (const [ch, yc] of [[node.left!, yl], [node.right!, yr]] as const) {
      const cg = topOf.get(ch.id), on = cg && cg.group_id === selectedId
      lines.push(<line key={`h${ch.id}`} x1={xn} x2={x(ch.height)} y1={yc} y2={yc} stroke={on ? INK : ABOVE} strokeWidth={on ? 1.8 : 1.2} />)
    }
    return (yl + yr) / 2
  }
  yOf(root)
  const instGroup = selectedInstance ? c.groupOf[selectedInstance] : undefined
  return (
    <div className="sg-ov" ref={boxRef} data-testid="group-overview">
      <svg width={W} height={height} role="group" aria-label={`${groups.length} structural groups covering ${total} RNA instances`}>
        {rows.map((g) => { const r = at.get(g.group_id)!; return g.group_id === selectedId && <rect key={`bg${g.group_id}`} className="sg-ov-sel" style={{ ['--g' as string]: g.color }} x={xBase - 2} y={r.y - 2} width={W - xBase + 1} height={r.h + 4} rx={6} /> })}
        <g className="sg-dendro-axis">
          <line x1={x(maxHeight)} x2={xZero} y1={TOP - 14} y2={TOP - 14} />
          <line x1={x(maxHeight)} x2={x(maxHeight)} y1={TOP - 17} y2={TOP - 11} />
          <line x1={xZero} x2={xZero} y1={TOP - 17} y2={TOP - 11} />
          <text x={x(maxHeight)} y={TOP - 21} textAnchor="start">{fmt(maxHeight)}</text>
          <text x={(x(maxHeight) + xZero) / 2} y={TOP - 21} textAnchor="middle">{quantity} at merge</text>
          <text x={xZero} y={TOP - 21} textAnchor="end">0</text>
        </g>
        {sides.filter((g) => g.group_id !== selectedId).map((g) => { const r = at.get(g.group_id)!; return <rect key={`cmp${g.group_id}`} className="sg-ov-cmp" style={{ ['--g' as string]: g.color }} x={xBase - 2} y={r.y - 2} width={W - xBase + 1} height={r.h + 4} rx={6} data-testid={`overview-compare-${g.group_id}`} /> })}
        <g className={route ? 'sg-dendro-dim' : undefined}>{wedges}{lines}</g>
        {route && <g className="sg-dendro-route" data-testid="overview-route" data-nodes={[...route].sort((a, b) => a - b).join(',')}>{routeLines}</g>}
        {floor != null && floor <= maxHeight && <line className="sg-dendro-floor" x1={x(floor)} x2={x(floor)} y1={TOP - 6} y2={bottom + 4} data-testid="overview-floor"><title>{`Duplicate floor ${fmt(floor)}`}</title></line>}
        <line className="sg-dendro-cut" x1={x(cutHeight)} x2={x(cutHeight)} y1={TOP - 6} y2={bottom + 4} />
        <text className="sg-ov-cut" x={x(cutHeight)} y={bottom + 17} textAnchor="middle" data-testid="overview-cut">cut at {fmt(cutHeight)}</text>
        {rows.map((g) => {
          const r = at.get(g.group_id)!, on = g.group_id === selectedId, n = g.n_instances, tall = r.h >= 58, ty = tall ? r.mid - 9 : r.mid
          return (
            <g key={g.group_id} className={'sg-ov-row' + (on ? ' on' : '')} role="button" tabIndex={0} aria-pressed={on} aria-label={`Select ${g.label}, ${n} RNA instances`}
              onClick={() => onSelectGroup(g.group_id)} onKeyDown={(e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); onSelectGroup(g.group_id) } }} data-testid={`overview-row-${g.group_id}`}>
              <title>{`${g.label}: ${n} of ${total} RNA instances. Representative ${instanceLabel(g.representative)}. Click to select.`}</title>
              <rect x={xStrip} y={r.y + 2} width={STRIP_W} height={r.h - 4} rx={3} fill={g.color} />
              <text className="sg-ov-name" x={xText} y={ty} dominantBaseline="central">{g.label}<tspan className="sg-ov-n" dx={8}>{n} {n === 1 ? 'RNA instance' : 'RNA instances'}</tspan></text>
              {tall && <text className="sg-ov-rep" x={W - 12} y={ty} textAnchor="end" dominantBaseline="central">rep. {instanceLabel(g.representative)}</text>}
              <rect className="sg-ov-track" x={xText} y={tall ? r.mid + 8 : r.y + r.h - 7} width={barMax} height={4} rx={2} />
              <rect x={xText} y={tall ? r.mid + 8 : r.y + r.h - 7} width={Math.max(3, barMax * (n / total))} height={4} rx={2} fill={g.color} />
              {compared.has(g.group_id) && <circle className="sg-dendro-mark" cx={xBase - 10} cy={r.mid} r={5} fill={g.color} data-testid={`overview-compare-mark-${g.group_id}`}><title>{`${g.label}: representative ${instanceLabel(g.representative)} compared`}</title></circle>}
              {instGroup === g.group_id && <circle className="sg-ov-inst" cx={xBase} cy={r.mid} r={4.5} data-testid="overview-instance-mark"><title>{`Selected instance ${instanceLabel(selectedInstance!)} is in this group`}</title></circle>}
              <rect className="sg-ov-hit" x={0} y={r.y - GAP / 2} width={W} height={r.h + GAP} />
            </g>
          )
        })}
      </svg>
    </div>
  )
}
