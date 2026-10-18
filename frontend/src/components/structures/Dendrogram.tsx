import type { TreeNode } from '../../lib/structureGroups'
import { fmt, pathBetween } from '../../lib/structureGroups'
import './structures.css'

/*
 * Instance-level dendrogram in the heatmap's left gutter, one leaf per matrix row. Ported from
 * antibody_antigen_interface/app-react/src/groups/Dendrogram.jsx: root on the left, linear merge heights, a short stub past 0, neutral
 * grey branches with the selected group's branch in ink, group colour only in the clickable strip. Added for RNA: the duplicate-floor line
 * (differences below it are not interpretable) and the selected instance's leaf mark.
 */
const STUB = 8, BAR_W = 16, BAR_GAP = 5, PAD_L = 6
const INK = '#15191f', QUIET = '#9aa2ad', ABOVE = '#c9ced6'

export interface DendroGroup { id: number; label: string; color: string; node: TreeNode; from: number; to: number; size: number }

export function Dendrogram({ root, groups, selectedId, onSelectGroup, cell, top, width, maxHeight, cutHeight, floor, quantity, markLeaf, compareLeaves = [] }: {
  root: TreeNode; groups: DendroGroup[]; selectedId: number | null; onSelectGroup: (id: number) => void; cell: number; top: number; width: number
  maxHeight: number; cutHeight: number | null; floor: number | null; quantity: string; markLeaf?: number | null
  /** representative comparison: leaves (matrix indices) in their group colours; the route is drawn when both are in this tree */
  compareLeaves?: { leaf: number; color: string; label: string }[]
}) {
  const n = root.leaves.length
  const height = top + n * cell
  const xLeaf = width - BAR_W - BAR_GAP, xZero = xLeaf - STUB
  const x = (h: number) => (maxHeight > 0 ? xZero - (Math.min(h, maxHeight) / maxHeight) * (xZero - PAD_L) : xZero)
  const pos = new Map(root.leaves.map((leaf, p) => [leaf, p]))
  const yLeaf = (leaf: number) => top + (pos.get(leaf)! + 0.5) * cell
  const nodeGroup = new Map<number, number>()
  const mark = (node: TreeNode, g: number) => { nodeGroup.set(node.id, g); if (node.left && node.right) { mark(node.left, g); mark(node.right, g) } }
  for (const g of groups) mark(g.node, g.id)
  const lines: React.ReactNode[] = [], routeLines: React.ReactNode[] = []
  const marked = compareLeaves.filter((m) => pos.has(m.leaf))
  const route = marked.length === 2 ? pathBetween(root, marked[0].leaf, marked[1].leaf) : null
  const stroke = (g: number | undefined) => (g == null ? ABOVE : g === selectedId ? INK : QUIET)
  const w = (g: number | undefined) => (cell < 4 ? 0.6 : g != null && g === selectedId ? 1.9 : 1.2)
  const yOf = (node: TreeNode): number => {
    if (!node.left || !node.right) return yLeaf(node.leaves[0])
    const yl = yOf(node.left), yr = yOf(node.right), xn = x(node.height), g = nodeGroup.get(node.id)
    if (route?.has(node.id)) for (const [ch, yc] of [[node.left, yl], [node.right, yr]] as const) {
      if (!route.has(ch.id)) continue
      routeLines.push(<line key={`rh${ch.id}`} x1={xn} x2={ch.left ? x(ch.height) : xLeaf} y1={yc} y2={yc} />)
      routeLines.push(<line key={`rv${ch.id}`} x1={xn} x2={xn} y1={yc} y2={(yl + yr) / 2} />)
    }
    lines.push(<line key={`v${node.id}`} x1={xn} x2={xn} y1={yl} y2={yr} stroke={stroke(g)} strokeWidth={w(g)} />)
    for (const [child, yc] of [[node.left, yl], [node.right, yr]] as const) {
      const gc = nodeGroup.get(child.id)
      lines.push(<line key={`h${child.id}`} x1={xn} x2={child.left ? x(child.height) : xLeaf} y1={yc} y2={yc} stroke={stroke(gc)} strokeWidth={w(gc)} />)
    }
    return (yl + yr) / 2
  }
  yOf(root)
  const xCut = cutHeight != null && cutHeight <= maxHeight ? x(cutHeight) : null
  const xFloor = floor != null && floor <= maxHeight ? x(floor) : null
  return (
    <svg className="sg-dendro" width={width} height={height} role="img" aria-label={`Dendrogram of ${n} RNA instances`} data-testid="instance-dendrogram">
      {top > 0 && (
        <g className="sg-dendro-axis">
          <line x1={x(maxHeight)} x2={xZero} y1={top - 8} y2={top - 8} />
          <line x1={x(maxHeight)} x2={x(maxHeight)} y1={top - 11} y2={top - 5} />
          <line x1={xZero} x2={xZero} y1={top - 11} y2={top - 5} />
          <text x={x(maxHeight)} y={top - 15} textAnchor="start">{fmt(maxHeight)}</text>
          <text x={(x(maxHeight) + xZero) / 2} y={top - 15} textAnchor="middle">{quantity}</text>
          <text x={xZero} y={top - 15} textAnchor="end">0</text>
        </g>
      )}
      <g className={route ? 'sg-dendro-dim' : undefined}>{lines}</g>
      {route && <g className="sg-dendro-route" data-testid="dendro-route" data-nodes={[...route].sort((a, b) => a - b).join(',')}>{routeLines}</g>}
      {marked.map((m) => <circle key={`cmp${m.leaf}`} className="sg-dendro-mark" cx={xLeaf - 8} cy={yLeaf(m.leaf)} r={4.5} fill={m.color} data-testid="dendro-compare-mark"><title>{`${m.label}, representative (compared)`}</title></circle>)}
      {markLeaf != null && pos.has(markLeaf) && <circle className="sg-dendro-mark" cx={xLeaf} cy={yLeaf(markLeaf)} r={4} fill="#1d4ed8" data-testid="dendro-instance-mark"><title>Selected instance</title></circle>}
      {xFloor != null && <line className="sg-dendro-floor" x1={xFloor} x2={xFloor} y1={Math.max(0, top - 2)} y2={height} data-testid="dendro-floor"><title>{`Duplicate floor ${fmt(floor!)}: differences below it are within the spread of duplicate structures`}</title></line>}
      {xCut != null && <line className="sg-dendro-cut" x1={xCut} x2={xCut} y1={Math.max(0, top - 2)} y2={height}><title>{`Cut at ${fmt(cutHeight!)}`}</title></line>}
      {groups.map((g) => {
        const span = (g.to - g.from + 1) * cell, inset = span >= 10 ? 2 : 0, y0 = top + g.from * cell + inset, h = span - 2 * inset, on = g.id === selectedId, cy = y0 + h / 2
        return (
          <g key={g.id} className={'sg-dendro-bar' + (on ? ' on' : '')} role="button" tabIndex={0} aria-pressed={on} aria-label={`Select ${g.label}, ${g.size} RNA instances`}
            onClick={() => onSelectGroup(g.id)} onKeyDown={(e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); onSelectGroup(g.id) } }} data-testid={`dendro-bar-${g.id}`}>
            <title>{`${g.label}: ${g.size} RNA instances. Click to select.`}</title>
            <rect x={xLeaf + BAR_GAP - 1} y={y0} width={BAR_W - 2} height={Math.max(2, h)} rx={3} fill={g.color} stroke={on ? INK : g.color} strokeWidth={on ? 2 : 1} />
            {h >= 46 && <text x={xLeaf + BAR_GAP + BAR_W / 2 - 1} y={cy} textAnchor="middle" dominantBaseline="central" transform={`rotate(-90 ${xLeaf + BAR_GAP + BAR_W / 2 - 1} ${cy})`}>{g.label.replace('Structural group', 'Group')}</text>}
          </g>
        )
      })}
    </svg>
  )
}
