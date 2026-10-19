import { memo, useMemo, useRef, useState, type ReactNode } from 'react'
import { arcPath, type CircleLayout } from '../../lib/interactions'
import { HoverCard } from '../InfoTip'

/*
 * R3DCID-inspired circular interaction diagram (SVG). Consensus positions sit clockwise on the circle (fixed layout from the Rfam span, so
 * every diagram of the family shares the same angles); interactions are circular arcs orthogonal to the circle, drawn back to front by
 * `order` (most important on top). Labels: the selected numbering's round numbers, structural-region band and names, and the selected or
 * looked-up positions. Arcs are memoised; hover and click are handled once on the SVG through invisible wide hit paths.
 */
export interface DiagramArc { id: string; i: number; j: number; color: string; width: number; dash?: string; opacity?: number; order: number
  /** several interactions between the same two nucleotides: drawn as one double-line arc */
  multi?: boolean
  /** category (or 'multiple') and comparison status, exposed as data attributes */
  category?: string; status?: string }
export interface DiagramLabel { short: string; text: string; fromReference: boolean }
const REGION_COLORS = ['#e7e5e4', '#d6d3d1']

export const CircularDiagram = memo(function CircularDiagram({ L, size, arcs, labelFor, regions, absent, selectedId, focus, onSelect, tooltip, testId, title, dimUnselected = false }: {
  L: CircleLayout; size: number; arcs: DiagramArc[]
  labelFor: (pos: number) => DiagramLabel
  regions: (string | null)[]
  /** positions without a residue (or coordinates) in the instance drawn by this diagram: hollow markers */
  absent?: Set<number>
  selectedId: string | null
  /** highlighted positions (selected interaction or lookup) */
  focus: number[]
  onSelect: (id: string) => void
  tooltip: (id: string) => ReactNode
  testId: string; title: string
  /** draw every other arc faded while one is selected */
  dimUnselected?: boolean
}) {
  const svgRef = useRef<SVGSVGElement>(null)
  const [hover, setHover] = useState<{ id: string; anchor: DOMRect } | null>(null)
  const positions = useMemo(() => Array.from({ length: L.last - L.first + 1 }, (_, k) => L.first + k), [L])

  // region band: runs of equal non-null region names
  const regionRuns = useMemo(() => {
    const runs: { name: string; from: number; to: number }[] = []
    for (const p of positions) {
      const r = regions[p - 1]; const last = runs[runs.length - 1]
      if (r && last && last.name === r && last.to === p - 1) last.to = p
      else if (r) runs.push({ name: r, from: p, to: p })
    }
    return runs
  }, [positions, regions])
  const ticks = useMemo(() => {
    const out: number[] = []
    for (const p of positions) {
      const l = labelFor(p)
      const round = l.fromReference ? /^-?\d+$/.test(l.short) && +l.short % 10 === 0 : p % 10 === 0
      if (round && out.every((q) => Math.abs(q - p) >= 4)) out.push(p)
    }
    return out
  }, [positions, labelFor])

  const sorted = useMemo(() => [...arcs].sort((a, b) => a.order - b.order), [arcs])
  const base = useMemo(() => (
    <g data-testid={`${testId}-arcs`}>
      {sorted.map((a) => <g key={a.id}><path d={arcPath(L, a.i, a.j)} fill="none" stroke={a.color} strokeWidth={a.multi ? a.width + 1.6 : a.width} strokeDasharray={a.dash} strokeOpacity={a.opacity ?? 0.9} strokeLinecap="round"
        data-testid={`${testId}-arc-${a.id}`} data-i={a.i} data-j={a.j} data-color={a.color} data-multi={a.multi ? 'true' : undefined} data-category={a.category} data-status={a.status} data-dash={a.dash ?? ''} />
        {a.multi && <path d={arcPath(L, a.i, a.j)} fill="none" stroke="white" strokeWidth={0.8} strokeOpacity={a.opacity ?? 0.9} />}</g>)}
    </g>), [sorted, L, testId])
  const hits = useMemo(() => (
    <g>{sorted.map((a) => <path key={a.id} d={arcPath(L, a.i, a.j)} fill="none" stroke="transparent" strokeWidth={8} data-hit={a.id} style={{ cursor: 'pointer' }} />)}</g>
  ), [sorted, L])
  const ring = useMemo(() => {
    const regionPath = (from: number, to: number, r: number) => {
      const a0 = L.angle(from - 0.45), a1 = L.angle(to + 0.45), large = a1 - a0 > Math.PI ? 1 : 0
      return `M${L.cx + r * Math.cos(a0)},${L.cy + r * Math.sin(a0)}A${r},${r} 0 ${large} 1 ${L.cx + r * Math.cos(a1)},${L.cy + r * Math.sin(a1)}`
    }
    return (
      <g>
        {regionRuns.map((run, k) => {
          const mid = (run.from + run.to) / 2, [tx, ty] = L.point(mid, L.r + 46)
          return <g key={`${run.name}-${run.from}`} data-testid={`${testId}-region-${run.from}`}>
            <path d={regionPath(run.from, run.to, L.r + 12)} stroke={REGION_COLORS[k % 2]} strokeWidth={7} fill="none" />
            <text x={tx} y={ty} fontSize={10} textAnchor="middle" dominantBaseline="middle" fill="#57534e">{run.name.replace(/^Helix /, '')}</text>
          </g>
        })}
        {positions.map((p) => { const [x, y] = L.point(p, L.r + 3); const off = absent?.has(p)
          return <circle key={p} cx={x} cy={y} r={1.8} fill={off ? 'white' : '#78716c'} stroke="#78716c" strokeWidth={off ? 0.7 : 0} data-testid={`${testId}-pos-${p}`} data-present={off ? 'false' : 'true'} /> })}
        {ticks.map((p) => { const [x, y] = L.point(p, L.r + 26), [x1, y1] = L.point(p, L.r + 6), [x2, y2] = L.point(p, L.r + 17)
          return <g key={p}><line x1={x1} y1={y1} x2={x2} y2={y2} stroke="#a8a29e" strokeWidth={0.8} />
            <text x={x} y={y} fontSize={9.5} textAnchor="middle" dominantBaseline="middle" fill="#0e7490" data-testid={`${testId}-tick-${p}`}>{labelFor(p).short}</text></g> })}
      </g>)
  }, [L, positions, regionRuns, ticks, labelFor, absent, testId])

  const sel = selectedId ? arcs.find((a) => a.id === selectedId) ?? null : null
  const hitAt = (e: React.MouseEvent) => (e.target as Element).getAttribute?.('data-hit')
  return (
    <div className="relative w-full max-w-[780px]">
      <svg ref={svgRef} viewBox={`0 0 ${size} ${size}`} className="block h-auto w-full select-none" role="img" aria-label={title} data-testid={testId}
        data-selected={selectedId ?? ''} data-arc-count={arcs.length}
        onMouseMove={(e) => { const id = hitAt(e); if (!id) { if (hover) setHover(null); return } if (hover?.id !== id) setHover({ id, anchor: new DOMRect(e.clientX - 4, e.clientY - 4, 8, 8) }) }}
        onMouseLeave={() => setHover(null)}
        onClick={(e) => { const id = hitAt(e); if (id) onSelect(id) }}>
        <circle cx={L.cx} cy={L.cy} r={L.r} fill="none" stroke="#d6d3d1" strokeWidth={0.6} />
        {ring}
        <g opacity={sel && dimUnselected ? 0.25 : 1}>{base}</g>
        {sel && <g data-testid={`${testId}-selected`} data-id={sel.id}>
          <path d={arcPath(L, sel.i, sel.j)} fill="none" stroke="#1c1917" strokeWidth={sel.width + 3.5} strokeOpacity={0.85} strokeLinecap="round" />
          <path d={arcPath(L, sel.i, sel.j)} fill="none" stroke={sel.color} strokeWidth={sel.width + 0.6} strokeDasharray={sel.dash} strokeLinecap="round" />
        </g>}
        {focus.map((p) => { const [x, y] = L.point(p, L.r + 3), [tx, ty] = L.point(p, L.r - 14)
          return <g key={`f${p}`} data-testid={`${testId}-focus-${p}`}><circle cx={x} cy={y} r={4} fill="#f59e0b" stroke="#1c1917" strokeWidth={0.8} />
            <text x={tx} y={ty} fontSize={10.5} fontWeight={600} textAnchor="middle" dominantBaseline="middle" fill="#1c1917" stroke="white" strokeWidth={2.5} paintOrder="stroke">{labelFor(p).text}</text></g> })}
        {hits}
      </svg>
      <HoverCard anchor={hover?.anchor ?? null} testId={`${testId}-hover`}>{hover ? tooltip(hover.id) : null}</HoverCard>
    </div>
  )
})
