import { useCallback, type MouseEvent } from 'react'
import type { ConsensusStats } from '../types/backend'
import type { PositionState } from '../lib/consensus'

const C = { seq: '#2a78d6', str: '#1baf7a', high: '#1baf7a', medium: '#eda100', low: '#e34948', sel: '#6250d6', ink: '#0b0b0b', ink2: '#52514e', grid: '#e6e5e1', del: '#c3c2b7' }
const CELL = 9, LEFT = 118, TOP = 18, ROW = 34, GAP = 6
interface Props {
  stats: ConsensusStats
  selected: ReadonlySet<number>
  onToggle: (pos: number, shiftKey: boolean) => void
  instanceStates?: PositionState[]
  instanceLabel?: string
  insertions?: Map<number, { insert: number; local_end: number }>
}
const ppColor = (v: number | null, t: { high: number; medium: number }) => (v == null ? C.grid : v >= t.high ? C.high : v >= t.medium ? C.medium : C.low)

/** SVG tracks on the Rfam consensus axis. Click toggles a position, shift-click extends from the last click. */
export function ConsensusTracks({ stats, selected, onToggle, instanceStates, instanceLabel, insertions }: Props) {
  const clen = stats.clen
  const width = LEFT + clen * CELL + 12
  const rows = instanceStates ? 4 : 3
  const height = TOP + rows * (ROW + GAP) + 40
  const x = (pos: number) => LEFT + (pos - 1) * CELL
  const onClick = useCallback((e: MouseEvent<SVGSVGElement>) => {
    const svg = e.currentTarget
    const pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY
    const m = svg.getScreenCTM(); if (!m) return
    const p = pt.matrixTransform(m.inverse())
    const pos = Math.floor((p.x - LEFT) / CELL) + 1
    if (pos >= 1 && pos <= clen && p.y >= TOP) onToggle(pos, e.shiftKey)
  }, [clen, onToggle])
  const sl = stats.sequence_level, st = stats.structure_level, thr = stats.pp_display_thresholds
  const track = (i: number, label: string, sub: string, bars: (pos: number) => { h: number; fill: string; title: string } | null) => {
    const y0 = TOP + i * (ROW + GAP)
    return (
      <g key={label}>
        <text x={4} y={y0 + 14} fontSize={11} fill={C.ink}>{label}</text>
        <text x={4} y={y0 + 27} fontSize={9} fill={C.ink2}>{sub}</text>
        <line x1={LEFT} x2={LEFT + clen * CELL} y1={y0 + ROW} y2={y0 + ROW} stroke={C.grid} />
        {Array.from({ length: clen }, (_, k) => k + 1).map((pos) => { const b = bars(pos); return b && <rect key={pos} x={x(pos) + 1} y={y0 + ROW - b.h} width={CELL - 2} height={Math.max(b.h, 0.5)} fill={b.fill}><title>{b.title}</title></rect> })}
      </g>
    )
  }
  return (
    <svg viewBox={`0 0 ${width} ${height}`} className="w-full select-none" style={{ maxHeight: height * 1.4 }} onClick={onClick} data-testid="consensus-tracks" role="img" aria-label="Rfam consensus tracks">
      <rect x={0} y={0} width={width} height={height} fill="#fcfcfb" />
      {/* selection backdrop */}
      {[...selected].map((pos) => <rect key={pos} x={x(pos)} y={TOP - 4} width={CELL} height={height - TOP - 18} fill={C.sel} opacity={0.14} />)}
      {/* axis */}
      {Array.from({ length: Math.floor(clen / 10) + 1 }, (_, k) => k * 10 || 1).map((pos) => (
        <g key={pos}><line x1={x(pos) + CELL / 2} x2={x(pos) + CELL / 2} y1={height - 30} y2={height - 26} stroke={C.ink2} /><text x={x(pos) + CELL / 2} y={height - 15} fontSize={9} textAnchor="middle" fill={C.ink2}>{pos}</text></g>))}
      <text x={LEFT} y={height - 3} fontSize={9} fill={C.ink2}>Rfam {stats.rfam_acc} consensus position (1–{clen}); click = toggle, shift-click = range</text>
      {track(0, 'Sequence coverage', `${sl.n_unique_sequences} unique seqs`, (pos) => { const c = sl.coverage[pos - 1]; return c == null ? null : { h: c * (ROW - 4), fill: C.seq, title: `Consensus ${pos}: sequence coverage ${sl.n_match[pos - 1]} / ${sl.n_unique_sequences} unique sequences = ${(c * 100).toFixed(1)} % (deleted in ${sl.n_deletion[pos - 1]}, uncovered in ${sl.n_uncovered[pos - 1]})` } })}
      {track(1, 'Structural coverage', `${st.n_instances} instances, coords`, (pos) => { const c = st.structural_coverage[pos - 1]; return c == null ? null : { h: c * (ROW - 4), fill: C.str, title: `Consensus ${pos}: structural coverage ${st.n_match_with_coordinates[pos - 1]} / ${st.n_instances} instances with observed coordinates = ${(c * 100).toFixed(1)} % (${st.n_match[pos - 1]} instances have the residue in sequence)` } })}
      {track(2, 'Alignment confidence', 'mean PP (seq level)', (pos) => { const v = sl.pp_mean[pos - 1]; return v == null ? null : { h: v * (ROW - 4), fill: ppColor(v, thr), title: `Consensus ${pos}: mean Infernal posterior probability ${v.toFixed(3)} over ${sl.pp_n[pos - 1]} aligned unique sequences; category ${v >= thr.high ? 'high (green)' : v >= thr.medium ? 'medium (amber)' : 'low (red)'} (provisional bands: high ≥ ${thr.high}, medium ≥ ${thr.medium})` } })}
      {instanceStates && (() => {
        const y0 = TOP + 3 * (ROW + GAP)
        return (
          <g>
            <text x={4} y={y0 + 14} fontSize={11} fill={C.ink}>{instanceLabel ?? 'Instance'}</text>
            <text x={4} y={y0 + 27} fontSize={9} fill={C.ink2}>state per position</text>
            {instanceStates.map((s, k) => {
              const pos = k + 1, cx = x(pos) + 1, w = CELL - 2
              if (s === 'match_observed') return <rect key={pos} x={cx} y={y0 + 8} width={w} height={18} fill={C.str}><title>{`pos ${pos}: residue with coordinates`}</title></rect>
              if (s === 'match_unobserved') return <rect key={pos} x={cx + 0.5} y={y0 + 8.5} width={w - 1} height={17} fill="none" stroke={C.str} strokeDasharray="2 1.5"><title>{`pos ${pos}: residue present in sequence, no coordinates`}</title></rect>
              if (s === 'deletion') return <rect key={pos} x={cx} y={y0 + 15} width={w} height={4} fill={C.del}><title>{`pos ${pos}: deleted in this chain`}</title></rect>
              return <rect key={pos} x={cx} y={y0 + 8} width={w} height={18} fill="url(#hatch)"><title>{`pos ${pos}: outside the aligned model span (uncovered)`}</title></rect>
            })}
            {insertions && [...insertions.entries()].map(([anchor, n]) => (
              <g key={anchor}><polygon points={`${x(anchor) + CELL},${y0 + 4} ${x(anchor) + CELL - 4},${y0 - 2} ${x(anchor) + CELL + 4},${y0 - 2}`} fill={n.local_end ? C.low : C.seq}><title>{`after consensus ${anchor}: ${n.insert} inserted, ${n.local_end} local-end residues (no consensus coordinate)`}</title></polygon></g>))}
          </g>
        )
      })()}
      <defs><pattern id="hatch" width={4} height={4} patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line x1={0} y1={0} x2={0} y2={4} stroke={C.del} strokeWidth={1.2} /></pattern></defs>
    </svg>
  )
}
