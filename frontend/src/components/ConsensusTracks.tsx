import { useCallback, useEffect, useId, useLayoutEffect, useRef, useState, type MouseEvent, type ReactNode } from 'react'
import type { ConsensusStats } from '../types/backend'
import type { PositionState } from '../lib/consensus'
import { HoverCard, InfoTip } from './InfoTip'
import { loopTooltip, type InstanceLoop, type LoopTrack } from '../lib/loops'
import { boundaryText, markerWidth, type InsertionMarker } from '../lib/insertions'
import { authorResidueId, chainAt, referenceTicks, residueValue, type Numbering } from '../lib/numbering'

const C = { seq: '#2a78d6', str: '#1baf7a', high: '#1baf7a', medium: '#eda100', low: '#e34948', sel: '#6250d6', ink: '#0b0b0b', ink2: '#52514e', grid: '#e6e5e1', del: '#c3c2b7' }
const LEFT = 140, TOP = 4, GAP = 5, AXIS = 18
/** Shared horizontal scale. Every strip is an SVG in "track units" drawn at a fixed S px per unit (text and row heights never stretch); the
 *  width of one consensus position, `cell` (units), is the only thing zoom changes, so every track uses the same x(pos) = LEFT + (pos-1) * cell.
 *  Fit all: the whole model fits the visible width; zoom levels are px per position. */
const S = 1.2
const ZOOM_PX = [1, 1.5, 2, 3, 4, 6, 8, 11, 16, 24, 32]
/** level-of-detail thresholds (px per position): letters, +N insertion markers, inner gaps between columns */
const LOD = { letters: 7, markers: 4 }
/** axis label spacing: the smallest step whose labels are at least this far apart (px) */
const TICK_MIN_PX = 60
const TICK_STEPS = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
/** primary rows (track units): sequence and structural coverage bars, selected-chain state */
const ROWS = [16, 16, 10]
const rowY = (i: number) => TOP + ROWS.slice(0, i).reduce((a, b) => a + b + GAP, 0)
/** reference rows (inside "Reference and alignment details"): mean PP bars, RF and SS_cons characters */
const REF_ROWS = [10, 9, 9]
const refY = (i: number) => TOP + REF_ROWS.slice(0, i).reduce((a, b) => a + b + GAP, 0)

export interface ConsensusReference { rfamAcc: string; rfamId?: string | null; rf: string; ssCons: string; clanMembers: string[] }
interface Props {
  stats: ConsensusStats
  selected: ReadonlySet<number>
  /** click: select the position (or clear it when already selected); no range selection */
  onToggle: (pos: number) => void
  instanceStates?: PositionState[]
  instanceLabel?: string
  /** Rfam insertion markers of the selected chain (between consensus columns; never extra columns) */
  insertionMarkers?: InsertionMarker[]
  /** marker whose details are open; run shown in 2D/3D */
  openInsertion?: string | null
  selectedInsertionRun?: string | null
  onInsertionClick?: (markerId: string, rect: DOMRect) => void
  /** Rfam consensus model (RF, SS_cons), shown with the alignment confidence in the collapsed reference section */
  reference?: ConsensusReference
  /** per-instance BGSU loops of the selected chain (null until loaded); the section is collapsed by default and loads on first opening */
  loops?: InstanceLoop[] | null
  loopsState?: 'idle' | 'loading' | 'ready' | 'error'
  onLoopsOpen?: () => void
  selectedLoop?: string | null
  onSelectLoop?: (loopId: string) => void
  /** Position numbering reference (display only): axis labels and tooltip rows use its author residue numbers; null = Rfam consensus numbering */
  numbering?: Numbering | null
  /** the "Position numbering" selector, shown in the section heading */
  numberingControl?: ReactNode
  /** selected-chain tooltip (Rfam position, author residues, status); defaults to the status-only tooltip */
  chainTipFor?: (pos: number, s: PositionState) => Tip
}
const LOOP_COLOR: Record<LoopTrack, string> = { HL: '#2f8a7e', IL: '#7a6ab3', JL: '#b8792e' }
const LOOP_LABEL: Record<LoopTrack, string> = { HL: 'Hairpin loops (HL)', IL: 'Internal loops (IL)', JL: 'Junction loops (JL)' }
const LANE = 8, LANE_GAP = 1.5
/** selected loop: solid red inset outline on every block of the loop id; hover: thin dashed dark outline (never red) */
const SEL_RED = '#d1202f', SEL_STROKE = 1.5, HOVER_STROKE = 0.9
/** Rfam insertion markers: restrained amber (never the red selection outline unless shown in 2D/3D) */
const INS = { fill: '#fef3c7', openFill: '#b45309', stroke: '#b45309', text: '#92400e' }
/** per-position tooltip: heading, main value, optional muted detail, optional key-value rows (rows replace main when main is empty) */
export interface Tip { title: string; main: string; secondary?: string; rows?: { key: string; label: string; value: string }[] }
/** deposited author identifiers of an insertion run: "U51" or "U51–A53" (a run is consecutive in the polymer sequence); "unavailable" when absent */
const runAuthors = (r: { residues: Parameters<typeof authorResidueId>[0][] }) => {
  const ids = r.residues.map(authorResidueId)
  if (!ids.length) return ''
  const a = ids[0] ?? 'unavailable', b = ids[ids.length - 1] ?? 'unavailable'
  return ids.length === 1 ? a : `${a}–${b}`
}
/** percentage with at most one decimal: 1 -> 100%, 0.952 -> 95.2% */
const pct = (f: number) => `${(f * 100).toFixed(1).replace(/\.0$/, '')}%`
/** selected-chain status at one consensus match position (deleted / outside the aligned region are never called "no coordinates") */
const chainTipDefault = (pos: number, s: PositionState): Tip => ({ title: `Selected chain · Position ${pos}`,
  ...(s === 'match_observed' ? { main: 'Coordinates available' }
    : s === 'match_unobserved' ? { main: 'Coordinates unavailable', secondary: 'Nucleotide present in the deposited sequence' }
    : s === 'deletion' ? { main: 'Deleted in chain', secondary: 'No nucleotide aligned to this position' }
    : s === 'uncovered' ? { main: 'Outside aligned region' }
    : { main: 'Not available in the mapping' }) })
const ppColor = (v: number | null, t: { high: number; medium: number }) => (v == null ? C.grid : v >= t.high ? C.high : v >= t.medium ? C.medium : C.low)

/**
 * Tracks on the Rfam consensus axis. Three primary rows (sequence coverage, structural coverage, selected RNA chain), then a collapsed
 * "Reference and alignment details" section (alignment confidence, RF, SS_cons), then one shared position axis. Every part is its own SVG
 * with the same viewBox width and scaling, so all columns align; a click on any of them selects that position (click again to clear).
 */
export function ConsensusTracks({ stats, selected, onToggle, instanceStates, instanceLabel, insertionMarkers = [], openInsertion = null, selectedInsertionRun = null, onInsertionClick, reference, loops = null, loopsState = 'idle', onLoopsOpen, selectedLoop = null, onSelectLoop, numbering = null, numberingControl, chainTipFor }: Props) {
  const clen = stats.clen   // from the loaded family (never hard-coded)
  // --- shared horizontal scale (zoom + fit) and the single scroll container ---
  const scroller = useRef<HTMLDivElement>(null)
  const [viewW, setViewW] = useState(0)                       // visible width of the scroll container (px)
  const [zoomPx, setZoomPx] = useState<number | null>(null)   // px per position; null = Fit all
  const [sectionOpen, setSectionOpen] = useState(true)        // whole "Sequence and annotations" section (expanded by default)
  useLayoutEffect(() => {
    const el = scroller.current; if (!el) return
    setViewW(el.clientWidth)
    if (typeof ResizeObserver === 'undefined') return
    const ro = new ResizeObserver(() => setViewW(el.clientWidth)); ro.observe(el); return () => ro.disconnect()
  }, [sectionOpen])
  const fitPx = viewW > 0 ? Math.max(0.25, (viewW / S - LEFT - 12) / clen * S) : 9 * S
  const cellPx = zoomPx == null ? fitPx : Math.max(zoomPx, fitPx)
  const CELL = cellPx / S                                     // units per consensus position (shared by every strip)
  const width = LEFT + clen * CELL + 12
  const x = (pos: number) => LEFT + (pos - 1) * CELL
  const pad = cellPx >= 5 ? 1 : cellPx >= 2.5 ? 0.5 / S : 0 // inner gap between adjacent columns (units); none at overview scale
  const showLetters = cellPx >= LOD.letters, showMarkers = cellPx >= LOD.markers
  const tickStep = TICK_STEPS.find((t) => t * cellPx >= TICK_MIN_PX) ?? TICK_STEPS[TICK_STEPS.length - 1]
  const zoomIn = ZOOM_PX.find((z) => z > cellPx + 0.01)
  const zoomOutTo = [...ZOOM_PX].reverse().find((z) => z < cellPx - 0.01)
  const atFit = zoomPx == null || cellPx <= fitPx + 0.01
  /** keep a consensus position at the same place in the viewport across a zoom change (the selection when there is one, else the centre) */
  const anchor = useRef<{ pos: number; px: number } | null>(null)
  const setZoom = (z: number | null) => {
    const el = scroller.current
    if (el) {
      const sel = [...selected].sort((a, b) => a - b)
      const vc = (LEFT * S + el.clientWidth) / 2                                 // centre of the visible data area (px from the scroller's left edge)
      const pos = sel.length ? sel[Math.floor(sel.length / 2)] : ((el.scrollLeft + vc) / S - LEFT) / CELL + 1
      anchor.current = { pos, px: vc }
    }
    setZoomPx(z)
  }
  useLayoutEffect(() => {
    const el = scroller.current, a = anchor.current
    if (!el || !a) return
    anchor.current = null
    el.scrollLeft = Math.max(0, (LEFT + (a.pos - 1) * CELL) * S - a.px)
  }, [CELL])  // eslint-disable-line react-hooks/exhaustive-deps
  // a selection made elsewhere (2D, 3D, loop) that lies outside the visible part is scrolled into view; never re-zooms
  const selKey = [...selected].sort((a, b) => a - b).join(',')
  useEffect(() => {
    const el = scroller.current; if (!el || !selKey) return
    const ps = selKey.split(',').map(Number), lo = Math.min(...ps), hi = Math.max(...ps)
    const a = (LEFT + (lo - 1) * CELL) * S, b = (LEFT + hi * CELL) * S
    const v0 = el.scrollLeft + LEFT * S, v1 = el.scrollLeft + el.clientWidth
    if (a >= v0 && b <= v1) return
    const mid = (a + b) / 2, half = (el.clientWidth - LEFT * S) / 2
    el.scrollLeft = Math.max(0, b - a <= 2 * half ? mid - LEFT * S - half : a - LEFT * S - 8)
  }, [selKey])  // eslint-disable-line react-hooks/exhaustive-deps
  const [open, setOpen] = useState(false)
  const bodyId = useId()
  const [loopsOpen, setLoopsOpen] = useState(false)          // kept across instance changes (this component stays mounted)
  const loopsId = useId()
  const [hoverLoop, setHoverLoop] = useState<{ id: string; rect: DOMRect } | null>(null)
  const [hoverIns, setHoverIns] = useState<{ id: string; rect: DOMRect } | null>(null)
  const [trackInfo, setTrackInfo] = useState<DOMRect | null>(null)
  // per-position tooltip (structured: heading, main value, optional muted detail), shared by every track row
  const [tip, setTip] = useState<(Tip & { rect: DOMRect }) | null>(null)
  /** with a numbering reference, position tooltips (titled "· Rfam {pos}") add the reference residue (display only) */
  const withRef = (t: Tip, pos: number): Tip => (!numbering ? t : { ...t,
    rows: [...(t.rows ?? []), ...(t.rows?.some((r) => r.key === 'reference' || r.key === 'rfam') ? [] :   // the chain tooltip handles its own rows
      [{ key: 'reference', label: `${numbering.label} reference residue`, value: residueValue(chainAt(numbering.idx, pos)) }])] })
  const tipProps = (t0: Tip, pos?: number) => { const t = pos == null ? t0 : withRef(t0, pos); return { onMouseEnter: (e: MouseEvent<SVGElement>) => setTip({ ...t, rect: e.currentTarget.getBoundingClientRect() }), onMouseLeave: () => setTip(null) } }
  useEffect(() => { setHoverIns(null); setTip(null) }, [insertionMarkers, instanceStates])          // markers of the previous instance never leave a tooltip behind
  const onClick = useCallback((e: MouseEvent<SVGSVGElement>) => {
    const svg = e.currentTarget
    const pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY
    const m = svg.getScreenCTM(); if (!m) return
    const p = pt.matrixTransform(m.inverse())
    const pos = Math.floor((p.x - LEFT) / CELL) + 1
    if (pos >= 1 && pos <= clen && p.y >= 0) onToggle(pos)
  }, [clen, onToggle, CELL])
  const sl = stats.sequence_level, st = stats.structure_level, thr = stats.pp_display_thresholds
  /** one strip on the shared grid: the data SVG (whole width, same x(pos) for every strip) under a sticky label column that stays in view
   *  while the tracks scroll horizontally together */
  const strip = (height: number, testId: string, label: string, labels: ReactNode, children: ReactNode, clickable = true) => (
    <div className="flex" style={{ width: width * S }}>
      <svg width={LEFT * S} height={height * S} viewBox={`0 0 ${LEFT} ${height}`} className="sticky left-0 z-[3] block shrink-0 select-none bg-[#fcfcfb]"
        style={{ marginRight: -LEFT * S }} data-testid={`${testId}-labels`}>{labels}</svg>
      <svg width={width * S} height={height * S} viewBox={`0 0 ${width} ${height}`} className="block shrink-0 select-none" onClick={clickable ? onClick : undefined}
        data-testid={testId} data-cell={CELL} data-left={LEFT} data-cell-px={cellPx} role="img" aria-label={label}>
        {[...selected].map((pos) => <rect key={`s${pos}`} x={x(pos)} y={0} width={CELL} height={height} fill={C.sel} opacity={0.14} />)}
        {children}
      </svg>
    </div>
  )
  const label = (y0: number, row: number, name: string, sub: string) => (row >= 16
    ? <><text x={4} y={y0 + 8} fontSize={9.5} fill={C.ink}>{name}</text><text x={4} y={y0 + 16} fontSize={7.5} fill={C.ink2}>{sub}</text></>
    : <text x={4} y={y0 + 8} fontSize={8.5} fill={C.ink}>{name}<tspan dx={3} fontSize={7} fill={C.ink2}>{sub ? `· ${sub}` : ''}</tspan></text>)
  const bars = (y0: number, row: number, key: string, bar: (pos: number) => { h: number; fill: string; tip: Tip } | null) => (
    <g key={key}>
      <line x1={LEFT} x2={LEFT + clen * CELL} y1={y0 + row} y2={y0 + row} stroke={C.grid} />
      {Array.from({ length: clen }, (_, k) => k + 1).map((pos) => { const b = bar(pos); return b && <rect key={pos} x={x(pos) + pad} y={y0 + row - b.h * row} width={CELL - 2 * pad} height={Math.max(b.h * row, 0.5)} fill={b.fill} data-row={key} data-tip-pos={pos} {...tipProps(b.tip, pos)} /> })}
    </g>
  )
  const nPrimary = instanceStates ? 3 : 2
  const mainH = rowY(nPrimary) - GAP + 2
  const refH = refY(reference ? 3 : 1) - GAP + 2
  return (
    <div className="rounded border border-stone-200 bg-[#fcfcfb]" data-testid="consensus-track-area" data-open={String(sectionOpen)}>
      {/* section heading: collapse the whole track area; zoom the shared Rfam axis */}
      <div className="flex flex-wrap items-center gap-x-2 gap-y-1 px-1 py-0.5" data-testid="tracks-heading">
        <button type="button" aria-expanded={sectionOpen} aria-controls={`${bodyId}-tracks`} onClick={() => { setSectionOpen((o) => !o); setTip(null); setHoverIns(null); setHoverLoop(null) }} data-testid="tracks-toggle"
          className="flex items-center gap-1.5 rounded px-1 py-0.5 text-left text-[11px] font-semibold text-stone-800 hover:bg-stone-100 focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500">
          <span aria-hidden="true" className={`inline-block text-[9px] text-stone-500 transition-transform ${sectionOpen ? 'rotate-90' : ''}`}>▶</span>Sequence and annotations
        </button>
        {sectionOpen && numberingControl}
        {sectionOpen && (
          <div className="ml-auto flex items-center gap-1 text-[11px]" role="group" aria-label="Zoom the Rfam position axis" data-testid="zoom-controls" data-cell-px={cellPx.toFixed(2)} data-fit={String(atFit)}>
            <button type="button" onClick={() => setZoom(zoomOutTo != null && zoomOutTo > fitPx + 0.01 ? zoomOutTo : null)} disabled={atFit} aria-label="Zoom out" title="Zoom out" data-testid="zoom-out"
              className="h-5 w-5 rounded border border-stone-300 bg-white leading-none text-stone-700 hover:bg-stone-50 disabled:opacity-40">−</button>
            <button type="button" onClick={() => zoomIn != null && setZoom(zoomIn)} disabled={zoomIn == null} aria-label="Zoom in" title="Zoom in" data-testid="zoom-in"
              className="h-5 w-5 rounded border border-stone-300 bg-white leading-none text-stone-700 hover:bg-stone-50 disabled:opacity-40">+</button>
            <button type="button" onClick={() => setZoom(null)} disabled={atFit} title="Fit all Rfam positions" data-testid="zoom-fit"
              className="rounded border border-stone-300 bg-white px-1.5 py-px text-stone-700 hover:bg-stone-50 disabled:opacity-40">Fit all</button>
          </div>
        )}
      </div>
      {sectionOpen && <div ref={scroller} id={`${bodyId}-tracks`} className="overflow-x-auto overflow-y-hidden border-t border-stone-200" data-testid="tracks-scroll">
      {strip(mainH, 'consensus-tracks', 'Rfam consensus tracks', <>
        {label(rowY(0), ROWS[0], 'Sequence coverage', `${sl.n_unique_sequences} unique sequences`)}
        {label(rowY(1), ROWS[1], 'Structural coverage', `${st.n_instances} instances`)}
        {instanceStates && <>
          {label(rowY(2), ROWS[2], instanceLabel ?? 'Instance', 'state per position')}
          {/* the coordinate-system note of this track */}
          <g className="cursor-help" tabIndex={0} role="button" aria-label="About the selected RNA chain track" data-testid="instance-track-info"
            onMouseEnter={(e) => setTrackInfo(e.currentTarget.getBoundingClientRect())} onMouseLeave={() => setTrackInfo(null)} onFocus={(e) => setTrackInfo(e.currentTarget.getBoundingClientRect())} onBlur={() => setTrackInfo(null)}>
            <circle cx={LEFT - 9} cy={rowY(2) + ROWS[2] / 2} r={3.6} fill="#fff" stroke={C.ink2} strokeWidth={0.7} />
            <text x={LEFT - 9} y={rowY(2) + ROWS[2] / 2 + 2} fontSize={5.5} textAnchor="middle" fill={C.ink2} fontWeight={600}>i</text>
          </g>
        </>}
      </>, <>
        {bars(rowY(0), ROWS[0], 'seq', (pos) => { const c = sl.coverage[pos - 1]; if (c == null) return null
          const del = sl.n_deletion[pos - 1], unc = sl.n_uncovered[pos - 1]
          return { h: c, fill: C.seq, tip: { title: `Sequence coverage · Rfam ${pos}`, main: `${sl.n_match[pos - 1]}/${sl.n_unique_sequences} unique sequences (${pct(c)})`,
            secondary: del || unc ? [del ? `deleted in ${del}` : '', unc ? `not covered by the alignment in ${unc}` : ''].filter(Boolean).join(' · ').replace(/^./, (m) => m.toUpperCase()) : undefined } } })}
        {bars(rowY(1), ROWS[1], 'str', (pos) => { const c = st.structural_coverage[pos - 1]; if (c == null) return null
          const obs = st.n_match_with_coordinates[pos - 1], inSeq = st.n_match[pos - 1]
          return { h: c, fill: C.str, tip: { title: `Structural coverage · Rfam ${pos}`, main: `${obs}/${st.n_instances} RNA instances (${pct(c)})`,
            secondary: inSeq !== obs ? `${inSeq} instance${inSeq === 1 ? '' : 's'} contain the nucleotide in sequence` : undefined } } })}
        {instanceStates && (() => {
          const y0 = rowY(2), ROW = ROWS[2]
          return (
            <g data-testid="instance-track" data-markers-shown={String(showMarkers)}>
              {instanceStates.map((s, k) => {
                const pos = k + 1, cx = x(pos) + pad, w = CELL - 2 * pad
                const chainTip = (p: number, st: PositionState) => (chainTipFor ? chainTipFor(p, st) : chainTipDefault(p, st))
                if (s === 'match_observed') return <rect key={pos} x={cx} y={y0 + 1} width={w} height={ROW - 2} fill={C.str} data-row="chain" data-tip-pos={pos} {...tipProps(chainTip(pos, s), pos)} />
                if (s === 'match_unobserved') return <rect key={pos} x={cx + Math.min(0.5, w / 4)} y={y0 + 1.5} width={Math.max(w - 1, w / 2)} height={ROW - 3} fill="none" pointerEvents="all" stroke={C.str} strokeDasharray="2 1.5" data-row="chain" data-tip-pos={pos} {...tipProps(chainTip(pos, s), pos)} />
                if (s === 'deletion') return <rect key={pos} x={cx} y={y0 + ROW / 2 - 1.5} width={w} height={3} fill={C.del} data-row="chain" data-tip-pos={pos} {...tipProps(chainTip(pos, s), pos)} />
                return <rect key={pos} x={cx} y={y0 + 1} width={w} height={ROW - 2} fill="url(#hatch)" data-row="chain" data-tip-pos={pos} {...tipProps(chainTip(pos, s), pos)} />
              })}
              {/* Rfam insertions: "+N" markers ON the boundary between two consensus columns (never an extra column; the grid does not move) */}
              {showMarkers && insertionMarkers.map((mk) => {
                const bx = LEFT + (mk.boundaries.reduce((a, b) => a + b, 0) / mk.boundaries.length) * CELL
                const w = markerWidth(mk.length), h = 9, py = y0 - 5.5
                const open = mk.id === openInsertion, chosen = mk.runs.some((r) => r.id === selectedInsertionRun)
                const localOnly = mk.runs.every((r) => r.kind === 'local_end')
                const label = `${mk.runs.length > 1 ? `${mk.runs.length} insertion runs` : localOnly ? 'Local-end residues' : 'Insertion'}, ${mk.length} nucleotides, ${mk.runs.map((r) => boundaryText(r).toLowerCase()).join('; ')}`
                return (
                  <g key={mk.id} data-testid={`insertion-marker-${mk.id}`} data-insertion-marker="" data-boundaries={mk.boundaries.join(',')} data-length={mk.length} data-open={String(open)} data-selected={String(chosen)}
                    role="button" tabIndex={0} aria-label={label} aria-expanded={open} className="cursor-pointer focus:outline-none"
                    onMouseEnter={(e) => setHoverIns({ id: mk.id, rect: e.currentTarget.getBoundingClientRect() })} onMouseLeave={() => setHoverIns(null)}
                    onFocus={(e) => setHoverIns({ id: mk.id, rect: e.currentTarget.getBoundingClientRect() })} onBlur={() => setHoverIns(null)}
                    onClick={(e) => { e.stopPropagation(); setHoverIns(null); onInsertionClick?.(mk.id, e.currentTarget.getBoundingClientRect()) }}
                    onKeyDown={(e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); setHoverIns(null); onInsertionClick?.(mk.id, e.currentTarget.getBoundingClientRect()) } }}>
                    {/* hit target larger than the visible label (gap above the row + the row itself) */}
                    <rect x={bx - Math.max(w, 8) / 2} y={py} width={Math.max(w, 8)} height={ROW + 5} fill="transparent" />
                    <line x1={bx} x2={bx} y1={py + h} y2={y0 + ROW} stroke={INS.stroke} strokeWidth={1.1} />
                    <rect x={bx - w / 2} y={py} width={w} height={h} rx={2} fill={open ? INS.openFill : INS.fill} stroke={chosen ? SEL_RED : INS.stroke} strokeWidth={chosen ? 1.2 : 0.7} strokeDasharray={localOnly ? '1.5 1' : undefined} />
                    <text x={bx} y={py + 6.9} fontSize={7.5} textAnchor="middle" fontWeight={600} fill={open ? '#fff' : INS.text} pointerEvents="none">+{mk.length}</text>
                  </g>)
              })}
            </g>
          )
        })()}
        <defs><pattern id="hatch" width={4} height={4} patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line x1={0} y1={0} x2={0} y2={4} stroke={C.del} strokeWidth={1.2} /></pattern></defs>
      </>)}
      {tip && <HoverCard anchor={tip.rect} testId="track-tooltip">
        <div className="font-semibold" data-testid="track-tooltip-title">{tip.title}</div>
        {tip.main && <div data-testid="track-tooltip-main">{tip.main}</div>}
        {tip.secondary && <div className="text-[11px] text-stone-500" data-testid="track-tooltip-secondary">{tip.secondary}</div>}
        {tip.rows && <dl className="mt-0.5 grid grid-cols-[auto_auto] gap-x-2 text-[11px] leading-4">{tip.rows.map((r) => <div key={r.key} className="contents" data-testid={`track-tooltip-row-${r.key}`}>
          <dt className="text-stone-500">{r.label}</dt><dd className={`${r.key === 'status' ? '' : 'font-mono '}${r.key === 'reference' ? 'text-violet-800' : 'text-stone-900'}`} data-testid={`track-tooltip-value-${r.key}`}>{r.value}</dd></div>)}</dl>}
      </HoverCard>}
      {trackInfo && <HoverCard anchor={trackInfo} testId="instance-track-info-card"><div className="font-semibold">Selected RNA chain</div><div data-testid="instance-track-info-text">Shows nucleotide availability at Rfam match positions. +N markers indicate insertions between positions.</div><div className="text-[11px] text-stone-500">Symbols are explained in Legend &amp; help.</div></HoverCard>}
      {(() => { const mk = hoverIns && insertionMarkers.find((m) => m.id === hoverIns.id); if (!mk || mk.id === openInsertion) return null
        const one = mk.runs.length === 1 ? mk.runs[0] : null
        return <HoverCard anchor={hoverIns!.rect} testId="insertion-tooltip">
          {one ? <>
            <div className="font-semibold" data-testid="insertion-tooltip-title">{one.kind === 'insert' ? 'Insertion' : 'Local-end residues'} · +{one.residues.length} nt</div>
            <div data-testid="insertion-tooltip-boundary">{boundaryText(one)}</div>
            <div data-testid="insertion-tooltip-sequence">Sequence: <span className="font-mono">{one.sequence}</span></div>
            <div data-testid="insertion-tooltip-author">{instanceLabel ? `${instanceLabel} author residues` : 'Author residues'}: <span className="font-mono">{runAuthors(one)}</span></div>
          </> : <>
            <div className="font-semibold" data-testid="insertion-tooltip-title">Insertions · +{mk.length} nt in {mk.runs.length} runs</div>
            {mk.runs.map((r) => <div key={r.id} data-testid="insertion-tooltip-run">{boundaryText(r)}: <span className="font-mono">{r.sequence}</span> (+{r.residues.length}; <span className="font-mono">{runAuthors(r)}</span>)</div>)}
          </>}
          <div className="mt-0.5 text-[10px] text-stone-500">Click for residue details</div>
        </HoverCard> })()}
      {/* RNA loop annotations of the selected chain: collapsed by default; loads loops.json on first opening */}
      <div className="sticky left-0 flex items-center gap-1 border-t border-stone-200 px-1" style={{ width: viewW || undefined }} data-testid="loops-section" data-open={String(loopsOpen)}>
        <button type="button" aria-expanded={loopsOpen} aria-controls={loopsId} onClick={() => { if (!loopsOpen) onLoopsOpen?.(); setLoopsOpen(!loopsOpen) }} data-testid="loops-toggle"
          className="flex items-center gap-1.5 rounded px-1 py-0.5 text-left text-[11px] font-medium text-stone-800 hover:bg-stone-100 focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500">
          <span aria-hidden="true" className={`inline-block text-[9px] text-stone-500 transition-transform ${loopsOpen ? 'rotate-90' : ''}`}>▶</span>RNA loop annotations
        </button>
        <InfoTip label="About the RNA loop annotations" title="RNA loop annotations" text="BGSU hairpin, internal and junction loops in this RNA chain. Select a block to highlight the full loop, including discontinuous segments." />
      </div>
      {loopsOpen && (
        <div id={loopsId} data-testid="loops-body">
          {loopsState === 'loading' && <p className="sticky left-0 px-2 text-[10px] text-stone-500">Loading loop annotations…</p>}
          {loopsState === 'error' && <p className="sticky left-0 px-2 text-[10px] text-amber-700">Loop annotations unavailable.</p>}
          {loops && (() => {
            const tracks: LoopTrack[] = ['HL', 'IL', 'JL']
            const lanes = (t: LoopTrack) => Math.max(1, ...loops.filter((l) => l.track === t && l.blocks.length).map((l) => l.lane + 1))
            const trackH = (t: LoopTrack) => lanes(t) * (LANE + LANE_GAP) - LANE_GAP
            const ys: number[] = []; let yy = TOP
            for (const t of tracks) { ys.push(yy); yy += trackH(t) + GAP }
            const H = yy - GAP + 2
            const offGrid = loops.filter((l) => !l.blocks.length)
            return <>
              {strip(H, 'loop-tracks', 'RNA loop annotations of the selected chain', <>
                {tracks.map((t, i) => <text key={t} x={4} y={ys[i] + 7} fontSize={8.5} fill={C.ink}>{LOOP_LABEL[t]}</text>)}
              </>, <>
                {tracks.map((t, i) => (
                  <g key={t} data-testid={`loop-track-${t}`} data-lanes={lanes(t)}>
                    <line x1={LEFT} x2={LEFT + clen * CELL} y1={ys[i] + trackH(t) + 0.5} y2={ys[i] + trackH(t) + 0.5} stroke={C.grid} />
                    {loops.filter((l) => l.track === t).flatMap((l) => l.blocks.map((b) => {
                      const sel = l.id === selectedLoop, hov = l.id === hoverLoop?.id
                      // one block per original segment, over that segment's own columns only (never merged into one interval); inset by a small gap
                      const g = Math.min(0.75, pad), bx = x(b.from) + g, by = ys[i] + l.lane * (LANE + LANE_GAP), bw = Math.max((b.to - b.from + 1) * CELL - 2 * g, 0.4)
                      const tip = loopTooltip(l)
                      // outlines are drawn INSIDE the block (inset by half the stroke): block geometry and the gaps to adjacent blocks never change
                      const sw = sel ? SEL_STROKE : HOVER_STROKE
                      return (
                        <g key={`${l.id}-${b.strand}`}>
                          <rect x={bx} y={by} width={bw} height={LANE} rx={Math.min(2.5, bw / 2)} fill={LOOP_COLOR[t]} opacity={hoverLoop && !hov && !sel ? 0.45 : 1}
                            className="cursor-pointer" data-testid={`loop-block-${l.id}-${b.strand}`} data-loop-id={l.id} data-from={b.from} data-to={b.to} data-lane={l.lane} data-selected={String(sel)} data-hover={String(hov)}
                            aria-label={`${tip.title} ${tip.id}, ${tip.counts}`}
                            onMouseEnter={(e) => setHoverLoop({ id: l.id, rect: e.currentTarget.getBoundingClientRect() })} onMouseLeave={() => setHoverLoop(null)} onClick={(e) => { e.stopPropagation(); onSelectLoop?.(l.id) }} />
                          {(sel || hov) && <rect x={bx + Math.min(sw / 2, bw / 4)} y={by + sw / 2} width={Math.max(bw - sw, bw / 2)} height={LANE - sw} rx={Math.min(2, bw / 2)} fill="none" pointerEvents="none"
                            stroke={sel ? SEL_RED : C.ink} strokeWidth={sw} strokeDasharray={sel ? undefined : '2 1.2'} data-testid={sel ? `loop-outline-${l.id}-${b.strand}` : `loop-hover-${l.id}-${b.strand}`} />}
                        </g>)
                    }))}
                  </g>))}
              </>, false)}
              {(() => { const l = hoverLoop && loops.find((x) => x.id === hoverLoop.id); if (!l) return null; const tip = loopTooltip(l)
                return <HoverCard anchor={hoverLoop!.rect} testId="loop-tooltip">
                  <div className="font-semibold" data-testid="loop-tooltip-title">{tip.title}</div>
                  <div className="font-mono text-[11px]" data-testid="loop-tooltip-id">{tip.id}</div>
                  <div className="text-stone-600" data-testid="loop-tooltip-counts">{tip.counts}{tip.otherChain ? ` · ${tip.otherChain}` : ''}</div>
                  {tip.merged && <div className="text-stone-600" data-testid="loop-tooltip-merged">{tip.merged}</div>}
                  {tip.motif && <div className="mt-0.5" data-testid="loop-tooltip-motif">{tip.motif}{tip.motifNote && <div className="text-[10px] text-stone-500">{tip.motifNote}</div>}</div>}
                </HoverCard> })()}
              {offGrid.length > 0 && <p className="sticky left-0 px-2 text-[10px] text-stone-500" style={{ width: viewW || undefined }} data-testid="loops-off-grid">Not on the consensus grid: {offGrid.map((l) => <button key={l.id} type="button" aria-pressed={l.id === selectedLoop} className={`mr-1 rounded px-0.5 underline ${l.id === selectedLoop ? 'outline outline-[1.5px] outline-[#d1202f]' : ''}`} onClick={() => onSelectLoop?.(l.id)}>{l.id}</button>)}</p>}
            </>
          })()}
        </div>
      )}
      {/* collapsed by default: heading and chevron only */}
      <div className="sticky left-0 flex items-center gap-1 border-t border-stone-200 px-1" style={{ width: viewW || undefined }} data-testid="reference-section" data-open={String(open)}>
        <button type="button" aria-expanded={open} aria-controls={bodyId} onClick={() => setOpen(!open)} data-testid="reference-toggle"
          className="flex items-center gap-1.5 rounded px-1 py-0.5 text-left text-[11px] font-medium text-stone-800 hover:bg-stone-100 focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500">
          <span aria-hidden="true" className={`inline-block text-[9px] text-stone-500 transition-transform ${open ? 'rotate-90' : ''}`}>▶</span>Reference and alignment details
        </button>
        <InfoTip label="About the reference and alignment details" title="Reference and alignment details" text="Rfam consensus sequence, reference secondary structure and alignment confidence." />
      </div>
      {open && (
        <div id={bodyId} data-testid="reference-body">
          {strip(refH, 'reference-tracks', 'Alignment confidence and Rfam consensus model', <>
            {label(refY(0), REF_ROWS[0], 'Alignment confidence', 'mean PP')}
            {reference && <>{label(refY(1), REF_ROWS[1], 'Consensus sequence (RF)', showLetters ? '' : 'zoom in to read')}{label(refY(2), REF_ROWS[2], 'Secondary structure (SS_cons)', showLetters ? '' : 'zoom in to read')}</>}
          </>, <>
            {bars(refY(0), REF_ROWS[0], 'pp', (pos) => { const v = sl.pp_mean[pos - 1]; return v == null ? null : { h: v, fill: ppColor(v, thr), tip: { title: `Alignment confidence · Rfam ${pos}`, main: `Mean cmalign PP ${v.toFixed(3)} (${v >= thr.high ? 'high' : v >= thr.medium ? 'medium' : 'low'})`, secondary: `${sl.pp_n[pos - 1]} aligned unique sequences` } } })}
            {reference && ([['rf', 'Consensus sequence (RF)', reference.rf], ['ss', 'Secondary structure (SS_cons)', reference.ssCons]] as const).map(([key, name, line], i) => {
              const y0 = refY(i + 1)
              return (
                <g key={key} data-testid={`reference-row-${key}`} data-letters={String(showLetters)}>
                  {/* letters only where a column is wide enough to read them; otherwise a hover target per position (zoom in to read) */}
                  {[...line].slice(0, clen).map((ch, k) => showLetters
                    ? <text key={k} x={x(k + 1) + CELL / 2} y={y0 + 7.5} fontSize={7.5} fontFamily="ui-monospace, Menlo, monospace" textAnchor="middle" fill={ch === ch.toUpperCase() && /[A-Z]/.test(ch) ? C.ink : C.ink2} data-row={key} data-tip-pos={k + 1} {...tipProps({ title: `${name} · Position ${k + 1}`, main: ch }, k + 1)}>{ch}</text>
                    : <rect key={k} x={x(k + 1)} y={y0} width={CELL} height={REF_ROWS[i + 1]} fill="transparent" data-row={key} data-tip-pos={k + 1} {...tipProps({ title: `${name} · Position ${k + 1}`, main: ch }, k + 1)} />)}
                </g>
              )
            })}
          </>)}
        </div>
      )}
      {/* one shared Rfam position axis, below the tracks and the optional reference section */}
      {strip(AXIS, 'consensus-axis', numbering ? `${numbering.label} author residue numbering axis` : `Rfam ${stats.rfam_acc} position axis`, <>
        {numbering
          ? <text x={4} y={AXIS - 4} fontSize={8} fill="#5b21b6" data-testid="axis-convention">{numbering.label} · author numbering</text>
          : <text x={4} y={AXIS - 4} fontSize={8} fill={C.ink2} data-testid="axis-convention">Rfam {stats.rfam_acc} position</text>}
      </>, numbering ? <>
        {/* reference numbering: labels only at positions where the reference chain has a match residue with an author number (no offsets, no
            invented numbers); the columns are the Rfam columns, unchanged */}
        {referenceTicks(numbering.numbers, cellPx).map((t) => (
          <g key={t.pos} data-ref-tick={t.pos}><line x1={x(t.pos) + CELL / 2} x2={x(t.pos) + CELL / 2} y1={2} y2={6} stroke="#5b21b6" /><text x={x(t.pos) + CELL / 2} y={AXIS - 4} fontSize={8} textAnchor="middle" fill="#5b21b6">{t.text}</text></g>))}
      </> : <>
        {/* tick spacing adapts to the zoom so labels never overlap; numbering is always the Rfam consensus numbering */}
        {Array.from({ length: Math.floor(clen / tickStep) + 1 }, (_, k) => k * tickStep || 1).filter((pos) => pos <= clen).map((pos) => (
          <g key={pos}><line x1={x(pos) + CELL / 2} x2={x(pos) + CELL / 2} y1={2} y2={6} stroke={C.ink2} /><text x={x(pos) + CELL / 2} y={AXIS - 4} fontSize={8} textAnchor="middle" fill={C.ink2}>{pos}</text></g>))}
      </>)}
      </div>}
    </div>
  )
}
