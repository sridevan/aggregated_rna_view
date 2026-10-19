import { useCallback, useEffect, useMemo, useRef, useState, type ReactNode } from 'react'
import { loadInteractions, loadStructureClustering } from '../api/data'
import { useFamily } from '../state/family'
import { useMapping } from '../state/useMapping'
import { buildNumbering, chainAt, type Numbering } from '../lib/numbering'
import { indexMapping, type MappingIndex } from '../lib/consensus'
import {
  CATEGORIES, CATEGORY_CODE, CATEGORY_STYLE, CATEGORY_TEXT, categoryAvailability, categoryColor, circleLayout, compareInstances, COMPARE_TEXT, COMPARE_TIP, DEFAULT_FILTER,
  classBreakdown, defaultComparePair, defaultFamilyOrder, DEFAULT_PAIR_MAX_RESOLUTION, helixCoverage, medianHelixCoverage, dominantClass, groupByEndpoints, groupId, groupOccurrence, indexInstance, matchesPositions, MULTI_RAMP, multiColor, pairIsCanonical,
  pairLabel, parseSearch, passes, positionLabel, showCircular, CIRCULAR_MAX_CLEN, STATE_TEXT, statusIn, supportingInstances, type CompareRow, type CompareStatus,
  type InstanceStatus, type InteractionFilter, type InteractionGroup, type SearchQuery,
} from '../lib/interactions'
import type { InstanceRow, InteractionCategory, InteractionPair, InteractionsDataset } from '../types/backend'
import { CircularDiagram, type DiagramArc } from './interactions/CircularDiagram'
import { InteractionPairViewer, type FitResult } from './interactions/InteractionPairViewer'
import { PdbeMolstarViewer } from './PdbeMolstarViewer'
import { InfoTip } from './InfoTip'

/*
 * Interactions tab: RNA-RNA interactions annotated by BGSU FR3D (base pairing Bp, base stacking Bs, base-ribose Br, base-phosphate Bph) on
 * Rfam consensus positions (interactions.json).
 * Family conservation (default): one circular diagram of every consensus-mapped interaction observed in the family (families up to
 * CIRCULAR_MAX_CLEN consensus positions), its table with observed / eligible chains, and the status of each interaction in the DISPLAYED 3D
 * structure (Present / Absent / Not comparable).
 * Compare two: two diagrams on the same layout, the comparison table and both chains in 3D with an interaction-local superposition; the
 * default pair has the highest valid geometric discrepancy.
 * Interactions between the same two nucleotides are shown together as "Multiple interactions" (display grouping; every member keeps its
 * own id, class and frequency). Two separate choices: the numbering reference (labels only) and the displayed structure (3D and status).
 */
const SIZE = 560, RADIUS = 204, ARC_WIDTH = 1.8
const COLOR_A = '#0072B2', COLOR_B = '#E69F00'
const COMPARE_COLORS: Record<CompareStatus, { color: string; order: number }> = {
  shared: { color: '#a8a29e', order: 0 }, not_comparable: { color: '#c9b98f', order: 1 },
  reclassified: { color: '#CC79A7', order: 2 }, only_a: { color: '#D55E00', order: 3 }, only_b: { color: '#D55E00', order: 3 },
}
const COMPARE_ORDER: CompareStatus[] = ['only_a', 'only_b', 'reclassified', 'not_comparable', 'shared']
const STATE_STYLE: Record<string, string> = { present: 'text-emerald-800', absent: 'text-stone-700', not_comparable: 'text-amber-800' }
const STATE_ICON: Record<string, string> = { present: '●', absent: '○', not_comparable: '–' }
const STATE_RANK = { present: 0, absent: 1, not_comparable: 2 }
const NO_POS: number[] = []
/** rows (top-level interactions) per table page; 0 = all */
const PAGE_SIZES = [15, 50, 0]
interface Paging { pageSize: number; setPageSize: (n: number) => void }
const pct = (k: number, n: number) => (n ? `${Math.round((100 * k) / n)}%` : '–')
const instLabel = (r: { pdb_id: string; auth_asym_id: string }) => `${r.pdb_id.toUpperCase()} ${r.auth_asym_id}`
const classText = (c: Record<string, number>) => Object.entries(c).map(([k, v]) => `${k} ${v}`).join(' · ')
/** structural region(s) of an interaction; '–' when the family has no validated region annotation at either position */
const regionText = (p: Pick<InteractionPair, 'region_1' | 'region_2'>) => (!p.region_1 && !p.region_2 ? '–' : p.region_1 && p.region_1 === p.region_2 ? p.region_1 : `${p.region_1 ?? '–'} × ${p.region_2 ?? '–'}`)
/** directional categories read base → contacted backbone nucleotide */
const interactionLabel = (p: { category: InteractionCategory; pos_1: number; pos_2: number }, N: Numbering | null) =>
  p.category === 'baseribose' || p.category === 'basephosphate' ? `${positionLabel(p.pos_1, N).text}→${positionLabel(p.pos_2, N).text}` : pairLabel(p.pos_1, p.pos_2, N)
const groupStatus = (ss: (InstanceStatus | null | undefined)[]): InstanceStatus['state'] | null =>
  ss.some((s) => s?.state === 'present') ? 'present' : ss.some((s) => s?.state === 'absent') ? 'absent' : ss.some((s) => s) ? 'not_comparable' : null

const TIP = {
  diagram: 'Interactions between corresponding nucleotides across this RNA family. Arc colour indicates the fraction of eligible instances where each interaction is observed.',
  frequency: 'Frequency is calculated from structures where both nucleotides are sufficiently resolved for interaction annotation.',
  klass: 'Interaction classifications follow BGSU/FR3D annotations. Codes describe the interaction geometry or contact type.\nBp: Leontis–Westhof family (cis/trans, edges W/H/S).\nBs: stacking face orientation (s35, s53, s33, s55).\nBr, Bph: base–ribose / base–phosphate class 0–9 (position of the base donor).\nn prefix: near interaction (not counted).',
  numbering: 'Residue labels use the selected reference structure. Equivalent nucleotides across RNA structures are identified using Rfam consensus positions. Changing the numbering reference changes labels only.',
  displayed: 'The structure shown in 3D. The “In displayed structure” column refers to it; it does not change the numbering or the family-wide frequencies.',
  multiple: 'More than one interaction annotation between the same two nucleotides. Each keeps its own class and frequency; “any” counts chains with at least one of them among chains where at least one can be evaluated.',
}

/** label_seq_id of the instance's residue at a consensus position when it has coordinates (what can be highlighted in 3D) */
const labelSeqAt = (idx: MappingIndex | null, pos: number) => { if (!idx) return null; const a = chainAt(idx, pos); return a.state === 'match_observed' && a.residue ? a.residue.polymer_seq_position : null }

/** Compare-two pair per family for this page session: kept across mode and tab switches; `manual` once the user picked a chain. */
const compareMemory = new Map<string, { a: string | null; b: string | null; auto: boolean }>()

export function InteractionsTab() {
  const { acc } = useFamily()
  const [ds, setDs] = useState<InteractionsDataset | null | 'loading'>('loading')
  useEffect(() => { let alive = true; setDs('loading'); loadInteractions(acc).then((d) => { if (alive) setDs(d) }); return () => { alive = false } }, [acc])
  if (ds === 'loading') return <p className="text-sm text-stone-600" data-testid="interactions-loading">Loading interactions…</p>
  if (!ds) return <p className="text-sm text-stone-600" data-testid="interactions-unavailable">No RNA–RNA interaction data have been prepared for this family.</p>
  return <InteractionsView ds={ds} />
}

function useNumberingOf(seqName: string | null): Numbering | null {
  const { acc, family, instances } = useFamily()
  const inst = seqName ? instances.find((i) => i.seq_name === seqName) ?? null : null
  const ms = useMapping(acc, inst)
  return useMemo(() => (inst && ms.state === 'ready' ? buildNumbering(inst.seq_name, inst.pdb_id, inst.auth_asym_id, ms.mapping, family.clen) : null), [inst, ms, family.clen])
}

function InteractionsView({ ds }: { ds: InteractionsDataset }) {
  const { acc, family, instances, orientation } = useFamily()
  const rowOf = useMemo(() => new Map(instances.map((i) => [i.seq_name, i])), [instances])
  const defaultRef = ds.default_numbering_reference.seq_name && rowOf.has(ds.default_numbering_reference.seq_name) ? ds.default_numbering_reference.seq_name : null
  const [mode, setMode] = useState<'family' | 'compare'>('family')
  const [refSeq, setRefSeq] = useState<string | null>(defaultRef)
  const [displayed, setDisplayed] = useState<string | null>(defaultRef ?? ds.instances[0]?.seq_name ?? null)
  const [cmp, setCmpState] = useState(() => compareMemory.get(acc) ?? { a: null, b: null, auto: false })
  const setCmp = useCallback((v: { a: string | null; b: string | null; auto: boolean }) => { compareMemory.set(acc, v); setCmpState(v) }, [acc])
  const [filter, setFilter] = useState<InteractionFilter>(DEFAULT_FILTER)
  const [selectedId, setSelectedId] = useState<string | null>(null)
  const [query, setQuery] = useState('')
  const [pageSize, setPageSize] = useState<number>(PAGE_SIZES[0])

  // default comparison: highest valid geometric discrepancy, only when no pair has been chosen in this session
  useEffect(() => {
    if (mode !== 'compare' || cmp.a || cmp.b) return
    let alive = true
    loadStructureClustering(acc).then((sc) => { if (!alive) return; const p = defaultComparePair(ds, sc, (sn) => rowOf.get(sn)?.resolution ?? null); setCmp(p ? { a: p.a, b: p.b, auto: true } : { a: displayed, b: null, auto: false }) })
    return () => { alive = false }
  }, [mode, cmp.a, cmp.b, acc, ds, displayed, setCmp, rowOf])

  const numbering = useNumberingOf(refSeq)
  const numberingReady = !refSeq || !!numbering
  const N = refSeq ? numbering : null
  const labelFor = useCallback((pos: number) => positionLabel(pos, N), [N])
  const pairById = useMemo(() => new Map(ds.pairs.map((p) => [p.id, p])), [ds])
  const canonical = useMemo(() => new Map(ds.pairs.map((p) => [p.id, pairIsCanonical(p)])), [ds])
  const circular = showCircular(ds.clen)
  const L = useMemo(() => {
    let first = ds.clen, last = 1
    for (const i of ds.instances) for (let k = 0; k < ds.clen; k++) if (i.eligibility[k] !== '-') { first = Math.min(first, k + 1); last = Math.max(last, k + 1) }
    return circleLayout(Math.min(first, last), Math.max(first, last), SIZE, RADIUS)
  }, [ds])
  const visiblePairs = useMemo(() => ds.pairs.filter((p) => passes(filter, p, canonical.get(p.id)!)), [ds, filter, canonical])
  const selected = selectedId ? pairById.get(selectedId) ?? null : null
  const select = useCallback((id: string | null) => setSelectedId(id), [])
  // ---- unified table search (filters the tables only; selection changes only on an explicit row choice)
  const terms = useMemo(() => [...new Set(ds.position_regions.filter((r): r is string => !!r)), ...CATEGORIES.map((c) => CATEGORY_TEXT[c])], [ds])
  const search = useMemo(() => parseSearch(query, N, ds.clen, terms), [query, N, ds.clen, terms])
  const refRow = refSeq ? rowOf.get(refSeq) : null

  const header = (
    <div className="flex flex-wrap items-center gap-x-4 gap-y-2 rounded border border-stone-200 bg-stone-50 px-3 py-2 text-xs" data-testid="interactions-controls">
      <div className="flex overflow-hidden rounded border border-stone-300" role="group" aria-label="Interactions view">
        {(['family', 'compare'] as const).map((m) => <button key={m} type="button" data-testid={`imode-${m}`} aria-pressed={mode === m} onClick={() => setMode(m)}
          className={`px-2.5 py-1 ${mode === m ? 'bg-blue-600 text-white' : 'bg-white text-stone-700 hover:bg-stone-100'}`}>{m === 'family' ? 'Family conservation' : 'Compare two'}</button>)}
      </div>
      <label className="flex items-center gap-1 text-stone-700" data-testid="inumbering-control" data-state={numberingReady ? 'ready' : 'loading'}>Residue numbering
        <select data-testid="inumbering-select" value={refSeq ?? ''} onChange={(e) => setRefSeq(e.target.value || null)} className="max-w-[15rem] rounded border border-stone-300 bg-white px-1 py-px"
          aria-describedby="inumbering-active">
          {defaultRef && <option value={defaultRef}>{instLabel(rowOf.get(defaultRef)!)} (default reference)</option>}
          <option value="">Rfam consensus positions</option>
          {instances.filter((i) => i.seq_name !== defaultRef).map((i) => <option key={i.seq_name} value={i.seq_name}>{instLabel(i)}</option>)}
        </select>
        <InfoTip text={TIP.numbering} label="About residue numbering" />
        <span id="inumbering-active" className="sr-only" data-testid="inumbering-active">{refRow ? `Labels: ${instLabel(refRow)} author numbering` : 'Labels: Rfam consensus positions'}</span>
        {!numberingReady && <span className="text-stone-500">loading…</span>}
      </label>
      <TypeDropdown ds={ds} filter={filter} setFilter={setFilter} />
      <MoreFilters filter={filter} setFilter={setFilter} regions={ds.position_regions} />
    </div>)

  const visibleIds = useMemo(() => new Set(visiblePairs.map((p) => p.id)), [visiblePairs])
  const renderSearch = (onEnter: () => void) => <SearchBox query={query} setQuery={setQuery} search={search} N={N} ds={ds} visibleIds={visibleIds} onEnter={onEnter} />
  const noTypes = !CATEGORIES.some((c) => filter[c])
  const paging = { pageSize, setPageSize }
  const common = { ds, L, circular, labelFor, N, filter, visiblePairs, pairById, canonical, selected, selectedId, select, search, renderSearch, noTypes, paging, rowOf, orientation }
  return (
    <div className="space-y-3" data-testid="interactions-tab" data-mode={mode} data-circular={circular ? 'true' : 'false'}>
      {header}
      {!circular && <p className="rounded border border-stone-200 bg-white px-3 py-2 text-xs text-stone-600" data-testid="icircular-off">
        Circular diagrams are currently available for RNA families with up to {CIRCULAR_MAX_CLEN} Rfam consensus positions. Use the interaction table and 3D viewer to explore this family.</p>}
      {mode === 'family'
        ? <FamilyMode {...common} displayed={displayed} setDisplayed={setDisplayed} defaultRef={defaultRef} />
        : <CompareMode {...common} a={cmp.a} b={cmp.b} auto={cmp.auto} setPair={(a, b) => setCmp({ a, b, auto: false })} onBack={() => setMode('family')} />}
      <p className="text-[11px] text-stone-500" data-testid="interactions-provenance">
        BGSU RNA 3D Hub FR3D annotations mapped through the validated Rfam {String(ds.provenance.rfam_release ?? '')} mappings ·
        {' '}{ds.summary.n_intra_chain_interactions as number} interactions in {ds.instances.length} chains of {ds.summary.n_entries as number} entries · Rfam {family.rfam_acc}.
      </p>
    </div>
  )
}

/** "Interaction types" multiselect (checkbox popover); a category without data is marked and cannot hide data that do not exist. */
function TypeDropdown({ ds, filter, setFilter }: { ds: InteractionsDataset; filter: InteractionFilter; setFilter: (f: InteractionFilter) => void }) {
  const [open, setOpen] = useState(false)
  const box = useRef<HTMLDivElement>(null)
  useEffect(() => {
    if (!open) return
    const onDoc = (e: MouseEvent) => { if (!box.current?.contains(e.target as Node)) setOpen(false) }
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') setOpen(false) }
    document.addEventListener('mousedown', onDoc); document.addEventListener('keydown', onKey)
    return () => { document.removeEventListener('mousedown', onDoc); document.removeEventListener('keydown', onKey) }
  }, [open])
  const on = CATEGORIES.filter((c) => filter[c])
  const selectable = CATEGORIES.filter((c) => categoryAvailability(ds, c) !== 'none')
  const summary = !on.length ? 'none selected' : on.length === CATEGORIES.length ? 'All types' : on.map((c) => CATEGORY_CODE[c]).join(', ')
  return (
    <div ref={box} className="relative" data-testid="itypes">
      <button type="button" aria-haspopup="true" aria-expanded={open} onClick={() => setOpen((o) => !o)} data-testid="itypes-toggle"
        className="rounded border border-stone-300 bg-white px-2 py-0.5 hover:bg-stone-100">Interaction types: <b data-testid="itypes-summary" className={on.length ? '' : 'text-amber-800'}>{summary}</b> ▾</button>
      {open && <div className="absolute left-0 z-30 mt-1 w-60 rounded border border-stone-300 bg-white p-2 shadow" role="group" aria-label="Interaction types">
        {CATEGORIES.map((c) => { const av = categoryAvailability(ds, c)
          return <label key={c} className={`flex items-center gap-1.5 py-0.5 ${av === 'none' ? 'text-stone-400' : ''}`} data-testid={`itype-${c}-label`} data-availability={av}>
            <input type="checkbox" checked={filter[c]} disabled={av === 'none'} onChange={(e) => setFilter({ ...filter, [c]: e.target.checked })} data-testid={`ifilter-${c}`} />
            <LineSample cat={c} where="d" />{CATEGORY_TEXT[c]} ({CATEGORY_CODE[c]}){av === 'none' && <span className="ml-auto text-[10px]">no data</span>}{av === 'empty' && <span className="ml-auto text-[10px] text-stone-500">none observed</span>}</label> })}
        <div className="mt-1 flex justify-between border-t border-stone-200 pt-1">
          <button type="button" className="text-blue-700 hover:underline" onClick={() => setFilter({ ...filter, ...Object.fromEntries(CATEGORIES.map((c) => [c, selectable.includes(c)])) })} data-testid="itypes-all">Select all</button>
          <button type="button" className="text-blue-700 hover:underline" onClick={() => setFilter({ ...filter, basepair: true, stacking: false, baseribose: false, basephosphate: false })} data-testid="itypes-reset">Reset to Bp</button>
        </div>
      </div>}
    </div>)
}

function MoreFilters({ filter, setFilter, regions }: { filter: InteractionFilter; setFilter: (f: InteractionFilter) => void; regions: (string | null)[] }) {
  const names = useMemo(() => [...new Set(regions.filter((r): r is string => !!r))], [regions])
  const active = (filter.canonical !== 'all' ? 1 : 0) + (filter.structure !== 'all' ? 1 : 0) + (filter.region !== null ? 1 : 0)
  const [open, setOpen] = useState(false)
  const sel = 'rounded border border-stone-300 bg-white px-1 py-px'
  return (
    <div className="relative" data-testid="ifilter-more">
      <button type="button" aria-expanded={open} onClick={() => setOpen((o) => !o)} className="text-stone-600 hover:text-stone-900" data-testid="ifilter-more-toggle">{open ? '▾' : '▸'} More filters{active ? ` (${active})` : ''}</button>
      {open && <div className="absolute left-0 z-20 mt-1 flex w-64 flex-col gap-1.5 rounded border border-stone-300 bg-white p-2 shadow">
        <label className="flex items-center justify-between gap-2">Base-pair type
          <select className={sel} value={filter.canonical} onChange={(e) => setFilter({ ...filter, canonical: e.target.value as InteractionFilter['canonical'] })} data-testid="ifilter-canonical">
            <option value="all">All</option><option value="canonical">Canonical (cWW AU/GC/GU)</option><option value="noncanonical">Noncanonical</option></select></label>
        <label className="flex items-center justify-between gap-2">Context
          <select className={sel} value={filter.structure} onChange={(e) => setFilter({ ...filter, structure: e.target.value as InteractionFilter['structure'] })} data-testid="ifilter-structure">
            <option value="all">All</option><option value="secondary">Secondary (Rfam SS_cons)</option><option value="tertiary">Tertiary (between elements)</option></select></label>
        <p className="text-[10px] text-stone-500">Type and context apply to base pairs only; other categories are hidden while they are set.</p>
        <label className="flex items-center justify-between gap-2">Region
          <select className={sel} value={filter.region ?? ''} onChange={(e) => setFilter({ ...filter, region: e.target.value || null })} data-testid="ifilter-region">
            <option value="">All</option>{names.map((n) => <option key={n} value={n}>{n}</option>)}</select></label>
        <div className="flex justify-between"><button type="button" className="text-stone-500 hover:text-stone-800" onClick={() => setFilter({ ...filter, canonical: 'all', structure: 'all', region: null })} data-testid="ifilter-reset">Reset</button>
          <button type="button" className="text-stone-500 hover:text-stone-800" onClick={() => setOpen(false)}>Close</button></div>
      </div>}
    </div>)
}

interface CommonProps {
  ds: InteractionsDataset; L: ReturnType<typeof circleLayout>; circular: boolean; labelFor: (pos: number) => ReturnType<typeof positionLabel>; N: Numbering | null
  filter: InteractionFilter; visiblePairs: InteractionPair[]; pairById: Map<string, InteractionPair>; canonical: Map<string, boolean>
  selected: InteractionPair | null; selectedId: string | null; select: (id: string | null) => void; search: SearchQuery; renderSearch: (onEnter: () => void) => ReactNode; noTypes: boolean; paging: Paging
  rowOf: Map<string, InstanceRow>; orientation: ReturnType<typeof useFamily>['orientation']
}

function StatusBadge({ s, testId, compact = false }: { s: InstanceStatus; testId?: string; compact?: boolean }) {
  const title = s.state === 'not_comparable' ? s.reason ?? '' : s.state === 'absent' ? `Resolved and eligible; not annotated.${s.other.length ? ` Other: ${s.other.join(', ')}` : ''}` : `Annotated: ${s.classes.join(', ')}`
  return <span className={`whitespace-nowrap ${STATE_STYLE[s.state]}`} title={title} data-testid={testId} data-state={s.state}><span aria-hidden className="mr-1">{STATE_ICON[s.state]}</span>{STATE_TEXT[s.state]}{s.state === 'present' && !compact && <span className="ml-1 font-mono text-[10px] text-stone-500">{s.classes.join(',')}</span>}</span>
}
const TypeCell = ({ cats }: { cats: InteractionCategory[] }) => cats.length === 1
  ? <span title={CATEGORY_TEXT[cats[0]]}>{CATEGORY_CODE[cats[0]]}</span>
  : <span className="rounded bg-violet-50 px-1 text-violet-900" title={`Multiple interactions: ${cats.map((c) => CATEGORY_TEXT[c]).join(', ')}`} data-testid="itype-multiple">Multiple <span className="whitespace-nowrap">({[...new Set(cats)].map((c) => CATEGORY_CODE[c]).join('+')})</span></span>

// ===================================================================================================== family conservation
interface FamilyGroup { group: InteractionGroup<InteractionPair>; id: string; members: InteractionPair[]; frequency: number | null; k: number; n: number; entries: number; state: InstanceStatus['state'] | null }

function FamilyMode({ ds, L, circular, labelFor, N, visiblePairs, selected, selectedId, select, search, renderSearch, noTypes, paging, rowOf, orientation, displayed, setDisplayed, defaultRef }: CommonProps & {
  displayed: string | null; setDisplayed: (sn: string) => void; defaultRef: string | null }) {
  const { acc } = useFamily()
  const dIx = displayed ? indexInstance(ds, displayed) : null
  const dRow = displayed ? rowOf.get(displayed) ?? null : null
  const dMs = useMapping(acc, dRow)
  const dIdx = useMemo(() => (dMs.state === 'ready' ? indexMapping(dMs.mapping) : null), [dMs])
  const statusOf = useMemo(() => new Map(ds.pairs.map((p) => [p.id, dIx ? statusIn(dIx, p.category, p.pos_1, p.pos_2) : null])), [ds, dIx])
  const rows = visiblePairs
  const groups = useMemo<FamilyGroup[]>(() => groupByEndpoints(rows).map((g) => {
    const members = [...g.members].sort((a, b) => defaultFamilyOrder({ entries: a.n_entries_observed, frequency: a.frequency, n: a.n_eligible, i: a.pos_1, j: a.pos_2 }, { entries: b.n_entries_observed, frequency: b.frequency, n: b.n_eligible, i: b.pos_1, j: b.pos_2 }) || CATEGORIES.indexOf(a.category) - CATEGORIES.indexOf(b.category))
    const occ = members.length === 1 ? { k: members[0].n_observed, n: members[0].n_eligible, frequency: members[0].frequency } : groupOccurrence(ds, members)
    return { group: g, id: groupId(g), members, frequency: occ.frequency, k: occ.k, n: occ.n, entries: Math.max(...members.map((m) => m.n_entries_observed)), state: null }
  }), [rows, ds])
  const groupsWithStatus = useMemo(() => groups.map((g) => ({ ...g, state: groupStatus(g.members.map((m) => statusOf.get(m.id))) })), [groups, statusOf])
  const selGroup = selectedId ? groupsWithStatus.find((g) => g.members.some((m) => m.id === selectedId)) ?? null : null
  const selectGroup = useCallback((id: string) => { const g = groups.find((x) => x.id === id); select(g ? g.members[0].id : id) }, [groups, select])
  const arcs = useMemo<DiagramArc[]>(() => groups.map((g) => ({ id: g.id, i: g.members.length === 1 ? g.members[0].pos_1 : g.group.i, j: g.members.length === 1 ? g.members[0].pos_2 : g.group.j,
    ...(g.members.length === 1 ? { color: categoryColor(g.members[0].category, g.frequency), dash: CATEGORY_STYLE[g.members[0].category].dash, category: g.members[0].category } : { color: multiColor(g.frequency), multi: true, category: 'multiple' }),
    width: ARC_WIDTH, order: (g.frequency ?? 0) + (g.members.some((m) => m.category === 'basepair') ? 0.001 : 0) })), [groups])
  const tooltip = useCallback((id: string) => { const g = groups.find((x) => x.id === id); if (!g) return null; const s = (m: InteractionPair) => statusOf.get(m.id)
    return <div data-testid="iarc-tooltip" className="space-y-0.5">
      <p className="font-semibold">{pairLabel(g.group.i, g.group.j, N)}{g.members.length > 1 ? ' · Multiple interactions' : ''}</p>
      <p className="text-stone-600">Rfam {g.group.i}–{g.group.j}</p>
      {g.members.map((m) => <p key={m.id}><span className="font-mono">{CATEGORY_CODE[m.category]} {dominantClass(m)}</span> · {m.n_observed}/{m.n_eligible} ({pct(m.n_observed, m.n_eligible)}){dRow && s(m) ? ` · ${STATE_TEXT[s(m)!.state]} in ${instLabel(dRow)}` : ''}</p>)}
      {g.members.length > 1 && <p className="text-stone-600">Any: {g.k}/{g.n} ({pct(g.k, g.n)})</p>}</div> }, [groups, statusOf, N, dRow])
  const focus = useMemo(() => (selected ? [selected.pos_1, selected.pos_2] : NO_POS), [selected])
  const viewerSel = focus.map((p) => labelSeqAt(dIdx, p)).filter((x): x is number => x != null)
  const selStatus = selected ? statusOf.get(selected.id) ?? null : null
  const details = <SelectedDetails ds={ds} p={selected} group={selGroup} statusOf={statusOf} N={N} rowOf={rowOf} displayed={displayed} defaultRef={defaultRef} setDisplayed={setDisplayed} select={select} />
  return (
    <>
      {circular ? <div className="grid gap-3 lg:grid-cols-[minmax(0,1fr)_minmax(0,20rem)]">
        <section className="rounded border border-stone-200 bg-white p-2" data-testid="ifamily-diagram-card">
          <h3 className="px-1 text-sm font-semibold text-stone-800">Interaction conservation<InfoTip text={TIP.diagram} label="About the conservation diagram" /></h3>
          {noTypes ? <NoTypes /> : <div className="flex justify-center"><CircularDiagram L={L} size={SIZE} arcs={arcs} labelFor={labelFor} regions={ds.position_regions} selectedId={selGroup?.id ?? null} focus={focus} onSelect={selectGroup} tooltip={tooltip} testId="idiagram" title="Family-wide interaction conservation diagram" /></div>}
        </section>
        <div className="flex flex-col gap-3"><FamilyLegend />{details}</div>
      </div> : details}
      <div className="grid gap-x-3 gap-y-4 lg:grid-cols-2 lg:gap-y-1.5" data-testid="ipanel-grid">
        <FamilyTable ds={ds} groups={groupsWithStatus} statusOf={statusOf} N={N} selectedId={selectedId} select={select} displayedLabel={dRow ? instLabel(dRow) : ''} search={search} renderSearch={renderSearch} noTypes={noTypes} paging={paging} />
        <section className="grid min-w-0 content-start gap-y-1.5 lg:row-span-3 lg:grid-rows-subgrid" data-testid="ifamily-viewer">
          <div className="flex flex-wrap items-center gap-2 self-end text-xs">
            <label className="flex items-center gap-1 font-medium text-stone-800">Displayed structure
              <select value={displayed ?? ''} onChange={(e) => setDisplayed(e.target.value)} className="max-w-[15rem] rounded border border-stone-300 bg-white px-1 py-px font-normal" data-testid="idisplayed-select">
                {[...rowOf.values()].map((i) => <option key={i.seq_name} value={i.seq_name}>{instLabel(i)}{i.resolution != null ? ` · ${i.resolution} Å` : ''}{i.seq_name === defaultRef ? ' · default' : ''}</option>)}
              </select></label>
            <InfoTip text={TIP.displayed} label="About the displayed structure" />
          </div>
          {dRow && <PdbeMolstarViewer pdbId={dRow.pdb_id} labelAsymId={dRow.label_asym_id} authAsymId={dRow.auth_asym_id} region={{ start: dRow.map_start, end: dRow.map_end }}
            orientation={orientation?.transforms[dRow.seq_name] ?? null} selection={viewerSel.map((label_seq_id) => ({ label_seq_id }))} focusOnSelection structureColor={0x64748b}
            onClearSelection={() => select(null)} canClear={!!selected} className="h-[460px]" />}
          <ViewerCaption p={selected} status={selStatus} N={N} dRow={dRow} highlighted={viewerSel} />
        </section>
      </div>
    </>
  )
}

/** line sample of a category: its pattern, drawn with its light-to-dark conservation ramp */
function LineSample({ cat, w = 30, where = 'l' }: { cat: InteractionCategory; w?: number; where?: string }) {
  const st = CATEGORY_STYLE[cat], id = `ils-${where}-${cat}`
  return <svg width={w} height="8" aria-hidden className="shrink-0" data-testid={`${where === 'l' ? 'ilegend' : 'itype-sample'}-${cat}`} data-dash={st.dash ?? ''}>
    <defs><linearGradient id={id} gradientUnits="userSpaceOnUse" x1="2" y1="0" x2={w - 2} y2="0">{st.ramp.map((c, k) => <stop key={c} offset={k / (st.ramp.length - 1)} stopColor={c} />)}</linearGradient></defs>
    <line x1="2" y1="4" x2={w - 2} y2="4" stroke={`url(#${id})`} strokeWidth={2.2} strokeDasharray={st.dash} strokeLinecap="round" /></svg>
}
/** small rare -> conserved bar of one colour ramp */
const RampBar = ({ ramp, testId }: { ramp: string[]; testId: string }) =>
  <span className="inline-block h-2 w-14 shrink-0 rounded-sm" style={{ background: `linear-gradient(to right, ${ramp.join(', ')})` }} data-testid={testId} data-ramp={ramp.join(',')} aria-hidden />
function FamilyLegend() {
  const row = 'grid grid-cols-[30px_2rem_1fr_3.5rem] items-center gap-x-1.5'
  return (
    <section className="rounded border border-stone-200 bg-white px-3 py-2 text-[11px] text-stone-700" data-testid="ilegend">
      <div className={`${row} text-[10px] text-stone-500`}><span /><span /><span /><span className="flex justify-between" data-testid="ilegend-scale-head"><span>0%</span><span>100%</span></span></div>
      <ul className="flex flex-col gap-0.5">
        {CATEGORIES.map((c) => <li key={c} className={row}><LineSample cat={c} /><span className="font-medium">{CATEGORY_CODE[c]}</span><span className="text-stone-500">{CATEGORY_TEXT[c]}</span><RampBar ramp={CATEGORY_STYLE[c].ramp} testId={`ilegend-ramp-${c}`} /></li>)}
        <li className={row}><svg width="30" height="8" aria-hidden className="shrink-0" data-testid="ilegend-multiple"><line x1="2" y1="4" x2="28" y2="4" stroke={MULTI_RAMP[2]} strokeWidth={3.4} strokeLinecap="round" /><line x1="2" y1="4" x2="28" y2="4" stroke="white" strokeWidth={0.9} /></svg>
          <span /><span className="text-stone-500" title="More than one interaction type between the same two nucleotides">Multiple interactions</span><RampBar ramp={MULTI_RAMP} testId="ilegend-ramp-multiple" /></li>
        <li className={row}><svg width="30" height="10" aria-hidden><circle cx="15" cy="5" r="3.5" fill="#f59e0b" stroke="#1c1917" strokeWidth="0.8" /></svg><span /><span className="text-stone-500">Selected nucleotide</span><span /></li>
      </ul>
      <p className="mt-1.5 text-stone-500" data-testid="ilegend-scale">Each type has its own scale: colour intensity = observed / eligible chain instances (lighter = rarer, darker = more conserved). Multiple interactions use the “any” occurrence.</p>
    </section>)
}
function NoTypes() {
  return <p className="px-2 py-6 text-center text-xs text-stone-600" data-testid="ino-types">No interaction types selected. Choose at least one in <b>Interaction types</b>.</p>
}

/** The search field above the tables: residues (selected numbering or Rfam), residue pairs, class codes or words. Typing filters the table;
 *  Enter selects the first result (an explicit choice); nothing else changes the selection. */
function SearchBox({ query, setQuery, search, N, ds, visibleIds, onEnter }: { query: string; setQuery: (q: string) => void; search: SearchQuery; N: Numbering | null; ds: InteractionsDataset
  visibleIds: Set<string>; onEnter: () => void }) {
  const hidden = search.kind === 'residue' ? ds.pairs.filter((p) => matchesPositions(p, search.positions) && !visibleIds.has(p.id)).length : 0
  const pos = search.kind === 'residue' ? search.positions : []
  return (
    <div className="ml-auto flex min-w-[14rem] flex-col items-end gap-0.5" data-testid="isearch-box">
      <input value={query} onChange={(e) => setQuery(e.target.value)} onKeyDown={(e) => { if (e.key === 'Enter') { e.preventDefault(); onEnter() } if (e.key === 'Escape') setQuery('') }}
        placeholder="Search interactions or residues" aria-label="Search interactions or residues" title="G37, 37, G37-C61, Rfam 41, Rfam 41-66, cWW, s35, 9BR, J2/3 … Enter selects the first result"
        className="w-56 max-w-full rounded border border-stone-300 px-1.5 py-0.5" data-testid="isearch" />
      {search.kind === 'error' && <p className="text-[11px] text-amber-800" data-testid="isearch-error">{search.message}</p>}
      {search.kind === 'residue' && <p className="whitespace-nowrap text-[11px] text-stone-600" data-testid="isearch-match">
        {pos.length === 2 ? <>Between {positionLabel(pos[0], N).text} and {positionLabel(pos[1], N).text} (Rfam {pos[0]}–{pos[1]})</> : <>Involving {positionLabel(pos[0], N).text} (Rfam {pos[0]})</>}
        {hidden > 0 && <span data-testid="isearch-hidden"> · {hidden} more in types not shown</span>}</p>}
    </div>)
}

function SelectedDetails({ ds, p, group, statusOf, N, rowOf, displayed, defaultRef, setDisplayed, select }: { ds: InteractionsDataset; p: InteractionPair | null; group: FamilyGroup | null
  statusOf: Map<string, InstanceStatus | null>; N: Numbering | null; rowOf: Map<string, InstanceRow>; displayed: string | null; defaultRef: string | null; setDisplayed: (sn: string) => void; select: (id: string | null) => void }) {
  const support = useMemo(() => (p ? supportingInstances(ds, p, (sn) => rowOf.get(sn)?.resolution ?? null, defaultRef ? [defaultRef] : []) : []), [ds, p, rowOf, defaultRef])
  if (!p) return <p className="self-start rounded border border-dashed border-stone-300 bg-white px-3 py-2 text-xs text-stone-500" data-testid="idetails-empty">Select an interaction from the diagram or table to inspect its conservation and 3D geometry.</p>
  const status = statusOf.get(p.id) ?? null
  const dRow = displayed ? rowOf.get(displayed) : undefined
  const nc = p.n_missing + p.n_incomplete + p.n_unmapped + p.n_annotation_unavailable
  return (
    <section className="space-y-2 self-start rounded border border-stone-200 bg-white p-3 text-xs text-stone-800" data-testid="idetails" data-id={p.id}>
      <div className="flex items-start justify-between gap-2"><h3 className="text-sm font-semibold" data-testid="idetails-title">{interactionLabel(p, N)} · {CATEGORY_CODE[p.category]} {status?.state === 'present' ? status.classes.join(', ') : dominantClass(p)}</h3>
        <button type="button" className="text-stone-500 hover:text-stone-800" onClick={() => select(null)} aria-label="Clear selection" data-testid="idetails-clear">×</button></div>
      {group && group.members.length > 1 && <div data-testid="idetails-multiple"><p className="font-medium">Multiple interactions <InfoTip text={TIP.multiple} label="About multiple interactions" /> <span className="font-normal text-stone-500">any: {group.k}/{group.n} ({pct(group.k, group.n)})</span></p>
        <ul className="mt-0.5 space-y-0.5">{group.members.map((m) => { const s = statusOf.get(m.id)
          return <li key={m.id}><button type="button" onClick={() => select(m.id)} aria-pressed={m.id === p.id} data-testid={`idetails-member-${m.id}`}
            className={`w-full rounded px-1 text-left ${m.id === p.id ? 'bg-amber-100' : 'hover:bg-stone-100'}`}><span className="font-mono">{CATEGORY_CODE[m.category]} {Object.keys(m.classes).join('/')}</span> · {m.n_observed}/{m.n_eligible} ({pct(m.n_observed, m.n_eligible)}){s ? ` · ${STATE_TEXT[s.state]}` : ''}</button></li> })}</ul></div>}
      <dl className="grid grid-cols-[7rem_1fr] gap-x-2 gap-y-0.5">
        <dt className="text-stone-500">Rfam positions</dt><dd data-testid="idetails-rfam">{p.pos_1}{p.directional ? '→' : '–'}{p.pos_2}{p.directional ? ' (base → backbone)' : ''}</dd>
        <dt className="text-stone-500">Type</dt><dd>{CATEGORY_TEXT[p.category]} ({CATEGORY_CODE[p.category]}){p.structure_class && p.structure_class !== 'unclassified' ? ` · ${p.structure_class}` : ''}</dd>
        <dt className="text-stone-500">Interaction class</dt><dd className="font-mono" data-testid="idetails-classes">{classText(p.classes)}</dd>
        <dt className="text-stone-500">Conservation</dt><dd data-testid="idetails-conservation">{p.n_observed}/{p.n_eligible} eligible chains ({pct(p.n_observed, p.n_eligible)})</dd>
        <dt className="text-stone-500">PDB entries</dt><dd>{p.n_entries_observed} of {p.n_entries_eligible} eligible</dd>
        <dt className="text-stone-500">Region</dt><dd>{regionText(p)}</dd>
        <dt className="text-stone-500">Not evaluated</dt><dd title="Chains where this interaction cannot be evaluated; never counted as absent">{nc} chains{p.n_near_only ? ` · ${p.n_near_only} near only` : ''}</dd>
      </dl>
      {status && dRow && <div className={`rounded border px-2 py-1.5 ${status.state === 'present' ? 'border-emerald-200 bg-emerald-50' : status.state === 'absent' ? 'border-stone-300 bg-stone-50' : 'border-amber-200 bg-amber-50'}`} data-testid="idetails-status" data-state={status.state}>
        <p><b>{STATE_TEXT[status.state]}</b> in displayed structure ({instLabel(dRow)}){status.state === 'present' ? `: ${status.classes.join(', ')}` : ''}</p>
        {status.state === 'absent' && <p>Resolved and eligible; not annotated.{status.other.length > 0 && <> Other: <span className="font-mono" data-testid="idetails-other">{status.other.join(', ')}</span></>}</p>}
        {status.state === 'not_comparable' && <p data-testid="idetails-reason">{status.reason}</p>}
      </div>}
      {status && status.state !== 'present' && <label className="flex flex-col gap-0.5" data-testid="isupport">
        <span className="font-medium">Inspect in another structure <span className="font-normal text-stone-500">({support.length} with this interaction)</span></span>
        <select className="rounded border border-stone-300 bg-white px-1 py-0.5" value="" onChange={(e) => e.target.value && setDisplayed(e.target.value)} data-testid="isupport-select">
          <option value="">Choose a chain…</option>
          {support.map((s) => { const r = rowOf.get(s.seqName); return <option key={s.seqName} value={s.seqName}>{r ? instLabel(r) : s.seqName}{r?.resolution != null ? ` · ${r.resolution} Å` : ''} · {s.classes.join(',')}{s.seqName === defaultRef ? ' · reference' : ''}</option> })}
        </select></label>}
    </section>)
}

function ViewerCaption({ p, status, N, dRow, highlighted }: { p: InteractionPair | null; status: InstanceStatus | null; N: Numbering | null; dRow: InstanceRow | null; highlighted: number[] }) {
  if (!dRow) return null
  const who = instLabel(dRow)
  let text: ReactNode = <>Select an interaction to highlight its nucleotides in {who}.</>
  if (p && status) {
    const lab = interactionLabel(p, N)
    text = status.state === 'present' ? <>{lab}: present in {who} ({status.classes.join(', ')}).</>
      : status.state === 'absent' ? <>{lab}: <b>not annotated</b> in {who}. Both nucleotides are highlighted; no interaction is implied.</>
        : <>{lab}: not comparable in {who} ({status.reason}).{highlighted.length ? ` ${highlighted.length === 1 ? 'The modelled nucleotide is' : 'The modelled nucleotides are'} highlighted.` : ''}</>
  }
  return <p className="text-xs text-stone-600" data-testid="iviewer-caption" data-highlighted={highlighted.join(',')}>{text}</p>
}

type SortKey = 'default' | 'pos' | 'type' | 'class' | 'cons' | 'entries' | 'region' | 'status'
function useSort<T>(rows: T[], key: (r: T, k: SortKey) => string | number, byDefault: (a: T, b: T) => number) {
  const [sort, setSort] = useState<{ k: SortKey; dir: 1 | -1 }>({ k: 'default', dir: 1 })
  const sorted = useMemo(() => [...rows].sort((a, b) => { if (sort.k === 'default') return byDefault(a, b); const x = key(a, sort.k), y = key(b, sort.k); return (x < y ? -1 : x > y ? 1 : 0) * sort.dir || byDefault(a, b) }), [rows, sort, key, byDefault])
  const th = (k: SortKey, label: string, cls = '', extra: ReactNode = null) => <th className={`px-1.5 py-1 text-left font-medium ${cls}`} aria-sort={sort.k === k ? (sort.dir === 1 ? 'ascending' : 'descending') : 'none'}>
    <span className="inline-flex items-center"><button type="button" className="select-none text-left" onClick={() => setSort((s) => ({ k, dir: s.k === k ? (s.dir === 1 ? -1 : 1) : (k === 'cons' || k === 'entries' ? -1 : 1) }))} data-testid={`isort-${k}`}>{label}{sort.k === k ? (sort.dir === 1 ? '\u00a0▲' : '\u00a0▼') : ''}</button>{extra}</span></th>
  const reset = sort.k !== 'default' ? <button type="button" className="text-stone-500 underline hover:text-stone-800" onClick={() => setSort({ k: 'default', dir: 1 })} data-testid="isort-reset">Default order</button> : null
  return { sorted, th, reset, sortKey: `${sort.k}:${sort.dir}` }
}
function useScrollToSelected(selectedId: string | null) {
  const box = useRef<HTMLDivElement>(null)
  useEffect(() => { if (!selectedId) return; const el = box.current?.querySelector(`[data-row-id="${CSS.escape(selectedId)}"]`) as HTMLElement | null; el?.scrollIntoView?.({ block: 'nearest' }) }, [selectedId])
  return box
}
/** expanded multiple-interaction rows; the group of the selected interaction is always open */
/** Table pages of top-level rows (an expanded group's members stay with their group). A new sort, search, filter or page size starts at
 *  page 1 (`resetKey`), so the new order is read from the top; a selection made elsewhere (diagram arc, Enter in the search) shows its page;
 *  paging never changes the selection. */
function usePaged<T extends { members: { id: string }[] }>(rows: T[], selectedId: string | null, { pageSize, setPageSize }: Paging, resetKey: string) {
  const [page, setPage] = useState(1)
  const pages = pageSize ? Math.max(1, Math.ceil(rows.length / pageSize)) : 1
  const idx = selectedId ? rows.findIndex((g) => g.members.some((m) => m.id === selectedId)) : -1
  const pageOf = (i: number) => (pageSize && i >= 0 ? Math.floor(i / pageSize) + 1 : 1)
  const first = useRef(true)
  useEffect(() => { if (first.current) { first.current = false; setPage(pageOf(idx)); return } setPage(1) }, [resetKey, pageSize])   // eslint-disable-line react-hooks/exhaustive-deps
  // a selection made outside the visible page (diagram, search Enter) shows its page; a click on a visible row keeps the page
  useEffect(() => { if (idx >= 0 && pageOf(idx) !== page) setPage(pageOf(idx)) }, [selectedId])   // eslint-disable-line react-hooks/exhaustive-deps
  const cur = Math.min(page, pages)
  const visible = pageSize ? rows.slice((cur - 1) * pageSize, cur * pageSize) : rows
  const btn = 'rounded border border-stone-300 bg-white px-2 py-0.5 disabled:opacity-40'
  const pager = (
    <div className="flex flex-wrap items-center gap-2 text-xs text-stone-700" data-testid="ipager">
      <button type="button" className={btn} disabled={cur <= 1} onClick={() => setPage(1)} aria-label="First page" title="First page" data-testid="ipage-first">« First</button>
      <button type="button" className={btn} disabled={cur <= 1} onClick={() => setPage(Math.max(1, cur - 1))} data-testid="ipage-prev">Previous</button>
      <span className="flex items-center gap-1 whitespace-nowrap" data-testid="ipage-info" aria-label={`Page ${cur} of ${pages}`}>Page
        <input key={cur} type="number" min={1} max={pages} defaultValue={cur} aria-label="Go to page" data-testid="ipage-input"
          className="w-12 rounded border border-stone-300 px-1 py-px text-center [appearance:textfield] [&::-webkit-inner-spin-button]:appearance-none"
          onKeyDown={(e) => { if (e.key === 'Enter') { const v = Math.round(+e.currentTarget.value); if (Number.isFinite(v)) setPage(Math.min(pages, Math.max(1, v))) } }}
          onBlur={(e) => { const v = Math.round(+e.currentTarget.value); if (Number.isFinite(v) && v !== cur) setPage(Math.min(pages, Math.max(1, v))) }} />
        of <span data-testid="ipage-count">{pages}</span></span>
      <button type="button" className={btn} disabled={cur >= pages} onClick={() => setPage(Math.min(pages, cur + 1))} data-testid="ipage-next">Next</button>
      <button type="button" className={btn} disabled={cur >= pages} onClick={() => setPage(pages)} aria-label="Last page" title="Last page" data-testid="ipage-last">Last »</button>
      <label className="ml-auto flex items-center gap-1 text-stone-600">Rows per page
        <select value={pageSize} onChange={(e) => setPageSize(+e.target.value)} className="rounded border border-stone-300 bg-white px-1 py-px" data-testid="ipage-size">
          {PAGE_SIZES.map((n) => <option key={n} value={n}>{n || 'All'}</option>)}</select></label>
    </div>)
  /** a full page of the default size fills the box (rows stretch); a short last page keeps natural row heights */
  const fill = pageSize > 0 && visible.length === pageSize
  return { visible, pager, fill }
}
function useExpanded(selectedId: string | null, groups: { id: string; members: { id: string }[] }[]) {
  const [open, setOpen] = useState<Set<string>>(new Set())
  const selGroup = selectedId ? groups.find((g) => g.members.length > 1 && g.members.some((m) => m.id === selectedId))?.id : undefined
  const isOpen = (id: string) => open.has(id) || id === selGroup
  const toggle = (id: string) => setOpen((s) => { const n = new Set(s); if (n.has(id)) n.delete(id); else n.add(id); return n })
  return { isOpen, toggle }
}

const famDefault = (a: FamilyGroup, b: FamilyGroup) => defaultFamilyOrder({ entries: a.entries, frequency: a.frequency, n: a.n, i: a.group.i, j: a.group.j }, { entries: b.entries, frequency: b.frequency, n: b.n, i: b.group.i, j: b.group.j })
function FamilyTable({ ds, groups, statusOf, N, selectedId, select, displayedLabel, search, renderSearch, noTypes, paging }: { groups: FamilyGroup[]; statusOf: Map<string, InstanceStatus | null>; N: Numbering | null
  selectedId: string | null; select: (id: string) => void; displayedLabel: string; search: SearchQuery; renderSearch: (onEnter: () => void) => ReactNode; noTypes: boolean; paging: Paging
  ds: InteractionsDataset }) {
  const pdbOf = useMemo(() => { const m = new Map(ds.instances.map((i) => [i.seq_name, i.pdb_id])); return (sn: string) => m.get(sn) ?? sn }, [ds])
  // interactions annotated with different FR3D classes in different chains: per-class breakdown rows on demand
  const [classOpen, setClassOpen] = useState<Set<string>>(new Set())
  const toggleClasses = (id: string) => setClassOpen((o) => { const n = new Set(o); if (n.has(id)) n.delete(id); else n.add(id); return n })
  const text = useCallback((g: FamilyGroup) => g.members.map((p) => `${interactionLabel(p, N)} rfam ${p.pos_1}-${p.pos_2} ${CATEGORY_TEXT[p.category]} ${CATEGORY_CODE[p.category]} ${Object.keys(p.classes).join(' ')} ${regionText(p)}`).join(' ').toLowerCase(), [N])
  const shown = useMemo(() => (search.kind === 'none' ? groups : search.kind === 'error' ? [] : search.kind === 'residue' ? groups.filter((g) => matchesPositions({ pos_1: g.group.i, pos_2: g.group.j }, search.positions))
    : groups.filter((g) => text(g).includes(search.text))), [groups, search, text])
  const key = useCallback((g: FamilyGroup, k: SortKey): string | number => { const p = g.members[0]
    return k === 'pos' ? g.group.i * 10000 + g.group.j : k === 'type' ? (g.members.length > 1 ? 'zz' : CATEGORY_CODE[p.category]) : k === 'class' ? dominantClass(p) : k === 'cons' ? (g.frequency ?? -1) * 1e6 + g.n
      : k === 'entries' ? g.entries : k === 'region' ? regionText(p) : STATE_RANK[g.state ?? 'not_comparable'] }, [])
  const { sorted, th, reset, sortKey } = useSort(shown, key, famDefault)
  const { visible, pager, fill } = usePaged(sorted, selectedId, paging, `${sortKey}|${JSON.stringify(search)}|${groups.length}`)
  const box = useScrollToSelected(selectedId)
  const { isOpen, toggle } = useExpanded(selectedId, groups)
  const memberRow = (p: InteractionPair, sub: boolean) => { const s = statusOf.get(p.id), on = p.id === selectedId
    return <tr key={p.id} data-row-id={p.id} data-testid={`irow-${p.id}`} data-status={s?.state ?? ''} data-frequency={p.frequency ?? ''} aria-selected={on} onClick={() => select(p.id)}
      className={`cursor-pointer border-t border-stone-100 ${on ? 'bg-amber-100' : sub ? 'bg-stone-50/60 hover:bg-stone-100' : 'hover:bg-stone-50'}`}>
      <td className={`py-1 pr-1.5 whitespace-nowrap ${sub ? 'pl-5' : 'pl-1.5'}`} title={`Rfam ${p.pos_1}${p.directional ? '→' : '–'}${p.pos_2} · ${regionText(p)}`}><span className="font-medium">{interactionLabel(p, N)}</span>{N && <span className="ml-1 text-[10px] text-stone-400">{p.pos_1}{p.directional ? '→' : '–'}{p.pos_2}</span>}</td>
      <td className="px-1.5 py-1"><TypeCell cats={[p.category]} /></td>
      <td className="px-1.5 py-1 font-mono whitespace-nowrap" title={classText(p.classes)}>{Object.keys(p.classes).length > 1
        ? <button type="button" className="text-left hover:text-blue-700" aria-expanded={classOpen.has(p.id)} aria-label={`${classOpen.has(p.id) ? 'Hide' : 'Show'} the per-class breakdown`} data-testid={`iclass-toggle-${p.id}`}
            onClick={(e) => { e.stopPropagation(); toggleClasses(p.id) }}><span className="mr-0.5 font-sans text-stone-500">{classOpen.has(p.id) ? '▾' : '▸'}</span>{Object.keys(p.classes).slice(0, 2).join(', ')}{Object.keys(p.classes).length > 2 ? ` +${Object.keys(p.classes).length - 2}` : ''}</button>
        : Object.keys(p.classes)[0]}</td>
      <ConsCell k={p.n_observed} n={p.n_eligible} f={p.frequency} testId={`icons-${p.id}`} cat={p.category} />
      <td className="px-1.5 py-1 text-right">{p.n_entries_observed}</td>
      <td className="px-1.5 py-1">{s ? <StatusBadge s={s} testId={`istatus-${p.id}`} compact /> : '–'}</td>
    </tr> }
  /** one row per FR3D class of an interaction: chains with that class / eligible chains, its PDB entries, and the displayed chain's class */
  const classRows = (p: InteractionPair, sub: boolean) => !classOpen.has(p.id) ? [] : classBreakdown(ds, p, pdbOf).map((c) => { const s = statusOf.get(p.id)
    const here = s?.state === 'present' ? (s.classes.includes(c.cls) ? 'present' : 'other') : s?.state ?? null
    return <tr key={`${p.id}#${c.cls}`} className="border-t border-stone-100 bg-stone-50 text-stone-700" data-testid={`iclass-row-${p.id}-${c.cls}`}>
      <td className={`py-0.5 pr-1.5 text-[11px] text-stone-500 ${sub ? 'pl-9' : 'pl-5'}`}>class</td><td />
      <td className="px-1.5 py-0.5 font-mono">{c.cls}</td>
      <td className="px-1.5 py-0.5 whitespace-nowrap" data-testid={`iclass-cons-${p.id}-${c.cls}`}>{c.chains}/{p.n_eligible} <span className="text-stone-500">{pct(c.chains, p.n_eligible)}</span></td>
      <td className="px-1.5 py-0.5 text-right">{c.entries}</td>
      <td className="px-1.5 py-0.5 text-[11px]">{here === 'present' ? <span className="text-emerald-800">● this class</span> : here === 'other' ? <span className="text-stone-500">other class ({s!.classes.join(', ')})</span> : here ? <span className="text-stone-500">{STATE_TEXT[here]}</span> : ''}</td>
    </tr> })
  return (
    <section className="grid min-w-0 content-start gap-y-1.5 lg:row-span-3 lg:grid-rows-subgrid" data-testid="ifamily-table">
      <div className="flex flex-wrap items-center gap-2 self-end text-xs">
        <h3 className="text-sm font-semibold text-stone-800">Interactions</h3><span className="text-stone-500" data-testid="ifamily-count">{shown.length} of {groups.length}</span>{reset}
        {renderSearch(() => { const g = sorted[0]; if (g) select(g.members[0].id) })}
      </div>
      <div ref={box} className="h-[460px] overflow-auto rounded border border-stone-200 bg-white" data-testid="itable-box">
        <table className={`w-full text-xs ${fill ? 'h-full' : ''}`} data-fill={fill ? 'true' : 'false'}>
          <thead className="sticky top-0 z-10 bg-stone-100 text-stone-700"><tr className="h-8">
            {th('pos', 'Interaction')}{th('type', 'Type')}{th('class', 'Interaction class', '', <InfoTip text={TIP.klass} label="About interaction classes" />)}{th('cons', 'Conservation', '', <InfoTip text={TIP.frequency} label="About conservation" />)}{th('entries', 'Entries', 'text-right')}{th('status', 'In displayed structure')}
          </tr></thead>
          <tbody data-testid="ifamily-tbody">
            {visible.flatMap((g) => g.members.length === 1 ? [memberRow(g.members[0], false), ...classRows(g.members[0], false)] : [
              <tr key={g.id} data-testid={`igroup-${g.id}`} data-row-id={g.id} data-frequency={g.frequency ?? ''} data-status={g.state ?? ''} aria-expanded={isOpen(g.id)}
                onClick={() => { if (!g.members.some((m) => m.id === selectedId)) select(g.members[0].id); else toggle(g.id) }}
                className={`cursor-pointer border-t border-stone-100 ${g.members.some((m) => m.id === selectedId) ? 'bg-amber-50' : 'hover:bg-stone-50'}`}>
                <td className="px-1.5 py-1 whitespace-nowrap" title={`Rfam ${g.group.i}–${g.group.j} · ${regionText(g.members[0])}`}><button type="button" className="mr-1 text-stone-500" aria-label={isOpen(g.id) ? 'Collapse' : 'Expand'} onClick={(e) => { e.stopPropagation(); toggle(g.id) }} data-testid={`igroup-toggle-${g.id}`}>{isOpen(g.id) ? '▾' : '▸'}</button>
                  <span className="font-medium">{pairLabel(g.group.i, g.group.j, N)}</span>{N && <span className="ml-1 text-[10px] text-stone-400">{g.group.i}–{g.group.j}</span>}</td>
                <td className="px-1.5 py-1"><TypeCell cats={g.members.map((m) => m.category)} /></td>
                <td className="px-1.5 py-1 font-mono text-[11px]">{g.members.map((m) => `${CATEGORY_CODE[m.category]} ${dominantClass(m)}`).join(' · ')}</td>
                <ConsCell k={g.k} n={g.n} f={g.frequency} testId={`icons-${g.id}`} cat="multiple" />
                <td className="px-1.5 py-1 text-right">{g.entries}</td>
                <td className="px-1.5 py-1 whitespace-nowrap">{g.state && <span className={STATE_STYLE[g.state]} data-state={g.state}><span aria-hidden className="mr-1">{STATE_ICON[g.state]}</span>{g.state === 'present' ? 'Present (any)' : STATE_TEXT[g.state]}</span>}</td>
              </tr>,
              ...(isOpen(g.id) ? g.members.flatMap((m) => [memberRow(m, true), ...classRows(m, true)]) : [])])}
            {!sorted.length && <tr><td colSpan={6} className="px-2 py-3 text-center text-stone-500">{noTypes ? 'No interaction types selected.' : search.kind !== 'none' ? 'No interaction matches the search in the selected types.' : 'No interactions match the filters.'}</td></tr>}
          </tbody>
        </table>
      </div>
      <div className="flex flex-col gap-1">{pager}<p className="text-[11px] text-stone-500">Conservation: chains with the interaction / eligible chains. “In displayed structure”: {displayedLabel}.</p></div>
    </section>)
}
function ConsCell({ k, n, f, testId, cat }: { k: number; n: number; f: number | null; testId: string; cat: InteractionCategory | 'multiple' }) {
  const any = cat === 'multiple', color = any ? multiColor(f) : categoryColor(cat, f)
  return <td className="px-1.5 py-1 whitespace-nowrap" data-testid={testId} data-color={color}>
    <span className="mr-1 inline-block h-2 w-2 rounded-sm align-middle" style={{ background: color }} aria-hidden />{any && <span className="text-stone-500">any </span>}{k}/{n} <span className="text-stone-500">{pct(k, n)}</span></td>
}

// ===================================================================================================== compare two
interface CompareGroup { group: InteractionGroup<CompareRow>; id: string; members: CompareRow[]; status: CompareStatus }
const worst = (rs: CompareRow[]) => COMPARE_ORDER.find((s) => rs.some((r) => r.status === s)) ?? 'shared'
const cmpDefault = (a: CompareGroup, b: CompareGroup) => COMPARE_ORDER.indexOf(a.status) - COMPARE_ORDER.indexOf(b.status) || a.group.i - b.group.i || a.group.j - b.group.j

function CompareMode({ ds, L, circular, labelFor, N, visiblePairs, pairById, selected, selectedId, select, search, renderSearch, noTypes, paging, rowOf, orientation, a, b, auto, setPair, onBack }: CommonProps & {
  a: string | null; b: string | null; auto: boolean; setPair: (a: string | null, b: string | null) => void; onBack: () => void }) {
  const { acc } = useFamily()
  const [emph, setEmph] = useState<'all' | 'diff'>('diff')
  const [statuses, setStatuses] = useState<Set<CompareStatus>>(new Set(COMPARE_ORDER))
  const [fit, setFit] = useState<FitResult>({ state: 'none' })
  const ixA = a ? indexInstance(ds, a) : null, ixB = b ? indexInstance(ds, b) : null
  const rowA = a ? rowOf.get(a) ?? null : null, rowB = b ? rowOf.get(b) ?? null : null
  const msA = useMapping(acc, rowA), msB = useMapping(acc, rowB)
  const idxA = useMemo(() => (msA.state === 'ready' ? indexMapping(msA.mapping) : null), [msA]), idxB = useMemo(() => (msB.state === 'ready' ? indexMapping(msB.mapping) : null), [msB])
  const visible = useMemo(() => new Set(visiblePairs.map((p) => p.id)), [visiblePairs])
  const all = useMemo(() => (ixA && ixB ? compareInstances(ixA, ixB) : []), [ixA, ixB])
  const rows = useMemo(() => all.filter((r) => visible.has(r.id)), [all, visible])
  const groups = useMemo<CompareGroup[]>(() => groupByEndpoints(rows).map((g) => ({ group: g, id: groupId(g), members: g.members, status: worst(g.members) })), [rows])
  const counts = useMemo(() => { const c: Record<string, number> = {}; for (const r of rows) c[r.status] = (c[r.status] ?? 0) + 1; return c }, [rows])
  const arcsFor = useCallback((side: 'a' | 'b') => groups.flatMap<DiagramArc>((g) => {
    const here = g.members.filter((r) => r[side].state === 'present'); if (!here.length) return []
    const status = worst(here), st = COMPARE_COLORS[status], muted = status === 'shared' || status === 'not_comparable', one = g.members.length === 1
    // colour = comparison status; line pattern = category (single interactions); double line = multiple interactions
    return [{ id: g.id, i: one ? g.members[0].pos_1 : g.group.i, j: one ? g.members[0].pos_2 : g.group.j, color: st.color, dash: one ? CATEGORY_STYLE[g.members[0].category].dash : undefined,
      category: one ? g.members[0].category : 'multiple', status, width: muted ? 1.3 : 2, opacity: emph === 'diff' && muted ? 0.25 : 0.95, order: st.order, multi: !one }] }), [groups, emph])
  const arcsA = useMemo(() => arcsFor('a'), [arcsFor]), arcsB = useMemo(() => arcsFor('b'), [arcsFor])
  const absentOf = (sn: string | null) => { const i = sn ? ds.instances.find((x) => x.seq_name === sn) : null; return new Set(i ? [...i.eligibility].flatMap((c, k) => (c === 'E' ? [] : [k + 1])) : []) }
  const absentA = useMemo(() => absentOf(a), [a, ds]), absentB = useMemo(() => absentOf(b), [b, ds])   // eslint-disable-line react-hooks/exhaustive-deps
  const tooltip = useCallback((id: string) => { const g = groups.find((x) => x.id === id); if (!g) return null
    const cls = (s: InstanceStatus) => (s.state === 'present' ? s.classes.join(', ') : STATE_TEXT[s.state])
    return <div data-testid="icmp-tooltip" className="space-y-0.5"><p className="font-semibold">{pairLabel(g.group.i, g.group.j, N)}{g.members.length > 1 ? ' · Multiple interactions' : ''}</p><p className="text-stone-600">Rfam {g.group.i}–{g.group.j}</p>
      {g.members.map((r) => <p key={r.id}><span className="font-mono">{CATEGORY_CODE[r.category]}</span> <b>{COMPARE_TEXT[r.status]}</b> · A: {cls(r.a)} · B: {cls(r.b)}</p>)}</div> }, [groups, N])
  const selGroup = selectedId ? groups.find((g) => g.members.some((m) => m.id === selectedId)) ?? null : null
  const selectGroup = useCallback((id: string) => { const g = groups.find((x) => x.id === id); select(g ? g.members[0].id : id) }, [groups, select])
  const selRow = selectedId ? all.find((r) => r.id === selectedId) ?? null : null
  const focus = useMemo(() => (selected ? [selected.pos_1, selected.pos_2] : selRow ? [selRow.pos_1, selRow.pos_2] : NO_POS), [selected, selRow])
  // 3D: the selected interaction's nucleotides in both chains; local fit on the first nucleotide (the second when the first is not modelled in both)
  const { viewerFocus, fitPos } = useMemo(() => {
    if (!a || !b || focus.length === 0) return { viewerFocus: null, fitPos: null }
    const residues = { [a]: focus.map((p) => labelSeqAt(idxA, p)).filter((x): x is number => x != null), [b]: focus.map((p) => labelSeqAt(idxB, p)).filter((x): x is number => x != null) }
    const fp = focus.length === 2 ? focus.find((p) => labelSeqAt(idxA, p) != null && labelSeqAt(idxB, p) != null) ?? null : null
    return { viewerFocus: { residues, fit: fp != null ? { [a]: labelSeqAt(idxA, fp)!, [b]: labelSeqAt(idxB, fp)! } : null }, fitPos: fp }
  }, [a, b, idxA, idxB, focus.join(',')])   // eslint-disable-line react-hooks/exhaustive-deps
  const partner = fitPos != null && focus.length === 2 ? focus.find((p) => p !== fitPos) ?? null : null
  const entries = useMemo(() => [rowA, rowB].flatMap((r, k) => (r ? [{ seqName: r.seq_name, label: instLabel(r), color: k === 0 ? COLOR_A : COLOR_B, pdbId: r.pdb_id, labelAsymId: r.label_asym_id, transform: orientation?.transforms[r.seq_name] ?? null }] : [])), [rowA, rowB, orientation])
  const options = [...rowOf.values()]
  const sel = (v: string | null, set: (sn: string) => void, id: string, other: string | null) => <select value={v ?? ''} onChange={(e) => set(e.target.value)} data-testid={id} className="max-w-[15rem] rounded border border-stone-300 bg-white px-1 py-px">
    {!v && <option value="">Choose a chain…</option>}{options.map((i) => <option key={i.seq_name} value={i.seq_name} disabled={i.seq_name === other}>{instLabel(i)}{i.resolution != null ? ` · ${i.resolution} Å` : ''}</option>)}</select>
  const fitDetails = fit.state === 'local' && fitPos != null ? `Chain B is superposed on chain A by a least-squares fit (quaternion method) of the base heavy atoms of ${labelFor(fitPos).text} (Rfam ${fitPos}) only: ${fit.nAtoms} corresponding atoms (by nucleotide identity), RMSD ${fit.rmsd?.toFixed(2)} Å. The partner nucleotide is not fitted, so its position relative to that base is preserved.` : ''
  return (
    <>
      <div className="flex flex-wrap items-center gap-x-4 gap-y-2 text-xs" data-testid="icompare-controls">
        <button type="button" onClick={onBack} className="rounded border border-stone-300 bg-white px-2 py-0.5 hover:bg-stone-100" data-testid="icompare-back">← Back to family conservation</button>
        <label className="flex items-center gap-1"><span className="inline-block h-2.5 w-2.5 rounded-full" style={{ background: COLOR_A }} />Structure A {sel(a, (v) => setPair(v, b), 'icompare-a', b)}</label>
        <label className="flex items-center gap-1"><span className="inline-block h-2.5 w-2.5 rounded-full" style={{ background: COLOR_B }} />Structure B {sel(b, (v) => setPair(a, v), 'icompare-b', a)}</label>
        <div className="flex overflow-hidden rounded border border-stone-300" role="group" aria-label="Arc emphasis">
          {(['diff', 'all'] as const).map((m) => <button key={m} type="button" aria-pressed={emph === m} onClick={() => setEmph(m)} data-testid={`iemph-${m}`}
            className={`px-2 py-0.5 ${emph === m ? 'bg-stone-700 text-white' : 'bg-white text-stone-700 hover:bg-stone-100'}`}>{m === 'all' ? 'All interactions' : 'Show differences'}</button>)}
        </div>
        {auto && a && b && <p className="basis-full text-stone-600" data-testid="icompare-auto">Showing the pair with the highest geometric discrepancy among eligible comparisons (both chains ≤ {DEFAULT_PAIR_MAX_RESOLUTION.toFixed(1)} Å).</p>}
        {a && b && <CoverageNote ds={ds} a={a} b={b} rowOf={rowOf} />}
      </div>
      {!(rowA && rowB) ? <p className="rounded border border-dashed border-stone-300 bg-white p-4 text-sm text-stone-600" data-testid="icompare-choose">Choose two RNA chains to compare their interactions.</p> : <>
        {circular && <>
          <div className="grid gap-3 md:grid-cols-2" data-testid="icompare-diagrams">
            {([['a', rowA, arcsA, absentA], ['b', rowB, arcsB, absentB]] as const).map(([side, r, arcs, absent]) => (
              <section key={side} className="rounded border border-stone-200 bg-white p-2">
                <h3 className="px-1 text-sm font-semibold text-stone-800"><span className="mr-1 inline-block h-2.5 w-2.5 rounded-full" style={{ background: side === 'a' ? COLOR_A : COLOR_B }} />{side.toUpperCase()} · {instLabel(r)} <span className="font-normal text-stone-500">{r.resolution != null ? `${r.resolution} Å · ` : ''}{arcs.length} shown</span></h3>
                <div className="flex justify-center"><CircularDiagram L={L} size={SIZE} arcs={arcs} labelFor={labelFor} regions={ds.position_regions} absent={absent} selectedId={selGroup?.id ?? null} focus={focus} onSelect={selectGroup} tooltip={tooltip} testId={`idiagram-${side}`} title={`Interactions of ${instLabel(r)}`} /></div>
              </section>))}
          </div>
          <div className="flex flex-wrap items-center gap-x-3 gap-y-1 px-1 text-[11px] text-stone-600" data-testid="icompare-legend">
            <span className="font-medium text-stone-700">Colour = status:</span>
            {([['shared', 'Shared'], ['only_a', 'Only in this structure'], ['reclassified', 'Different class'], ['not_comparable', 'Not comparable in the other']] as const).map(([k, t]) => { const c = COMPARE_COLORS[k]
              return <span key={k} className="flex items-center gap-1" title={COMPARE_TIP[k]} data-testid={`icmp-legend-${k}`}><svg width="22" height="8" aria-hidden><line x1="1" y1="4" x2="21" y2="4" stroke={c.color} strokeWidth={k === 'shared' || k === 'not_comparable' ? 1.4 : 2.2} /></svg>{t}</span> })}
            <span className="font-medium text-stone-700">Pattern = type:</span>
            {CATEGORIES.map((c) => <span key={c} className="flex items-center gap-1"><svg width="24" height="8" aria-hidden><line x1="2" y1="4" x2="22" y2="4" stroke="#57534e" strokeWidth={1.8} strokeDasharray={CATEGORY_STYLE[c].dash} strokeLinecap="round" /></svg>{CATEGORY_CODE[c]}</span>)}
            <span className="flex items-center gap-1"><svg width="22" height="8" aria-hidden><line x1="1" y1="4" x2="21" y2="4" stroke="#78716c" strokeWidth={3} /><line x1="1" y1="4" x2="21" y2="4" stroke="white" strokeWidth={0.9} /></svg>multiple interactions</span>
            <span className="flex items-center gap-1"><svg width="10" height="10" aria-hidden><circle cx="5" cy="5" r="2.5" fill="white" stroke="#78716c" /></svg>nucleotide missing, incomplete or unmapped</span>
          </div>
        </>}
        <div className="grid gap-x-3 gap-y-4 lg:grid-cols-2 lg:gap-y-1.5" data-testid="ipanel-grid">
          <CompareTable groups={groups} counts={counts} statuses={statuses} setStatuses={setStatuses} N={N} pairById={pairById} selectedId={selectedId} select={select} nOff={[ixA?.nOffConsensus ?? 0, ixB?.nOffConsensus ?? 0]} search={search} renderSearch={renderSearch} noTypes={noTypes} paging={paging} />
          <section className="grid min-w-0 content-start gap-y-1.5 lg:row-span-3 lg:grid-rows-subgrid" data-testid="icompare-viewer">
            <h3 className="self-end text-sm font-semibold text-stone-800">3D superposition</h3>
            <InteractionPairViewer entries={entries} focus={viewerFocus} onFit={setFit} />
            <p className="text-xs text-stone-600" data-testid="icompare-fit" data-state={fit.state}>
              {fit.state === 'local' && fitPos != null ? <>Locally aligned on {labelFor(fitPos).text}{partner != null ? <> to compare the relative position of {labelFor(partner).text}</> : null}.<InfoTip text={fitDetails} label="Local superposition details" /></>
                : fit.state === 'kept' ? <>No local alignment for this selection ({fit.message}); previous alignment kept.</>
                  : fit.state === 'failed' ? <>No local alignment possible ({fit.message}); whole-RNA superposition.</>
                    : <>Whole-RNA superposition (as in the Structures tab). Select an interaction to align locally.</>}
              {selRow && (selRow.a.state === 'not_comparable' || selRow.b.state === 'not_comparable') && <span className="block text-amber-800" data-testid="icompare-nc">Not comparable: {[selRow.a.state === 'not_comparable' && `A — ${selRow.a.reason}`, selRow.b.state === 'not_comparable' && `B — ${selRow.b.reason}`].filter(Boolean).join('; ')}</span>}
            </p>
          </section>
        </div>
      </>}
    </>)
}

/** FR3D annotation coverage of the two compared chains (consensus helix pairs annotated / eligible), with the family median for reference;
 *  a chain clearly below the family (< 80 %) is highlighted, because its "only in the other" differences may be annotation gaps. */
function CoverageNote({ ds, a, b, rowOf }: { ds: InteractionsDataset; a: string; b: string; rowOf: Map<string, InstanceRow> }) {
  const med = useMemo(() => medianHelixCoverage(ds), [ds])
  const cov = useMemo(() => [a, b].map((sn) => ({ sn, ...helixCoverage(ds, sn) })), [ds, a, b])
  return (
    <p className="basis-full text-stone-600" data-testid="icompare-coverage">
      Helix base pairs annotated by FR3D:
      {cov.map((c, k) => { const low = c.fraction != null && c.fraction < 0.8; const r = rowOf.get(c.sn)
        return <span key={c.sn} className={`ml-2 ${low ? 'font-medium text-amber-800' : ''}`} data-testid={`icompare-coverage-${k === 0 ? 'a' : 'b'}`} data-fraction={c.fraction ?? ''}>
          {k === 0 ? 'A' : 'B'} {r ? instLabel(r) : c.sn}: {c.k} of {c.n} eligible ({pct(c.k, c.n)}){low ? ' — sparse annotation; “only in” differences may be annotation gaps' : ''}{k === 0 ? ' ·' : ''}</span> })}
      {med != null && <span className="ml-2 text-stone-500">family median {Math.round(100 * med)}%</span>}
      <InfoTip text="Of the Rfam consensus (SS_cons) helix base pairs that can be evaluated in the chain (both bases resolved and complete), how many BGSU/FR3D annotates. Distorted or moderate-resolution models can have paired nucleotides that FR3D does not annotate; BGSU exports no near base pairs." label="About annotation coverage" />
    </p>)
}

function CompareTable({ groups, counts, statuses, setStatuses, N, pairById, selectedId, select, nOff, search, renderSearch, noTypes, paging }: { groups: CompareGroup[]; counts: Record<string, number>; statuses: Set<CompareStatus>; setStatuses: (s: Set<CompareStatus>) => void
  N: Numbering | null; pairById: Map<string, InteractionPair>; selectedId: string | null; select: (id: string) => void; nOff: [number, number]; search: SearchQuery; renderSearch: (onEnter: () => void) => ReactNode; noTypes: boolean; paging: Paging }) {
  const hit = useCallback((r: CompareRow) => { const p = pairById.get(r.id)
    return search.kind === 'none' ? true : search.kind === 'error' ? false : search.kind === 'residue' ? matchesPositions(r, search.positions)
      : `${interactionLabel(r, N)} rfam ${r.pos_1}-${r.pos_2} ${CATEGORY_CODE[r.category]} ${CATEGORY_TEXT[r.category]} ${r.a.classes.join(' ')} ${r.b.classes.join(' ')} ${COMPARE_TEXT[r.status]} ${p ? regionText(p) : ''}`.toLowerCase().includes(search.text) }, [search, N, pairById])
  const shown = useMemo(() => groups.map((g) => ({ ...g, members: g.members.filter((r) => statuses.has(r.status) && hit(r)) })).filter((g) => g.members.length).map((g) => ({ ...g, status: worst(g.members) })), [groups, statuses, hit])
  const key = useCallback((g: CompareGroup, k: SortKey): string | number => { const r = g.members[0]
    return k === 'pos' ? g.group.i * 10000 + g.group.j : k === 'type' ? (g.members.length > 1 ? 'zz' : CATEGORY_CODE[r.category]) : k === 'status' ? COMPARE_ORDER.indexOf(g.status) : k === 'region' ? regionText(pairById.get(r.id) ?? { region_1: null, region_2: null }) : r.a.classes.join() }, [pairById])
  const { sorted, th, reset, sortKey } = useSort(shown, key, cmpDefault)
  const { visible, pager, fill } = usePaged(sorted, selectedId, paging, `${sortKey}|${JSON.stringify(search)}|${groups.length}|${[...statuses].join()}`)
  const box = useScrollToSelected(selectedId)
  const { isOpen, toggle } = useExpanded(selectedId, shown)
  const toggleStatus = (s: CompareStatus) => { const n = new Set(statuses); if (n.has(s)) n.delete(s); else n.add(s); setStatuses(n) }
  const cls = (s: InstanceStatus) => (s.state === 'present' ? <span className="font-mono">{s.classes.join(', ')}</span> : <span className={STATE_STYLE[s.state]} title={s.reason ?? (s.other.length ? `Other: ${s.other.join(', ')}` : 'not annotated')}>{s.state === 'absent' ? '—' : 'n/c'}</span>)
  const statusCell = (st: CompareStatus) => <><span className="mr-1 inline-block h-2 w-2 rounded-full align-middle" style={{ background: COMPARE_COLORS[st].color }} />{COMPARE_TEXT[st]}</>
  const memberRow = (r: CompareRow, sub: boolean) => { const on = r.id === selectedId, p = pairById.get(r.id)
    return <tr key={r.id} data-row-id={r.id} data-testid={`icmp-row-${r.id}`} data-status={r.status} aria-selected={on} onClick={() => select(r.id)} className={`cursor-pointer border-t border-stone-100 ${on ? 'bg-amber-100' : sub ? 'bg-stone-50/60 hover:bg-stone-100' : 'hover:bg-stone-50'}`}>
      <td className={`py-1 pr-1.5 whitespace-nowrap ${sub ? 'pl-5' : 'pl-1.5'}`} title={`Rfam ${r.pos_1}–${r.pos_2}`}><span className="font-medium">{interactionLabel(r, N)}</span>{N && <span className="ml-1 text-[10px] text-stone-400">{r.pos_1}{p?.directional ? '→' : '–'}{r.pos_2}</span>}</td>
      <td className="px-1.5 py-1"><TypeCell cats={[r.category]} /></td>
      <td className="px-1.5 py-1 whitespace-nowrap">{statusCell(r.status)}</td>
      <td className="px-1.5 py-1">{cls(r.a)}</td><td className="px-1.5 py-1">{cls(r.b)}</td>
      <td className="px-1.5 py-1 whitespace-nowrap text-stone-600">{p ? regionText(p) : '–'}</td></tr> }
  return (
    <section className="grid min-w-0 content-start gap-y-1.5 lg:row-span-3 lg:grid-rows-subgrid" data-testid="icompare-table">
      <div className="flex flex-wrap items-center gap-1.5 self-end text-xs">
        <h3 className="mr-1 text-sm font-semibold text-stone-800">Comparison</h3>
        {COMPARE_ORDER.map((s) => <button key={s} type="button" aria-pressed={statuses.has(s)} onClick={() => toggleStatus(s)} title={COMPARE_TIP[s]} data-testid={`icmp-filter-${s}`}
          className={`rounded-full border px-2 py-px ${statuses.has(s) ? 'border-stone-500 bg-stone-700 text-white' : 'border-stone-300 bg-white text-stone-500'}`}>{COMPARE_TEXT[s]} (<span data-testid={`icmp-count-${s}`}>{counts[s] ?? 0}</span>)</button>)}
        {reset}
        {renderSearch(() => { const g = sorted[0]; if (g) select(g.members[0].id) })}
      </div>
      <div ref={box} className="h-[460px] overflow-auto rounded border border-stone-200 bg-white" data-testid="itable-box">
        <table className={`w-full text-xs ${fill ? 'h-full' : ''}`} data-fill={fill ? 'true' : 'false'}>
          <thead className="sticky top-0 z-10 bg-stone-100 text-stone-700"><tr className="h-8">{th('pos', 'Interaction')}{th('type', 'Type')}{th('status', 'Status')}{th('class', 'Class in A', '', <InfoTip text={TIP.klass} label="About interaction classes" />)}<th className="px-1.5 py-1 text-left font-medium">Class in B</th>{th('region', 'Region')}</tr></thead>
          <tbody data-testid="icompare-tbody">
            {visible.flatMap((g) => g.members.length === 1 ? [memberRow(g.members[0], false)] : [
              <tr key={g.id} data-testid={`icmp-group-${g.id}`} data-row-id={g.id} data-status={g.status} aria-expanded={isOpen(g.id)} onClick={() => { if (!g.members.some((m) => m.id === selectedId)) select(g.members[0].id); else toggle(g.id) }}
                className={`cursor-pointer border-t border-stone-100 ${g.members.some((m) => m.id === selectedId) ? 'bg-amber-50' : 'hover:bg-stone-50'}`}>
                <td className="px-1.5 py-1 whitespace-nowrap"><button type="button" className="mr-1 text-stone-500" aria-label={isOpen(g.id) ? 'Collapse' : 'Expand'} onClick={(e) => { e.stopPropagation(); toggle(g.id) }}>{isOpen(g.id) ? '▾' : '▸'}</button><span className="font-medium">{pairLabel(g.group.i, g.group.j, N)}</span></td>
                <td className="px-1.5 py-1"><TypeCell cats={g.members.map((m) => m.category)} /></td>
                <td className="px-1.5 py-1 whitespace-nowrap">{statusCell(g.status)}</td>
                <td className="px-1.5 py-1 font-mono text-[11px]" colSpan={2}>{g.members.map((m) => `${CATEGORY_CODE[m.category]}: ${COMPARE_TEXT[m.status]}`).join(' · ')}</td>
                <td className="px-1.5 py-1 whitespace-nowrap text-stone-600">{(() => { const p = pairById.get(g.members[0].id); return p ? regionText(p) : '–' })()}</td></tr>,
              ...(isOpen(g.id) ? g.members.map((m) => memberRow(m, true)) : [])])}
            {!sorted.length && <tr><td colSpan={6} className="px-2 py-3 text-center text-stone-500">{noTypes ? 'No interaction types selected.' : search.kind !== 'none' ? 'No interaction matches the search.' : 'No interactions match the filters.'}</td></tr>}
          </tbody>
        </table>
      </div>
      <div className="flex flex-col gap-1">{pager}<p className="text-[11px] text-stone-500">“—”: eligible, not annotated. “n/c”: not comparable. Interactions with an inserted or flanking nucleotide are not compared ({nOff[0]} in A, {nOff[1]} in B).</p></div>
    </section>)
}
