import type { InteractionCategory, InteractionInstance, InteractionPair, InteractionsDataset, StructureClusteringDataset } from '../types/backend'
import type { Numbering } from './numbering'

/*
 * Interactions tab logic (pure). The family-wide numbers (observed / eligible, entries, class distributions) are precomputed by
 * `rna_explorer interactions` (interactions.json); here only per-instance status and two-instance comparisons are derived from the same
 * per-instance observations and eligibility codes, with the same definitions: an interaction is evaluated in an instance only when its chain
 * is FR3D-annotated and both positions are match residues with coordinates and complete bases; otherwise it is "not comparable", never absent.
 * Identifiers are Rfam-position based (`${category}:${pos_1}-${pos_2}`); displayed labels never identify anything.
 */

export const pairKey = (cat: InteractionCategory, i: number, j: number) => `${cat}:${i}-${j}`
export const CATEGORIES: InteractionCategory[] = ['basepair', 'stacking', 'baseribose', 'basephosphate']
export const CATEGORY_TEXT: Record<InteractionCategory, string> = { basepair: 'Base pairing', stacking: 'Base stacking', baseribose: 'Base–ribose', basephosphate: 'Base–phosphate' }
/** BGSU short codes */
export const CATEGORY_CODE: Record<InteractionCategory, string> = { basepair: 'Bp', stacking: 'Bs', baseribose: 'Br', basephosphate: 'Bph' }
export const isDirectional = (c: InteractionCategory) => c === 'baseribose' || c === 'basephosphate'
/** Unordered endpoint key: interactions of any category between the same two nucleotides share it (display grouping only). */
export const endpointKey = (i: number, j: number) => (i < j ? `${i}-${j}` : `${j}-${i}`)

/** Circular diagrams only for families whose Rfam consensus model has at most this many positions (UX decision; table and 3D otherwise). */
export const CIRCULAR_MAX_CLEN = 300
export const showCircular = (clen: number) => clen <= CIRCULAR_MAX_CLEN

/** Availability of a category in a dataset: 'none' = no source files (or a dataset built before the category existed), 'empty' = annotated
 *  but no consensus interaction observed, 'ok' otherwise. */
export function categoryAvailability(ds: InteractionsDataset, c: InteractionCategory): 'ok' | 'empty' | 'none' {
  const by = (ds.summary as { by_category?: Record<string, { n_entries_with_file?: number; n_consensus_pairs?: number }> }).by_category
  const has = by ? (by[c]?.n_entries_with_file ?? 0) > 0 : c === 'basepair' || c === 'stacking'
  if (!has) return 'none'
  return ds.pairs.some((p) => p.category === c) ? 'ok' : 'empty'
}

/** One instance's interactions on consensus position pairs. */
export interface InstanceIndex {
  inst: InteractionInstance
  /** pairKey -> non-near FR3D classes (sorted) */
  byKey: Map<string, string[]>
  /** endpointKey -> every annotation between the two positions ("Bs s35", "Bp ncWW", "Br 9BR"), near ones included */
  byPos: Map<string, string[]>
  /** interactions with an inserted or flanking (non-consensus) endpoint: kept in the data, not compared on consensus positions */
  nOffConsensus: number
}
const indexCache = new WeakMap<InteractionsDataset, Map<string, InstanceIndex>>()
export function indexInstance(ds: InteractionsDataset, seqName: string): InstanceIndex | null {
  let m = indexCache.get(ds); if (!m) { m = new Map(); indexCache.set(ds, m) }
  const hit = m.get(seqName); if (hit) return hit
  const inst = ds.instances.find((x) => x.seq_name === seqName); if (!inst) return null
  const byKey = new Map<string, string[]>(), byPos = new Map<string, string[]>(); let off = 0
  // schema 1.2+: compact per-chain index (observed pair indices by class + near annotations); older files: every observation inline
  const rows: [InteractionCategory, string, number | string | null, number | string | null][] = ds.observations
    ? (ds.observations[seqName] ?? []).map((o) => [o[0], o[1], o[4], o[5]])
    : [...Object.entries(inst.observed ?? {}).flatMap(([cls, ks]) => ks.map((k) => [ds.pairs[k].category, cls, ds.pairs[k].pos_1, ds.pairs[k].pos_2] as [InteractionCategory, string, number, number])),
       ...(inst.near ?? [])]
  if (!ds.observations) off = inst.n_off_consensus ?? 0
  for (const [cat, cls, p1, p2] of rows) {
    if (typeof p1 !== 'number' || typeof p2 !== 'number') { off++; continue }
    const pk = endpointKey(p1, p2); (byPos.get(pk) ?? byPos.set(pk, []).get(pk)!).push(`${CATEGORY_CODE[cat]} ${cls}`)
    if (cls.startsWith('n')) continue
    const k = pairKey(cat, p1, p2), l = byKey.get(k) ?? byKey.set(k, []).get(k)!
    if (!l.includes(cls)) { l.push(cls); l.sort() }
  }
  const ix = { inst, byKey, byPos, nOffConsensus: off }; m.set(seqName, ix); return ix
}

type InteractionIndexLike = Pick<InstanceIndex, 'inst' | 'byKey' | 'byPos'>
export type InstanceState = 'present' | 'absent' | 'not_comparable'
export const STATE_TEXT: Record<InstanceState, string> = { present: 'Present', absent: 'Absent', not_comparable: 'Not comparable' }
export interface InstanceStatus {
  state: InstanceState
  /** FR3D classes of the interaction in this instance (present only) */
  classes: string[]
  /** other FR3D annotations between the same two positions (other category, or near), shown with an absent interaction */
  other: string[]
  /** why the interaction cannot be evaluated (not comparable only) */
  reason: string | null
}
const CODE_REASON: Record<string, string> = { M: 'nucleotide not modelled', '-': 'no residue at this Rfam position (deletion or outside the aligned region)' }
const INCOMPLETE: Record<'base' | 'sugar' | 'phosphate', string> = { base: 'incomplete base', sugar: 'incomplete ribose', phosphate: 'incomplete or absent phosphate' }
/**
 * Status of one interaction in one chain. Category-specific eligibility: base pairs and stacking need both bases complete; base-ribose needs
 * the base of pos_1 and the ribose of pos_2; base-phosphate the base of pos_1 and the phosphate of pos_2 (backbone codes from the dataset).
 */
export function statusIn(ix: InteractionIndexLike, cat: InteractionCategory, i: number, j: number): InstanceStatus {
  if (ix.inst.annotation !== 'available') return { state: 'not_comparable', classes: [], other: [], reason: `FR3D annotation unavailable (${ix.inst.annotation_note ?? 'no records'})` }
  const gap = (ix.inst.annotation_gaps ?? []).find(([a, b]) => (i >= a && i <= b) || (j >= a && j <= b))
  if (gap) return { state: 'not_comparable', classes: [], other: [], reason: `BGSU annotation incomplete for this chain (no annotation at Rfam ${gap[0]}–${gap[1]})` }
  const second = cat === 'baseribose' ? ix.inst.sugar_eligibility : cat === 'basephosphate' ? ix.inst.phosphate_eligibility : ix.inst.eligibility
  if (second == null || ix.inst.unavailable_categories?.includes(cat)) return { state: 'not_comparable', classes: [], other: [], reason: `${CATEGORY_TEXT[cat]} annotations unavailable for this chain` }
  const kind2 = cat === 'baseribose' ? 'sugar' : cat === 'basephosphate' ? 'phosphate' : 'base'
  const bad = ([[i, ix.inst.eligibility[i - 1] ?? '-', 'base'], [j, second[j - 1] ?? '-', kind2]] as const).filter(([, c]) => c !== 'E')
  if (bad.length) return { state: 'not_comparable', classes: [], other: [], reason: bad.map(([p, c, k]) => `Rfam ${p}: ${c === 'I' ? INCOMPLETE[k] : CODE_REASON[c] ?? c}`).join('; ') }
  const classes = ix.byKey.get(pairKey(cat, i, j)) ?? []
  const own = (s: string) => s.startsWith(`${CATEGORY_CODE[cat]} `) && classes.includes(s.slice(CATEGORY_CODE[cat].length + 1))
  const other = (ix.byPos.get(endpointKey(i, j)) ?? []).filter((s) => !own(s))
  return classes.length ? { state: 'present', classes, other, reason: null } : { state: 'absent', classes: [], other, reason: null }
}

/** canonical Watson-Crick or wobble pair: cWW with AU, GC or GU */
const CANONICAL = new Set(['AU', 'UA', 'GC', 'CG', 'GU', 'UG'])
export const isCanonical = (cls: string, nts: string) => cls === 'cWW' && CANONICAL.has(nts.toUpperCase())
/** family-wide: canonical when at least half of the observed instances are canonical cWW */
export const pairIsCanonical = (p: InteractionPair) => p.category === 'basepair' && p.n_observed > 0 && 2 * p.n_canonical >= p.n_observed
export const dominantClass = (p: InteractionPair) => Object.keys(p.classes)[0] ?? ''

/** Interaction filters (shared by both modes): one switch per BGSU category (Bp only by default), and base-pair-only refinements. */
export interface InteractionFilter { basepair: boolean; stacking: boolean; baseribose: boolean; basephosphate: boolean; canonical: 'all' | 'canonical' | 'noncanonical'; structure: 'all' | 'secondary' | 'tertiary'; region: string | null }
export const DEFAULT_FILTER: InteractionFilter = { basepair: true, stacking: false, baseribose: false, basephosphate: false, canonical: 'all', structure: 'all', region: null }
export function passes(f: InteractionFilter, p: Pick<InteractionPair, 'category' | 'structure_class' | 'region_1' | 'region_2'>, canonical: boolean): boolean {
  if (!f[p.category]) return false
  if (p.category === 'basepair') {
    if (f.canonical === 'canonical' && !canonical) return false
    if (f.canonical === 'noncanonical' && canonical) return false
    if (f.structure !== 'all' && p.structure_class !== f.structure) return false
  } else if (f.canonical !== 'all' || f.structure !== 'all') return false   // canonical / secondary / tertiary describe base pairs only (other categories are hidden while they are set)
  if (f.region && p.region_1 !== f.region && p.region_2 !== f.region) return false
  return true
}

/** Two-instance comparison of corresponding interactions. */
export type CompareStatus = 'shared' | 'only_a' | 'only_b' | 'reclassified' | 'not_comparable'
export const COMPARE_TEXT: Record<CompareStatus, string> = { shared: 'Shared', only_a: 'Only in A', only_b: 'Only in B', reclassified: 'Different class', not_comparable: 'Not comparable' }
/** short explanations (chip and legend tooltips); 'reclassified' is the internal key of "Different class" */
export const COMPARE_TIP: Record<CompareStatus, string> = {
  shared: 'Same interaction type and FR3D class in both structures', only_a: 'Annotated in A; both nucleotides eligible in B but not annotated there',
  only_b: 'Annotated in B; both nucleotides eligible in A but not annotated there', reclassified: 'Same interaction type, different FR3D class (e.g. cWW in A, tWW in B)',
  not_comparable: 'Cannot be evaluated in one structure (missing, incomplete or unmapped nucleotide, or no annotation)' }
export interface CompareRow { id: string; category: InteractionCategory; pos_1: number; pos_2: number; status: CompareStatus; a: InstanceStatus; b: InstanceStatus }
export function compareInstances(a: InstanceIndex, b: InstanceIndex): CompareRow[] {
  const keys = new Set([...a.byKey.keys(), ...b.byKey.keys()])
  const rows: CompareRow[] = []
  for (const k of keys) {
    const [cat, rest] = k.split(':') as [InteractionCategory, string]; const [i, j] = rest.split('-').map(Number)
    const sa = statusIn(a, cat, i, j), sb = statusIn(b, cat, i, j)
    const status: CompareStatus = sa.state === 'not_comparable' || sb.state === 'not_comparable' ? 'not_comparable'
      : sa.state === 'present' && sb.state === 'present' ? (sa.classes.join() === sb.classes.join() ? 'shared' : 'reclassified')
        : sa.state === 'present' ? 'only_a' : 'only_b'
    rows.push({ id: k, category: cat, pos_1: i, pos_2: j, status, a: sa, b: sb })
  }
  return rows.sort((x, y) => x.pos_1 - y.pos_1 || x.pos_2 - y.pos_2 || x.category.localeCompare(y.category))
}

/** Instances in which a family-wide interaction is observed (eligible, non-near), best supported first: the given preferred instances
 *  (e.g. the numbering reference), then by resolution. */
export function supportingInstances(ds: InteractionsDataset, p: InteractionPair, resolution: (sn: string) => number | null, preferred: string[] = []): { seqName: string; classes: string[] }[] {
  const out: { seqName: string; classes: string[] }[] = []
  for (const inst of ds.instances) {
    const ix = indexInstance(ds, inst.seq_name)!; const s = statusIn(ix, p.category, p.pos_1, p.pos_2)
    if (s.state === 'present') out.push({ seqName: inst.seq_name, classes: s.classes })
  }
  const rank = (sn: string) => { const k = preferred.indexOf(sn); return k < 0 ? preferred.length : k }
  return out.sort((x, y) => rank(x.seqName) - rank(y.seqName) || (resolution(x.seqName) ?? 99) - (resolution(y.seqName) ?? 99) || x.seqName.localeCompare(y.seqName))
}

// ------------------------------------------------------------------------------------------------ numbering labels and lookup
/** Label of a consensus position in the selected numbering: the reference residue ("G37", "PSU516", "A306B") or, when the reference has
 *  no match residue there (or the numbering is Rfam consensus), the Rfam position ("Rfam 41"). `short` is for the circle ("37" / "41"). */
export function positionLabel(pos: number, numbering: Numbering | null): { text: string; short: string; fromReference: boolean } {
  if (numbering) {
    const n = numbering.numbers[pos - 1]
    if (n != null) {
      const c = numbering.idx.byConsensus.get(pos), r = c?.polymer != null ? numbering.idx.byPolymer.get(c.polymer) : undefined
      return { text: `${r?.residue_name ?? ''}${n}`, short: n, fromReference: true }
    }
  }
  return { text: `Rfam ${pos}`, short: String(pos), fromReference: false }
}
export const pairLabel = (i: number, j: number, numbering: Numbering | null) => `${positionLabel(i, numbering).text}–${positionLabel(j, numbering).text}`

export type LookupResult = { kind: 'position'; positions: [number] } | { kind: 'pair'; positions: [number, number] } | { kind: 'error'; message: string }
const TOKEN = /^([A-Za-z]{0,3}?)(-?\d+)([A-Za-z]?)$/
/**
 * Residue / interaction lookup: "G37", "37", "G37-C61", "Rfam 41", "41-66". Numbers are read in the selected numbering (reference author
 * numbers incl. insertion code; Rfam positions when the numbering is Rfam consensus or the query starts with "Rfam"). An optional residue name
 * must agree with the reference residue. Inserted reference residues and numbers matching several positions are reported, never guessed.
 */
export function resolveLookup(query: string, numbering: Numbering | null, clen: number): LookupResult | null {
  let q = query.trim(); if (!q) return null
  const rfam = /^rfam\b/i.test(q) || !numbering
  q = q.replace(/^rfam\s*/i, '')
  const parts = q.split(/(?<=\d[A-Za-z]?)\s*(?:[-–—:/]|\s)\s*(?=[A-Za-z]{0,3}-?\d)/)
  if (parts.length > 2) return { kind: 'error', message: 'Use one residue (G37) or one pair (G37-C61).' }
  const out: number[] = []
  for (const t of parts) {
    const m = TOKEN.exec(t.trim()); if (!m) return { kind: 'error', message: `“${t.trim()}” is not a residue (e.g. G37, 37, Rfam 41).` }
    const [, name, num, ins] = m
    if (rfam) {
      const p = +num
      if (ins || p < 1 || p > clen) return { kind: 'error', message: `Rfam positions are 1–${clen}.` }
      out.push(p); continue
    }
    const nb = numbering!
    const want = `${num}${ins.toUpperCase()}`
    const hits = nb.numbers.flatMap((n, k) => (n === want ? [k + 1] : []))
    if (hits.length > 1) return { kind: 'error', message: `${want} is ambiguous in ${nb.label} (Rfam ${hits.join(', ')}).` }
    if (!hits.length) {
      const r = [...nb.idx.byPolymer.values()].find((x) => `${x.author_residue_number}${x.author_insertion_code ?? ''}` === want)
      return { kind: 'error', message: r ? `${r.residue_name}${want} of ${nb.label} is ${r.alignment_state === 'insert' ? 'an inserted residue' : 'outside the aligned region'}: it has no Rfam consensus position.` : `No residue ${want} in ${nb.label}.` }
    }
    const label = positionLabel(hits[0], nb).text
    if (name && !label.toUpperCase().startsWith(name.toUpperCase())) return { kind: 'error', message: `${want} of ${nb.label} is ${label}, not ${name.toUpperCase()}${want}.` }
    out.push(hits[0])
  }
  if (out.length === 2) { const [i, j] = out[0] < out[1] ? out : [out[1], out[0]]; return i === j ? { kind: 'position', positions: [i] } : { kind: 'pair', positions: [i, j] } }
  return { kind: 'position', positions: [out[0]] }
}

// ------------------------------------------------------------------------------------------------ circular layout
export interface CircleLayout { first: number; last: number; cx: number; cy: number; r: number; angle: (pos: number) => number; point: (pos: number, r?: number) => [number, number] }
/** Consensus positions clockwise from the top, 5' end just right of 12 o'clock, a fixed gap before it; identical for every diagram that
 *  shares (first, last, size). */
export function circleLayout(first: number, last: number, size: number, r: number, gapDeg = 10): CircleLayout {
  const n = Math.max(1, last - first), span = ((360 - gapDeg) * Math.PI) / 180, g = (gapDeg * Math.PI) / 360
  const angle = (pos: number) => -Math.PI / 2 + g + ((pos - first) / n) * span
  const cx = size / 2, cy = size / 2
  return { first, last, cx, cy, r, angle, point: (pos, rr = r) => [cx + rr * Math.cos(angle(pos)), cy + rr * Math.sin(angle(pos))] }
}
/** SVG path of an interaction: the circular arc orthogonal to the layout circle (bowing towards the centre); a straight chord for
 *  diametrically opposite positions. */
export function arcPath(L: CircleLayout, i: number, j: number): string {
  const [x1, y1] = L.point(i), [x2, y2] = L.point(j)
  let d = Math.abs(L.angle(j) - L.angle(i)); if (d > Math.PI) d = 2 * Math.PI - d
  if (Math.abs(d - Math.PI) < 0.02 || d < 1e-6) return `M${x1.toFixed(1)},${y1.toFixed(1)}L${x2.toFixed(1)},${y2.toFixed(1)}`
  const rr = L.r * Math.tan(d / 2)
  // centre of the orthogonal circle: along the bisector, outside the layout circle
  const mx = (x1 + x2) / 2 - L.cx, my = (y1 + y2) / 2 - L.cy, ml = Math.hypot(mx, my), dist = L.r / Math.cos(d / 2)
  const ccx = L.cx + (mx / ml) * dist, ccy = L.cy + (my / ml) * dist
  const dx = x2 - x1, dy = y2 - y1
  const sweep = (ccx - x1) * -dy + (ccy - y1) * dx > 0 ? 1 : 0
  return `M${x1.toFixed(1)},${y1.toFixed(1)}A${rr.toFixed(1)},${rr.toFixed(1)} 0 0 ${sweep} ${x2.toFixed(1)},${y2.toFixed(1)}`
}

/**
 * Visual encoding of the four BGSU categories (shared by every diagram and legend): a colour family AND a line pattern, so the category is
 * never conveyed by colour alone. Within a category, colour intensity (light -> dark) encodes conservation (observed / eligible chains);
 * the lightest step keeps >= 2:1 contrast on white so rare interactions stay visible. Mid-tones (#2171b5, #238b45, #d94801, #6a51a3) pass
 * the dataviz palette validator (CVD separation, contrast). Dash arrays are drawn with round caps (dotted = round dots).
 */
export const CATEGORY_STYLE: Record<InteractionCategory, { family: string; ramp: string[]; dash?: string; pattern: string }> = {
  basepair: { family: 'blue', ramp: ['#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b'], pattern: 'solid' },
  stacking: { family: 'green', ramp: ['#74c476', '#41ab5d', '#238b45', '#006d2c', '#00441b'], dash: '6 3.5', pattern: 'dashed' },
  baseribose: { family: 'orange', ramp: ['#fd8d3c', '#f16913', '#d94801', '#a63603', '#7f2704'], dash: '0.1 3.4', pattern: 'dotted' },
  basephosphate: { family: 'purple', ramp: ['#9e9ac8', '#807dba', '#6a51a3', '#54278f', '#3f007d'], dash: '7 2.6 0.1 2.6', pattern: 'dash-dot' },
}
/** "Multiple interactions" (display grouping): a neutral grey ramp (no category colour is borrowed), drawn as a double line */
export const MULTI_RAMP = ['#a8a29e', '#78716c', '#57534e', '#44403c', '#292524']
function interpolate(ramp: string[], f: number | null): string {
  const t = Math.min(1, Math.max(0, f ?? 0)) * (ramp.length - 1), k = Math.min(ramp.length - 2, Math.floor(t)), u = t - k
  const rgb = (h: string) => [1, 3, 5].map((o) => parseInt(h.slice(o, o + 2), 16))
  const [a, b] = [rgb(ramp[k]), rgb(ramp[k + 1])]
  return `#${a.map((x, n) => Math.round(x + (b[n] - x) * u).toString(16).padStart(2, '0')).join('')}`
}
/** colour of an interaction of `cat` observed in a fraction `f` of eligible chains */
export const categoryColor = (cat: InteractionCategory, f: number | null) => interpolate(CATEGORY_STYLE[cat].ramp, f)
/** colour of a "Multiple interactions" group from its explicit "any" occurrence */
export const multiColor = (f: number | null) => interpolate(MULTI_RAMP, f)

/**
 * Unified table search. Residue-like queries ("G37", "37", "G37-C61", "Rfam 41", "Rfam 41-66") are resolved with the validated lookup
 * (selected numbering, Rfam correspondence); FR3D class codes ("cWW", "s35", "9BR", "0BPh") and other words are matched as text.
 */
const CLASS_CODE = /^n?([ct][WHS][WHS]|s[35][35]|\d{1,2}BR|\d{1,2}BPh)$/i
const RESIDUE_QUERY = /^(rfam\s*)?[A-Za-z]{0,3}-?\d+[A-Za-z]?(\s*[-–—:/\s]\s*[A-Za-z]{0,3}-?\d+[A-Za-z]?)?$/i
export type SearchQuery = { kind: 'none' } | { kind: 'text'; text: string } | { kind: 'residue'; positions: number[] } | { kind: 'error'; message: string }
/** single-letter residue names that can start a residue query (others, e.g. the region names "J2/3", "L2", "P1", are text) */
const ONE_LETTER = /^[ACGUTN]$/i
export function parseSearch(q: string, numbering: Numbering | null, clen: number, textTerms: string[] = []): SearchQuery {
  const t = q.trim(); if (!t) return { kind: 'none' }
  const lower = t.toLowerCase()
  if (CLASS_CODE.test(t) || !RESIDUE_QUERY.test(t) || (/[a-z]/i.test(t) && textTerms.some((x) => x.toLowerCase().includes(lower)))) return { kind: 'text', text: lower }
  const names = t.replace(/^rfam\s*/i, '').split(/[^A-Za-z0-9-]+|(?<=\d[A-Za-z]?)-(?=[A-Za-z]{0,3}\d)/).map((x) => /^([A-Za-z]*)/.exec(x)![1])
  if (names.some((n) => n.length === 1 && !ONE_LETTER.test(n))) return { kind: 'text', text: lower }
  const r = resolveLookup(t, numbering, clen)
  if (!r) return { kind: 'none' }
  if (r.kind === 'error') return { kind: 'error', message: r.message }
  return { kind: 'residue', positions: r.positions }
}
/** does an interaction (any direction) match a residue query: one position = involves it; two = between exactly these two */
export const matchesPositions = (p: { pos_1: number; pos_2: number }, pos: number[]) =>
  pos.length === 1 ? p.pos_1 === pos[0] || p.pos_2 === pos[0] : Math.min(p.pos_1, p.pos_2) === pos[0] && Math.max(p.pos_1, p.pos_2) === pos[1]

/** Display grouping of interactions that share their two nucleotides ("Multiple interactions"); members keep their own ids and frequencies. */
export interface InteractionGroup<T extends { id: string; pos_1: number; pos_2: number }> { key: string; i: number; j: number; members: T[] }
export function groupByEndpoints<T extends { id: string; pos_1: number; pos_2: number }>(items: T[]): InteractionGroup<T>[] {
  const m = new Map<string, InteractionGroup<T>>()
  for (const x of items) { const k = endpointKey(x.pos_1, x.pos_2); const g = m.get(k) ?? m.set(k, { key: k, i: Math.min(x.pos_1, x.pos_2), j: Math.max(x.pos_1, x.pos_2), members: [] }).get(k)!; g.members.push(x) }
  return [...m.values()]
}
/** Arc / row id of a group: the member id when alone, `multi:i-j` otherwise */
export const groupId = <T extends { id: string; pos_1: number; pos_2: number }>(g: InteractionGroup<T>) => (g.members.length === 1 ? g.members[0].id : `multi:${g.key}`)
/**
 * Pair-level occurrence of a group (explicit, never an average or sum of member frequencies): denominator = chains eligible for at least one
 * member interaction; numerator = those chains in which at least one member interaction is observed.
 */
export function groupOccurrence(ds: InteractionsDataset, members: InteractionPair[]): { k: number; n: number; frequency: number | null } {
  let k = 0, n = 0
  for (const inst of ds.instances) {
    const ix = indexInstance(ds, inst.seq_name)!; let el = false, pr = false
    for (const p of members) { const s = statusIn(ix, p.category, p.pos_1, p.pos_2); if (s.state !== 'not_comparable') el = true; if (s.state === 'present') pr = true }
    if (el) { n++; if (pr) k++ }
  }
  return { k, n, frequency: n ? k / n : null }
}

/** Default family-table order (user decision 2026-10-19): supporting PDB entries desc, then frequency desc, eligible chains desc, Rfam
 *  positions asc (numeric values, never formatted strings). */
export function defaultFamilyOrder(a: { entries: number; frequency: number | null; n: number; i: number; j: number }, b: typeof a): number {
  return b.entries - a.entries || (b.frequency ?? -1) - (a.frequency ?? -1) || b.n - a.n || a.i - b.i || a.j - b.j
}

/** Default comparison pair: both chains at this resolution or better (Å). Moderate-resolution models can have distorted geometry that FR3D does
 *  not annotate (e.g. RF00167 9V50 A, 3.25 Å), and their large discrepancy then partly reflects model quality. User decision, 2026-10-19. */
export const DEFAULT_PAIR_MAX_RESOLUTION = 3.0
/**
 * Default pair of Compare two: the highest geometric discrepancy (structure_clustering.json, not recomputed) among chains of the clustered
 * set (the Structures tab's inclusion and coverage rule), with a valid (non-null) value, both chains at <= DEFAULT_PAIR_MAX_RESOLUTION Å
 * (chains without a resolution are excluded), both FR3D-annotated, and at least one interaction comparable in both.
 * Ties: lexicographically smallest (A, B) seq_names; A is the smaller id.
 */
export function defaultComparePair(ds: InteractionsDataset, sc: StructureClusteringDataset | null, resolution: (sn: string) => number | null = () => null,
  maxResolution: number | null = DEFAULT_PAIR_MAX_RESOLUTION): { a: string; b: string; value: number } | null {
  const okRes = (sn: string) => { if (maxResolution == null) return true; const r = resolution(sn); return r != null && r <= maxResolution }
  const m = sc?.metrics?.discrepancy as { values_upper?: (number | null)[] } | undefined
  const ids = sc?.dataset?.instances as string[] | undefined
  if (!m?.values_upper || !ids) return null
  const n = ids.length, cand: { a: string; b: string; value: number }[] = []
  let k = 0
  for (let r = 0; r < n; r++) for (let c = r + 1; c < n; c++, k++) {
    const v = m.values_upper[k]; if (v == null || !Number.isFinite(v)) continue
    const [a, b] = ids[r] < ids[c] ? [ids[r], ids[c]] : [ids[c], ids[r]]
    cand.push({ a, b, value: v })
  }
  cand.sort((x, y) => y.value - x.value || x.a.localeCompare(y.a) || x.b.localeCompare(y.b))
  for (const p of cand) {
    const A = indexInstance(ds, p.a), B = indexInstance(ds, p.b)
    if (!okRes(p.a) || !okRes(p.b)) continue
    if (!A || !B || A.inst.annotation !== 'available' || B.inst.annotation !== 'available') continue
    if (compareInstances(A, B).some((r) => r.status !== 'not_comparable')) return p
  }
  return null
}

// ------------------------------------------------------------------------------------------------ local superposition
type V3 = [number, number, number]
const PURINE_RING = ['N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4'], PYR_RING = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6']
const EXOCYCLIC: Record<string, string[]> = { A: ['N6'], G: ['O6', 'N2'], C: ['O2', 'N4'], U: ['O2', 'O4'] }
/**
 * Corresponding base heavy atoms of two nucleotides (atom names of A, atom names of B): same parent nucleotide -> ring + exocyclic atoms
 * present in both; two purines or two pyrimidines -> the ring atoms present in both; purine vs pyrimidine -> the glycosidic frame
 * (N9/N1, C4/C2, C8/C6) plus C1'. Modified residues are matched through their parent nucleotide and the atoms they actually have.
 */
export function baseCorrespondence(a: { names: Set<string>; parent: string }, b: { names: Set<string>; parent: string }): [string, string][] {
  const pa = a.names.has('N9'), pb = b.names.has('N9')
  if (pa === pb) {
    const ring = pa ? PURINE_RING : PYR_RING
    const extra = a.parent === b.parent ? EXOCYCLIC[a.parent.toUpperCase()] ?? [] : []
    return [...ring, ...extra].filter((n) => a.names.has(n) && b.names.has(n)).map((n) => [n, n])
  }
  const map: [string, string][] = [['N9', 'N1'], ['C4', 'C2'], ['C8', 'C6'], ["C1'", "C1'"]]
  return (pa ? map : map.map(([x, y]) => [y, x] as [string, string])).filter(([x, y]) => a.names.has(x) && b.names.has(y))
}

/** Symmetric eigen-decomposition (cyclic Jacobi); eigenvectors in columns. */
function jacobi(A: number[][]): { values: number[]; vectors: number[][] } {
  const n = A.length, a = A.map((r) => [...r]), V: number[][] = a.map((_, i) => a.map((__, j) => (i === j ? 1 : 0)))
  for (let sweep = 0; sweep < 50; sweep++) {
    let off = 0; for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += a[p][q] ** 2
    if (off < 1e-20) break
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(a[p][q]) < 1e-15) continue
      const th = (a[q][q] - a[p][p]) / (2 * a[p][q]), t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c
      for (let k = 0; k < n; k++) { const kp = a[k][p], kq = a[k][q]; a[k][p] = c * kp - s * kq; a[k][q] = s * kp + c * kq }
      for (let k = 0; k < n; k++) { const pk = a[p][k], qk = a[q][k]; a[p][k] = c * pk - s * qk; a[q][k] = s * pk + c * qk }
      for (let k = 0; k < n; k++) { const kp = V[k][p], kq = V[k][q]; V[k][p] = c * kp - s * kq; V[k][q] = s * kp + c * kq }
    }
  }
  return { values: a.map((r, i) => r[i]), vectors: V }
}
/** Least-squares rigid fit (Horn's quaternion method) of `mobile` onto `target` (corresponding points): x' = R x + t, with the RMSD. */
export function kabsch(mobile: V3[], target: V3[]): { R: number[][]; t: V3; rmsd: number } | null {
  const n = mobile.length; if (n < 3 || target.length !== n) return null
  const cm = [0, 1, 2].map((k) => mobile.reduce((s, p) => s + p[k], 0) / n), ct = [0, 1, 2].map((k) => target.reduce((s, p) => s + p[k], 0) / n)
  const S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
  for (let k = 0; k < n; k++) for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) S[i][j] += (mobile[k][i] - cm[i]) * (target[k][j] - ct[j])
  const [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = S
  const N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx], [yz - zy, xx - yy - zz, xy + yx, zx + xz], [zx - xz, xy + yx, -xx + yy - zz, yz + zy], [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]
  const { values, vectors } = jacobi(N); const k = values.indexOf(Math.max(...values))
  const [w, x, y, z] = [0, 1, 2, 3].map((i) => vectors[i][k])
  const R = [[w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)], [2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)], [2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z]]
  const t = [0, 1, 2].map((i) => ct[i] - (R[i][0] * cm[0] + R[i][1] * cm[1] + R[i][2] * cm[2])) as V3
  const apply = (p: V3) => [0, 1, 2].map((i) => R[i][0] * p[0] + R[i][1] * p[1] + R[i][2] * p[2] + t[i])
  const rmsd = Math.sqrt(mobile.reduce((s, p, i) => { const q = apply(p); return s + (q[0] - target[i][0]) ** 2 + (q[1] - target[i][1]) ** 2 + (q[2] - target[i][2]) ** 2 }, 0) / n)
  return { R, t, rmsd }
}

/**
 * Annotation coverage of a chain: of the Rfam consensus helix base pairs (SS_cons position pairs observed in the family) that can be evaluated
 * in the chain (both bases eligible), how many FR3D annotates (any base-pair class). Low coverage flags chains whose models FR3D annotates
 * sparsely (e.g. distorted geometry at moderate resolution), so that "only in A/B" differences are read with care.
 */
export function helixCoverage(ds: InteractionsDataset, seqName: string): { k: number; n: number; fraction: number | null } {
  const ix = indexInstance(ds, seqName); if (!ix) return { k: 0, n: 0, fraction: null }
  let k = 0, n = 0
  for (const p of ds.pairs) {
    if (p.category !== 'basepair' || !p.in_ss_cons) continue
    const st = statusIn(ix, 'basepair', p.pos_1, p.pos_2)
    if (st.state === 'not_comparable') continue
    n++; if (st.state === 'present') k++
  }
  return { k, n, fraction: n ? k / n : null }
}
/** median helix coverage over the family's annotated chains (reference for one chain's value) */
export function medianHelixCoverage(ds: InteractionsDataset): number | null {
  const f = ds.instances.map((i) => helixCoverage(ds, i.seq_name).fraction).filter((x): x is number => x != null).sort((a, b) => a - b)
  if (!f.length) return null
  const m = Math.floor(f.length / 2); return f.length % 2 ? f[m] : (f[m - 1] + f[m]) / 2
}

/**
 * Per-class breakdown of one family interaction whose chains carry different FR3D classes: for each class, the chains (and their PDB entries)
 * in which it is annotated, out of the interaction's eligible chains. A chain with two annotated classes (e.g. alternate locations) counts in
 * both, so the class counts can sum to more than the observed count.
 */
export function classBreakdown(ds: InteractionsDataset, p: InteractionPair, pdbOf: (sn: string) => string): { cls: string; chains: number; entries: number }[] {
  const by = new Map<string, { chains: number; entries: Set<string> }>()
  for (const inst of ds.instances) {
    const s = statusIn(indexInstance(ds, inst.seq_name)!, p.category, p.pos_1, p.pos_2)
    if (s.state !== 'present') continue
    for (const c of s.classes) { const e = by.get(c) ?? by.set(c, { chains: 0, entries: new Set() }).get(c)!; e.chains++; e.entries.add(pdbOf(inst.seq_name)) }
  }
  return [...by.entries()].map(([cls, v]) => ({ cls, chains: v.chains, entries: v.entries.size })).sort((a, b) => b.chains - a.chains || a.cls.localeCompare(b.cls))
}

/** Construct rule: chains whose consensus-position sequences differ only at a few positions near their ends are the same RNA in
 *  different constructs (e.g. RF00059 thiM 2GDI vs 2HOJ: P1 ends). */
export const CONSTRUCT_MAX_END_DIFFS = 3, CONSTRUCT_END_WINDOW = 5
export type SequenceRelation = { kind: 'identical' | 'construct' | 'different'; shared: number; diffs: number[] }
/**
 * Sequence relation of two compared chains: 'identical' = the same chain sequence (sequence_key); 'construct' = a different chain sequence
 * but the same nucleotides on every shared consensus position, or differing only at <= CONSTRUCT_MAX_END_DIFFS positions all within
 * CONSTRUCT_END_WINDOW consensus positions of either chain's first or last matched position; 'different' otherwise.
 */
export function sequenceRelation(ds: InteractionsDataset, a: { seq_name: string; sequence_key: string }, b: { seq_name: string; sequence_key: string }): SequenceRelation {
  const na = ds.instances.find((i) => i.seq_name === a.seq_name)?.nucleotides ?? '', nb = ds.instances.find((i) => i.seq_name === b.seq_name)?.nucleotides ?? ''
  const shared: number[] = [], diffs: number[] = []
  for (let k = 0; k < Math.min(na.length, nb.length); k++) if (na[k] !== '-' && nb[k] !== '-') { shared.push(k + 1); if (na[k] !== nb[k]) diffs.push(k + 1) }
  if (a.sequence_key === b.sequence_key) return { kind: 'identical', shared: shared.length, diffs }
  const ends = (s: string) => { const m = [...s].flatMap((c, k) => (c !== '-' ? [k + 1] : [])); return m.length ? [m[0], m[m.length - 1]] : [] }
  const e = [...ends(na), ...ends(nb)]
  const nearEnd = (p: number) => e.some((x) => Math.abs(p - x) <= CONSTRUCT_END_WINDOW)
  const construct = shared.length > 0 && (diffs.length === 0 || (diffs.length <= CONSTRUCT_MAX_END_DIFFS && diffs.every(nearEnd)))
  return { kind: construct ? 'construct' : 'different', shared: shared.length, diffs }
}
