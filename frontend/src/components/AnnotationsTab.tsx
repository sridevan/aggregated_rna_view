import { useEffect, useLayoutEffect, useMemo, useRef, useState, type ReactNode } from 'react'
import { useFamily } from '../state/family'
import { useAnnotations } from '../state/useAnnotations'
import { familyFileUrl } from '../api/data'
import { SectionCard } from './SectionCard'
import { InfoTip } from './InfoTip'
import { NumberingSelect, useNumbering } from '../state/numbering'
import { useNavigation } from '../state/navigation'
import { chainAt, referenceTicks, residueValue } from '../lib/numbering'
import type { LigandsDataset, LoopTrackPosition, ModificationTracksDataset, PocketTracksDataset, SasaTracksDataset, TracksDataset } from '../types/backend'

/*
 * PDBe-KB-style residue annotation tracks on the Rfam consensus axis.
 * Geometry: every row is an HTML flex row with a fixed label column and an SVG stretched to the available width (preserveAspectRatio "none")
 * with a fixed pixel height. Zoom (px per consensus position, Fit all by default) and the horizontal scroll offset are one shared window:
 * every row's viewBox and the axis show the same visible span of positions, so labels stay fixed, categories keep their own vertical scroll,
 * and all rows stay aligned by construction. Marks are fixed size; the value of a position is encoded by colour through the track's `scale`.
 * Scales are per track (frequency today; quantitative scales such as median RSA later), so legends are generated from the scales in use.
 */
const CELL = 10, MARK_H = 16, ROW_PX = 30, LABEL_W = 176
const SEL = '#6250d6', UNEVAL = '#e4e2dc', GRID = '#ececea'
/** zoom levels (px per consensus position) and axis label spacing, as in the Family Overview tracks */
const ZOOM_PX = [1, 1.5, 2, 3, 4, 6, 8, 11, 16, 24, 32]
const TICK_MIN_PX = 60
const TICK_STEPS = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
/** px per position at "Fit all" for a plot width, and the axis label step for a scale (smallest step with labels >= TICK_MIN_PX apart) */
export const fitCellPx = (viewW: number, clen: number) => (viewW > 0 ? viewW / clen : 9)
export const tickStepFor = (cellPx: number) => TICK_STEPS.find((t) => t * cellPx >= TICK_MIN_PX) ?? TICK_STEPS[TICK_STEPS.length - 1]
/** consensus positions (1-based, inclusive) visible in the window [scrollX, scrollX + viewW) at cellPx */
export function visibleRange(scrollX: number, viewW: number, cellPx: number, clen: number): [number, number] {
  return [Math.max(1, Math.floor(scrollX / cellPx) + 1), Math.min(clen, Math.ceil((scrollX + viewW) / cellPx))]
}

/** Colour scale of a track. `frequency`: discretised light->dark bins of a fraction 0..1. `quantity`: continuous light->dark over [domainMin, domainMax] in `unit`.
 *  `presence`: one fixed colour (`dark`) wherever the track has a mark, whatever the value (for sparse data; the exact value is in the tooltip). */
export interface ScaleDef { id: string; kind: 'frequency' | 'quantity' | 'presence'; label: string; light: string; dark: string; bins: number; binLabels: string[]; unit?: string; domainMin?: number; domainMax?: number }
/** `summary`: rows of the pinned details panel (measured quantity, result with numerator / denominator / unit, absence wording, what was not
 *  counted and why). `details`: the longer evidence and provenance lines shown under "Evidence details". */
interface Mark { value: number | null; tip: { title: string; lines: string[] }; summary: [string, string][]; details: string[] }
interface TrackDef { id: string; label: string; sub: string; scale: ScaleDef; at: (pos: number) => Mark | 'unevaluable' | null }

const FREQ_BINS = ['<20%', '20–40%', '40–60%', '60–80%', '≥80%']
/** Shared occurrence-frequency scale (teal) of RNA loops, ligand binding sites and predicted pockets: one colour per frequency bin whatever the source;
 *  numerators/denominators stay source-specific (see each track). SASA keeps its own quantity scale; nonstandard nucleotides a presence colour. */
const FREQ_SCALE: ScaleDef = { id: 'occurrence-frequency', kind: 'frequency', label: 'Occurrence frequency (0–100%)', light: '#d2ecef', dark: '#0a4f5c', bins: 5, binLabels: FREQ_BINS }
const LOOP_SCALE = FREQ_SCALE, POCKET_SCALE = FREQ_SCALE, LIGAND_SCALE = FREQ_SCALE
const MOD_SCALE: ScaleDef = { id: 'nonstandard-presence', kind: 'presence', label: 'Nonstandard nucleotide observed', light: '#c2185b', dark: '#c2185b', bins: 1, binLabels: ['observed'] }

const hex = (c: string) => [1, 3, 5].map((i) => parseInt(c.slice(i, i + 2), 16))
const mix = (s: ScaleDef, t: number) => { const a = hex(s.light), b = hex(s.dark); return `rgb(${a.map((x, i) => Math.round(x + (b[i] - x) * t)).join(',')})` }
/** frequency: discretised light→dark colour for a value in [0, 1]; quantity: continuous colour over the scale domain (values outside are clamped) */
export function scaleColor(s: ScaleDef, v: number): string {
  if (s.kind === 'presence') return s.dark
  if (s.kind === 'quantity') {
    const lo = s.domainMin ?? 0, hi = s.domainMax ?? 1
    return mix(s, Math.min(1, Math.max(0, (v - lo) / Math.max(hi - lo, 1e-9))))
  }
  const bin = Math.min(s.bins - 1, Math.max(0, Math.floor(v * s.bins)))
  return mix(s, s.bins === 1 ? 1 : bin / (s.bins - 1))
}
/** "82.5%", "100%", "76.3%": one decimal at most, as used in the residue tooltips (display only; values are never recomputed from it) */
export const pct = (n: number, d: number) => `${d ? Math.round((1000 * n) / d) / 10 : 0}%`
const plural = (n: number, one: string, many: string) => `${n} ${n === 1 ? one : many}`
/** user-facing wording for the backend's pocket missing-data reasons (keys stay unchanged in pocket_tracks.json) */
const MISSING_REASON: Record<string, string> = { no_coordinates: 'no coordinates', failed_run: 'prediction failed', no_mapping: 'not mapped' }
const DETERMINABLE = 'Frequencies are calculated from instances where annotation presence or absence can be determined.'
const titleCase = (name: string | null) => (name ? (name === name.toUpperCase() ? name[0] + name.slice(1).toLowerCase() : name) : '')

/** Category-level explanations (concise by design; the full methodology lives in the downloadable JSON and the milestone documents). */
const POCKET_DEFAULT_FILTER = 'all_score_ge_0'   // fpocketR's published default (all pockets with fpocket Score >= 0); the other variants stay in the downloadable pocket_tracks.json
const CATEGORY_INFO = {
  pockets: 'RNA pockets predicted using fpocketR and mapped to Rfam positions.\n\nColours indicate pocket occurrence frequency. No prediction does not necessarily mean a binding pocket is absent.',
  sasa: 'Per-nucleotide SASA calculated for isolated experimental RNA chains and aggregated by Rfam position.\n\nColours indicate median SASA (Å²), reflecting nucleotide exposure.',
  modifications: 'Nonstandard nucleotide components in deposited RNA polymer sequences, mapped to Rfam residue positions.\n\nMarks indicate where nonstandard components have been observed. These do not necessarily represent naturally occurring RNA modifications.',
  loops: 'Hairpin (HL), internal (IL) and junction (JL) loop annotations from the BGSU RNA 3D Motif Atlas, mapped to Rfam residue positions. Loop classifications may overlap and vary due to evolutionary or conformational differences, experimental limitations, and methodological criteria.',
  ligands: 'Ligand-binding residues from PDBe annotations, mapped to Rfam residue positions and aggregated across PDB entries. Colours indicate contact frequency. Annotations are entity-level and may not represent contacts in every RNA chain.',
}

/** "4 RNA instances: no observed nucleotide here; 1 RNA instance: source unavailable" (only non-zero parts), or "None" */
const notCounted = (parts: [number, string][], unit: [string, string]) => parts.filter(([n]) => n > 0).map(([n, why]) => `${plural(n, unit[0], unit[1])}: ${why}`).join('; ') || 'None'
const INSTANCE: [string, string] = ['RNA instance', 'RNA instances'], ENTRY: [string, string] = ['PDB entry', 'PDB entries']

function loopTrack(id: 'HL' | 'IL' | 'JL', typeName: string, tracks: TracksDataset): TrackDef {
  const pick = (t: TracksDataset['positions'][number]): LoopTrackPosition | TracksDataset['positions'][number]['junction_loops_union'] => id === 'HL' ? t.hairpin_loops : id === 'IL' ? t.internal_loops : t.junction_loops_union
  return {
    id, label: `${id} · ${typeName.toLowerCase()}s`, sub: 'Motif Atlas', scale: LOOP_SCALE,
    at: (pos) => {
      const t = tracks.positions[pos - 1]; const p = pick(t); const d = p.denominators
      if (!d.n_observable_and_evaluated) return 'unevaluable'
      if (!p.n_instances) return null
      const k = d.n_positive_observable_and_evaluated, n = d.n_observable_and_evaluated
      const details = [`${p.n_loop_observations} loop observations on ${plural(p.n_instances, 'RNA instance', 'RNA instances')} in ${plural(p.n_entries, 'PDB entry', 'PDB entries')}`,
        `Not annotated as ${id} by the source in ${d.n_not_annotated_by_source} of ${n} RNA instances; source unavailable for ${d.n_source_unavailable}${d.n_positive_not_observable ? `; ${d.n_positive_not_observable} annotated instances without coordinates here are not counted` : ''}`]
      if (id === 'JL') { const bt = t.junction_loops_union.by_type_n_instances ?? { J3: t.junction_loops_J3.n_instances, J4: t.junction_loops_J4.n_instances }
        details.push(`${Object.entries(bt).map(([k, v]) => `${k}: ${v} instances`).join(', ')} (union, never summed)`) }
      if ('by_motif_group' in p) details.push(`Motif groups: ${Object.entries(p.by_motif_group).map(([g, v]) => `${g} (${v.n_loops})`).join(', ') || 'none'}`)
      details.push(`${DETERMINABLE} BGSU loop classifications can overlap.`)
      // RNA instances that are not in n: no observed nucleotide at this position (deleted, outside the aligned region or no coordinates), or
      // an observed nucleotide whose instance has no record from the source; both counts come from the unchanged backend denominators
      const noObs = tracks.n_instances - d.n_observable_instances, unavailable = d.n_observable_instances - n
      const summary: [string, string][] = [
        ['Annotation', `${typeName} (${id}), BGSU RNA 3D Motif Atlas`],
        ['Observed', `Observed in ${k} of ${n} RNA instances (${pct(k, n)})`],
        ['Not annotated', `${plural(n - k, 'RNA instance', 'RNA instances')} (not annotated by the source; not a confirmed absence)`],
        ['Not counted', notCounted([[noObs, 'no observed nucleotide at this position'], [unavailable, 'source unavailable']], INSTANCE)],
        ['Counting unit', 'RNA chain instance'],
      ]
      return { value: p.frequency ?? 0, tip: { title: `Rfam residue ${pos}`, lines: [`${typeName} (${id})`, `Observed in ${k} of ${n} RNA instances (${pct(k, n)})`] }, summary, details }
    },
  }
}

/** Solvent accessibility: one quantitative track (median SASA, Å²). Colour domain = 0 .. max median rounded up to 50 Å² (from the dataset, documented in sasa_tracks.json). */
function sasaTrack(sasa: SasaTracksDataset): TrackDef {
  const measured = sasa.summary.n_instances_measured ?? 0
  const medians = sasa.positions.map((p) => p.median_sasa_A2).filter((v): v is number => v != null)
  const max = sasa.colour_domain_suggestion_A2?.max ?? (medians.length ? Math.ceil(Math.max(...medians) / 50) * 50 : 1)
  const scale: ScaleDef = { id: 'sasa-median', kind: 'quantity', label: 'Median SASA (Å²)', light: '#fdf0d5', dark: '#8a2a0b', bins: 5, binLabels: [], unit: 'Å²', domainMin: 0, domainMax: max }
  return {
    id: 'sasa-median', label: 'Median SASA (Å²)', sub: `${sasa.summary.n_instances_measured ?? ''} RNA instances`, scale,
    at: (pos) => {
      const p = sasa.positions[pos - 1]
      if (!p.n_mapped_instances_with_coordinates) return 'unevaluable'          // no observed nucleotide at this position in any instance
      if (p.median_sasa_A2 == null) return 'unevaluable'                        // observed but no complete nucleotide qualified
      const n = p.n_contributing_instances
      const summary: [string, string][] = [
        ['Measurement', 'Solvent-accessible surface area of the nucleotide in the isolated RNA chain (Å²)'],
        ['Median', `${p.median_sasa_A2.toFixed(1)} Å² (interquartile range ${p.q1_sasa_A2?.toFixed(1)}–${p.q3_sasa_A2?.toFixed(1)} Å²) over ${plural(n, 'RNA instance', 'RNA instances')}`],
        ['Not counted', notCounted([[measured - p.n_mapped_instances_with_coordinates, 'no observed nucleotide at this position'], [p.n_excluded_incomplete, 'incomplete nucleotide'], [p.n_excluded_unsupported, 'unsupported nucleotide']], INSTANCE)],
        ['Counting unit', 'RNA chain instance'],
      ]
      return { value: p.median_sasa_A2, tip: { title: `Rfam residue ${pos}`, lines: [`Median SASA: ${p.median_sasa_A2.toFixed(1)} Å²`, `Based on ${n} RNA ${n === 1 ? 'instance' : 'instances'}`] }, summary,
        details: [`Mean ${p.mean_sasa_A2?.toFixed(1)} Å², min ${p.min_sasa_A2?.toFixed(1)}, Q1 ${p.q1_sasa_A2?.toFixed(1)}, Q3 ${p.q3_sasa_A2?.toFixed(1)}, max ${p.max_sasa_A2?.toFixed(1)}${p.sd_sasa_A2 != null ? `, sd ${p.sd_sasa_A2.toFixed(1)}` : ''}`,
                  `${p.n_mapped_instances_with_coordinates} instances with an observed nucleotide here; excluded: ${p.n_excluded_incomplete} incomplete, ${p.n_excluded_unsupported} unsupported`,
                  `Nucleotide identity among contributing instances: ${Object.entries(p.nucleotide_identity).map(([k, v]) => `${k} ${v}`).join(', ')}`] }
    },
  }
}

/** Predicted binding pockets: one frequency track (default filter). Grey = pocket presence/absence cannot be determined in any instance at the position; empty = determined, in no pocket. */
function pocketTrack(pt: PocketTracksDataset): TrackDef | null {
  const positions = pt.tracks[POCKET_DEFAULT_FILTER]
  if (!positions) return null
  return {
    id: 'pocket-occurrence', label: 'Pocket occurrence', sub: `fpocketR · ${pt.n_instances} RNA instances`, scale: POCKET_SCALE,
    at: (pos) => {
      const p = positions[pos - 1]
      if (!p.n_evaluable_instances) return 'unevaluable'
      if (!p.n_supporting_instances) return null
      const n = p.n_supporting_instances, d = p.n_evaluable_instances
      const reasons = Object.entries(p.unevaluable_reasons ?? {}).map(([k, v]) => [v, MISSING_REASON[k] ?? k.replace(/_/g, ' ')] as [number, string])
      const summary: [string, string][] = [
        ['Prediction', `Predicted binding pocket (fpocketR, ${pt.filters[POCKET_DEFAULT_FILTER] ?? POCKET_DEFAULT_FILTER})`],
        ['Observed', `In a predicted pocket in ${n} of ${d} RNA instances (${pct(n, d)})`],
        ['No pocket', `${plural(d - n, 'RNA instance', 'RNA instances')} (prediction ran, no pocket here; not a confirmed absence of binding)`],
        ['Not counted', notCounted(reasons.length ? reasons : [[p.n_unevaluable_instances ?? 0, 'prediction not available']], INSTANCE)],
        ['Counting unit', 'RNA chain instance'],
      ]
      return { value: p.frequency ?? n / d, tip: { title: `Rfam residue ${pos}`, lines: [`Pocket occurrence: ${n} of ${d} RNA instances (${pct(n, d)})`] }, summary,
        details: [`${p.n_supporting_entries} PDB entries; ${p.n_distinct_pockets} distinct predicted pockets${p.score_summary ? ` (fpocket Score ${p.score_summary.min.toFixed(2)}–${p.score_summary.max.toFixed(2)}, median ${p.score_summary.median.toFixed(2)})` : ''}`,
                  ...(p.n_unevaluable_instances ? [`Not counted: ${plural(p.n_unevaluable_instances, 'RNA instance', 'RNA instances')}${p.unevaluable_reasons && Object.keys(p.unevaluable_reasons).length ? ` (${Object.entries(p.unevaluable_reasons).map(([k, v]) => `${MISSING_REASON[k] ?? k.replace(/_/g, ' ')}: ${v}`).join(', ')})` : ''}`] : []),
                  `Pocket filter: ${pt.filters[POCKET_DEFAULT_FILTER] ?? POCKET_DEFAULT_FILTER}`,
                  `${DETERMINABLE} Successful predictions without a pocket here count as absent; failed predictions and missing residues are not counted.`] }
    },
  }
}

/** Nonstandard nucleotides: one presence row per CCD component that occurs on the consensus axis (like the per-ligand rows). Fixed colour wherever an
 *  instance carries that component at the position (sparse data); grey = no instance with a determined component there; empty = none observed.
 *  Frequency = instances with this component at the position / instances with a determined component there. Exact counts in the tooltip. */
const PARENT_SRC: Record<string, string> = { CCD: 'CCD', NAKB: 'NAKB' }
function modificationTracks(md: ModificationTracksDataset): TrackDef[] {
  const comps = Object.entries(md.components).filter(([, c]) => c.n_on_consensus_axis > 0)
    .sort((a, b) => b[1].n_instances - a[1].n_instances || a[0].localeCompare(b[0]))
  return comps.map(([ccd, ci]) => {
    const offAxis = ci.n_observations - ci.n_on_consensus_axis
    return {
      id: `mod-${ccd}`, label: ccd, sub: `${titleCase(ci.chem_comp_name)} · ${plural(ci.n_instances, 'instance', 'instances')}`, scale: MOD_SCALE,
      at: (pos) => {
        const p = md.positions[pos - 1]
        if (!p || !p.denominator) return 'unevaluable'
        const c = p.ccd.find((x) => x.ccd_id === ccd)
        if (!c) return null
        const k = c.n_instances, n = p.denominator
        const par = ci.parent ? `parent ${ci.parent}${ci.parent_source ? ` (${PARENT_SRC[ci.parent_source]})` : ''}` : 'parent unknown'
        const details = [`${ccd} — ${ci.chem_comp_name ?? 'name not available'}; ${par}`,
          `Observed in ${k} of ${n} RNA instances (${pct(k, n)}) in ${plural(c.n_entries, 'PDB entry', 'PDB entries')}`,
          ...(p.ccd.length > 1 ? [`All nonstandard components at this position: ${p.n_positive} of ${n} RNA instances (${p.ccd.map((x) => `${x.ccd_id} ${x.n_instances}`).join(', ')})`] : []),
          ...(p.n_denominator_without_coordinates ? [`${p.n_denominator_without_coordinates} of the ${n} RNA instances have no coordinates for this nucleotide; its identity comes from the deposited sequence`] : []),
          ...(p.n_missing_identity ? [`Not counted: ${plural(p.n_missing_identity, 'RNA instance', 'RNA instances')} with undetermined component identity`] : []),
          `${DETERMINABLE} Instances without a residue at this position are not counted.${offAxis ? ` ${plural(offAxis, `further ${ccd} residue lies`, `further ${ccd} residues lie`)} outside the mapped Rfam region (downloadable data only).` : ''}`]
        const others = p.n_positive - k
        const summary: [string, string][] = [
          ['Component', `${ccd} — ${titleCase(ci.chem_comp_name) || 'name not available'} (deposited polymer sequence; ${par})`],
          ['Observed', `Observed in ${k} of ${n} RNA instances (${pct(k, n)})`],
          ['Other identities', [`${plural(p.n_negative, 'RNA instance', 'RNA instances')} with a standard nucleotide`, ...(others > 0 ? [`${others} with another nonstandard component`] : [])].join('; ')],
          ['Not counted', notCounted([[p.n_no_residue, 'no residue at this position (deleted or outside the aligned region)'], [p.n_missing_identity, 'component identity undetermined'], [p.n_source_unavailable, 'source unavailable']], INSTANCE)],
          ['Counting unit', 'RNA chain instance'],
        ]
        return { value: k / n, tip: { title: `Rfam residue ${pos}`, lines: [`${ccd} — ${titleCase(ci.chem_comp_name)}`, `Observed in ${k} of ${n} RNA instances (${pct(k, n)})`] }, summary, details }
      },
    }
  })
}

function buildTracks(tracks: TracksDataset, ligands: LigandsDataset) {
  const structure: TrackDef[] = [loopTrack('HL', 'Hairpin loop', tracks), loopTrack('IL', 'Internal loop', tracks), loopTrack('JL', 'Junction loop', tracks)]
  const perCcdPos = new Map<string, Map<number, Set<string>>>()
  for (const bm of ligands.bound_molecules) {
    if (!bm.in_default_track) continue
    const m = perCcdPos.get(bm.ccd_id) ?? new Map<number, Set<string>>(); perCcdPos.set(bm.ccd_id, m)
    for (const pos of bm.consensus_positions) (m.get(pos) ?? m.set(pos, new Set()).get(pos)!).add(bm.pdb_id)
  }
  const ccds = Object.entries(ligands.by_ccd).filter(([, v]) => v.n_entries_in_default_track > 0).sort((a, b) => b[1].n_entries_in_default_track - a[1].n_entries_in_default_track || a[0].localeCompare(b[0]))
  const ligandTracks: TrackDef[] = ccds.map(([ccd, v]) => ({
    id: `lig-${ccd}`, label: ccd, sub: `${titleCase(v.name)} · ${v.n_entries_in_default_track} ${v.n_entries_in_default_track === 1 ? 'entry' : 'entries'}`, scale: LIGAND_SCALE,
    at: (pos) => {
      const lp = tracks.positions[pos - 1].ligand_binding_sites
      if (!lp.denominators.n_observable_and_evaluated) return 'unevaluable'
      const entries = perCcdPos.get(ccd)?.get(pos)
      const e = lp.by_ccd_entries[ccd]
      if (!entries?.size || !e) return null
      // counted per PDB entry among the entries with bound ccd that have an observed nucleotide here (backend by_ccd_entries); entries without one are not negatives
      if (!e.n_observable_and_evaluated) return 'unevaluable'
      const k = e.n_positive_observable_and_evaluated, n = e.n_observable_and_evaluated
      const summary: [string, string][] = [
        ['Annotation', `${ccd} — ${titleCase(v.name)} in contact (PDBe binding-site annotation, entity level)`],
        ['Recorded contact', `${k} of ${plural(n, 'observable PDB entry', 'observable PDB entries')} with bound ${ccd} (${pct(k, n)})`],
        ['No contact', `${plural(n - k, 'PDB entry', 'PDB entries')} with bound ${ccd} and an observed nucleotide: no contact at this position recorded by PDBe (not a confirmed absence)`],
        ['Not counted', notCounted([[e.n_not_observable, 'no observed nucleotide at this position'], [e.n_source_unavailable, 'ligand annotation unavailable for the chains observed here']], ENTRY)],
        ['Counting unit', 'PDB entry (chains of one entry counted once)'],
      ]
      return { summary, value: k / n,
        tip: { title: `Rfam residue ${pos}`, lines: [`${ccd} — ${titleCase(v.name)}`, `Recorded contact in ${k} of ${plural(n, 'observable PDB entry', 'observable PDB entries')} (${pct(k, n)})`] },
        details: [`Supporting entries: ${[...entries].sort().map((x) => x.toUpperCase()).join(', ')}`,
                  `Counted per PDB entry: entries with ${ccd} in contact at this position, of the ${n} (of ${v.n_entries_in_default_track}) PDB entries containing bound ${ccd} that have an observed nucleotide here. An observed-contact frequency, not a binding probability`,
                  ...(e.n_positive > k ? [`${plural(e.n_positive - k, 'PDB entry', 'PDB entries')} with a recorded contact but no observed nucleotide here: not counted`] : []),
                  `Ligand-contact presence or absence can be determined at this position for ${lp.denominators.n_observable_and_evaluated} PDB entry entities (all ligands)`] }
    },
  }))
  return { structure, ligandTracks }
}

/** Category order (user decision): Solvent accessibility, RNA loop annotations, Nonstandard nucleotides, Ligand binding sites, Predicted binding pockets.
 *  Compact PDBe-KB-style annotation viewer: fixed-size marks coloured by value, collapsible and independently scrollable categories, concise tooltips, click-to-pin details, download of the full evidence. */
export function AnnotationsTab() {
  const { acc, family } = useFamily()
  const st = useAnnotations(acc)
  const [inspect, setInspect] = useState<{ trackId: string; label: string; pos: number; title: string; summary: [string, string][]; details: string[] } | null>(null)
  const nav = useNavigation()
  const [tip, setTip] = useState<{ x: number; y: number; title: string; lines: string[] } | null>(null)
  const clen = family.clen
  const x = (pos: number) => (pos - 1) * CELL
  // Position numbering (display only, shared with the Family Overview): labels and tooltip lines; marks stay on Rfam positions
  const { numbering } = useNumbering()
  const refLine = (pos: number) => (numbering ? [`${numbering.label} reference residue: ${residueValue(chainAt(numbering.idx, pos))}`] : [])
  // --- shared zoom + horizontal window (measured on the axis plot area; every row's plot area has the same width) ---
  const plotRef = useRef<HTMLDivElement>(null), barRef = useRef<HTMLDivElement>(null)
  const [viewW, setViewW] = useState(0)
  const [zoomPx, setZoomPx] = useState<number | null>(null)   // px per position; null = Fit all
  const [scrollRaw, setScrollX] = useState(0)                 // px from the start of the consensus axis
  const ready = st.state === 'ready'
  useLayoutEffect(() => {
    const el = plotRef.current; if (!el) return
    setViewW(el.clientWidth)
    if (typeof ResizeObserver === 'undefined') return
    const ro = new ResizeObserver(() => setViewW(el.clientWidth)); ro.observe(el); return () => ro.disconnect()
  }, [ready])
  const fitPx = fitCellPx(viewW, clen)
  const cellPx = zoomPx == null ? fitPx : Math.max(zoomPx, fitPx)
  const totalW = clen * cellPx
  const scrollX = Math.min(Math.max(0, scrollRaw), Math.max(0, totalW - viewW))
  const [lo, hi] = visibleRange(scrollX, viewW || totalW, cellPx, clen)
  const vb = { x0: (scrollX / cellPx) * CELL, w: ((viewW || totalW) / cellPx) * CELL }   // visible window in track units
  const pad = Math.min(1, CELL / cellPx)                      // 1 px gap each side of a mark (10% of a cell at overview scale)
  const zoomIn = ZOOM_PX.find((z) => z > cellPx + 0.01)
  const zoomOutTo = [...ZOOM_PX].reverse().find((z) => z < cellPx - 0.01)
  const atFit = zoomPx == null || cellPx <= fitPx + 0.01
  /** zoom keeping the inspected position (else the window centre) at the same place */
  const setZoom = (z: number | null) => {
    const nextPx = z == null ? fitPx : Math.max(z, fitPx)
    const ip = inspect ? (inspect.pos - 0.5) * cellPx - scrollX : -1   // an inspected position in view stays put; one out of view is centred
    const keep = inspect ? { pos: inspect.pos - 0.5, px: ip >= 0 && ip <= viewW ? ip : viewW / 2 } : { pos: (scrollX + viewW / 2) / cellPx, px: viewW / 2 }
    setZoomPx(z); setScrollX(Math.max(0, keep.pos * nextPx - keep.px))
  }
  useLayoutEffect(() => { if (barRef.current && Math.abs(barRef.current.scrollLeft - scrollX) > 0.5) barRef.current.scrollLeft = scrollX }, [scrollX, totalW])
  /** horizontal trackpad / shift+wheel over the tracks pans the shared window (vertical wheel keeps scrolling the categories) */
  const onWheel = (e: React.WheelEvent) => {
    const dx = Math.abs(e.deltaX) > Math.abs(e.deltaY) ? e.deltaX : e.shiftKey ? e.deltaY : 0
    if (dx && totalW > viewW) setScrollX(Math.min(Math.max(0, scrollX + dx), totalW - viewW))
  }
  const built = useMemo(() => (st.state === 'ready' ? buildTracks(st.tracks, st.ligands) : null), [st])
  const sasa = useMemo(() => (st.state === 'ready' && st.sasa ? sasaTrack(st.sasa) : null), [st])
  const pockets = useMemo(() => (st.state === 'ready' && st.pockets ? pocketTrack(st.pockets) : null), [st])
  const mods = useMemo(() => (st.state === 'ready' && st.modifications ? modificationTracks(st.modifications) : null), [st])
  const scales = useMemo(() => { const m = new Map<string, ScaleDef>(); for (const t of [...(sasa ? [sasa] : []), ...(built?.structure ?? []), ...(mods ?? []), ...(built?.ligandTracks ?? []), ...(pockets ? [pockets] : [])]) m.set(t.scale.id, t.scale); return [...m.values()] }, [built, sasa, pockets, mods])
  const inspectedPos = inspect?.pos ?? null
  if (typeof window !== 'undefined') window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), annotationTracks: built ? { structure: built.structure.map((t) => t.id), ligands: built.ligandTracks.map((t) => t.id), pockets: pockets ? [pockets.id] : [], sasa: sasa ? [sasa.id] : [], modifications: mods ? mods.map((t) => t.id) : [] } : null }

  const marksSvg = (t: TrackDef) => (
    <svg viewBox={`${vb.x0} 0 ${vb.w} ${ROW_PX}`} preserveAspectRatio="none" className="block w-full select-none" style={{ height: ROW_PX }} role="img" aria-label={t.label}>
      {inspectedPos != null && <rect x={x(inspectedPos)} y={0} width={CELL} height={ROW_PX} fill={SEL} opacity={0.16} />}
      <line x1={0} x2={clen * CELL} y1={ROW_PX - 4} y2={ROW_PX - 4} stroke={GRID} vectorEffect="non-scaling-stroke" />
      {Array.from({ length: hi - lo + 1 }, (_, k) => lo + k).map((pos) => {
        const m = t.at(pos)
        if (!m) return null
        const y = (ROW_PX - MARK_H) / 2
        if (m === 'unevaluable') return <rect key={pos} x={x(pos) + pad} y={y} width={CELL - 2 * pad} height={MARK_H} fill={UNEVAL} data-pos={pos} data-state="unevaluable" />
        return <rect key={pos} x={x(pos) + pad} y={y} width={CELL - 2 * pad} height={MARK_H} fill={scaleColor(t.scale, m.value ?? 0)} data-pos={pos} data-value={(m.value ?? 0).toFixed(3)} className="cursor-pointer"
          onMouseEnter={(e) => setTip({ x: e.clientX, y: e.clientY, title: m.tip.title, lines: [...m.tip.lines, ...refLine(pos)] })} onMouseMove={(e) => setTip((s) => (s ? { ...s, x: e.clientX, y: e.clientY } : s))} onMouseLeave={() => setTip(null)}
          onClick={() => setInspect((cur) => (cur && cur.trackId === t.id && cur.pos === pos ? null : { trackId: t.id, label: t.label, pos, title: `${t.label} · ${m.tip.title}`, summary: m.summary, details: m.details }))} />
      })}
    </svg>
  )
  const row = (t: TrackDef) => (
    <div key={t.id} className="flex items-center" style={{ height: ROW_PX }} data-testid={`track-${t.id}`}>
      <div className="min-w-0 shrink-0 pr-2 leading-tight" style={{ width: LABEL_W }} title={`${t.label} — ${t.sub}`}>
        <div className="truncate text-[11px] text-stone-900">{t.label}</div><div className="truncate text-[9px] text-stone-500">{t.sub}</div>
      </div>
      <div className="min-w-0 flex-1">{marksSvg(t)}</div>
    </div>
  )
  const tickStep = tickStepFor(cellPx)
  const left = (pos: number) => (pos - 0.5) * cellPx - scrollX   // px of a position's centre in the plot area
  const axis = (   // px-[5px] = category border (1 px) + body padding (4 px) so the axis shares the rows' horizontal extent exactly
    <div className="px-[5px]">
      <div className="flex items-center" data-testid="annotation-axis" data-cell-px={cellPx.toFixed(2)} data-first-pos={lo} data-last-pos={hi}>
        {numbering
          ? <div className="shrink-0 pr-1 text-[10px] text-violet-800" style={{ width: LABEL_W }} data-testid="annotation-axis-convention">{numbering.label} · author numbering</div>
          : <div className="shrink-0 pr-1 text-[10px] text-stone-500" style={{ width: LABEL_W }} data-testid="annotation-axis-convention">Rfam {acc} consensus position</div>}
        <div ref={plotRef} className="relative h-6 min-w-0 flex-1 overflow-hidden text-[9px] text-stone-600" data-testid="annotation-axis-plot">
          {inspectedPos != null && <div className="absolute inset-y-0" style={{ left: left(inspectedPos) - cellPx / 2, width: cellPx, background: SEL, opacity: 0.16 }} />}
          <div className="absolute inset-x-0 bottom-1 border-t border-stone-500" />
          {/* label step from the zoom level (labels at least TICK_MIN_PX apart; labels cut by the window edge are dropped): 10 for riboswitch-sized models at Fit all, 100 for 16S */}
          {/* reference numbering: the reference chain's author numbers at the Rfam columns (only where it has a mapped residue; never invented) */}
          {(numbering ? referenceTicks(numbering.numbers, cellPx) : Array.from({ length: Math.floor(clen / tickStep) + 1 }, (_, k) => k * tickStep || 1).filter((p) => p <= clen).map((p) => ({ pos: p, text: String(p) })))
            .filter((t) => { const half = t.text.length * 2.6; return left(t.pos) >= half && left(t.pos) <= (viewW || totalW) - half }).map((t) => (
            <div key={t.pos} className={`absolute bottom-0 flex -translate-x-1/2 flex-col items-center ${numbering ? 'text-violet-800' : ''}`} style={{ left: left(t.pos) }} data-ref-tick={numbering ? t.pos : undefined}><span>{t.text}</span><span className={`h-1.5 w-px ${numbering ? 'bg-violet-700' : 'bg-stone-500'}`} /></div>))}
        </div>
      </div>
      {/* one horizontal scrollbar for all categories (only when zoomed in beyond the plot width) */}
      <div className="flex">
        <div className="shrink-0" style={{ width: LABEL_W }} />
        <div ref={barRef} className="annotations-hscroll min-w-0 flex-1 overflow-x-auto overflow-y-hidden" style={{ height: totalW > viewW + 0.5 ? 10 : 0 }} onScroll={(e) => setScrollX(e.currentTarget.scrollLeft)} data-testid="annotations-hscroll">
          <div style={{ width: totalW, height: 1 }} />
        </div>
      </div>
    </div>
  )
  const zoomControls = (
    <div className="flex flex-wrap items-center gap-1 text-[11px]">
    <NumberingSelect />
    <div className="ml-auto flex items-center gap-1" role="group" aria-label="Zoom the Rfam position axis" data-testid="annotations-zoom-controls" data-cell-px={cellPx.toFixed(2)} data-fit={String(atFit)}>
      <span className="mr-1 text-[10px] text-stone-500">Zoom</span>
      <button type="button" onClick={() => setZoom(zoomOutTo != null && zoomOutTo > fitPx + 0.01 ? zoomOutTo : null)} disabled={atFit} aria-label="Zoom out" title="Zoom out" data-testid="annotations-zoom-out"
        className="h-5 w-5 rounded border border-stone-300 bg-white leading-none text-stone-700 hover:bg-stone-50 disabled:opacity-40">−</button>
      <button type="button" onClick={() => zoomIn != null && setZoom(zoomIn)} disabled={zoomIn == null} aria-label="Zoom in" title="Zoom in" data-testid="annotations-zoom-in"
        className="h-5 w-5 rounded border border-stone-300 bg-white leading-none text-stone-700 hover:bg-stone-50 disabled:opacity-40">+</button>
      <button type="button" onClick={() => setZoom(null)} disabled={atFit} title="Fit the whole consensus model in view" data-testid="annotations-zoom-fit"
        className="rounded border border-stone-300 bg-white px-1.5 py-px text-stone-700 hover:bg-stone-50 disabled:opacity-40">Fit all</button>
    </div>
    </div>
  )
  /** One compact legend: the shared occurrence-frequency scale (loops, ligands, pockets), the SASA quantity scale, the nonstandard-nucleotide
   *  presence colour, and the two non-value states. Built from the scales in use, so a family without e.g. SASA shows no SASA entry. */
  const KIND_ORDER = { frequency: 0, quantity: 1, presence: 2 } as const
  const SW = 'inline-block h-2.5 w-2.5 rounded-[2px] ring-1 ring-inset ring-black/10'
  const legend = (
    <div className="flex flex-wrap items-center gap-x-4 gap-y-0.5 text-[10px] leading-4 text-stone-600" data-testid="scale-legend">
      {[...scales].sort((x, y) => KIND_ORDER[x.kind] - KIND_ORDER[y.kind]).map((s) => s.kind === 'quantity' ? (
        <div key={s.id} className="flex items-center gap-1" data-testid={`legend-${s.id}`}>
          <span>{s.label}:</span><span>{s.domainMin} {s.unit}</span>
          <span className="inline-block h-2.5 w-20 rounded-[2px]" style={{ background: `linear-gradient(to right, ${s.light}, ${s.dark})` }} data-testid={`legend-gradient-${s.id}`} />
          <span>≥ {s.domainMax} {s.unit}</span>
        </div>
      ) : s.kind === 'presence' ? (
        <div key={s.id} className="flex items-center gap-1" data-testid={`legend-${s.id}`}>
          <span className={SW} style={{ background: s.dark }} /><span>{s.label}</span>
        </div>
      ) : (
        <div key={s.id} className="flex items-center gap-1" data-testid={`legend-${s.id}`}>
          <span>{s.label}:</span>
          {Array.from({ length: s.bins }, (_, i) => <span key={i} className="inline-flex items-center gap-0.5"><span className={SW} style={{ background: scaleColor(s, (i + 0.5) / s.bins) }} /><span>{s.binLabels[i]}</span></span>)}
        </div>
      ))}
      <div className="flex items-center gap-1" data-testid="legend-no-data" title="No RNA instance can be assessed at this position: no observed nucleotide, or the source data are unavailable (SASA: no qualifying measurement)"><span className={SW} style={{ background: UNEVAL }} /> No data</div>
      <div className="flex items-center gap-1" data-testid="legend-no-annotation" title="Assessed, but nothing recorded at this position: no loop annotated by the source, no PDBe ligand contact, no predicted pocket, or a standard nucleotide. For loops, ligands and pockets this is not a confirmed absence."><span className={`${SW} bg-white ring-stone-300`} /> Not observed</div>
    </div>
  )
  // labels without "(JSON)" (the button says it); each description sits behind an info icon (checked against the files' contents)
  const datasets: DownloadItem[] = [
    { id: 'loops', label: 'RNA loop annotations', file: 'loops.json', detail: 'BGSU RNA 3D Motif Atlas loops (HL, IL, J3, J4) with motif-group assignments, nucleotide and segment membership and Rfam consensus mapping.' },
    { id: 'ligands', label: 'Ligand binding annotations', file: 'ligands.json', detail: 'Bound molecules (CCD id, name, class), binding-site nucleotides, supporting entries and chains, and provenance.' },
    { id: 'tracks', label: 'Aggregated annotation tracks', file: 'tracks.json', detail: 'Per-position annotation frequencies with numerators, denominators and missing-data counts.' },
    { id: 'pockets', label: 'Predicted pockets · Per structure', file: 'pockets.json', detail: 'fpocketR pockets per RNA chain with scores, volumes, descriptors and lining nucleotides mapped to Rfam positions.' },
    { id: 'pocket-tracks', label: 'Predicted pockets · Aggregated', file: 'pocket_tracks.json', detail: 'Per-position pocket frequencies for each of the three pocket filters.' },
    { id: 'sasa', label: 'Solvent accessibility · Per nucleotide', file: 'sasa.json', detail: 'Per-chain nucleotide solvent accessibility (Å²) with atom completeness and calculation details.' },
    { id: 'sasa-tracks', label: 'Solvent accessibility · Aggregated', file: 'sasa_tracks.json', detail: 'Per-position solvent accessibility summaries (median, quartiles) and contributing-chain counts.' },
    { id: 'modifications', label: 'Nonstandard nucleotides · Detailed', file: 'modifications.json', detail: 'Nonstandard nucleotide identities, locations, parent nucleotides, Rfam mappings and source information.' },
    { id: 'modification-tracks', label: 'Nonstandard nucleotides · Aggregated', file: 'modification_tracks.json', detail: 'Per-position nonstandard-nucleotide frequencies, counts and contributing components.' },
  ]
  return (
    <div className="space-y-3" data-testid="annotations-tab">
      <SectionCard title="Annotations" actions={st.state === 'ready' ? <DownloadMenu acc={acc} items={datasets} /> : undefined}
        description={<span data-testid="annotations-summary">Residue-level annotations across {family.counts.instances} RNA instances from {family.counts.pdb_entries} PDB entries.</span>}>
        {st.state === 'loading' && <p className="text-sm text-stone-500" data-testid="annotations-loading">Loading annotation tracks…</p>}
        {st.state === 'error' && <p className="text-sm text-red-700" data-testid="annotations-error">Annotation data not available: {st.message}</p>}
        {st.state === 'ready' && built && (
          <div className="space-y-2">
            {legend}
            {inspect && (
              <div className="rounded border border-violet-200 bg-violet-50 px-2 py-1.5 text-xs text-stone-800" data-testid="annotation-inspect">
                <div className="flex items-start justify-between gap-2"><p className="font-semibold">{inspect.title}</p><button type="button" className="rounded border border-stone-300 bg-white px-1.5 text-xs" onClick={() => setInspect(null)} data-testid="annotations-clear">Close</button></div>
                {/* measured quantity, numerator / denominator / unit, absence wording and what was not counted (all from the unchanged backend fields) */}
                <dl className="mt-0.5 grid grid-cols-[auto_1fr] gap-x-3 gap-y-0.5 leading-4" data-testid="annotation-summary">
                  {[...inspect.summary, ...(numbering ? [['Reference residue', `${numbering.label} reference residue: ${residueValue(chainAt(numbering.idx, inspect.pos))}`] as [string, string]] : [])].map(([k, v]) => (
                    <div key={k} className="contents" data-testid={`annotation-summary-${k.toLowerCase().replace(/[^a-z]+/g, '-')}`}><dt className="text-stone-500">{k}</dt><dd>{v}</dd></div>))}
                </dl>
                {nav && (
                  <div className="mt-1 flex flex-wrap items-center gap-x-2 gap-y-0.5">
                    <button type="button" className="rounded border border-violet-300 bg-white px-1.5 py-px text-[11px] text-violet-800 hover:bg-violet-100" onClick={() => nav.viewInOverview(inspect.pos, `${inspect.label} at Rfam position ${inspect.pos}`)} data-testid="view-in-overview"
                      title="Select this Rfam position for the RNA instance shown in the Family Overview">View in Family Overview</button>
                    <span className="text-[10px] text-stone-500">Highlights Rfam position {inspect.pos} in the selected RNA instance; the annotation may not be present in that chain.</span>
                  </div>)}
                <details className="mt-1"><summary className="cursor-pointer text-[11px] text-stone-600">Evidence details</summary>{inspect.details.map((l, i) => <p key={i} className="mt-0.5 text-[11px] leading-4 text-stone-700">{l}</p>)}</details>
              </div>
            )}
            {zoomControls}
            {axis}
            <div className="space-y-2 overscroll-x-contain" onWheel={onWheel} data-testid="annotation-categories">
            {sasa && <Category id="sasa" title="Solvent accessibility" count={1} info={CATEGORY_INFO.sasa}>{row(sasa)}</Category>}
            <Category id="loops" title="RNA loop annotations" count={built.structure.length} info={CATEGORY_INFO.loops}>{built.structure.map(row)}</Category>
            {mods && mods.length > 0 && <Category id="modifications" title="Nonstandard nucleotides" count={mods.length} info={CATEGORY_INFO.modifications}>{mods.map(row)}</Category>}
            <Category id="ligands" title="Ligand binding sites" count={built.ligandTracks.length} info={CATEGORY_INFO.ligands}>{built.ligandTracks.map(row)}</Category>
            {pockets && <Category id="pockets" title="Predicted binding pockets" count={1} info={CATEGORY_INFO.pockets}>{row(pockets)}</Category>}
            </div>
            {/* categories without data for this family are named, never silently dropped */}
            {(() => { const missing = [!sasa && 'Solvent accessibility', !mods && 'Nonstandard nucleotides', !pockets && 'Predicted binding pockets'].filter(Boolean)
              return <>
                {missing.length > 0 && <p className="text-[11px] text-stone-500" data-testid="categories-unavailable">Not available for {acc}: {missing.join(', ')}.</p>}
                {mods && mods.length === 0 && <p className="text-[11px] text-stone-500" data-testid="modifications-none">Nonstandard nucleotides: none observed at Rfam positions.</p>}
              </> })()}
          </div>
        )}
      </SectionCard>
      {tip && <ResidueTooltip {...tip} />}
    </div>
  )
}

/** Residue hover tooltip: below-right of the pointer, flipped above the pointer when it would leave the viewport, clamped horizontally. */
const TIP_W = 230, TIP_PAD = 8
function ResidueTooltip({ x, y, title, lines }: { x: number; y: number; title: string; lines: string[] }) {
  const ref = useRef<HTMLDivElement>(null)
  const [h, setH] = useState(0)
  useLayoutEffect(() => { if (ref.current) setH(ref.current.offsetHeight) }, [title, lines])
  const below = y + 14, above = y - 10 - h
  const top = h && below + h > window.innerHeight - TIP_PAD ? Math.max(TIP_PAD, above) : below
  const left = Math.max(TIP_PAD, Math.min(x + 12, window.innerWidth - TIP_W - TIP_PAD))
  return (
    <div ref={ref} role="tooltip" data-testid="annotation-tooltip" style={{ left, top, width: TIP_W }} className="pointer-events-none fixed z-50 rounded border border-stone-300 bg-white px-2 py-1.5 text-xs text-stone-800 shadow">
      <p className="font-semibold">{title}</p>
      {lines.map((l, i) => <p key={i} className="leading-4">{l}</p>)}
    </div>
  )
}

/** Collapsible category: header (toggle button + ⓘ popover, kept outside the body) above a body that scrolls vertically beyond MAX_BODY px (natural height below). */
const MAX_BODY = 220
function Category({ id, title, count, info, children }: { id: string; title: string; count: number; info: string; children: ReactNode }) {
  const [open, setOpen] = useState(true)
  return (
    <section className="rounded border border-stone-200" data-testid={`category-${id}`} data-open={String(open)}>
      <div className="flex items-center gap-1 bg-sky-50 pr-2" data-testid={`category-header-${id}`}>
        <button type="button" aria-expanded={open} onClick={() => setOpen((o) => !o)} className="flex min-w-0 flex-1 items-center gap-2 px-3 py-1 text-left text-sm font-medium text-stone-900 hover:bg-sky-100" data-testid={`category-toggle-${id}`}>
          <span aria-hidden="true" className={`inline-block text-xs text-stone-500 transition-transform ${open ? 'rotate-90' : ''}`}>▶</span>{title} <span className="text-xs font-normal text-stone-500">({count} {count === 1 ? 'track' : 'tracks'})</span>
        </button>
        <InfoTip label={`About ${title}`} text={info} />
      </div>
      {open && <div className="overflow-y-auto overflow-x-hidden px-1 py-0.5" style={{ maxHeight: MAX_BODY }} data-testid={`category-body-${id}`}>{children}</div>}
    </section>
  )
}

/** Download dropdown: icon + label + chevron; items are plain download links of the family's static JSON files (no new endpoints). Unavailable files are
 *  detected with a HEAD request when the menu opens and shown disabled; Escape / outside click closes; arrow keys move between items. */
interface DownloadItem { id: string; label: string; file: string; detail: string }
function DownloadMenu({ acc, items }: { acc: string; items: DownloadItem[] }) {
  const [open, setOpen] = useState(false)
  const [avail, setAvail] = useState<Record<string, boolean | null>>({})
  const root = useRef<HTMLDivElement>(null)
  useEffect(() => {
    if (!open) return
    const onDoc = (e: MouseEvent) => { if (!root.current?.contains(e.target as Node)) setOpen(false) }
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') setOpen(false) }
    document.addEventListener('mousedown', onDoc); document.addEventListener('keydown', onKey)
    for (const it of items) if (avail[it.id] === undefined) {
      setAvail((a) => ({ ...a, [it.id]: null }))
      // available = served as JSON (an SPA fallback answers a missing file with HTML and status 200)
      fetch(familyFileUrl(acc, it.file), { method: 'HEAD' }).then((r) => setAvail((a) => ({ ...a, [it.id]: r.ok && (r.headers.get('content-type') ?? '').includes('json') }))).catch(() => setAvail((a) => ({ ...a, [it.id]: false })))
    }
    return () => { document.removeEventListener('mousedown', onDoc); document.removeEventListener('keydown', onKey) }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [open, acc, items])
  const move = (e: React.KeyboardEvent<HTMLDivElement>) => {
    if (e.key !== 'ArrowDown' && e.key !== 'ArrowUp') return
    const els = [...(root.current?.querySelectorAll<HTMLElement>('[role="menuitem"]:not([aria-disabled="true"])') ?? [])]
    const i = els.indexOf(document.activeElement as HTMLElement); const n = els.length
    if (n) { els[(i + (e.key === 'ArrowDown' ? 1 : n - 1)) % n].focus(); e.preventDefault() }
  }
  return (
    <div ref={root} className="relative" onKeyDown={move} data-testid="download-menu">
      <button type="button" aria-haspopup="menu" aria-expanded={open} onClick={() => setOpen((o) => !o)} data-testid="download-button"
        className="inline-flex items-center gap-1.5 rounded border border-stone-300 bg-white px-2 py-1 text-xs text-stone-800 hover:bg-stone-50 focus:outline-none focus:ring-2 focus:ring-blue-500">
        <svg width="12" height="12" viewBox="0 0 16 16" aria-hidden="true" fill="none" stroke="currentColor" strokeWidth="1.6"><path d="M8 2v8m0 0 3-3m-3 3L5 7M3 12v1.5h10V12" strokeLinecap="round" strokeLinejoin="round" /></svg>
        Download JSON
        <svg width="10" height="10" viewBox="0 0 16 16" aria-hidden="true" fill="none" stroke="currentColor" strokeWidth="1.8"><path d="m4 6 4 4 4-4" strokeLinecap="round" strokeLinejoin="round" /></svg>
      </button>
      {open && (
        <div role="menu" aria-label="Download annotation data" className="absolute right-0 z-40 mt-1 max-h-[70vh] w-72 max-w-[calc(100vw-2rem)] overflow-y-auto rounded border border-stone-300 bg-white p-1 text-xs shadow" data-testid="download-items">
          {items.map((it) => {
            const ok = avail[it.id]
            // one row: the download (label, left) and an info icon (right) that only shows the description; the icon is outside the link
            return (
              <div key={it.id} className="flex items-center gap-1 rounded hover:bg-stone-100" data-testid={`download-row-${it.id}`}>
                {ok === false
                  ? <div role="menuitem" aria-disabled="true" className="flex-1 px-2 py-1.5 text-stone-400" data-testid={`download-${it.id}`}>{it.label} <span className="text-[10px]">· not available for {acc}</span></div>
                  : <a role="menuitem" href={familyFileUrl(acc, it.file)} download={`${acc}_${it.file}`} onClick={() => setOpen(false)} data-testid={`download-${it.id}`} title={`${acc}_${it.file}`}
                      className="flex-1 rounded px-2 py-1.5 font-medium text-stone-800 focus:bg-stone-100 focus:outline-none">{it.label}{ok === null ? ' …' : ''}</a>}
                <span className="pr-1" data-testid={`download-info-${it.id}`}><InfoTip label={`About ${it.label}`} text={it.detail} /></span>
              </div>
            )
          })}
        </div>
      )}
    </div>
  )
}
