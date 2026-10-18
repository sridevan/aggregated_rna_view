import { useId, useMemo, useState, type ReactNode } from 'react'
import { useStructures } from '../../state/structures'
import type { Clustering, Group } from '../../lib/structureGroups'
import type { InstanceRow, LigandsDataset, ModificationsDataset } from '../../types/backend'
import { shortMethod } from '../../lib/format'
import { InfoTip } from '../InfoTip'

/*
 * Associated features of the structural groups: neutral, counts-first comparison tables (adapted in layout from the reference demo's
 * AssociatedFeatures, without its enrichment chips, arrows, colours or 15-point display threshold).
 * Observation unit everywhere: the clustered RNA chain instance of the active metric (not-clustered instances are never counted).
 * Each instance counts at most once per feature. Denominators are the instances whose feature status is known for that source
 * (ligands: PDBe entry-entity evaluated; nonstandard nucleotides: deposited sequence read; species/method: value present); unknown
 * status is never counted as a negative. Provenance is kept: ligands are entity-level PDBe annotations, nonstandard nucleotides are
 * chain-specific and split by location (inside / outside the Rfam-mapped region).
 */
export interface FeatureColumn { id: string; label: string; members: string[]; color: string | null; selected: boolean }
export interface CountCell { k: number; n: number }
export interface FeatureRow { key: string; cells: CountCell[]; total: number }

/** Categorical comparison: per column, k = members with the feature (deduplicated by instance), n = members with known status. */
export function categoricalRows(columns: FeatureColumn[], known: (sn: string) => boolean, featuresOf: (sn: string) => Set<string>): FeatureRow[] {
  const keys = new Set<string>()
  for (const col of columns) for (const sn of col.members) if (known(sn)) for (const f of featuresOf(sn)) keys.add(f)
  const rows = [...keys].map((key) => {
    const cells = columns.map((col) => { const kn = col.members.filter(known); return { k: kn.filter((sn) => featuresOf(sn).has(key)).length, n: kn.length } })
    return { key, cells, total: cells.reduce((s, x) => s + x.k, 0) }
  })
  return rows.sort((a, b) => b.total - a.total || a.key.localeCompare(b.key))   // predictable order: most observed overall first
}
/** Median over members with a value; n = number of contributing instances (missing values excluded, never zero-filled). */
export function medianOf(members: string[], get: (sn: string) => number | null | undefined): { median: number | null; n: number } {
  const v = members.map(get).filter((x): x is number => x != null && Number.isFinite(x)).sort((a, b) => a - b)
  if (!v.length) return { median: null, n: 0 }
  const m = Math.floor(v.length / 2)
  return { median: v.length % 2 ? v[m] : (v[m - 1] + v[m]) / 2, n: v.length }
}
/** Columns follow the viewing mode of the Structures tab: the selected group alone, or, while representatives are compared, the two compared
 *  groups (selected first, as in the comparison controls). Derived from the existing comparison state; never every other instance. */
export function featureColumns(selected: Group, compared: Group | null): FeatureColumn[] {
  const col = (g: Group, sel: boolean): FeatureColumn => ({ id: `g${g.group_id}`, label: g.label, members: g.members, color: g.color, selected: sel })
  return compared && compared.group_id !== selected.group_id ? [col(selected, true), col(compared, false)] : [col(selected, true)]
}

const ROWS_SHOWN = 5   // per category by default; each category expands independently
const pct = (k: number, n: number) => (n ? `${Math.round((100 * k) / n)}%` : '–')

export function AssociatedFeatures({ c, group, compared = null, instances, ligands, mods, clen }: { c: Clustering; group: Group; compared?: Group | null; instances: InstanceRow[]; ligands: LigandsDataset | null; mods: ModificationsDataset | null; clen: number }) {
  const byId = useMemo(() => new Map(instances.map((r) => [r.seq_name, r])), [instances])
  const { featuresOpen: open, setFeaturesOpen: setOpen } = useStructures()
  const bodyId = useId()
  const cols = featureColumns(group, compared)
  const comparing = cols.length === 2
  const lig = useMemo(() => {
    const of = new Map<string, Set<string>>(), names = new Map<string, string | null>()
    for (const b of ligands?.bound_molecules ?? []) {
      if (!b.in_default_track) continue
      names.set(b.ccd_id, b.name)
      for (const sn of (b as unknown as { entity_chains?: string[] }).entity_chains ?? []) (of.get(sn) ?? of.set(sn, new Set()).get(sn)!).add(b.ccd_id)
    }
    const src = (ligands as unknown as { entry_entity_sources?: Record<string, { status: string }> } | null)?.entry_entity_sources ?? {}
    const known = (sn: string) => { const r = byId.get(sn); return !!r && src[`${r.pdb_id}:${r.entity_id}`]?.status === 'evaluated' }
    return { of, names, known }
  }, [ligands, byId])
  const mod = useMemo(() => {
    const inside = new Map<string, Set<string>>(), outside = new Map<string, Set<string>>(), names = new Map<string, string | null>()
    for (const o of mods?.observations ?? []) {
      names.set(o.ccd_id, o.chem_comp_name)
      const m = o.outside_mapped_region ? outside : inside
      ;(m.get(o.seq_name) ?? m.set(o.seq_name, new Set()).get(o.seq_name)!).add(o.ccd_id)
    }
    const ok = new Set((mods?.instances ?? []).filter((r) => r.status === 'ok').map((r) => r.seq_name))
    return { inside, outside, names, known: (sn: string) => ok.has(sn) }
  }, [mods])
  const species = (sn: string) => byId.get(sn)?.organism_species ?? byId.get(sn)?.organism ?? null
  const method = (sn: string) => byId.get(sn)?.experimental_method ?? null
  const sections: { id: string; title: string; tip: string; available: boolean; rows: FeatureRow[]; label: (key: string) => ReactNode; title_: (key: string) => string; empty: string }[] = [
    { id: 'ligands', title: 'Ligand associations', available: !!ligands, tip: 'Ligand associations are derived from PDBe entity-level annotations and do not necessarily indicate direct contact with each RNA chain. A count of 0 means not annotated by PDBe, not confirmed absence.',
      rows: categoricalRows(cols, lig.known, (sn) => lig.of.get(sn) ?? new Set()), label: (k) => <><span className="font-mono">{k}</span>{lig.names.get(k) ? <span className="text-stone-500"> — {titleCase(lig.names.get(k)!)}</span> : null}</>, title_: (k) => `${k}${lig.names.get(k) ? ` — ${lig.names.get(k)}` : ''}`, empty: `No ligand associations annotated for ${comparing ? 'these groups' : 'this group'}.` },
    { id: 'species', title: 'Species', available: true, tip: 'Counts indicate RNA instances associated with each organism.',
      rows: categoricalRows(cols, (sn) => species(sn) != null, (sn) => new Set([species(sn)!])), label: (k) => <span className="italic">{k}</span>, title_: (k) => k, empty: 'No organism information.' },
    { id: 'nonstandard', title: 'Nonstandard nucleotides (Rfam-mapped region)', available: !!mods, tip: 'Counts indicate RNA instances containing the specified nonstandard nucleotide component within the Rfam-mapped region. Components are not necessarily post-transcriptional modifications.',
      rows: categoricalRows(cols, mod.known, (sn) => mod.inside.get(sn) ?? new Set()), label: (k) => <><span className="font-mono">{k}</span>{mod.names.get(k) ? <span className="text-stone-500"> — {short(mod.names.get(k)!)}</span> : null}</>, title_: (k) => `${k}${mod.names.get(k) ? ` — ${mod.names.get(k)}` : ''}`, empty: `None in the mapped region for ${comparing ? 'these groups' : 'this group'}.` },
    { id: 'nonstandard-outside', title: 'Nonstandard nucleotides (outside the mapped region)', available: !!mods, tip: 'Counts indicate RNA instances containing the specified nonstandard component outside the Rfam-mapped region (for example terminal nucleotides of the construct).',
      rows: categoricalRows(cols, mod.known, (sn) => mod.outside.get(sn) ?? new Set()), label: (k) => <><span className="font-mono">{k}</span>{mod.names.get(k) ? <span className="text-stone-500"> — {short(mod.names.get(k)!)}</span> : null}</>, title_: (k) => `${k}${mod.names.get(k) ? ` — ${mod.names.get(k)}` : ''}`, empty: `None outside the mapped region for ${comparing ? 'these groups' : 'this group'}.` },
    { id: 'method', title: 'Experimental method', available: true, tip: 'Counts indicate RNA instances determined using each experimental method.',
      rows: categoricalRows(cols, (sn) => method(sn) != null, (sn) => new Set([method(sn)!])), label: (k) => <span>{shortMethod(k)}</span>, title_: (k) => k, empty: 'No method information.' },
  ]
  const nums = [
    { id: 'resolution', label: 'Median resolution (Å)', get: (sn: string) => byId.get(sn)?.resolution, digits: 2 },
    { id: 'strcov', label: 'Median structural coverage (%)', get: (sn: string) => { const r = byId.get(sn); return r ? (r.n_match_with_coordinates / clen) * 100 : null }, digits: 1 },
    { id: 'length', label: 'Median chain length (nt)', get: (sn: string) => byId.get(sn)?.chain_length, digits: 0 },
  ]
  const provisional = c.metric.clustering.validation === 'provisional'
  const head = (
    <thead>
      <tr className="text-[11px] text-stone-600">
        <th className="py-1 pr-2 text-left font-medium" />
        {cols.map((col) => (
          <th key={col.id} className={`w-[9.5rem] px-2 py-1 text-right align-bottom font-medium leading-tight ${comparing && col.selected ? 'bg-stone-100' : ''}`} data-testid={`feature-col-${col.id}`} data-selected={col.selected}>
            <span className="inline-flex items-center justify-end gap-1">{col.color && <span className="inline-block h-2 w-2 shrink-0 rounded-full" style={{ background: col.color }} />}{col.label}</span>
            <span className="block font-normal text-stone-500">(n = {col.members.length} RNA instances)</span>
          </th>))}
      </tr>
    </thead>
  )
  return (
    <section className="rounded border border-stone-200 bg-white" data-testid="associated-features" data-open={String(open)}>
      {/* collapsed by default: heading, selected-group badge and chevron only (state kept in the shared Structures state) */}
      <h2 className="text-base font-semibold text-stone-900">
        <button type="button" aria-expanded={open} aria-controls={bodyId} onClick={() => setOpen(!open)} data-testid="features-toggle"
          className="flex w-full items-center gap-2 rounded px-3 py-2 text-left hover:bg-stone-50 focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500">
          <span aria-hidden="true" className={`inline-block text-xs text-stone-500 transition-transform ${open ? 'rotate-90' : ''}`}>▶</span>
          Associated features
          {cols.map((col, k) => <span key={col.id} className="inline-flex items-center gap-1">{k > 0 && <span className="text-xs font-normal text-stone-500">vs</span>}<span className="inline-flex items-center gap-1 rounded-full border px-2 py-0.5 text-xs font-normal" style={{ borderColor: col.color ?? undefined }} data-testid={k === 0 ? 'features-group-badge' : 'features-compare-badge'}><span className="inline-block h-2 w-2 rounded-full" style={{ background: col.color ?? undefined }} />{col.label}</span></span>)}
        </button>
      </h2>
      {open && <div id={bodyId} className="px-3 pb-3" data-testid="features-body">
      <p className="flex flex-wrap items-center gap-x-1 gap-y-1 text-xs text-stone-600">
        <span data-testid="features-summary">Comparisons are descriptive; group sizes and repeated observations should be considered when interpreting differences.</span>
        <InfoTip label="About associated features" text="Features are summarised across clustered RNA instances. These are descriptive comparisons, not statistical enrichment tests. RNA instances from the same PDB entry are not independent observations." />
        {provisional && <span className="ml-1 rounded border border-amber-300 bg-amber-50 px-1 text-[10px] font-normal text-amber-800" title={c.metric.clustering.rationale}>provisional groups</span>}
      </p>
      <div className="mt-2 grid min-w-0 gap-x-6 gap-y-3 lg:grid-cols-2 [&>*]:min-w-0">
        {sections.map((sec) => <FeatureSection key={sec.id} id={sec.id} title={sec.title} tip={sec.tip} head={head} cols={cols} rows={sec.rows} available={sec.available} label={sec.label} rowTitle={sec.title_} empty={sec.empty} />)}
        <div data-testid="features-numeric">
          <p className="mb-1 flex items-center text-xs font-semibold text-stone-700">Continuous properties<InfoTip label="About continuous properties" text="Medians are calculated from RNA instances with available measurements." /></p>
          <div className="overflow-x-auto"><table className="w-full text-xs">
            {head}
            <tbody>
              {nums.map((x) => (
                <tr key={x.id} className="border-t border-stone-100">
                  <td className="py-0.5 pr-2 text-stone-700">{x.label}</td>
                  {cols.map((col) => { const m = medianOf(col.members, x.get); return <td key={col.id} className={`px-2 py-0.5 text-right tabular-nums ${comparing && col.selected ? 'bg-stone-50' : ''}`} title={`Median of ${m.n} RNA instances with a value`} data-testid={`median-${x.id}-${col.id}`} data-n={m.n}>{m.median == null ? '–' : m.median.toFixed(x.digits)}</td> })}
                </tr>))}
            </tbody>
          </table></div>
        </div>
      </div>
      </div>}
    </section>
  )
}

function FeatureSection({ id, title, tip, head, cols, rows, available, label, rowTitle, empty }: { id: string; title: string; tip: string; head: ReactNode; cols: FeatureColumn[]; rows: FeatureRow[]; available: boolean; label: (k: string) => ReactNode; rowTitle: (k: string) => string; empty: string }) {
  const [all, setAll] = useState(false)
  const shown = all ? rows : rows.slice(0, ROWS_SHOWN)
  const denomNote = cols.some((col, i) => rows[0] && rows[0].cells[i].n !== col.members.length)
  return (
    <div data-testid={`features-${id}`}>
      <p className="mb-1 flex items-center text-xs font-semibold text-stone-700">{title}<InfoTip label={`About ${title}`} text={tip} /></p>
      {!available ? <p className="text-xs text-stone-500">Source unavailable.</p> : rows.length === 0 ? <p className="text-xs text-stone-500" data-testid={`features-${id}-empty`}>{empty}</p> : (
        <>
          <div className="overflow-x-auto"><table className="w-full text-xs">
            {head}
            <tbody>
              {shown.map((r) => (
                <tr key={r.key} className="border-t border-stone-100" data-testid={`feature-row-${id}-${r.key}`}>
                  <td className="max-w-[13rem] truncate py-0.5 pr-2 text-stone-800" title={rowTitle(r.key)}>{label(r.key)}</td>
                  {r.cells.map((cell, i) => (
                    <td key={cols[i].id} className={`px-2 py-0.5 text-right tabular-nums ${cols.length > 1 && cols[i].selected ? 'bg-stone-50' : ''} ${cell.k === 0 ? 'text-stone-400' : 'text-stone-900'}`} title={`${pct(cell.k, cell.n)} of ${cell.n} RNA instances`} data-testid={`cell-${id}-${r.key}-${cols[i].id}`}>{cell.k}/{cell.n}</td>
                  ))}
                </tr>))}
            </tbody>
          </table></div>
          {rows.length > ROWS_SHOWN && <button type="button" className="mt-0.5 text-[11px] text-blue-700 underline" aria-expanded={all} onClick={() => setAll(!all)} data-testid={`features-${id}-toggle`}>{all ? 'Show fewer' : `Show all ${rows.length}`}</button>}
          {denomNote && <p className="mt-0.5 text-[10px] text-stone-500">Denominators count only RNA instances whose status is known for this source.</p>}
        </>
      )}
    </div>
  )
}
const titleCase = (s: string) => (s === s.toUpperCase() ? s[0] + s.slice(1).toLowerCase() : s)
const short = (s: string) => { const t = titleCase(s); return t.length > 34 ? `${t.slice(0, 33)}…` : t }
