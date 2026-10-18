import { useEffect, useMemo, useRef, useState, type ReactNode } from 'react'
import type { InstanceRow } from '../types/backend'
import { displayOrganism, shortMethod } from '../lib/format'
import { InfoTip } from './InfoTip'

const PAGE = 10
type Col = { key: string; label: string; get: (r: InstanceRow) => string | number | null; num?: boolean; tip?: string; render?: (r: InstanceRow) => React.ReactNode }

/** Optional structural-group information (Structures tab): badge per row, dendrogram-order sorting, selected-group row marking. */
export interface TableGroups {
  of: (seqName: string) => { kind: 'group'; id: number; label: string; color: string; representative: boolean; order: number } | { kind: 'excluded'; reason: string } | null
  selectedId: number | null
  onSelectGroup: (id: number) => void
}
const COMPACT = ['pdb_id', 'chain', 'group', 'title', 'method', 'resolution']

export function InstanceTable({ rows, selected, onSelect, clen, groups, compact = false, toolbar, restrict = null }: { rows: InstanceRow[]; selected: string | null; onSelect: (seqName: string) => void; clen: number; groups?: TableGroups; compact?: boolean; toolbar?: ReactNode; restrict?: ((r: InstanceRow) => boolean) | null }) {
  // `selected` is the complete instance identity (seq_name = pdb_entity_authChain_labelAsym), never the PDB id alone
  const [q, setQ] = useState('')
  // null = the default order: by structural group (active group first, its representative at the top) when groups exist, else by PDB ID
  const [userSort, setUserSort] = useState<{ key: string; dir: 1 | -1 } | null>(null)
  const sort = userSort ?? (groups ? { key: 'group', dir: 1 as const } : { key: 'pdb_id', dir: 1 as const })
  // sorting the Structural group column ascending again is the default order itself (page 1, representative on top)
  const backToDefault = useRef(false)
  const setSort = (f: (s: { key: string; dir: 1 | -1 }) => { key: string; dir: 1 | -1 }) => { const n = f(sort); const d = !!groups && n.key === 'group' && n.dir === 1; backToDefault.current = d; setUserSort(d ? null : n) }
  const defaultOrder = !userSort && !!groups
  const [page, setPage] = useState(1)
  const allCols = useMemo<Col[]>(() => [
    { key: 'pdb_id', label: 'PDB ID', get: (r) => r.pdb_id.toUpperCase(),
      render: (r) => <ExtLink href={`https://pdbe.org/${r.pdb_id.toLowerCase()}`} title={`Open PDB entry ${r.pdb_id.toUpperCase()} at PDBe`} testId={`link-pdb-${r.seq_name}`}>{r.pdb_id.toUpperCase()}</ExtLink> },
    { key: 'chain', label: 'Chain', get: (r) => r.auth_asym_id, tip: 'Author chain identifier. The label (mmCIF) asym id is kept internally for Mol* selections.',
      render: (r) => <ExtLink href={`https://www.ebi.ac.uk/pdbe/entry/pdb/${r.pdb_id.toLowerCase()}?activeTab=macromolecules&entity=${r.entity_id}`} title={`Open molecule (entity ${r.entity_id}) of ${r.pdb_id.toUpperCase()} at PDBe`} testId={`link-chain-${r.seq_name}`}>{r.auth_asym_id}</ExtLink> },
    ...(groups ? [{ key: 'group', label: 'Structural group', tip: 'Structural group of the RNA instance under the selected comparison metric (descriptive; not a biological state). Default order: the selected group first with its representative (medoid) at the top, then dendrogram order. Click a badge to select its group. Instances without enough observed core positions are not clustered.',
      get: (r: InstanceRow) => groupRank(groups, r.seq_name), num: true,   // numeric rank (a string sort put 10 before 9)
      render: (r: InstanceRow) => { const g = groups.of(r.seq_name)
        if (!g) return '–'
        if (g.kind === 'excluded') return <span className="rounded-full border border-stone-300 bg-stone-50 px-1.5 py-px text-[10px] text-stone-500" title={`Not clustered: ${g.reason}`} data-testid={`not-clustered-${r.seq_name}`}>Not clustered</span>
        const on = g.id === groups.selectedId
        return <span className="inline-flex items-center gap-1"><button type="button" onClick={(e) => { e.stopPropagation(); groups.onSelectGroup(g.id) }} data-testid={`group-badge-${r.seq_name}`} data-group={g.id}
          className={`inline-flex items-center gap-1 rounded-full border px-1.5 py-px text-[10px] ${on ? 'font-semibold' : ''}`} style={{ borderColor: g.color, background: on ? `color-mix(in srgb, ${g.color} 14%, white)` : undefined }}>
          <span className="inline-block h-2 w-2 rounded-full" style={{ background: g.color }} />{g.label.replace('Structural group', 'Group')}</button>{g.representative && <span className="rounded bg-stone-100 px-1 text-[9px] text-stone-600" title="Representative (medoid) of its structural group">medoid</span>}</span> } }] : []),
    { key: 'title', label: 'Title', get: (r) => r.title, render: (r) => <span className={`block truncate ${groups ? 'max-w-[11rem]' : 'max-w-[16rem]'}`} title={r.title ?? ''} data-testid="cell-title">{r.title ?? '–'}</span> },
    { key: 'organism', label: 'Organism', get: (r) => displayOrganism(r).name, render: (r) => { const o = displayOrganism(r); return o.name ? <span className={o.isSpecies ? 'italic' : ''} title={r.organism ?? ''} data-testid="cell-organism">{o.name}</span> : '–' } },
    { key: 'method', label: 'Method', get: (r) => shortMethod(r.experimental_method), render: (r) => <span title={r.experimental_method ?? ''} data-testid="cell-method">{shortMethod(r.experimental_method) ?? '–'}</span> },
    { key: 'resolution', label: 'Res. (Å)', get: (r) => r.resolution, num: true, tip: 'Resolution in ångström (PDBe); lower is better. Not applicable to some cryo-EM entries.' },
    { key: 'length', label: 'Length (nt)', get: (r) => r.mapped_length, num: true, tip: 'Length in nucleotides of the PDBe/Rfam-mapped RNA region (hover a value for its polymer label_seq_id range).',
      render: (r) => <span title={`RNA region: ${r.map_start}–${r.map_end}`} data-testid="cell-length">{r.mapped_length}</span> },
    { key: 'seqcov', label: 'Seq. cov. %', get: (r) => +(r.fraction_clen_matched * 100).toFixed(1), num: true, tip: `Consensus positions of the Rfam model to which this chain contributes a match-state residue, divided by the model length (${clen}). Insertions and local-end residues do not count.` },
    { key: 'strcov', label: 'Struct. cov. %', get: (r) => +((r.n_match_with_coordinates / clen) * 100).toFixed(1), num: true, tip: `Match-state residues that have experimental coordinates in this chain, divided by the model length (${clen}). Lower than sequence coverage when parts of the chain are unmodelled.` },
    { key: 'avgpp', label: 'Avg PP', get: (r) => r.cmalign?.avg_pp ?? null, num: true, tip: 'Average Infernal (cmalign) posterior probability of the aligned residues of the mapped region, as reported by cmalign; 1.0 means every residue placement is near-certain under the model.' },
    { key: 'bits', label: 'Rfam bits', get: (r) => r.rfam_bit_score, num: true, tip: 'Rfam cmsearch bit score of this chain against the family covariance model (from the Rfam structures endpoint). Higher = closer to the family model; it is not a resolution or quality measure.' },
  ], [clen, groups])
  const cols = useMemo(() => (compact ? allCols.filter((c) => COMPACT.includes(c.key)) : allCols), [allCols, compact])
  const filtered = useMemo(() => {
    const needle = q.trim().toLowerCase()
    const base = restrict ? rows.filter(restrict) : rows
    const f = needle ? base.filter((r) => [r.pdb_id, r.auth_asym_id, r.label_asym_id, r.organism, r.organism_species, r.experimental_method, r.title, r.molecule_name].join(' ').toLowerCase().includes(needle)) : base.slice()
    const col = cols.find((c) => c.key === sort.key) ?? allCols.find((c) => c.key === sort.key) ?? cols[0]
    return f.sort((a, b) => {
      const va = col.get(a), vb = col.get(b)
      if (va == null && vb == null) return 0
      if (va == null) return 1
      if (vb == null) return -1
      return (col.num ? (va as number) - (vb as number) : String(va).localeCompare(String(vb))) * sort.dir
    })
  }, [rows, q, sort, cols, allCols, restrict])
  const pages = Math.max(1, Math.ceil(filtered.length / PAGE))
  const selectedIndex = selected ? filtered.findIndex((r) => r.seq_name === selected) : -1
  const selectedHidden = !!selected && rows.some((r) => r.seq_name === selected) && selectedIndex < 0
  // first visit: page 1 in the default group order (the active group's representative on top); otherwise the selected instance's page
  const mounted = useRef(false)
  // search/sort changes by the user: go to the page holding the selected instance when it is still listed, else page 1
  useEffect(() => { if (!mounted.current) return; if (backToDefault.current) { backToDefault.current = false; setPage(1); return } setPage(selectedIndex >= 0 ? Math.floor(selectedIndex / PAGE) + 1 : 1) }, [q, userSort, restrict])   // eslint-disable-line react-hooks/exhaustive-deps
  // another active group in the default order: its representative is now first -> page 1 (a selection made in the same step wins below)
  useEffect(() => { if (mounted.current && defaultOrder) setPage(1) }, [groups?.selectedId])   // eslint-disable-line react-hooks/exhaustive-deps
  // selection changes ("Use as selected instance", heatmap/dendrogram, Annotations selector): navigate to its page; user paging never alters the selection
  useEffect(() => { if (!mounted.current) return; if (selectedIndex >= 0) setPage(Math.floor(selectedIndex / PAGE) + 1) }, [selected])   // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => { if (!groups && selectedIndex >= 0) setPage(Math.floor(selectedIndex / PAGE) + 1); mounted.current = true }, [])   // eslint-disable-line react-hooks/exhaustive-deps
  const cur = Math.min(page, pages)
  const shown = filtered.slice((cur - 1) * PAGE, cur * PAGE)
  const from = filtered.length ? (cur - 1) * PAGE + 1 : 0, to = Math.min(cur * PAGE, filtered.length)
  return (
    <div className="rounded border border-stone-200 bg-white">
      <div className="flex flex-wrap items-center gap-2 border-b border-stone-200 p-2">
        <input data-testid="instance-search" value={q} onChange={(e) => setQ(e.target.value)} placeholder="Search PDB id, chain, organism, method, title…" className="w-full max-w-md flex-1 rounded border border-stone-300 px-2 py-1 text-sm" />
        {toolbar}
      </div>
      <div className="overflow-x-auto">
        <table className="w-full text-left text-xs" data-testid="instance-table">
          <thead className="bg-stone-100 text-xs uppercase text-stone-600">
            <tr>{cols.map((c) => (
              <th key={c.key} scope="col" className="whitespace-nowrap px-1.5 py-1">
                <button type="button" className="cursor-pointer select-none" onClick={() => setSort((s) => ({ key: c.key, dir: s.key === c.key ? (s.dir === 1 ? -1 : 1) : 1 }))} aria-sort={sort.key === c.key ? (sort.dir === 1 ? 'ascending' : 'descending') : undefined}>
                  {c.label}{sort.key === c.key ? (sort.dir === 1 ? ' ▲' : ' ▼') : ''}
                </button>
                {c.tip && <InfoTip text={c.tip} label={`About ${c.label}`} />}
              </th>))}</tr>
          </thead>
          <tbody>
            {shown.map((r) => (
              <tr key={r.seq_name} data-testid={`row-${r.seq_name}`} data-selected={r.seq_name === selected ? 'true' : 'false'} onClick={() => onSelect(r.seq_name)} aria-selected={r.seq_name === selected} tabIndex={0}
                data-in-group={(() => { const g = groups?.of(r.seq_name); return g?.kind === 'group' && g.id === groups?.selectedId ? 'true' : undefined })()}
                style={(() => { const g = groups?.of(r.seq_name); return g?.kind === 'group' && g.id === groups?.selectedId ? { boxShadow: `inset 3px 0 0 ${g.color}` } : undefined })()}
                onKeyDown={(e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); onSelect(r.seq_name) } }}
                className={`cursor-pointer border-t border-stone-100 hover:bg-blue-50 focus:bg-blue-50 focus:outline-none ${r.seq_name === selected ? 'bg-blue-100 font-medium ring-1 ring-inset ring-blue-300' : ''}`}>
                {cols.map((c) => <td key={c.key} className="whitespace-nowrap px-1.5 py-1">{c.render ? c.render(r) : (c.get(r) ?? '–')}</td>)}
              </tr>))}
            {shown.length === 0 && <tr><td colSpan={cols.length} className="px-2 py-3 text-stone-500">No instances match the search.</td></tr>}
          </tbody>
        </table>
      </div>
      <div className="flex flex-wrap items-center justify-between gap-x-4 gap-y-2 border-t border-stone-200 p-2 text-sm" data-testid="table-footer">
        <span className="text-xs text-stone-500" data-testid="table-count">Showing {from}–{to} of {filtered.length} instances{filtered.length !== rows.length ? ` (filtered from ${rows.length})` : ''}</span>
        {selectedHidden && <span className="rounded bg-amber-50 px-2 py-0.5 text-xs text-amber-800" data-testid="selected-hidden">Selected instance {selected!.split('_')[0].toUpperCase()} chain {selected!.split('_')[2]} is hidden by the current search or filter. <button className="underline" onClick={() => setQ('')}>Clear search</button></span>}
        <nav className="ml-auto flex items-center gap-3" aria-label="Pagination" data-testid="pagination">
          <button data-testid="page-prev" disabled={cur <= 1} onClick={() => setPage((p) => Math.max(1, p - 1))} className="rounded border border-stone-300 px-2 py-0.5 disabled:opacity-40">Previous</button>
          <span data-testid="page-info" className="whitespace-nowrap">Page {cur} of {pages}</span>
          <button data-testid="page-next" disabled={cur >= pages} onClick={() => setPage((p) => Math.min(pages, p + 1))} className="rounded border border-stone-300 px-2 py-0.5 disabled:opacity-40">Next</button>
        </nav>
      </div>
    </div>
  )
}

/** External link inside a selectable row: opens in a new tab and never selects the row (click and keyboard). */
function ExtLink({ href, title, testId, children }: { href: string; title: string; testId: string; children: ReactNode }) {
  return <a href={href} target="_blank" rel="noreferrer" title={title} data-testid={testId} className="text-blue-700 underline-offset-2 hover:underline"
    onClick={(e) => e.stopPropagation()} onKeyDown={(e) => e.stopPropagation()}>{children}</a>
}

/** Default sort rank of the Structural group column: the active group first with its representative on top, then every group in
 *  dendrogram order; not-clustered instances (null) sort last. */
export function groupRank(groups: Pick<TableGroups, 'of' | 'selectedId'>, seqName: string): number | null {
  const g = groups.of(seqName)
  if (g?.kind !== 'group') return null
  const active = g.id === groups.selectedId
  return (active ? 0 : 1_000_000) + (active && g.representative ? -1 : g.order)
}
