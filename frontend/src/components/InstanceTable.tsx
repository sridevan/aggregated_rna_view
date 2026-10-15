import { useEffect, useMemo, useState } from 'react'
import type { InstanceRow } from '../types/backend'
import { displayOrganism } from '../lib/format'
import { InfoTip } from './InfoTip'

const PAGE = 10
type Col = { key: string; label: string; get: (r: InstanceRow) => string | number | null; num?: boolean; tip?: string; render?: (r: InstanceRow) => React.ReactNode }

export function InstanceTable({ rows, selected, onSelect, clen }: { rows: InstanceRow[]; selected: string | null; onSelect: (seqName: string) => void; clen: number }) {
  // `selected` is the complete instance identity (seq_name = pdb_entity_authChain_labelAsym), never the PDB id alone
  const [q, setQ] = useState('')
  const [sort, setSort] = useState<{ key: string; dir: 1 | -1 }>({ key: 'pdb_id', dir: 1 })
  const [page, setPage] = useState(1)
  const cols = useMemo<Col[]>(() => [
    { key: 'pdb_id', label: 'PDB ID', get: (r) => r.pdb_id.toUpperCase() },
    { key: 'chain', label: 'Chain', get: (r) => r.auth_asym_id, tip: 'Author chain identifier. The label (mmCIF) asym id is kept internally for Mol* selections.' },
    { key: 'title', label: 'Title', get: (r) => r.title, render: (r) => <span className="block max-w-[16rem] truncate" title={r.title ?? ''} data-testid="cell-title">{r.title ?? '–'}</span> },
    { key: 'organism', label: 'Organism', get: (r) => displayOrganism(r).name, render: (r) => { const o = displayOrganism(r); return o.name ? <span className={o.isSpecies ? 'italic' : ''} title={r.organism ?? ''} data-testid="cell-organism">{o.name}</span> : '–' } },
    { key: 'method', label: 'Method', get: (r) => r.experimental_method },
    { key: 'resolution', label: 'Res. (Å)', get: (r) => r.resolution, num: true, tip: 'Resolution in ångström (PDBe); lower is better. Not applicable to some cryo-EM entries.' },
    { key: 'region', label: 'Region', get: (r) => `${r.map_start}–${r.map_end} (${r.mapped_length} nt)`, tip: 'Polymer (label_seq_id) positions of the PDBe/Rfam-mapped segment and its length in nucleotides.' },
    { key: 'seqcov', label: 'Seq. cov. %', get: (r) => +(r.fraction_clen_matched * 100).toFixed(1), num: true, tip: `Consensus positions of the Rfam model to which this chain contributes a match-state residue, divided by the model length (${clen}). Insertions and local-end residues do not count.` },
    { key: 'strcov', label: 'Struct. cov. %', get: (r) => +((r.n_match_with_coordinates / clen) * 100).toFixed(1), num: true, tip: `Match-state residues that have experimental coordinates in this chain, divided by the model length (${clen}). Lower than sequence coverage when parts of the chain are unmodelled.` },
    { key: 'avgpp', label: 'Avg PP', get: (r) => r.cmalign?.avg_pp ?? null, num: true, tip: 'Average Infernal (cmalign) posterior probability of the aligned residues of the mapped region, as reported by cmalign; 1.0 means every residue placement is near-certain under the model.' },
    { key: 'bits', label: 'Rfam bits', get: (r) => r.rfam_bit_score, num: true, tip: 'Rfam cmsearch bit score of this chain against the family covariance model (from the Rfam structures endpoint). Higher = closer to the family model; it is not a resolution or quality measure.' },
  ], [clen])
  const filtered = useMemo(() => {
    const needle = q.trim().toLowerCase()
    const f = needle ? rows.filter((r) => [r.pdb_id, r.auth_asym_id, r.label_asym_id, r.organism, r.organism_species, r.experimental_method, r.title, r.molecule_name].join(' ').toLowerCase().includes(needle)) : rows.slice()
    const col = cols.find((c) => c.key === sort.key)!
    return f.sort((a, b) => {
      const va = col.get(a), vb = col.get(b)
      if (va == null && vb == null) return 0
      if (va == null) return 1
      if (vb == null) return -1
      return (col.num ? (va as number) - (vb as number) : String(va).localeCompare(String(vb))) * sort.dir
    })
  }, [rows, q, sort, cols])
  const pages = Math.max(1, Math.ceil(filtered.length / PAGE))
  const selectedIndex = selected ? filtered.findIndex((r) => r.seq_name === selected) : -1
  const selectedHidden = !!selected && rows.some((r) => r.seq_name === selected) && selectedIndex < 0
  // search/sort changes: go to the page holding the selected instance when it is still listed, else page 1
  useEffect(() => { setPage(selectedIndex >= 0 ? Math.floor(selectedIndex / PAGE) + 1 : 1) }, [q, sort])   // eslint-disable-line react-hooks/exhaustive-deps
  // selection changes (first visit, "Use as selected instance", Annotations selector): navigate to its page; user paging never alters the selection
  useEffect(() => { if (selectedIndex >= 0) setPage(Math.floor(selectedIndex / PAGE) + 1) }, [selected])   // eslint-disable-line react-hooks/exhaustive-deps
  const cur = Math.min(page, pages)
  const shown = filtered.slice((cur - 1) * PAGE, cur * PAGE)
  const from = filtered.length ? (cur - 1) * PAGE + 1 : 0, to = Math.min(cur * PAGE, filtered.length)
  return (
    <div className="rounded border border-stone-200 bg-white">
      <div className="border-b border-stone-200 p-2">
        <input data-testid="instance-search" value={q} onChange={(e) => setQ(e.target.value)} placeholder="Search PDB id, chain, organism, method, title…" className="w-full max-w-md rounded border border-stone-300 px-2 py-1 text-sm" />
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
        {selectedHidden && <span className="rounded bg-amber-50 px-2 py-0.5 text-xs text-amber-800" data-testid="selected-hidden">Selected instance {selected!.split('_')[0].toUpperCase()} chain {selected!.split('_')[2]} is hidden by the current search (still shown in 3D). <button className="underline" onClick={() => setQ('')}>Clear search</button></span>}
        <nav className="ml-auto flex items-center gap-3" aria-label="Pagination" data-testid="pagination">
          <button data-testid="page-prev" disabled={cur <= 1} onClick={() => setPage((p) => Math.max(1, p - 1))} className="rounded border border-stone-300 px-2 py-0.5 disabled:opacity-40">Previous</button>
          <span data-testid="page-info" className="whitespace-nowrap">Page {cur} of {pages}</span>
          <button data-testid="page-next" disabled={cur >= pages} onClick={() => setPage((p) => Math.min(pages, p + 1))} className="rounded border border-stone-300 px-2 py-0.5 disabled:opacity-40">Next</button>
        </nav>
      </div>
    </div>
  )
}
