import { lazy, Suspense, useMemo, useState } from 'react'
/** the Mol* viewer is loaded only when a user first opens an ideal 3D structure */
const LigandIdealViewer = lazy(() => import('./molstar/LigandIdealViewer').then((m) => ({ default: m.LigandIdealViewer })))
import { useFamily } from '../state/family'
import { useAnnotations } from '../state/useAnnotations'
import { SectionCard } from './SectionCard'
import { InfoTip } from './InfoTip'
import { defaultLigandCcds, DownloadMenu, ligandTracks, pocketTrack, titleCase, TrackBrowser, type CategorySpec, type DownloadItem } from './AnnotationsTab'
import type { LigandsDataset } from '../types/backend'

/*
 * Ligands tab: (A) gallery of RNA-contacting ligands (the default ligand set of ligands.json: PDBe binding-site annotations of the family's
 * RNA entities, class "ligand", not flagged), (B) the per-CCD ligand-contact residue tracks and (C) the fpocketR predicted-pocket track,
 * both in the shared TrackBrowser (same axis, numbering, zoom, legend, tooltips and pinned details as the Annotations tab). Nothing is
 * recomputed here: counts are by_ccd.n_entries_in_default_track (distinct PDB entries) and the tracks read tracks.json / pocket_tracks.json.
 */

/** 2D depictions by CCD id: PDBe chemistry images (as on PDBe-KB), RCSB CCD images as fallback; both serve every default-set CCD (checked 2026-10-22). */
export const depictionUrls = (ccd: string) => [
  `https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/${encodeURIComponent(ccd)}_200.svg`,
  `https://cdn.rcsb.org/images/ccd/unlabeled/${encodeURIComponent(ccd[0])}/${encodeURIComponent(ccd)}.svg`,
]
export const ligandPageUrl = (ccd: string) => `https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/${encodeURIComponent(ccd)}`

export const TEXT = {
  gallery: 'Ligands shown here have contacts reported by PDBe for RNA entities belonging to this Rfam family. Contacts are derived from Arpeggio and may not represent every RNA chain or biological assembly in an entry.',
  residues: 'Ligand-contact residues are obtained from PDBe binding-site annotations and mapped to Rfam consensus positions. The annotations are aggregated at RNA entity level rather than individual chain level.',
  pockets: 'Computationally predicted RNA binding pockets identified using fpocketR. These annotations are independent of the PDBe ligand-contact data.',
  frequency: 'Contact frequency is the proportion of eligible PDB entries containing a ligand in which PDBe reports a contact at this position. Entries with unobserved residues are excluded from the denominator. Pocket colours show the proportion of RNA instances with a predicted pocket.',
  numbering: 'Residue numbering is taken from the selected RNA chain. Ligand-contact annotations remain aggregated at RNA entity level and are not necessarily specific to the selected chain. Switching numbering systems changes residue labels only; the underlying ligand contacts and frequencies remain unchanged.',
}
/** gallery page size: 3 columns x 3 rows on desktop */
export const PAGE_SIZE = 9
/** page numbers to show: all when they fit, else first, last and the current page with its neighbours ('…' marks a gap) */
export function pageList(page: number, pages: number): (number | '…')[] {
  if (pages <= 7) return Array.from({ length: pages }, (_, i) => i + 1)
  const keep = [...new Set([1, pages, page - 1, page, page + 1].filter((n) => n >= 1 && n <= pages))].sort((a, b) => a - b)
  return keep.flatMap((n, i) => (i && n - keep[i - 1] > 1 ? ['…' as const, n] : [n]))
}

/** Default-set ligands in track order, with the PDB entries in which PDBe records an RNA contact; search by CCD id or name (case-insensitive). */
export function galleryItems(ligands: LigandsDataset, query = '') {
  const entries = new Map<string, Set<string>>()
  for (const b of ligands.bound_molecules) if (b.in_default_track) (entries.get(b.ccd_id) ?? entries.set(b.ccd_id, new Set()).get(b.ccd_id)!).add(b.pdb_id)
  const q = query.trim().toLowerCase()
  return defaultLigandCcds(ligands)
    .map(([ccd, v]) => ({ ccd, name: titleCase(v.name), nEntries: v.n_entries_in_default_track, entries: [...(entries.get(ccd) ?? [])].sort() }))
    .filter((it) => !q || it.ccd.toLowerCase().includes(q) || it.name.toLowerCase().includes(q))
}

/** PDBe depiction, then RCSB, then a plain placeholder; a light skeleton while loading. Images are never generated locally. */
function Depiction({ ccd }: { ccd: string }) {
  const urls = depictionUrls(ccd)
  const [k, setK] = useState(0)
  const [loaded, setLoaded] = useState(false)
  if (k >= urls.length) return <div className="flex h-full w-full items-center justify-center text-center text-[10px] leading-3 text-stone-400" data-testid={`depiction-missing-${ccd}`}>No depiction available</div>
  return (
    <>
      {!loaded && <div className="absolute inset-1 animate-pulse rounded bg-stone-100" aria-hidden="true" />}
      <img src={urls[k]} alt={`2D structure of ${ccd}`} loading="lazy" className={`h-full w-full object-contain ${loaded ? '' : 'opacity-0'}`} data-testid={`depiction-${ccd}`} data-source={k === 0 ? 'pdbe' : 'rcsb'}
        onLoad={() => setLoaded(true)} onError={() => { setLoaded(false); setK((i) => i + 1) }} />
    </>
  )
}

export function LigandsTab() {
  const { acc, family } = useFamily()
  const st = useAnnotations(acc)
  const [query, setQuery] = useState('')
  const [page, setPage] = useState(1)
  const [viewing, setViewing] = useState<{ ccd: string; name: string } | null>(null)
  const [selected, setSelected] = useState<string | null>(null)
  const [revealKey, setRevealKey] = useState(0)
  const ligTracks = useMemo(() => (st.state === 'ready' ? ligandTracks(st.tracks, st.ligands) : null), [st])
  const pockets = useMemo(() => (st.state === 'ready' && st.pockets ? pocketTrack(st.pockets) : null), [st])
  const all = useMemo(() => (st.state === 'ready' ? galleryItems(st.ligands) : []), [st])
  const items = useMemo(() => (st.state === 'ready' ? galleryItems(st.ligands, query) : []), [st, query])
  // pagination over the filtered list; a new query starts again at page 1
  const pages = Math.max(1, Math.ceil(items.length / PAGE_SIZE))
  const cur = Math.min(page, pages)
  const shown = items.slice((cur - 1) * PAGE_SIZE, cur * PAGE_SIZE)
  const onQuery = (q: string) => { setQuery(q); setPage(1) }
  if (typeof window !== 'undefined') window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), annotationTracks: ligTracks ? { ligands: ligTracks.map((t) => t.id), pockets: pockets ? [pockets.id] : [] } : null }
  const categories = useMemo<CategorySpec[]>(() => !ligTracks ? [] : [
    { id: 'ligands', title: 'Ligand-contact residues', info: TEXT.residues, tracks: ligTracks, empty: 'No ligand-contact residue tracks for this family.' },
    ...(pockets ? [{ id: 'pockets', title: 'Predicted binding pockets', info: TEXT.pockets, tracks: [pockets] }] : []),
  ], [ligTracks, pockets])
  const choose = (ccd: string) => { setSelected((cur) => (cur === ccd ? null : ccd)); setRevealKey((k) => k + 1) }
  const datasets: DownloadItem[] = [
    { id: 'ligands', label: 'Ligand binding annotations', file: 'ligands.json', detail: 'Bound molecules (CCD id, name, class), binding-site nucleotides, supporting entries and chains, and provenance.' },
    { id: 'tracks', label: 'Aggregated annotation tracks', file: 'tracks.json', detail: 'Per-position annotation frequencies with numerators, denominators and missing-data counts.' },
    { id: 'pockets', label: 'Predicted pockets · Per structure', file: 'pockets.json', detail: 'fpocketR pockets per RNA chain with scores, volumes, descriptors and lining nucleotides mapped to Rfam positions.' },
    { id: 'pocket-tracks', label: 'Predicted pockets · Aggregated', file: 'pocket_tracks.json', detail: 'Per-position pocket frequencies for each of the three pocket filters.' },
  ]
  return (
    <div className="space-y-3" data-testid="ligands-tab">
      {st.state === 'loading' && <p className="text-sm text-stone-500" data-testid="ligands-loading">Loading ligand annotations…</p>}
      {st.state === 'error' && <p className="text-sm text-red-700" data-testid="ligands-error">Ligand data not available: {st.message}</p>}
      {st.state === 'ready' && ligTracks && <>
        {/* A. gallery */}
        <SectionCard testId="ligand-gallery-section" actions={<DownloadMenu acc={acc} items={datasets} />}
          title={<span className="inline-flex items-center gap-1.5">RNA-contacting ligands <span className="text-sm font-normal text-stone-500" data-testid="ligand-gallery-count">({all.length})</span><InfoTip label="About RNA-contacting ligands" text={TEXT.gallery} /></span>}
          description={<span data-testid="ligand-gallery-summary">Ligands with PDBe-reported contacts to RNA chains of this family, with the number of PDB entries in which a contact is recorded. Metal ions, crystallisation additives and flagged artefacts are not shown. Select a structure image or “Locate track” to highlight its residue track; the CCD id opens the PDBe chemical component page.</span>}>
          {all.length === 0 ? <p className="text-xs text-stone-500" data-testid="ligand-gallery-empty">No ligands with PDBe-reported RNA contacts in the default ligand set for {acc}.</p> : <>
            <div className="mb-2 flex flex-wrap items-center gap-x-4 gap-y-1">
              <label className="flex min-w-0 max-w-sm flex-1 items-center gap-2 text-[11px] text-stone-700">Filter
                <input type="search" value={query} onChange={(e) => onQuery(e.target.value)} placeholder="Name or CCD id, e.g. adenine, HPA" data-testid="ligand-search"
                  className="min-w-0 flex-1 rounded border border-stone-300 bg-white px-1.5 py-0.5 text-xs placeholder:text-stone-400" />
              </label>
              {/* the highlighted track stays highlighted on every gallery page and through searches until cleared */}
              {selected && (
                <span className="inline-flex items-center gap-1 rounded border border-violet-300 bg-violet-50 px-1.5 py-px text-[11px] text-violet-900" data-testid="ligand-selected-chip">
                  Highlighted track: <span className="font-mono">{selected}</span>
                  <button type="button" onClick={() => setSelected(null)} aria-label={`Clear the highlighted ${selected} track`} className="ml-0.5 rounded px-1 hover:bg-violet-100" data-testid="ligand-selected-clear">✕</button>
                </span>)}
            </div>
            {items.length === 0 && <p className="text-xs text-stone-500" data-testid="ligand-gallery-no-match">No ligand matches “{query}”.</p>}
            <ul className="grid grid-cols-1 gap-3 sm:grid-cols-2 lg:grid-cols-3" data-testid="ligand-gallery" data-page={cur} data-pages={pages}>
              {shown.map((it) => (
                <li key={it.ccd} className={`flex h-full gap-3 rounded border p-2 ${selected === it.ccd ? 'border-violet-400 bg-violet-50 ring-1 ring-violet-300' : 'border-stone-200 bg-white'}`} data-testid={`ligand-card-${it.ccd}`} data-selected={String(selected === it.ccd)}>
                  <button type="button" onClick={() => choose(it.ccd)} aria-pressed={selected === it.ccd} aria-label={selected === it.ccd ? `Clear the highlighted ${it.ccd} track` : `Locate the ${it.ccd} residue track`}
                    title={selected === it.ccd ? 'Clear the highlighted track' : `Locate the ${it.ccd} residue track`}
                    className="relative h-32 w-32 shrink-0 rounded border border-stone-200 bg-white p-1 hover:border-violet-400 focus:outline-none focus:ring-2 focus:ring-violet-500" data-testid={`ligand-select-${it.ccd}`}>
                    <Depiction ccd={it.ccd} />
                  </button>
                  <div className="flex min-w-0 flex-1 flex-col text-xs leading-4">
                    <a href={ligandPageUrl(it.ccd)} target="_blank" rel="noreferrer" title={`PDBe chemical component page for ${it.ccd} (opens in a new tab)`} data-testid={`ligand-ccd-${it.ccd}`}
                      className="self-start rounded font-mono text-sm font-semibold text-blue-700 underline decoration-blue-300 underline-offset-2 hover:decoration-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-500">{it.ccd}<span aria-hidden="true" className="ml-0.5 text-[11px] no-underline">↗</span><span className="sr-only"> (opens the PDBe page in a new tab)</span></a>
                    <p className="mt-0.5 break-words text-stone-800" data-testid={`ligand-name-${it.ccd}`}>{it.name}</p>
                    <p className="mt-0.5 text-stone-600" title={`RNA contact recorded by PDBe in: ${it.entries.map((e) => e.toUpperCase()).join(', ')}`} data-testid={`ligand-entries-${it.ccd}`}>
                      Found in {it.nEntries} PDB {it.nEntries === 1 ? 'entry' : 'entries'}</p>
                    <div className="mt-auto flex flex-wrap gap-1.5 pt-1.5">
                      <button type="button" onClick={() => choose(it.ccd)} aria-pressed={selected === it.ccd} data-testid={`ligand-locate-${it.ccd}`}
                        className="rounded border border-stone-300 bg-white px-1.5 py-px text-[11px] text-stone-700 hover:bg-stone-50 focus:outline-none focus:ring-2 focus:ring-violet-500">{selected === it.ccd ? 'Clear track highlight' : 'Locate track'}</button>
                      <button type="button" onClick={() => setViewing({ ccd: it.ccd, name: it.name })} data-testid={`ligand-view3d-${it.ccd}`} title="Ideal 3D structure from the PDB Chemical Component Dictionary"
                        className="rounded border border-stone-300 bg-white px-1.5 py-px text-[11px] text-stone-700 hover:bg-stone-50 focus:outline-none focus:ring-2 focus:ring-blue-500">View ideal 3D</button>
                    </div>
                  </div>
                </li>))}
            </ul>
            {items.length > 0 && (
              <nav className="mt-2 flex flex-wrap items-center gap-x-3 gap-y-1 text-[11px] text-stone-600" aria-label="Ligand gallery pages" data-testid="ligand-pagination">
                <span data-testid="ligand-page-summary">Showing {(cur - 1) * PAGE_SIZE + 1}–{(cur - 1) * PAGE_SIZE + shown.length} of {items.length} {items.length === 1 ? 'ligand' : 'ligands'}{query.trim() ? ` matching “${query.trim()}”` : ''}</span>
                {pages > 1 && (
                  <span className="ml-auto flex items-center gap-1">
                    <button type="button" onClick={() => setPage(cur - 1)} disabled={cur === 1} aria-label="Previous page" data-testid="ligand-page-prev"
                      className="rounded border border-stone-300 bg-white px-1.5 py-px text-stone-700 hover:bg-stone-50 disabled:opacity-40">‹ Previous</button>
                    {pageList(cur, pages).map((n, i) => n === '…' ? <span key={`gap-${i}`} className="px-0.5">…</span> : (
                      <button key={n} type="button" onClick={() => setPage(n)} aria-current={n === cur ? 'page' : undefined} aria-label={`Page ${n}`} data-testid={`ligand-page-${n}`}
                        className={`min-w-[1.5rem] rounded border px-1 py-px ${n === cur ? 'border-blue-600 bg-blue-50 font-medium text-blue-800' : 'border-stone-300 bg-white text-stone-700 hover:bg-stone-50'}`}>{n}</button>))}
                    <button type="button" onClick={() => setPage(cur + 1)} disabled={cur === pages} aria-label="Next page" data-testid="ligand-page-next"
                      className="rounded border border-stone-300 bg-white px-1.5 py-px text-stone-700 hover:bg-stone-50 disabled:opacity-40">Next ›</button>
                    <span className="ml-1" data-testid="ligand-page-indicator">Page {cur} of {pages}</span>
                  </span>)}
              </nav>)}
          </>}
        </SectionCard>

        {/* B. ligand-contact residues and C. predicted pockets, on one shared axis */}
        <SectionCard testId="ligand-tracks-section" title="Ligand-contact residues and predicted pockets"
          description={<span data-testid="ligand-tracks-summary">Rfam positions with PDBe-reported ligand contacts (one track per ligand) and fpocketR-predicted pockets, across {family.counts.instances} RNA instances from {family.counts.pdb_entries} PDB entries.</span>}>
          <TrackBrowser categories={categories} origin="ligands" numberingInfo={TEXT.numbering} legendInfo={TEXT.frequency} highlight={selected ? `lig-${selected}` : null} revealKey={revealKey} />
          {!pockets && <p className="mt-2 text-[11px] text-stone-500" data-testid="categories-unavailable">Not available for {acc}: Predicted binding pockets.</p>}
        </SectionCard>

        {/* ideal 3D structure of one chemical component, on demand (one Mol* plugin while open, disposed on close) */}
        {viewing && <Suspense fallback={<div className="fixed inset-0 z-50 flex items-center justify-center bg-stone-900/40 text-sm text-white" data-testid="ligand-viewer-loading-code">Loading 3D viewer…</div>}>
          <LigandIdealViewer ccd={viewing.ccd} name={viewing.name} pageUrl={ligandPageUrl(viewing.ccd)} onClose={() => setViewing(null)} /></Suspense>}

        {/* methodology, collapsed by default */}
        <SectionCard testId="ligand-methods" title="About these annotations" headingLevel={3} collapsible defaultOpen={false}>
          <dl className="grid grid-cols-1 gap-x-4 gap-y-1 text-xs leading-5 text-stone-700 sm:grid-cols-[12rem_1fr]" data-testid="ligand-methods-list">
            {[['Source', 'PDBe binding-site annotations (PDBe API v2 binding_sites), one record per PDB entry, RNA entity and ligand.'],
              ['Contact identification', 'Interactions computed by PDBe with Arpeggio. These are computed contacts in the deposited structures, not proof of an experimentally validated binding event.'],
              ['Annotation level', 'RNA entity: contacts of all PDBe-processed copies of a ligand are pooled on the RNA sequence, so they are not always attributable to one RNA chain.'],
              ['Coverage', "Ligand copies processed by PDBe in each entry's preferred biological assembly. Copies in other assemblies of the same entry are not annotated."],
              ['Coordinate mapping', 'Contacting nucleotides are mapped to Rfam consensus positions through the validated chain-to-Rfam alignment.'],
              ['Alternative numbering', 'Author residue numbering of a selected RNA chain, used as labels only.'],
              ['Frequency unit', 'Distinct PDB entries containing the ligand, counted once whatever the number of RNA chains or ligand copies; entries without an observed nucleotide at a position are not counted there.'],
              ['Limitation', 'Contact attribution to individual RNA chains and ligand copies is not always available; ligand-copy counts are therefore not shown.'],
              ['Predicted pockets', 'fpocketR predictions on each RNA chain, independent of the PDBe contact data and never combined with it.'],
            ].map(([k, v]) => <div key={k} className="contents"><dt className="font-medium text-stone-900">{k}</dt><dd>{v}</dd></div>)}
          </dl>
        </SectionCard>
      </>}
    </div>
  )
}
