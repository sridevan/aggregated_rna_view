import { useMemo, useState } from 'react'
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
/** gallery cards shown before "Show all" (a search always shows every match) */
export const GALLERY_INITIAL = 12

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
  const [showAll, setShowAll] = useState(false)
  const [selected, setSelected] = useState<string | null>(null)
  const [revealKey, setRevealKey] = useState(0)
  const ligTracks = useMemo(() => (st.state === 'ready' ? ligandTracks(st.tracks, st.ligands) : null), [st])
  const pockets = useMemo(() => (st.state === 'ready' && st.pockets ? pocketTrack(st.pockets) : null), [st])
  const all = useMemo(() => (st.state === 'ready' ? galleryItems(st.ligands) : []), [st])
  const items = useMemo(() => (st.state === 'ready' ? galleryItems(st.ligands, query) : []), [st, query])
  const shown = query.trim() || showAll ? items : items.slice(0, GALLERY_INITIAL)
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
          description={<span data-testid="ligand-gallery-summary">Ligands with PDBe-reported contacts to RNA chains of this family, with the number of PDB entries in which a contact is recorded. Metal ions, crystallisation additives and flagged artefacts are not shown. Select a ligand to locate its residue track.</span>}>
          {all.length === 0 ? <p className="text-xs text-stone-500" data-testid="ligand-gallery-empty">No ligands with PDBe-reported RNA contacts in the default ligand set for {acc}.</p> : <>
            <label className="mb-2 flex max-w-sm items-center gap-2 text-[11px] text-stone-700">Filter
              <input type="search" value={query} onChange={(e) => setQuery(e.target.value)} placeholder="Name or CCD id, e.g. adenine, HPA" data-testid="ligand-search"
                className="min-w-0 flex-1 rounded border border-stone-300 bg-white px-1.5 py-0.5 text-xs placeholder:text-stone-400" />
            </label>
            {items.length === 0 && <p className="text-xs text-stone-500" data-testid="ligand-gallery-no-match">No ligand matches “{query}”.</p>}
            <ul className="grid grid-cols-1 gap-2 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4" data-testid="ligand-gallery">
              {shown.map((it) => (
                <li key={it.ccd} className={`flex gap-2 rounded border p-1.5 ${selected === it.ccd ? 'border-violet-400 bg-violet-50 ring-1 ring-violet-300' : 'border-stone-200 bg-white'}`} data-testid={`ligand-card-${it.ccd}`} data-selected={String(selected === it.ccd)}>
                  <button type="button" onClick={() => choose(it.ccd)} aria-pressed={selected === it.ccd} title={selected === it.ccd ? 'Clear the highlighted track' : `Locate the ${it.ccd} residue track`}
                    className="relative h-20 w-20 shrink-0 rounded border border-stone-200 bg-white p-0.5 hover:border-violet-400 focus:outline-none focus:ring-2 focus:ring-violet-500" data-testid={`ligand-select-${it.ccd}`}>
                    <Depiction ccd={it.ccd} />
                  </button>
                  <div className="min-w-0 flex-1 text-xs leading-4">
                    <button type="button" onClick={() => choose(it.ccd)} className="font-mono font-semibold text-stone-900 hover:text-violet-800 hover:underline" data-testid={`ligand-ccd-${it.ccd}`}>{it.ccd}</button>
                    <p className="line-clamp-2 break-words text-stone-700" title={it.name} data-testid={`ligand-name-${it.ccd}`}>{it.name}</p>
                    <p className="mt-0.5 text-stone-600" title={`RNA contact recorded by PDBe in: ${it.entries.map((e) => e.toUpperCase()).join(', ')}`} data-testid={`ligand-entries-${it.ccd}`}>
                      Found in {it.nEntries} PDB {it.nEntries === 1 ? 'entry' : 'entries'}</p>
                    <a href={ligandPageUrl(it.ccd)} target="_blank" rel="noreferrer" className="text-[11px] text-blue-700 hover:underline" data-testid={`ligand-link-${it.ccd}`}>PDBe ligand page ↗</a>
                  </div>
                </li>))}
            </ul>
            {!query.trim() && items.length > GALLERY_INITIAL && (
              <button type="button" onClick={() => setShowAll((v) => !v)} className="mt-2 self-start rounded border border-stone-300 bg-white px-2 py-0.5 text-xs text-stone-700 hover:bg-stone-50" data-testid="ligand-gallery-toggle">
                {showAll ? 'Show fewer' : `Show all ${items.length} ligands`}</button>)}
          </>}
        </SectionCard>

        {/* B. ligand-contact residues and C. predicted pockets, on one shared axis */}
        <SectionCard testId="ligand-tracks-section" title="Ligand-contact residues and predicted pockets"
          description={<span data-testid="ligand-tracks-summary">Rfam positions with PDBe-reported ligand contacts (one track per ligand) and fpocketR-predicted pockets, across {family.counts.instances} RNA instances from {family.counts.pdb_entries} PDB entries.</span>}>
          <TrackBrowser categories={categories} origin="ligands" numberingInfo={TEXT.numbering} legendInfo={TEXT.frequency} highlight={selected ? `lig-${selected}` : null} revealKey={revealKey} />
          {!pockets && <p className="mt-2 text-[11px] text-stone-500" data-testid="categories-unavailable">Not available for {acc}: Predicted binding pockets.</p>}
        </SectionCard>

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
