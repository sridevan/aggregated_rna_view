import { useEffect, useState } from 'react'
import { loadFamily } from '../api/data'
import { AVAILABLE_ACCS, FAMILY_BLURBS, FAMILY_NAMES } from '../lib/route'
import { shortMethod } from '../lib/format'
import type { FamilyOverview } from '../types/backend'

/** Home page ("/"): what the app is, what a family page contains, and one card per available family. Card numbers come from each family's
 *  family.json (small file; nothing else is loaded here), so they always match the family page. */
type CardState = { state: 'loading' } | { state: 'error' } | { state: 'ready'; family: FamilyOverview }

const PAGE_PARTS: { title: string; text: string }[] = [
  { title: 'Family Overview', text: 'Sequence and structural coverage along the Rfam consensus, the RNA chains of the family, and a 2D diagram and 3D structure kept in sync: select a position, loop or nucleotide once and it is shown in every view.' },
  { title: 'Structures', text: 'How the experimental structures differ in 3D: pairwise structural comparison, structural groups and their representative structures, superposed in a shared frame.' },
  { title: 'Annotations', text: 'Residue-level annotations aggregated per consensus position: RNA 3D loops, ligand contacts, solvent accessibility, nonstandard nucleotides and predicted pockets, each with its own counts and missing-data rules.' },
]

export function HomePage() {
  const [cards, setCards] = useState<Record<string, CardState>>(() => Object.fromEntries(AVAILABLE_ACCS.map((a) => [a, { state: 'loading' }])))
  useEffect(() => {
    let alive = true
    for (const acc of AVAILABLE_ACCS) loadFamily(acc)
      .then((family) => { if (alive) setCards((c) => ({ ...c, [acc]: { state: 'ready', family } })) })
      .catch(() => { if (alive) setCards((c) => ({ ...c, [acc]: { state: 'error' } })) })
    return () => { alive = false }
  }, [])
  const release = Object.values(cards).find((c): c is Extract<CardState, { state: 'ready' }> => c.state === 'ready')?.family.provenance.rfam_release
  return (
    <div className="mx-auto max-w-6xl px-4 py-6" data-testid="home-page">
      <header className="border-b border-stone-200 pb-4">
        <p className="text-xs font-medium uppercase tracking-wide text-stone-500">Aggregated RNA 3D View</p>
        <h1 className="mt-1 text-2xl font-semibold leading-tight text-stone-900 sm:text-3xl">Experimental RNA structures, aggregated by RNA family</h1>
        <p className="mt-2 text-sm leading-6 text-stone-700" data-testid="home-intro">
          Each page brings together experimentally determined RNA structures from an <a href="https://rfam.org" target="_blank" rel="noreferrer" className="text-blue-700 underline">Rfam</a> family.
          Individual RNA chains are mapped to a shared set of Rfam consensus positions, enabling comparisons of sequence coverage, structural variation
          and residue-level annotations across structures.
        </p>
      </header>

      <section className="mt-6" aria-labelledby="families-heading">
        <h2 id="families-heading" className="text-lg font-semibold text-stone-900">Families</h2>
        <div className="mt-3 grid gap-4 md:grid-cols-2" data-testid="family-cards">
          {AVAILABLE_ACCS.map((acc) => <FamilyCard key={acc} acc={acc} card={cards[acc]} />)}
        </div>
      </section>

      <section className="mt-8" aria-labelledby="pages-heading">
        <h2 id="pages-heading" className="text-lg font-semibold text-stone-900">What a family page shows</h2>
        <div className="mt-3 grid gap-3 md:grid-cols-3" data-testid="page-parts">
          {PAGE_PARTS.map((p) => (
            <div key={p.title} className="rounded border border-stone-200 bg-white p-3">
              <h3 className="text-sm font-semibold text-stone-900">{p.title}</h3>
              <p className="mt-1 text-xs leading-5 text-stone-700">{p.text}</p>
            </div>))}
        </div>
      </section>

      <section className="mt-8 text-xs leading-5 text-stone-600" aria-labelledby="notes-heading" data-testid="home-notes">
        <h2 id="notes-heading" className="text-sm font-semibold text-stone-800">Good to know</h2>
        <ul className="mt-1 list-disc space-y-1 pl-5">
          <li>Only nucleotides aligned to a consensus position are compared across structures. Insertions, nucleotides outside the aligned region and missing coordinates are reported separately, never counted as absent.</li>
          <li>Residue numbers can be shown in the author numbering of any chosen structure (Position numbering); the underlying comparison always uses the consensus positions.</li>
          <li>Structural groups are labelled stable or provisional depending on how robust the grouping is.</li>
        </ul>
      </section>

      <footer className="mt-8 border-t border-stone-200 pt-3 text-xs text-stone-500" data-testid="home-footer">
        Data: {release ? `Rfam ${release}` : 'Rfam'} · PDBe · BGSU RNA 3D Motif Atlas. 3D: PDBe Mol* and 3Dmol.js. 2D: PDBe RNA Viewer.
      </footer>
    </div>
  )
}

function FamilyCard({ acc, card }: { acc: string; card: CardState }) {
  const f = card.state === 'ready' ? card.family : null
  const curated = (f?.provenance as { selection?: { n_chains?: number } } | undefined)?.selection
  const rep = f?.display_representative ?? f?.representative ?? null
  const methods = f ? Object.entries(f.experimental_methods).sort((a, b) => b[1] - a[1]).map(([m, n]) => `${shortMethod(m)} ${n}`).join(' · ') : ''
  return (
    // the whole card is one link (Enter activates it); the accession is a plain badge, so there are no nested interactive elements
    <a href={`/${acc}`} aria-label={`${FAMILY_NAMES[acc] ?? acc} (${acc}): open family page`}
      className="group flex cursor-pointer flex-col rounded-lg border border-stone-200 bg-white p-4 shadow-sm transition-colors hover:border-blue-300 hover:bg-blue-50/30 hover:shadow focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-600 focus-visible:ring-offset-2"
      data-testid={`family-card-${acc}`} data-state={card.state}>
      <div className="flex items-start gap-3 sm:gap-4">
        <Thumbnail acc={acc} name={FAMILY_NAMES[acc] ?? acc} rep={rep} />
        <div className="min-w-0 flex-1">
          <div className="flex items-start justify-between gap-2">
            <h3 className="text-base font-semibold text-stone-900 group-hover:text-blue-800">{FAMILY_NAMES[acc] ?? acc}</h3>
            <span className="shrink-0 rounded border border-stone-300 bg-stone-100 px-1.5 py-0.5 font-mono text-[11px] text-stone-700">{acc}</span>
          </div>
          {f && <p className="text-xs text-stone-500" data-testid="card-description">Rfam: {f.description}</p>}
          <p className="mt-2 text-sm leading-5 text-stone-700">{FAMILY_BLURBS[acc]}</p>
        </div>
      </div>
      {card.state === 'loading' && <p className="mt-3 text-xs text-stone-500">Loading…</p>}
      {card.state === 'error' && <p className="mt-3 text-xs text-amber-700">Summary unavailable; the family page may still open.</p>}
      {f && (
        <dl className="mt-3 grid grid-cols-3 gap-2 text-center" data-testid="card-stats">
          <Stat label="PDB entries" value={f.counts.pdb_entries} testId="card-entries" />
          <Stat label="RNA chains" value={f.counts.instances} testId="card-chains" />
          <Stat label="consensus positions" value={f.clen} testId="card-clen" />
        </dl>
      )}
      {f && (
        <div className="mt-2 space-y-0.5 text-xs text-stone-600">
          <p data-testid="card-methods">{methods}</p>
          {rep && <p data-testid="card-representative">Shown by default: {rep.pdb_id.toUpperCase()} chain {rep.auth_asym_id}</p>}
          {curated && <p className="text-amber-800" data-testid="card-curated">Curated selection of {curated.n_chains} chains</p>}
        </div>
      )}
    </a>
  )
}

/** Static thumbnail of the default representative chain (public/thumbnails, generated by scripts/make-thumbnails.mjs from the real
 *  coordinates; no viewer on this page). A neutral placeholder replaces a missing image. Part of the card link: a click opens the family. */
function Thumbnail({ acc, name, rep }: { acc: string; name: string; rep: { pdb_id: string; auth_asym_id: string } | null }) {
  const [failed, setFailed] = useState(false)
  const box = 'h-[96px] w-[96px] shrink-0 rounded-md border border-stone-200 sm:h-[120px] sm:w-[120px]'
  if (failed) return <div className={`${box} flex items-center justify-center bg-stone-100 text-[10px] text-stone-400`} data-testid={`card-thumbnail-placeholder-${acc}`} aria-hidden="true">No image</div>
  return <img src={`/thumbnails/${acc}.webp`} width={120} height={120} loading="lazy" decoding="async" onError={() => setFailed(true)}
    alt={`Representative RNA structure of the ${name} family${rep ? `: ${rep.pdb_id.toUpperCase()} chain ${rep.auth_asym_id}` : ''}`}
    className={`${box} bg-white object-contain`} data-testid={`card-thumbnail-${acc}`} />
}

function Stat({ label, value, testId }: { label: string; value: number; testId: string }) {
  return (
    <div className="flex flex-col-reverse rounded bg-stone-50 px-1 py-1.5">
      <dt className="text-[10px] uppercase tracking-wide text-stone-500">{label}</dt>
      <dd className="text-lg font-semibold text-stone-900" data-testid={testId}>{value.toLocaleString('en-US')}</dd>
    </div>
  )
}
