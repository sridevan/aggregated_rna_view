import { useCallback, useEffect, useMemo, useState } from 'react'
import { FamilyProvider, useFamilyStatus } from './state/family'
import { SelectionProvider, useSelection } from './state/selection'
import { NumberingProvider } from './state/numbering'
import { StructuresProvider } from './state/structures'
import { NavigationProvider, type NavigationOrigin } from './state/navigation'
import { Tabs, type TabId } from './components/Tabs'
import { FamilyOverview } from './components/FamilyOverview'
import { StructuresTab } from './components/StructuresTab'
import { AnnotationsTab } from './components/AnnotationsTab'
import { LigandsTab } from './components/LigandsTab'
import { InteractionsTab } from './components/InteractionsTab'
import { HomePage } from './components/HomePage'

import { AVAILABLE_ACCS, canonicalPath, FAMILY_NAMES } from './lib/route'

/** /{rfam_acc} selects the family; "/" is the home page. Other paths are replaced (no reload, no history entry) by their canonical path
 *  (upper-case accession, or "/"). Returns the accession, or null for the home page. */
function familyFromLocation(): string | null {
  const canon = canonicalPath(window.location.pathname)
  if (canon !== window.location.pathname) window.history.replaceState(null, '', `${canon}${window.location.search}${window.location.hash}`)
  return canon === '/' ? null : canon.slice(1)
}

function Shell({ acc }: { acc: string }) {
  const status = useFamilyStatus()
  const sel = useSelection()
  const [tab, setTab] = useState<TabId>('overview')
  // explicit "View in Family Overview" from the Annotations or Ligands tab: select the Rfam position (current instance and numbering kept), open the tab
  const [origin, setOrigin] = useState<NavigationOrigin | null>(null)
  const { select } = sel
  const viewInOverview = useCallback((pos: number, label: string, from: NavigationOrigin['from'] = 'annotations') => { select(pos); setOrigin({ from, pos, label }); setTab('overview') }, [select])
  // the origin note belongs to that transferred selection only: any other selection (or clearing it) drops it
  const active = sel.active
  useEffect(() => { if (origin && !(active?.kind === 'position' && active.pos === origin.pos)) setOrigin(null) }, [active, origin])
  const nav = useMemo(() => ({ viewInOverview, origin }), [viewInOverview, origin])
  useEffect(() => {   // default instance = backend representative
    if (status.state === 'ready' && sel.instance === null) { const d = status.data.family.display_representative ?? status.data.family.representative; if (d) sel.setInstance(d.seq_name) }
  }, [status, sel])
  if (status.state === 'loading') return <p className="p-6 text-stone-600">Loading {acc}…</p>
  if (status.state === 'error') return (
    <div className="p-6 text-sm" data-testid="family-unavailable">
      {AVAILABLE_ACCS.includes(acc) ? <p className="text-red-700">Could not load family data: {status.message}</p>
        : <p className="text-stone-800">Rfam family <b>{acc}</b> is not available in this view.</p>}
      <p className="mt-2 text-stone-600">Available: {AVAILABLE_ACCS.map((a) => <a key={a} href={`/${a}`} className="mr-2 text-blue-700 underline">{a}</a>)}</p>
    </div>)
  const { family } = status.data
  return (
    <div className="mx-auto max-w-7xl px-4 py-4">
      <header className="mb-2 border-b border-stone-200 pb-3" data-testid="app-header">
        <p className="text-xs font-medium uppercase tracking-wide text-stone-500"><a href="/" className="hover:text-blue-700 hover:underline" data-testid="app-title" title="All families">Aggregated RNA 3D View</a></p>
        <div className="mt-1 flex flex-wrap items-center gap-x-3 gap-y-1">
          <h1 className="text-2xl font-semibold leading-tight text-stone-900 sm:text-3xl" data-testid="family-heading">{family.description}</h1>
          <a href={`https://rfam.org/family/${family.rfam_acc}`} target="_blank" rel="noreferrer" data-testid="family-badge"
            title={`Rfam family ${family.rfam_acc} (${family.rfam_id}); opens the Rfam entry`}
            className="rounded border border-stone-300 bg-stone-100 px-2 py-0.5 font-mono text-xs text-stone-700 hover:bg-stone-200">{family.rfam_acc}</a>
          {AVAILABLE_ACCS.length > 1 && (
            <label className="ml-auto flex items-center gap-1 text-xs text-stone-600">Family
              {/* each family has its own URL (/{rfam_acc}); switching loads that page */}
              <select data-testid="family-select" value={acc} onChange={(e) => { window.location.assign(`/${e.target.value}`) }} className="rounded border border-stone-300 bg-white px-1 py-0.5 text-xs text-stone-800">
                {AVAILABLE_ACCS.map((a) => <option key={a} value={a}>{a} · {FAMILY_NAMES[a] ?? a}</option>)}
              </select>
            </label>
          )}
        </div>
        <p className="mt-1 max-w-3xl text-sm text-stone-600">Explore aggregated RNA structures, conformational diversity, RNA–RNA interactions and residue-level annotations.</p>
      </header>
      <Tabs value={tab} onChange={setTab} />
      <main className="mt-4">
        {/* Position numbering is shared by the Family Overview, Ligands and Annotations tabs (kept across tab switches) */}
        <NumberingProvider>
        <NavigationProvider value={nav}>
        {tab === 'overview' && <FamilyOverview />}
        {tab === 'structures' && <StructuresTab />}
        {tab === 'ligands' && <LigandsTab />}
        {tab === 'interactions' && <InteractionsTab />}
        {tab === 'annotations' && <AnnotationsTab />}
        </NavigationProvider>
        </NumberingProvider>
      </main>
      <footer className="mt-8 border-t border-stone-200 pt-3 text-xs text-stone-500">
        Data: Rfam {family.provenance.rfam_release} · CM {family.provenance.cm_checksum} · {family.provenance.infernal_version} · generated {family.provenance.generated}.
        3D: PDBe Mol* (Apache-2.0). 2D: PDBe RNA Viewer 0.3.0 (Apache-2.0, vendored).
      </footer>
    </div>
  )
}

export default function App() {
  const [acc] = useState(familyFromLocation)
  if (acc === null) return <HomePage />
  return (
    <FamilyProvider acc={acc}>
      <SelectionProvider>
        <StructuresProvider>
        <Shell acc={acc} />
        </StructuresProvider>
      </SelectionProvider>
    </FamilyProvider>
  )
}
