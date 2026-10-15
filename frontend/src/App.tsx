import { useEffect, useState } from 'react'
import { FamilyProvider, useFamilyStatus } from './state/family'
import { SelectionProvider, useSelection } from './state/selection'
import { Tabs, type TabId } from './components/Tabs'
import { FamilyOverview } from './components/FamilyOverview'
import { StructuresTab } from './components/StructuresTab'
import { AnnotationsTab } from './components/AnnotationsTab'

const ACC = 'RF00167'

function Shell() {
  const status = useFamilyStatus()
  const sel = useSelection()
  const [tab, setTab] = useState<TabId>('overview')
  useEffect(() => {   // default instance = backend representative
    if (status.state === 'ready' && sel.instance === null) { const d = status.data.family.display_representative ?? status.data.family.representative; if (d) sel.setInstance(d.seq_name) }
  }, [status, sel])
  if (status.state === 'loading') return <p className="p-6 text-stone-600">Loading {ACC}…</p>
  if (status.state === 'error') return <p className="p-6 text-red-700">Could not load family data: {status.message}</p>
  const { family } = status.data
  return (
    <div className="mx-auto max-w-7xl px-4 py-4">
      <header className="mb-2 border-b border-stone-200 pb-3" data-testid="app-header">
        <p className="text-xs font-medium uppercase tracking-wide text-stone-500" data-testid="app-title">Aggregated RNA 3D View</p>
        <div className="mt-1 flex flex-wrap items-center gap-x-3 gap-y-1">
          <h1 className="text-2xl font-semibold leading-tight text-stone-900 sm:text-3xl" data-testid="family-heading">{family.description}</h1>
          <a href={family.rfam_url ?? `https://rfam.org/family/${family.rfam_acc}`} target="_blank" rel="noreferrer" data-testid="family-badge"
            title={`Rfam family ${family.rfam_acc} (${family.rfam_id}); opens the Rfam entry`}
            className="rounded border border-stone-300 bg-stone-100 px-2 py-0.5 font-mono text-xs text-stone-700 hover:bg-stone-200">{family.rfam_acc}</a>
        </div>
        <p className="mt-1 max-w-3xl text-sm text-stone-600">Explore aggregated RNA structures, conformational diversity and residue-level annotations.</p>
      </header>
      <Tabs value={tab} onChange={setTab} />
      <main className="mt-4">
        {tab === 'overview' && <FamilyOverview />}
        {tab === 'structures' && <StructuresTab />}
        {tab === 'annotations' && <AnnotationsTab />}
      </main>
      <footer className="mt-8 border-t border-stone-200 pt-3 text-xs text-stone-500">
        Data: Rfam {family.provenance.rfam_release} · CM {family.provenance.cm_checksum} · {family.provenance.infernal_version} · generated {family.provenance.generated}.
        3D: PDBe Mol* (Apache-2.0). 2D: PDBe RNA Viewer 0.3.0 (Apache-2.0, vendored).
      </footer>
    </div>
  )
}

export default function App() {
  return (
    <FamilyProvider acc={ACC}>
      <SelectionProvider>
        <Shell />
      </SelectionProvider>
    </FamilyProvider>
  )
}
