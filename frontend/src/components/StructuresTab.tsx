import { useFamily } from '../state/family'
import { useSelection } from '../state/selection'
import { displayOrganism } from '../lib/format'
import { InstanceTable } from './InstanceTable'
import { PdbeMolstarViewer } from './PdbeMolstarViewer'
import { PendingPanel } from './PendingPanel'

export function StructuresTab() {
  const { family, instances, orientation } = useFamily()
  const sel = useSelection()
  const current = instances.find((i) => i.seq_name === sel.instance) ?? family.display_representative ?? family.representative
  const o = current ? displayOrganism(current) : null
  return (
    <div className="space-y-4">
      <InstanceTable rows={instances} selected={current?.seq_name ?? null} onSelect={sel.setInstance} clen={family.clen} />
      <div className="grid gap-4 lg:grid-cols-3">
        <section className="lg:col-span-2">
          {current ? (
            <>
              <p className="mb-1 text-sm" data-testid="structures-current">Selected: <span className="font-mono">{current.pdb_id.toUpperCase()}</span> chain {current.auth_asym_id} · <span className={o?.isSpecies ? 'italic' : ''}>{o?.name}</span> · {current.experimental_method}{current.resolution != null ? ` ${current.resolution} Å` : ''} · <span title={current.title ?? ''}>{current.title}</span></p>
              <PdbeMolstarViewer pdbId={current.pdb_id} labelAsymId={current.label_asym_id} authAsymId={current.auth_asym_id} region={{ start: current.map_start, end: current.map_end }} orientation={orientation?.transforms[current.seq_name] ?? null} selection={[]} className="h-[520px]" />
            </>
          ) : <p>Select an instance.</p>}
        </section>
        <section className="space-y-3">
          <PendingPanel title="Pairwise structural similarity" detail="no similarity matrix has been generated yet (Milestone 3)." />
          <PendingPanel title="Hierarchical clustering / dendrogram" detail="no clustering has been computed yet (Milestone 3)." />
          <PendingPanel title="Structural groups and representatives" detail="depends on clustering; not available." />
        </section>
      </div>
    </div>
  )
}
