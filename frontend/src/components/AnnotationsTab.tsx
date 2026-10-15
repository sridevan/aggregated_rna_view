import { useCallback, useMemo, useState } from 'react'
import { useFamily } from '../state/family'
import { useSelection } from '../state/selection'
import { useMapping } from '../state/useMapping'
import { indexMapping, insertionsByAnchor, positionStates, residueToConsensus, translateConsensus, type ResidueTarget } from '../lib/consensus'
import { ConsensusTracks } from './ConsensusTracks'
import { PdbeMolstarViewer } from './PdbeMolstarViewer'
import { PdbRnaViewer2D } from './PdbRnaViewer2D'
import { TrackLegend } from './TrackLegend'
import { displayOrganism } from '../lib/format'

export function AnnotationsTab() {
  const { acc, family, instances, stats, orientation } = useFamily()
  const sel = useSelection()
  const current = instances.find((i) => i.seq_name === sel.instance) ?? family.display_representative ?? family.representative ?? null
  const ms = useMapping(acc, current)
  const mapping = ms.state === 'ready' ? ms.mapping : null
  const idx = useMemo(() => (mapping ? indexMapping(mapping) : null), [mapping])
  const translation = useMemo(() => (mapping && idx ? translateConsensus(mapping, sel.positions, idx) : null), [mapping, idx, sel.positions])
  const instanceTargets = useMemo<ResidueTarget[]>(() => (mapping && idx ? [...sel.instanceResidues].map((id) => residueToConsensus(mapping, id, idx)).filter((t): t is ResidueTarget => !!t) : []), [mapping, idx, sel.instanceResidues])
  const states = useMemo(() => (mapping && idx ? positionStates(mapping, family.clen, idx) : undefined), [mapping, idx, family.clen])
  const insertions = useMemo(() => (mapping ? insertionsByAnchor(mapping) : undefined), [mapping])
  // residues sent to the viewers: consensus-derived residues with coordinates + instance-level residues with coordinates
  const viewerIds = useMemo(() => [...(translation?.highlightable ?? []), ...instanceTargets.filter((t) => t.has_coordinates)].map((t) => t.label_seq_id), [translation, instanceTargets])
  const [goto, setGoto] = useState('')

  // click from either viewer -> consensus position when it has one, else instance-level residue
  const onResidueClick = useCallback((labelSeqId: number) => {
    if (!mapping || !idx) return
    const t = residueToConsensus(mapping, labelSeqId, idx)
    if (!t) return
    if (t.rfam_consensus_pos != null) sel.toggle(t.rfam_consensus_pos)
    else sel.toggleInstanceResidue(labelSeqId)
  }, [mapping, idx, sel])
  const onToggle = useCallback((pos: number, shift: boolean) => { if (shift && sel.anchor != null) sel.range(sel.anchor, pos); else sel.toggle(pos) }, [sel])

  if (typeof window !== 'undefined') window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), selection: { positions: [...sel.positions], instance: sel.instance, instanceResidues: [...sel.instanceResidues] }, translation, viewerIds }

  return (
    <div className="space-y-3">
      <div className="flex flex-wrap items-center gap-3 rounded border border-stone-200 bg-white p-2 text-sm">
        <label>Instance <select data-testid="instance-select" value={current?.seq_name ?? ''} onChange={(e) => sel.setInstance(e.target.value)} className="rounded border border-stone-300 px-1 py-0.5">
          {instances.map((i) => <option key={i.seq_name} value={i.seq_name}>{i.pdb_id.toUpperCase()} {i.auth_asym_id} · {displayOrganism(i).name} · {i.experimental_method?.split(' ')[0]} {i.resolution ?? ''}{i.pdbe_2d_available === false ? ' · no 2D' : ''}</option>)}
        </select></label>
        <form onSubmit={(e) => { e.preventDefault(); const n = +goto; if (n >= 1 && n <= family.clen) sel.toggle(n); setGoto('') }} className="flex items-center gap-1">
          <label>Consensus position <input data-testid="goto-input" value={goto} onChange={(e) => setGoto(e.target.value)} type="number" min={1} max={family.clen} className="w-20 rounded border border-stone-300 px-1 py-0.5" /></label>
          <button data-testid="goto-button" className="rounded bg-blue-600 px-2 py-0.5 text-white">Toggle</button>
        </form>
        <button data-testid="clear-selection" onClick={sel.clear} className="rounded border border-stone-300 px-2 py-0.5">Clear selection</button>
        <span className="text-stone-500">{sel.positions.size} consensus position(s) selected{sel.instanceResidues.size ? `, ${sel.instanceResidues.size} instance-only residue(s)` : ''}</span>
      </div>
      <ConsensusTracks stats={stats} selected={sel.positions} onToggle={onToggle} instanceStates={states} insertions={insertions} instanceLabel={current ? `${current.pdb_id.toUpperCase()} ${current.auth_asym_id}` : undefined} />
      <TrackLegend stats={stats} instanceLabel={current ? `${current.pdb_id.toUpperCase()} ${current.auth_asym_id}` : undefined} />
      <SelectionPanel translation={translation} instanceTargets={instanceTargets} mappingState={ms.state} />
      <div className="grid gap-3 lg:grid-cols-2">
        {current && <PdbRnaViewer2D pdbId={current.pdb_id} entityId={current.entity_id} authAsymId={current.auth_asym_id} selection={viewerIds} onResidueClick={onResidueClick} className="h-[560px]" />}
        {current && <PdbeMolstarViewer pdbId={current.pdb_id} labelAsymId={current.label_asym_id} authAsymId={current.auth_asym_id} region={{ start: current.map_start, end: current.map_end }} orientation={orientation?.transforms[current.seq_name] ?? null} selection={viewerIds.map((label_seq_id) => ({ label_seq_id }))} onResidueClick={onResidueClick} className="h-[560px]" />}
      </div>
      <p className="text-xs text-stone-500">2D: PDBe RNA Viewer (author chain {current?.auth_asym_id}; selection by label_seq_id). 3D: PDBe Mol* (label asym {current?.label_asym_id}, label_seq_id), showing the Rfam-mapped region by default; selection colours residues without moving the camera (use Focus selection). Only residues with coordinates are highlighted. Residues without a consensus coordinate (insertions, local ends, flanks) are selectable per instance and are cleared when the instance changes.</p>
    </div>
  )
}

function SelectionPanel({ translation, instanceTargets, mappingState }: { translation: ReturnType<typeof translateConsensus> | null; instanceTargets: ResidueTarget[]; mappingState: string }) {
  if (!translation) return <p className="text-sm text-stone-500" data-testid="selection-panel">{mappingState === 'loading' ? 'Loading chain mapping…' : 'No mapping loaded.'}</p>
  const fmt = (t: ResidueTarget) => `${t.rfam_consensus_pos ?? '–'} → ${t.nucleotide_parent}${t.label_seq_id} (auth ${t.author_residue_number}${t.author_insertion_code}${t.residue_name !== t.nucleotide_parent ? `, ${t.residue_name}` : ''}, PP ${t.posterior_prob_char || '?'})`
  return (
    <div className="grid gap-2 rounded border border-stone-200 bg-white p-2 text-xs md:grid-cols-5" data-testid="selection-panel">
      <div><span className="font-medium text-stone-700">Highlighted ({translation.highlightable.length})</span><ul data-testid="sel-highlightable">{translation.highlightable.map((t) => <li key={t.label_seq_id}>{fmt(t)}</li>)}</ul></div>
      <div><span className="font-medium text-stone-700">No coordinates ({translation.unobserved.length})</span><ul data-testid="sel-unobserved">{translation.unobserved.map((t) => <li key={t.label_seq_id}>{fmt(t)}</li>)}</ul></div>
      <div><span className="font-medium text-stone-700">Deleted in chain ({translation.deleted.length})</span><div data-testid="sel-deleted">{translation.deleted.join(', ')}</div></div>
      <div><span className="font-medium text-stone-700">Uncovered by alignment ({translation.uncovered.length})</span><div data-testid="sel-uncovered">{translation.uncovered.join(', ')}</div></div>
      <div><span className="font-medium text-stone-700">Instance-only residues ({instanceTargets.length})</span><ul data-testid="sel-instance">{instanceTargets.map((t) => <li key={t.label_seq_id}>{t.alignment_state}{t.insertion_anchor != null ? ` after ${t.insertion_anchor}.${t.insertion_offset}` : ''}: {t.nucleotide_parent}{t.label_seq_id} (auth {t.author_residue_number}){t.has_coordinates ? '' : ' no coords'}</li>)}</ul></div>
    </div>
  )
}
