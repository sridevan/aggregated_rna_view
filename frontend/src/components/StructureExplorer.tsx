import { useCallback, useEffect, useMemo, useState } from 'react'
import { loadInstanceLoops, loadLigandClasses } from '../api/data'
import { instanceLoops, loopSelection } from '../lib/loops'
import { authorLabel, boundaryText, insertionMarkers, insertionRuns, insertionSelection, type InsertionMarker, type InsertionRun } from '../lib/insertions'
import { HoverCard } from './InfoTip'
import type { ChainMapping, OrientationDataset, OrientationTransform } from '../types/backend'
import { isApplicable } from '../lib/orientation'
import type { LoopRecord } from '../types/backend'
import { useFamily } from '../state/family'
import { useSelection } from '../state/selection'
import { useMapping } from '../state/useMapping'
import { NumberingSelect, useNumbering } from '../state/numbering'
import { useNavigation } from '../state/navigation'
import { indexMapping, positionStates, residueToConsensus, translateConsensus, type ResidueTarget } from '../lib/consensus'
import { ConsensusTracks, type Tip } from './ConsensusTracks'
import { authorResidueId, chainAt, residueValue, STATUS_TEXT } from '../lib/numbering'
import type { PositionState } from '../lib/consensus'
import { PdbeMolstarViewer } from './PdbeMolstarViewer'
import { PdbRnaViewer2D } from './PdbRnaViewer2D'
import { TrackLegend } from './TrackLegend'
import { displayOrganism } from '../lib/format'

/**
 * Structural exploration of one Rfam-mapped chain instance: shared instance selector, consensus/coverage tracks with the instance state,
 * selection summary, PDBe RNA 2D viewer and PDBe Mol* 3D viewer, all synchronised through the consensus-position selection.
 * (Moved unchanged from the former Annotations tab into Family Overview.)
 */
/** one default RNA colour in the Family Overview Mol* viewer for every instance (slate grey-blue #64748B); never derived from PDB, chain or group */
const FAMILY_OVERVIEW_RNA_COLOR = 0x64748b

export function StructureExplorer() {
  const { acc, family, instances, stats, orientation } = useFamily()
  const sel = useSelection()
  const current = instances.find((i) => i.seq_name === sel.instance) ?? family.display_representative ?? family.representative ?? null
  const ms = useMapping(acc, current)
  const mapping = ms.state === 'ready' ? ms.mapping : null
  const idx = useMemo(() => (mapping ? indexMapping(mapping) : null), [mapping])
  const translation = useMemo(() => (mapping && idx ? translateConsensus(mapping, sel.positions, idx) : null), [mapping, idx, sel.positions])
  const instanceTargets = useMemo<ResidueTarget[]>(() => (mapping && idx ? [...sel.instanceResidues].map((id) => residueToConsensus(mapping, id, idx)).filter((t): t is ResidueTarget => !!t) : []), [mapping, idx, sel.instanceResidues])
  const states = useMemo(() => (mapping && idx ? positionStates(mapping, family.clen, idx) : undefined), [mapping, idx, family.clen])
  // Position numbering (display only): Rfam consensus by default, or the author numbering of a chosen chain, independent of the active instance.
  // Its mapping is the chain's validated mapping file (the active chain's own mapping is reused when it is the same chain).
  const { numbering } = useNumbering()
  const origin = useNavigation()?.origin ?? null
  const chainLabel = current ? `${current.pdb_id.toUpperCase()} ${current.auth_asym_id}` : ''
  /** selected-chain tooltip: Rfam position, the active chain's author residue, the reference residue (when another chain numbers the axis), status */
  const chainTipFor = useCallback((pos: number, st: PositionState): Tip => {
    const own = idx ? chainAt(idx, pos) : { state: st, residue: null }
    const other = numbering && numbering.seqName !== current?.seq_name ? numbering : null
    return { title: `Selected RNA chain · ${chainLabel}`, main: '', rows: [
      { key: 'rfam', label: 'Rfam position', value: String(pos) },
      { key: 'author', label: other ? `${chainLabel} author residue` : 'Author residue', value: residueValue(own) },
      ...(other ? [{ key: 'reference', label: `${other.label} reference residue`, value: residueValue(chainAt(other.idx, pos)) }] : []),
      { key: 'status', label: 'Status', value: STATUS_TEXT[st] },
    ] }
  }, [idx, numbering, current?.seq_name, chainLabel])
  const refSuffix = (pos: number) => (numbering && numbering.seqName !== current?.seq_name ? ` · ${numbering.label} reference ${residueValue(chainAt(numbering.idx, pos))}` : '')
  const numberingControl = <NumberingSelect />
  // Rfam insertion events of this chain (insert / local-end alignment states only) and their markers on the shared grid (CELL = 9 units)
  const insRuns = useMemo(() => (mapping ? insertionRuns(mapping, family.clen) : []), [mapping, family.clen])
  const insMarkers = useMemo(() => insertionMarkers(insRuns, 9), [insRuns])
  const [insOpen, setInsOpen] = useState<{ id: string; rect: DOMRect } | null>(null)
  useEffect(() => { setInsOpen(null) }, [current?.seq_name])           // details of the previous instance never survive a switch
  const activeIns = sel.active?.kind === 'group' && sel.active.type === 'INS' ? sel.active : null
  const activeInsRun = activeIns ? insRuns.find((r) => activeIns.id === `insertion:${current?.seq_name}:${r.id}`) ?? null : null
  const toggleInsertionView = useCallback((r: InsertionRun) => { if (current) sel.toggleGroup(insertionSelection(r, current.pdb_id, current.seq_name)) }, [current, sel])
  // marker click: open its details (never a consensus selection); the open marker clicked again closes them and clears its 2D/3D view if shown
  const onInsertionClick = useCallback((id: string, rect: DOMRect) => {
    if (insOpen?.id === id) {
      setInsOpen(null)
      const mk = insMarkers.find((m) => m.id === id)
      if (activeInsRun && mk?.runs.some((r) => r.id === activeInsRun.id)) toggleInsertionView(activeInsRun)
    } else setInsOpen({ id, rect })
  }, [insOpen, insMarkers, activeInsRun, toggleInsertionView])
  // per-instance BGSU loops: the selected chain's own loop file, loaded once the loop section has been opened (and again for each new instance)
  const [loopsWanted, setLoopsWanted] = useState(false)
  const [loopData, setLoopData] = useState<{ sn: string; records: LoopRecord[] } | null>(null)
  const [loopsState, setLoopsState] = useState<'idle' | 'loading' | 'ready' | 'error'>('idle')
  useEffect(() => {
    if (!loopsWanted || !current) return
    let alive = true
    setLoopsState('loading')
    loadInstanceLoops(acc, current.seq_name).then((r) => { if (alive) { setLoopData({ sn: current.seq_name, records: r }); setLoopsState('ready') } }).catch(() => { if (alive) setLoopsState('error') })
    return () => { alive = false }
  }, [loopsWanted, acc, current?.seq_name])   // eslint-disable-line react-hooks/exhaustive-deps
  const loopRecords = loopData && loopData.sn === current?.seq_name ? loopData.records : null
  const openLoops = useCallback(() => setLoopsWanted(true), [])
  const loops = useMemo(() => (loopRecords && current ? instanceLoops(loopRecords, current.seq_name) : null), [loopRecords, current])
  const activeLoop = sel.active?.kind === 'group' ? loops?.find((l) => l.id === (sel.active as { id: string }).id) ?? null : null
  const onSelectLoop = useCallback((id: string) => {
    const l = loops?.find((x) => x.id === id), rec = loopRecords?.find((r) => r.loop_id === id)
    // click on the active loop again clears it (and restores the default views); any other loop replaces it without an intermediate reset
    if (l && rec && current) sel.toggleGroup(loopSelection(l, current.pdb_id, current.seq_name, rec.strands))
  }, [loops, loopRecords, current, sel])
  // residues sent to the viewers. A loop: exactly its members in this chain that have coordinates (from the loop's own membership, never a
  // residue range). Otherwise: consensus-derived residues with coordinates + instance-level residues with coordinates.
  const viewerIds = useMemo(() => (activeLoop ? activeLoop.inChain.filter((n) => n.has_coordinates).map((n) => n.label_seq_id)
    : [...(translation?.highlightable ?? []), ...instanceTargets.filter((t) => t.has_coordinates)].map((t) => t.label_seq_id)), [activeLoop, translation, instanceTargets])

  // click from either viewer -> consensus position when it has one, else instance-level residue
  const onResidueClick = useCallback((labelSeqId: number) => {
    if (!mapping || !idx) return
    const t = residueToConsensus(mapping, labelSeqId, idx)
    if (!t) return
    if (t.rfam_consensus_pos != null) sel.toggle(t.rfam_consensus_pos)
    else sel.toggleInstanceResidue(labelSeqId)
  }, [mapping, idx, sel])
  const onToggle = useCallback((pos: number) => sel.toggle(pos), [sel])
  // concise availability of the active selection in the current instance
  const status = useMemo(() => {
    const a = sel.active
    if (!a) return null
    if (a.kind === 'position') {
      const ref = refSuffix(a.pos)
      if (!translation) return { ok: false, text: `Rfam position ${a.pos}: mapping not loaded${ref}` }
      const h = translation.highlightable[0], u = translation.unobserved[0]
      const id = (t: ResidueTarget) => authorResidueId({ residue_name: t.residue_name, author_residue_number: t.author_residue_number, author_insertion_code: t.author_insertion_code }) ?? 'author residue unavailable'
      if (h) return { ok: true, text: `Rfam position ${a.pos} → ${chainLabel} ${id(h)}${ref} · highlighted and centred in 2D/3D` }
      if (u) return { ok: false, text: `Rfam position ${a.pos} → ${chainLabel} ${id(u)} has no experimental coordinates: nothing highlighted${ref}` }
      if (translation.deleted.length) return { ok: false, text: `Rfam position ${a.pos} is deleted in ${chainLabel}: nothing highlighted${ref}` }
      return { ok: false, text: `Rfam position ${a.pos} is not covered by the alignment of ${chainLabel}: nothing highlighted${ref}` }
    }
    if (a.kind === 'residue') { const t = instanceTargets[0]; return t ? { ok: t.has_coordinates, text: `${t.alignment_state} residue ${t.nucleotide_parent}${t.label_seq_id} (auth ${t.author_residue_number}) of ${a.seq_name}${t.has_coordinates ? '' : ', no coordinates'}` } : { ok: false, text: `Residue of another instance (${a.seq_name})` } }
    if (activeLoop) {
      const extra = [activeLoop.nOtherChain ? `${activeLoop.nOtherChain} in another chain (not shown)` : '', activeLoop.nWithoutCoordinates ? `${activeLoop.nWithoutCoordinates} without coordinates` : ''].filter(Boolean).join(' · ')
      return { ok: true, text: `${activeLoop.typeName} · ${activeLoop.id} · ${activeLoop.nStrands} segment${activeLoop.nStrands === 1 ? '' : 's'} · ${activeLoop.nucleotides.length} nucleotides${extra ? ` · ${extra}` : ''}`,
        title: activeLoop.nucleotides.map((n) => n.unit_id).join(', ') }
    }
    if (a.type === 'INS') {
      if (!activeInsRun) return { ok: false, text: 'Insertion of another instance' }
      const shown = activeInsRun.residues.filter((x) => x.has_coordinates)
      return { ok: shown.length > 0, text: `Insertion +${activeInsRun.residues.length} · ${boundaryText(activeInsRun).toLowerCase()} · ${activeInsRun.residues.map((x) => `${x.nucleotide_parent}${x.label_seq_id}`).join(' ')} · ${shown.length} of ${activeInsRun.residues.length} with coordinates shown in 2D/3D (no consensus position)` }
    }
    return { ok: (translation?.highlightable.length ?? 0) > 0, text: `${a.type} ${a.id}: ${translation?.highlightable.length ?? 0} of ${a.members.length} member nucleotides highlighted in this instance` }
  }, [sel.active, translation, instanceTargets, current, activeLoop, activeInsRun, numbering, chainLabel])   // eslint-disable-line react-hooks/exhaustive-deps

  if (typeof window !== 'undefined') window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), selection: { positions: [...sel.positions], instance: sel.instance, instanceResidues: [...sel.instanceResidues], active: sel.active, viewReset: sel.viewReset }, translation, viewerIds,
    insertions: insMarkers.map((m) => ({ id: m.id, boundaries: m.boundaries, length: m.length, runs: m.runs.map((r) => ({ id: r.id, kind: r.kind, left: r.left, right: r.right, sequence: r.sequence, labelSeqIds: r.residues.map((x) => x.label_seq_id) })) })),
    numbering: numbering ? { seqName: numbering.seqName, label: numbering.label } : null,
    loops: loops ? loops.map((l) => ({ id: l.id, track: l.track, lane: l.lane, blocks: l.blocks, members: l.inChain.map((n) => n.label_seq_id) })) : null }

  return (
    <div className="space-y-2">
      <div className="flex flex-wrap items-center gap-x-3 gap-y-1 rounded border border-stone-200 bg-white px-2 py-1 text-sm" data-testid="instance-toolbar">
        <label className="flex min-w-0 max-w-full items-center gap-1">Instance <select data-testid="instance-select" value={current?.seq_name ?? ''} onChange={(e) => sel.setInstance(e.target.value)} className="min-w-0 max-w-[calc(100vw-7rem)] rounded border border-stone-300 px-1 py-0.5 sm:max-w-md">
          {instances.map((i) => <option key={i.seq_name} value={i.seq_name}>{i.pdb_id.toUpperCase()} {i.auth_asym_id} · {displayOrganism(i).name} · {i.experimental_method?.split(' ')[0]} {i.resolution ?? ''}{i.pdbe_2d_available === false ? ' · no 2D' : ''}</option>)}
        </select></label>
        <button data-testid="clear-selection" onClick={sel.clear} disabled={!sel.active} className="rounded border border-stone-300 px-2 py-0.5 disabled:cursor-default disabled:opacity-40">Clear selection</button>
        {/* one short hint while nothing is selected; the concise summary of the active selection replaces it */}
        {!status && <span className="text-xs text-stone-500" data-testid="selection-hint">Click a nucleotide or loop to highlight it in 2D and 3D.</span>}
        {status && <span className={`text-xs ${status.ok ? 'text-stone-600' : 'text-amber-700'}`} data-testid="selection-status" data-ok={String(status.ok)} title={'title' in status ? (status as { title?: string }).title : undefined}>{status.text}</span>}
        {/* a position transferred from the Annotations or Ligands tab: the evidence there is family-level, not an observation in this chain */}
        {origin && sel.active?.kind === 'position' && sel.active.pos === origin.pos && (
          <span className="basis-full text-[11px] text-stone-500" data-testid="selection-origin">Opened from {origin.from === 'ligands' ? 'Ligands' : 'Annotations'} ({origin.label}).{origin.seqName && origin.seqName === current?.seq_name ? ` Recorded for this RNA chain (${chainLabel}).` : ` Family-level evidence: not necessarily observed in ${chainLabel}.`}</span>)}
      </div>
      <ConsensusTracks stats={stats} selected={sel.positions} onToggle={onToggle} instanceStates={states} insertionMarkers={insMarkers} openInsertion={insOpen?.id ?? null} selectedInsertionRun={activeInsRun?.id ?? null} onInsertionClick={onInsertionClick} instanceLabel={current ? `${current.pdb_id.toUpperCase()} ${current.auth_asym_id}` : undefined}
        reference={{ rfamAcc: family.rfam_acc, rfamId: family.rfam_id, rf: family.consensus_sequence, ssCons: family.ss_cons, clanMembers: family.clan_members }}
        loops={loops} loopsState={loopsState} onLoopsOpen={openLoops} selectedLoop={activeLoop?.id ?? null} onSelectLoop={onSelectLoop}
        numbering={numbering} numberingControl={numberingControl} chainTipFor={chainTipFor} />
      {insOpen && mapping && (() => { const mk = insMarkers.find((m) => m.id === insOpen.id); return mk ? <InsertionDetails marker={mk} mapping={mapping} anchor={insOpen.rect} activeRun={activeInsRun?.id ?? null} onView={toggleInsertionView} onClose={() => setInsOpen(null)} /> : null })()}
      <TrackLegend stats={stats} instanceLabel={current ? `${current.pdb_id.toUpperCase()} ${current.auth_asym_id}` : undefined} reference={{ rfamAcc: family.rfam_acc, rfamId: family.rfam_id, rf: family.consensus_sequence, ssCons: family.ss_cons, clanMembers: family.clan_members }} />
      <SelectionPanel translation={translation} instanceTargets={instanceTargets} mappingState={ms.state} compact={!!activeLoop} />
      <div className="grid gap-3 lg:grid-cols-2">
        {current && <PdbRnaViewer2D pdbId={current.pdb_id} entityId={current.entity_id} authAsymId={current.auth_asym_id} selection={viewerIds} resetKey={sel.viewReset} onResidueClick={onResidueClick} className="h-[560px]" />}
        {current && <PdbeMolstarViewer pdbId={current.pdb_id} labelAsymId={current.label_asym_id} authAsymId={current.auth_asym_id} region={{ start: current.map_start, end: current.map_end }} orientation={orientation?.transforms[current.seq_name] ?? null} selection={viewerIds.map((label_seq_id) => ({ label_seq_id }))} focusOnSelection resetKey={sel.viewReset} structureColor={FAMILY_OVERVIEW_RNA_COLOR} onClearSelection={sel.clear} canClear={!!sel.active} onResidueClick={onResidueClick} ligandClasses={() => loadLigandClasses(acc)} className="h-[560px]" />}
      </div>
      <details className="rounded border border-stone-200 bg-white text-xs text-stone-700" data-testid="viewer-details">
        <summary className="cursor-pointer select-none px-2 py-1 font-medium text-stone-800" data-testid="viewer-details-toggle">Viewer details</summary>
        <div className="grid gap-x-6 gap-y-2 border-t border-stone-200 p-2 md:grid-cols-[3fr_2fr]" data-testid="viewer-details-body">
          <section data-testid="viewer-details-viewers">
            <h4 className="mb-0.5 font-semibold text-stone-800">Molecular viewers</h4>
            <p><b>2D viewer:</b> Displays the RNA secondary structure and highlights selected nucleotides or loops.</p>
            <p><b>3D viewer:</b> Displays the Rfam-mapped RNA region by default. Select <i>Complete structure</i> to view the entire deposited structure.</p>
            <p className="mt-0.5">Selecting a nucleotide or loop highlights it in both viewers and focuses it in 3D. Click the selection again, or use <i>Clear selection</i>, to restore the default views.</p>
          </section>
          <AlignmentDetails transform={current ? orientation?.transforms[current.seq_name] ?? null : null} reference={orientation?.reference ?? null} />
        </div>
      </details>
    </div>
  )
}

function SelectionPanel({ translation, instanceTargets, mappingState, compact = false }: { translation: ReturnType<typeof translateConsensus> | null; instanceTargets: ResidueTarget[]; mappingState: string; compact?: boolean }) {
  if (!translation) return <p className="text-sm text-stone-500" data-testid="selection-panel">{mappingState === 'loading' ? 'Loading chain mapping…' : 'No mapping loaded.'}</p>
  const fmt = (t: ResidueTarget) => `${t.rfam_consensus_pos ?? '–'} → ${t.nucleotide_parent}${t.label_seq_id} (auth ${t.author_residue_number}${t.author_insertion_code}${t.residue_name !== t.nucleotide_parent ? `, ${t.residue_name}` : ''}, PP ${t.posterior_prob_char || '?'})`
  return (
    <div className="grid gap-x-2 gap-y-0.5 rounded border border-stone-200 bg-white px-2 py-1 text-xs leading-4 md:grid-cols-5" data-testid="selection-panel">
      {compact   // a whole loop: the member list is folded (the concise loop summary is in the toolbar)
        ? <details data-testid="sel-highlightable-details"><summary className="cursor-pointer font-medium text-stone-700">Highlighted ({translation.highlightable.length})</summary><ul data-testid="sel-highlightable">{translation.highlightable.map((t) => <li key={t.label_seq_id}>{fmt(t)}</li>)}</ul></details>
        : <div><span className="font-medium text-stone-700">Highlighted ({translation.highlightable.length})</span><ul data-testid="sel-highlightable">{translation.highlightable.map((t) => <li key={t.label_seq_id}>{fmt(t)}</li>)}</ul></div>}
      <div><span className="font-medium text-stone-700">No coordinates ({translation.unobserved.length})</span><ul data-testid="sel-unobserved">{translation.unobserved.map((t) => <li key={t.label_seq_id}>{fmt(t)}</li>)}</ul></div>
      <div><span className="font-medium text-stone-700">Deleted in chain ({translation.deleted.length})</span><div data-testid="sel-deleted">{translation.deleted.join(', ')}</div></div>
      <div><span className="font-medium text-stone-700">Uncovered by alignment ({translation.uncovered.length})</span><div data-testid="sel-uncovered">{translation.uncovered.join(', ')}</div></div>
      <div><span className="font-medium text-stone-700">Instance-only residues ({instanceTargets.length})</span><ul data-testid="sel-instance">{instanceTargets.map((t) => <li key={t.label_seq_id}>{t.alignment_state}{t.insertion_anchor != null ? ` after ${t.insertion_anchor}.${t.insertion_offset}` : ''}: {t.nucleotide_parent}{t.label_seq_id} (auth {t.author_residue_number}){t.has_coordinates ? '' : ' no coords'}</li>)}</ul></div>
    </div>
  )
}

/** Click-open details of one insertion marker: every run with its residues (label_seq_id, author number + insertion code, coordinates). */
function InsertionDetails({ marker, mapping, anchor, activeRun, onView, onClose }: { marker: InsertionMarker; mapping: ChainMapping; anchor: DOMRect; activeRun: string | null; onView: (r: InsertionRun) => void; onClose: () => void }) {
  const pv = mapping.provenance
  return (
    <HoverCard anchor={anchor} interactive onDismiss={onClose} ignoreSelector="[data-insertion-marker]" testId="insertion-details" label="Insertion details">
      <div className="flex items-start justify-between gap-3">
        <span className="font-semibold">{marker.runs.length > 1 ? `${marker.runs.length} insertion runs · +${marker.length} nt` : `${marker.runs[0].kind === 'insert' ? 'Insertion' : 'Local-end residues'} · +${marker.length} nt`}</span>
        <button type="button" onClick={onClose} aria-label="Close insertion details" data-testid="insertion-details-close" className="rounded px-1 leading-none text-stone-500 hover:bg-stone-100">×</button>
      </div>
      {marker.runs.map((r) => (
        <div key={r.id} className="mt-1 border-t border-stone-100 pt-1" data-testid={`insertion-run-${r.id}`}>
          <div data-testid="insertion-run-boundary">{boundaryText(r)}</div>
          <div>Sequence: <span className="font-mono" data-testid="insertion-run-sequence">{r.sequence}</span> ({r.residues.length} nt, {mapping.pdb_id.toUpperCase()} chain {mapping.auth_asym_id})</div>
          <table className="mt-0.5 w-full text-[11px]" data-testid="insertion-run-residues">
            <thead className="text-stone-500"><tr><th className="pr-2 text-left font-normal">Nt</th><th className="pr-2 text-left font-normal">label_seq_id</th><th className="pr-2 text-left font-normal">Author no.</th><th className="text-left font-normal">Coordinates</th></tr></thead>
            <tbody>{r.residues.map((x) => <tr key={x.label_seq_id}><td className="pr-2 font-mono">{x.nucleotide_parent}{x.residue_name !== x.nucleotide_parent ? ` (${x.residue_name})` : ''}</td><td className="pr-2">{x.label_seq_id}</td><td className="pr-2">{authorLabel(x)}</td><td>{x.has_coordinates ? 'yes' : 'no'}</td></tr>)}</tbody>
          </table>
          {r.kind === 'local_end' && <div className="text-[10px] text-stone-500">Skipped by a local alignment end of the model, not an insert state.</div>}
          <button type="button" onClick={() => onView(r)} disabled={!r.residues.some((x) => x.has_coordinates)} data-testid={`insertion-view-${r.id}`} aria-pressed={activeRun === r.id}
            className={`mt-1 rounded border px-1.5 py-0.5 text-[11px] disabled:opacity-40 ${activeRun === r.id ? 'border-[#d1202f] text-[#d1202f]' : 'border-stone-300 text-stone-700'}`}>{activeRun === r.id ? 'Clear from 2D and 3D' : 'Show in 2D and 3D'}</button>
        </div>))}
      <div className="mt-1 border-t border-stone-100 pt-1 text-[10px] text-stone-500" data-testid="insertion-details-provenance">Rfam {pv.rfam_acc} {pv.rfam_release} covariance model, {pv.infernal_version.replace('INFERNAL', 'Infernal')} cmalign{pv.pipeline_version ? `, mapping pipeline ${pv.pipeline_version}` : ''}. Inserted nucleotides have no consensus position.</div>
    </HoverCard>
  )
}

/** "3D alignment" in Viewer details: the precomputed superposition of the selected instance onto the family reference (orientation.json). */
function AlignmentDetails({ transform: t, reference }: { transform: OrientationTransform | null; reference: OrientationDataset['reference'] | null }) {
  const atoms = t ? Object.entries(t.atoms_used).filter(([, n]) => n > 0).map(([a, n]) => `${n} ${a.replace("'", '′')}`).join(', ') : ''
  const rows: [string, string][] = !t ? [] : [
    ['Reference structure', reference ? `${reference.pdb_id.toUpperCase()} ${reference.auth_asym_id}${t.is_identity ? ' (this structure)' : ''}` : '–'],
    ['Superposition', isApplicable(t) ? 'Superposed on the reference' : 'Shown as deposited (fit insufficient)'],
    ...(t.rmsd != null ? [['RMSD', `${t.rmsd.toFixed(2)} Å`] as [string, string]] : []),
    ['Matched atom pairs', atoms ? `${t.n_pairs} (${atoms})` : String(t.n_pairs)],
    ['Alignment quality', t.fit_quality],
  ]
  return (
    <section data-testid="viewer-details-alignment" className="text-stone-600">
      <h4 className="mb-0.5 font-semibold text-stone-800">3D alignment</h4>
      {t ? <dl className="grid grid-cols-[auto_1fr] gap-x-3 gap-y-0">{rows.map(([k, v]) => <div key={k} className="contents"><dt className="text-stone-500">{k}</dt><dd data-testid={`alignment-${k.toLowerCase().replace(/ /g, '-')}`}>{v}</dd></div>)}</dl>
        : <p data-testid="alignment-none">No alignment information for this structure; shown as deposited.</p>}
    </section>
  )
}
