import { createContext, useCallback, useContext, useMemo, useReducer, type ReactNode } from 'react'

/**
 * Authoritative selection state shared by the consensus tracks, the 2D viewer and the 3D viewer of the Family Overview.
 * Exactly one selection is active at a time (or none); the selected RNA instance is kept separately and never changes the selection.
 *
 * Two selection kinds:
 *  - `position`: one Rfam consensus position (family coordinate, 1..CLEN). It is preserved across instance changes and re-resolved
 *    through each instance's mapping; when the instance lacks the residue (deleted, uncovered, no coordinates) nothing is highlighted.
 *  - `residue`: one nucleotide of a specific instance that has no consensus coordinate (insert, local end, flank). Instance-bound.
 *  - `group`: an annotation group of nucleotides (e.g. a BGSU RNA 3D Motif Atlas loop: HL/IL/J3/J4) with its exact, possibly
 *    discontinuous membership: consensus positions plus instance-bound residues without consensus positions. Prepared for the
 *    loop tracks of a later milestone; nothing creates groups yet.
 */
export interface GroupSelection {
  kind: 'group'
  /** unique identifier (e.g. BGSU loop id HL_1Y26_001) */
  id: string
  label: string
  /** e.g. 'HL' | 'IL' | 'J3' | 'J4' */
  type: string
  /** source PDB entry and chain instance the group was observed on */
  source: { pdb_id: string; seq_name: string } | null
  /** exact membership, in source order; original unit ids kept for provenance */
  members: { unit_id: string; rfam_consensus_pos: number | null; seq_name: string; label_seq_id: number }[]
  /** segment boundaries as [startIndex, endIndex] into `members` (discontinuous strands keep their order) */
  segments: [number, number][]
}
export type ActiveSelection = { kind: 'position'; pos: number } | { kind: 'residue'; seq_name: string; label_seq_id: number } | GroupSelection
export interface SelectionState {
  active: ActiveSelection | null; instance: string | null
  /** incremented on every explicit deselection (click on the active nucleotide or loop again, Clear selection): the viewers restore their default 2D/3D views */
  viewReset: number
}
export type Action = { type: 'togglePosition'; pos: number } | { type: 'selectPosition'; pos: number } | { type: 'clear' } | { type: 'instance'; seqName: string | null }
  | { type: 'toggleResidue'; seqName: string; labelSeqId: number } | { type: 'group'; group: GroupSelection } | { type: 'toggleGroup'; group: GroupSelection }

/** exported for unit tests */
export function reducer(s: SelectionState, a: Action): SelectionState {
  switch (a.type) {
    case 'togglePosition': return s.active?.kind === 'position' && s.active.pos === a.pos ? { ...s, active: null, viewReset: s.viewReset + 1 } : { ...s, active: { kind: 'position', pos: a.pos } }
    case 'selectPosition': return { ...s, active: { kind: 'position', pos: a.pos } }
    case 'clear': return { ...s, active: null, viewReset: s.viewReset + 1 }
    case 'instance': {
      if (a.seqName === s.instance) return s
      // a consensus position (or a group with consensus members) survives an instance change; an instance-bound residue does not
      // a loop (group) belongs to one instance and is cleared; a group without a source instance would be kept
      const keep = s.active && (s.active.kind === 'position' || (s.active.kind === 'group' && (!s.active.source || s.active.source.seq_name === a.seqName)))
      return { ...s, instance: a.seqName, active: keep ? s.active : null }
    }
    case 'toggleResidue': return s.active?.kind === 'residue' && s.active.seq_name === a.seqName && s.active.label_seq_id === a.labelSeqId ? { ...s, active: null, viewReset: s.viewReset + 1 } : { ...s, active: { kind: 'residue', seq_name: a.seqName, label_seq_id: a.labelSeqId } }
    case 'group': return { ...s, active: a.group }
    // click on the active group (same id, same source instance) again clears it and restores the default views; any other group replaces it
    case 'toggleGroup': return s.active?.kind === 'group' && s.active.id === a.group.id && s.active.source?.seq_name === a.group.source?.seq_name ? { ...s, active: null, viewReset: s.viewReset + 1 } : { ...s, active: a.group }
  }
}

export interface SelectionApi extends SelectionState {
  /** consensus positions of the active selection (0 or 1 for a position, the members' positions for a group) */
  positions: ReadonlySet<number>
  /** label_seq_ids of residues of the *current* instance that belong to the active selection but have no consensus coordinate */
  instanceResidues: ReadonlySet<number>
  /** click semantics: select the position, or clear it when it is already the selection */
  toggle: (pos: number) => void
  select: (pos: number) => void
  clear: () => void
  setInstance: (seqName: string | null) => void
  /** select / deselect one instance-bound residue without a consensus coordinate (insert, local end, flank) */
  toggleInstanceResidue: (labelSeqId: number) => void
  selectGroup: (group: GroupSelection) => void
  /** select a group, or clear it when it is already the active selection */
  toggleGroup: (group: GroupSelection) => void
}
const Ctx = createContext<SelectionApi | null>(null)
export function SelectionProvider({ children, initialInstance = null }: { children: ReactNode; initialInstance?: string | null }) {
  const [s, dispatch] = useReducer(reducer, { active: null, instance: initialInstance, viewReset: 0 })
  const toggle = useCallback((pos: number) => dispatch({ type: 'togglePosition', pos }), [])
  const select = useCallback((pos: number) => dispatch({ type: 'selectPosition', pos }), [])
  const clear = useCallback(() => dispatch({ type: 'clear' }), [])
  const setInstance = useCallback((seqName: string | null) => dispatch({ type: 'instance', seqName }), [])
  const selectGroup = useCallback((group: GroupSelection) => dispatch({ type: 'group', group }), [])
  const toggleGroup = useCallback((group: GroupSelection) => dispatch({ type: 'toggleGroup', group }), [])
  const toggleInstanceResidue = useCallback((labelSeqId: number) => { if (s.instance) dispatch({ type: 'toggleResidue', seqName: s.instance, labelSeqId }) }, [s.instance])
  const api = useMemo<SelectionApi>(() => {
    const a = s.active
    // group: only members of the current instance give consensus positions (never translate another chain's member onto this chain)
    const positions = new Set<number>(a?.kind === 'position' ? [a.pos] : a?.kind === 'group' ? a.members.filter((m) => m.seq_name === s.instance).map((m) => m.rfam_consensus_pos).filter((p): p is number => p != null) : [])
    const instanceResidues = new Set<number>(a?.kind === 'residue' && a.seq_name === s.instance ? [a.label_seq_id]
      : a?.kind === 'group' ? a.members.filter((m) => m.rfam_consensus_pos == null && m.seq_name === s.instance).map((m) => m.label_seq_id) : [])
    return { ...s, positions, instanceResidues, toggle, select, clear, setInstance, toggleInstanceResidue, selectGroup, toggleGroup }
  }, [s, toggle, select, clear, setInstance, toggleInstanceResidue, selectGroup, toggleGroup])
  return <Ctx.Provider value={api}>{children}</Ctx.Provider>
}
export function useSelection(): SelectionApi {
  const v = useContext(Ctx)
  if (!v) throw new Error('useSelection outside SelectionProvider')
  return v
}
