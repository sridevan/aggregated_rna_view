import { createContext, useCallback, useContext, useMemo, useReducer, type ReactNode } from 'react'

/**
 * Authoritative selection state shared by the consensus tracks, the 2D viewer and the 3D viewer (and a future 2D family viewer).
 * - `positions`: Rfam consensus positions (1..CLEN). Family coordinates; independent of the selected instance and preserved across instance changes.
 * - `instanceResidues`: label_seq_ids of residues of the *current* instance that have no consensus coordinate (insert, local end, flank).
 *   They are selectable at instance level only and are cleared whenever the instance changes.
 */
export interface SelectionState { positions: ReadonlySet<number>; anchor: number | null; instance: string | null; instanceResidues: ReadonlySet<number> }
type Action =
  | { type: 'toggle'; pos: number } | { type: 'range'; from: number; to: number } | { type: 'set'; positions: Iterable<number> }
  | { type: 'clear' } | { type: 'instance'; seqName: string | null } | { type: 'toggleInstanceResidue'; labelSeqId: number }

function reducer(s: SelectionState, a: Action): SelectionState {
  switch (a.type) {
    case 'toggle': { const p = new Set(s.positions); if (p.has(a.pos)) p.delete(a.pos); else p.add(a.pos); return { ...s, positions: p, anchor: a.pos } }
    case 'range': { const p = new Set(s.positions); const [lo, hi] = a.from <= a.to ? [a.from, a.to] : [a.to, a.from]; for (let i = lo; i <= hi; i++) p.add(i); return { ...s, positions: p, anchor: a.to } }
    case 'set': return { ...s, positions: new Set(a.positions), anchor: null }
    case 'clear': return { ...s, positions: new Set(), anchor: null, instanceResidues: new Set() }
    case 'instance': return a.seqName === s.instance ? s : { ...s, instance: a.seqName, instanceResidues: new Set() }   // positions deliberately untouched
    case 'toggleInstanceResidue': { const r = new Set(s.instanceResidues); if (r.has(a.labelSeqId)) r.delete(a.labelSeqId); else r.add(a.labelSeqId); return { ...s, instanceResidues: r } }
  }
}
export interface SelectionApi extends SelectionState {
  toggle: (pos: number) => void; range: (from: number, to: number) => void; set: (p: Iterable<number>) => void; clear: () => void
  setInstance: (seqName: string | null) => void; toggleInstanceResidue: (labelSeqId: number) => void
}
const Ctx = createContext<SelectionApi | null>(null)
export function SelectionProvider({ children, initialInstance = null }: { children: ReactNode; initialInstance?: string | null }) {
  const [s, dispatch] = useReducer(reducer, { positions: new Set<number>(), anchor: null, instance: initialInstance, instanceResidues: new Set<number>() })
  const toggle = useCallback((pos: number) => dispatch({ type: 'toggle', pos }), [])
  const range = useCallback((from: number, to: number) => dispatch({ type: 'range', from, to }), [])
  const set = useCallback((positions: Iterable<number>) => dispatch({ type: 'set', positions }), [])
  const clear = useCallback(() => dispatch({ type: 'clear' }), [])
  const setInstance = useCallback((seqName: string | null) => dispatch({ type: 'instance', seqName }), [])
  const toggleInstanceResidue = useCallback((labelSeqId: number) => dispatch({ type: 'toggleInstanceResidue', labelSeqId }), [])
  const api = useMemo(() => ({ ...s, toggle, range, set, clear, setInstance, toggleInstanceResidue }), [s, toggle, range, set, clear, setInstance, toggleInstanceResidue])
  return <Ctx.Provider value={api}>{children}</Ctx.Provider>
}
export function useSelection(): SelectionApi {
  const v = useContext(Ctx)
  if (!v) throw new Error('useSelection outside SelectionProvider')
  return v
}
