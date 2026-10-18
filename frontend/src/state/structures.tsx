import { createContext, useCallback, useContext, useMemo, useState, type ReactNode } from 'react'
import type { StructureMetricId } from '../types/backend'

/**
 * Shared state of the Structures tab (kept above the tabs so it survives switching tabs): the comparison metric, the selected
 * structural group PER METRIC (group ids are unrelated across metrics), the table mode, the representative comparison and the structures
 * shown in the 3D view. The selected RNA instance itself stays in state/selection.tsx (shared with the Family Overview).
 */
interface StructuresState {
  metric: StructureMetricId
  groupByMetric: Partial<Record<StructureMetricId, number>>
  tableMode: 'compact' | 'expanded'
  compareGroup: number | null
  displayed: string[] | null      // null = the group's default
  hidden: string[]
  featuresOpen: boolean          // Associated features panel; collapsed by default, kept across groups, metrics and tabs
}
interface Ctx extends StructuresState {
  setMetric: (m: StructureMetricId) => void
  selectGroup: (m: StructureMetricId, id: number) => void
  setTableMode: (t: 'compact' | 'expanded') => void
  setCompareGroup: (id: number | null) => void
  setDisplayed: (d: string[] | null) => void
  toggleHidden: (seqName: string) => void
  setFeaturesOpen: (open: boolean) => void
}
const StructuresContext = createContext<Ctx | null>(null)

export function StructuresProvider({ children }: { children: ReactNode }) {
  const [s, set] = useState<StructuresState>({ metric: 'discrepancy', groupByMetric: {}, tableMode: 'compact', compareGroup: null, displayed: null, hidden: [], featuresOpen: false })
  // switching metric keeps nothing that belonged to the other metric's groups: its own (or default largest) group, no comparison, default 3D content
  const setMetric = useCallback((m: StructureMetricId) => set((x) => (x.metric === m ? x : { ...x, metric: m, compareGroup: null, displayed: null, hidden: [] })), [])
  const selectGroup = useCallback((m: StructureMetricId, id: number) => set((x) => (x.groupByMetric[m] === id && x.metric === m ? x : { ...x, metric: m, groupByMetric: { ...x.groupByMetric, [m]: id }, compareGroup: null, displayed: null, hidden: [] })), [])
  const setTableMode = useCallback((t: 'compact' | 'expanded') => set((x) => ({ ...x, tableMode: t })), [])
  const setCompareGroup = useCallback((id: number | null) => set((x) => ({ ...x, compareGroup: id, displayed: null, hidden: [] })), [])
  const setDisplayed = useCallback((d: string[] | null) => set((x) => ({ ...x, displayed: d, hidden: x.hidden.filter((h) => d?.includes(h)) })), [])
  const toggleHidden = useCallback((sn: string) => set((x) => ({ ...x, hidden: x.hidden.includes(sn) ? x.hidden.filter((h) => h !== sn) : [...x.hidden, sn] })), [])
  const setFeaturesOpen = useCallback((open: boolean) => set((x) => ({ ...x, featuresOpen: open })), [])
  const value = useMemo(() => ({ ...s, setMetric, selectGroup, setTableMode, setCompareGroup, setDisplayed, toggleHidden, setFeaturesOpen }), [s, setMetric, selectGroup, setTableMode, setCompareGroup, setDisplayed, toggleHidden, setFeaturesOpen])
  return <StructuresContext.Provider value={value}>{children}</StructuresContext.Provider>
}

export function useStructures(): Ctx {
  const c = useContext(StructuresContext)
  if (!c) throw new Error('useStructures outside StructuresProvider')
  return c
}
