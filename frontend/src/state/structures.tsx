import { createContext, useCallback, useContext, useMemo, useState, type ReactNode } from 'react'
import type { StructureMetricId } from '../types/backend'

/**
 * Shared state of the Structures tab (kept above the tabs so it survives switching tabs): the comparison metric, the selected
 * structural group PER METRIC (group ids are unrelated across metrics), the table mode, the representative comparison and the structures
 * shown in the 3D view. The selected RNA instance itself stays in state/selection.tsx (shared with the Family Overview).
 */
export type StructuresView = 'clustering' | 'dissimilarity'
interface StructuresState {
  view: StructuresView             // Structural Dissimilarity (default) or Clustering mode; kept while the family page is open
  metric: StructureMetricId
  /** Structural Dissimilarity mode's own 3D content (separate from the Clustering overlay): null = the family representative alone, else one instance or a heatmap pair */
  dissimDisplayed: string[] | null
  groupByMetric: Partial<Record<StructureMetricId, number>>
  tableMode: 'compact' | 'expanded'
  compareGroup: number | null
  displayed: string[] | null      // Clustering mode: null = the selected group's representative, else the one instance shown instead
  featuresOpen: boolean          // Associated features panel; collapsed by default, kept across groups, metrics and tabs
}
interface Ctx extends StructuresState {
  setView: (v: StructuresView) => void
  setDissimDisplayed: (d: string[] | null) => void
  setMetric: (m: StructureMetricId) => void
  selectGroup: (m: StructureMetricId, id: number) => void
  setTableMode: (t: 'compact' | 'expanded') => void
  setCompareGroup: (id: number | null) => void
  setDisplayed: (d: string[] | null) => void
  setFeaturesOpen: (open: boolean) => void
}
const StructuresContext = createContext<Ctx | null>(null)

export function StructuresProvider({ children }: { children: ReactNode }) {
  const [s, set] = useState<StructuresState>({ view: 'dissimilarity', metric: 'discrepancy', dissimDisplayed: null, groupByMetric: {}, tableMode: 'compact', compareGroup: null, displayed: null, featuresOpen: false })
  // switching metric keeps nothing that belonged to the other metric's groups: its own (or default largest) group, no comparison, default 3D content
  const setView = useCallback((v: StructuresView) => set((x) => (x.view === v ? x : { ...x, view: v })), [])
  const setDissimDisplayed = useCallback((d: string[] | null) => set((x) => ({ ...x, dissimDisplayed: d && d.length ? d.slice(0, 2) : null })), [])   // at most two structures
  const setMetric = useCallback((m: StructureMetricId) => set((x) => (x.metric === m ? x : { ...x, metric: m, compareGroup: null, displayed: null })), [])
  const selectGroup = useCallback((m: StructureMetricId, id: number) => set((x) => (x.groupByMetric[m] === id && x.metric === m ? x : { ...x, metric: m, groupByMetric: { ...x.groupByMetric, [m]: id }, compareGroup: null, displayed: null })), [])
  const setTableMode = useCallback((t: 'compact' | 'expanded') => set((x) => ({ ...x, tableMode: t })), [])
  const setCompareGroup = useCallback((id: number | null) => set((x) => ({ ...x, compareGroup: id, displayed: null })), [])
  const setDisplayed = useCallback((d: string[] | null) => set((x) => ({ ...x, displayed: d && d.length ? d.slice(0, 1) : null })), [])   // one instance; a pair only via compareGroup
  const setFeaturesOpen = useCallback((open: boolean) => set((x) => ({ ...x, featuresOpen: open })), [])
  const value = useMemo(() => ({ ...s, setView, setDissimDisplayed, setMetric, selectGroup, setTableMode, setCompareGroup, setDisplayed, setFeaturesOpen }), [s, setView, setDissimDisplayed, setMetric, selectGroup, setTableMode, setCompareGroup, setDisplayed, setFeaturesOpen])
  return <StructuresContext.Provider value={value}>{children}</StructuresContext.Provider>
}

export function useStructures(): Ctx {
  const c = useContext(StructuresContext)
  if (!c) throw new Error('useStructures outside StructuresProvider')
  return c
}
