import { createContext, useContext, useEffect, useState, type ReactNode } from 'react'
import { loadFamily, loadInstances, loadOrientation, loadStats } from '../api/data'
import type { ConsensusStats, FamilyOverview, InstanceRow, OrientationDataset } from '../types/backend'

export interface FamilyData { acc: string; family: FamilyOverview; instances: InstanceRow[]; stats: ConsensusStats; orientation: OrientationDataset | null }
export type FamilyStatus = { state: 'loading' } | { state: 'error'; message: string } | { state: 'ready'; data: FamilyData }
const Ctx = createContext<FamilyStatus>({ state: 'loading' })

export function FamilyProvider({ acc, children }: { acc: string; children: ReactNode }) {
  const [status, setStatus] = useState<FamilyStatus>({ state: 'loading' })
  useEffect(() => {
    let alive = true
    Promise.all([loadFamily(acc), loadInstances(acc), loadStats(acc), loadOrientation(acc)])
      .then(([family, instances, stats, orientation]) => { if (alive) setStatus({ state: 'ready', data: { acc, family, instances, stats, orientation } }) })
      .catch((e: Error) => { if (alive) setStatus({ state: 'error', message: e.message }) })
    return () => { alive = false }
  }, [acc])
  return <Ctx.Provider value={status}>{children}</Ctx.Provider>
}
export const useFamilyStatus = () => useContext(Ctx)
export function useFamily(): FamilyData {
  const s = useContext(Ctx)
  if (s.state !== 'ready') throw new Error('family not loaded')
  return s.data
}
