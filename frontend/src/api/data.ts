import type { ChainMapping, ConsensusStats, FamilyOverview, InstanceRow, OrientationDataset } from '../types/backend'

const BASE = `${import.meta.env.BASE_URL.replace(/\/$/, '')}/families`

async function getJson<T>(url: string): Promise<T> {
  const r = await fetch(url)
  if (!r.ok) throw new Error(`${r.status} ${r.statusText} for ${url}`)
  return (await r.json()) as T
}
export const loadFamily = (acc: string) => getJson<FamilyOverview>(`${BASE}/${acc}/family.json`)
export const loadInstances = (acc: string) => getJson<InstanceRow[]>(`${BASE}/${acc}/instances.json`)
export const loadStats = (acc: string) => getJson<ConsensusStats>(`${BASE}/${acc}/consensus_stats.json`)
/** Orientation dataset is optional: a family without it is shown untransformed. */
export const loadOrientation = (acc: string) => getJson<OrientationDataset>(`${BASE}/${acc}/orientation.json`).catch(() => null)

/** Chain mappings are fetched only when a chain is selected, and memoised. Nothing is prefetched at start-up. */
const mappingCache = new Map<string, Promise<ChainMapping>>()
export function loadMapping(acc: string, mappingFile: string): Promise<ChainMapping> {
  const key = `${acc}/${mappingFile}`
  let p = mappingCache.get(key)
  if (!p) {
    p = getJson<ChainMapping>(`${BASE}/${acc}/${mappingFile}`)
    mappingCache.set(key, p)
    p.catch(() => mappingCache.delete(key))
  }
  return p
}
export const mappingCacheSize = () => mappingCache.size
