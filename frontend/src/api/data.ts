import type { ChainMapping, ConsensusStats, FamilyOverview, InstanceRow, LigandsDataset, OrientationDataset, PocketTracksDataset, SasaTracksDataset, TracksDataset, ModificationTracksDataset, SequenceDifferenceTracksDataset, StructureClusteringDataset, ModificationsDataset, LoopsDataset, InteractionsDataset } from '../types/backend'

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

/** Annotation datasets (Milestone 3B) are loaded only when the Annotations tab is opened, and memoised. */
const annotationCache = new Map<string, Promise<unknown>>()
function memo<T>(key: string, make: () => Promise<T>): Promise<T> {
  let p = annotationCache.get(key) as Promise<T> | undefined
  if (!p) { p = make(); annotationCache.set(key, p); p.catch(() => annotationCache.delete(key)) }
  return p
}
export const loadTracks = (acc: string) => memo(`${acc}/tracks`, () => getJson<TracksDataset>(`${BASE}/${acc}/tracks.json`))
export const loadLigands = (acc: string) => memo(`${acc}/ligands`, () => getJson<LigandsDataset>(`${BASE}/${acc}/ligands.json`))
/** Optional dataset: a family without SASA shows no Solvent accessibility category. */
export const loadSasaTracks = (acc: string) => memo(`${acc}/sasa_tracks`, () => getJson<SasaTracksDataset>(`${BASE}/${acc}/sasa_tracks.json`).catch(() => null))

/** URL of a file of the family data directory (used for download links). */
export const familyFileUrl = (acc: string, file: string) => `${BASE}/${acc}/${file}`
/** Optional dataset: a family without nonstandard-nucleotide data shows no Nonstandard nucleotides category. */
export const loadSeqDiffTracks = (acc: string) => memo(`${acc}/sequence_difference_tracks`, () => getJson<SequenceDifferenceTracksDataset>(`${BASE}/${acc}/sequence_difference_tracks.json`).catch(() => null))
export const loadModificationTracks = (acc: string) => memo(`${acc}/modification_tracks`, () => getJson<ModificationTracksDataset>(`${BASE}/${acc}/modification_tracks.json`).catch(() => null))
/** Optional dataset: a family without pocket predictions shows no Predicted binding pockets category. */
export const loadPocketTracks = (acc: string) => memo(`${acc}/pocket_tracks`, () => getJson<PocketTracksDataset>(`${BASE}/${acc}/pocket_tracks.json`).catch(() => null))

/** Optional datasets of the Structures tab: a family without structure_clustering.json shows the table without structural groups. */
export const loadStructureClustering = (acc: string) => memo(`${acc}/structure_clustering`, () => getJson<StructureClusteringDataset>(`${BASE}/${acc}/structure_clustering.json`).catch(() => null))
export const loadModifications = (acc: string) => memo(`${acc}/modifications`, () => getJson<ModificationsDataset>(`${BASE}/${acc}/modifications.json`).catch(() => null))
/** BGSU loops with complete membership (Family Overview loop tracks); loaded on first use. */
export const loadLoops = (acc: string) => memo(`${acc}/loops`, () => getJson<LoopsDataset>(`${BASE}/${acc}/loops.json`))
/** Loops touching one chain: its own small file (loops_by_instance/<seq_name>.json); families annotated before that file existed fall back to
 *  the family-wide loops.json. */
export const loadInstanceLoops = (acc: string, seqName: string) => memo(`${acc}/loops/${seqName}`, () =>
  getJson<{ loops: LoopsDataset['loops'] }>(`${BASE}/${acc}/loops_by_instance/${seqName}.json`).then((d) => d.loops).catch(() => loadLoops(acc).then((d) => d.loops)))
/** Optional dataset of the Interactions tab (interactions.json); loaded when the tab is opened. */
export const loadInteractions = (acc: string) => memo(`${acc}/interactions`, () => getJson<InteractionsDataset>(`${BASE}/${acc}/interactions.json`).catch(() => null))
/** CCD id -> ligand class of a family (ligands.json by_ccd), for the Mol* "Show ligands" filter; loaded only when first needed. */
export const loadLigandClasses = (acc: string) => loadLigands(acc).then((d) => Object.fromEntries(Object.entries(d.by_ccd).map(([k, v]) => [k, v.ligand_class])))
