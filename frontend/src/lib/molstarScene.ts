import type { OrientationTransform } from '../types/backend'

/*
 * Mol*-free helpers of the shared 3D superposition scene (components/molstar/MolstarSuperposition.tsx): entry and atom types, the residue
 * hover label, and the per-chain BinaryCIF cache (PDBe ModelServer, label_asym_id, model 1) with its fetch counter. Kept apart so modules
 * and unit tests that only need these never load Mol*.
 */
const MODEL_SERVER = 'https://www.ebi.ac.uk/pdbe/model-server/v1'

export interface SceneEntry { seqName: string; label: string; color: string; pdbId: string; labelAsymId: string; transform: OrientationTransform | null }
/** atoms of one residue in the family frame (orientation only, without any extra transform) */
export interface SceneAtom { name: string; x: number; y: number; z: number; seq: number; resn: string }
export interface SceneApi { residueAtoms: (seqName: string, labelSeqId: number) => SceneAtom[] }
export interface SceneHooks { loaded: string[]; fetches: number; transformed: boolean[]; probe: (seqName: string, i: number) => { x: number; y: number; label: string } | null }

/** Hover label of a residue: PDB id, author chain, residue name and author number (with insertion code), e.g. "5J7L.AA G230 (author)". */
export const residueLabel = (pdbId: string, a: { chain: string; resn: string; resi: string; icode: string }) => `${pdbId.toUpperCase()}.${a.chain} ${a.resn}${a.resi}${a.icode} (author)`

const bytesCache = new Map<string, Promise<Uint8Array<ArrayBuffer>>>()
let fetches = 0
export const chainUrl = (e: Pick<SceneEntry, 'pdbId' | 'labelAsymId'>) => `${MODEL_SERVER}/${e.pdbId.toLowerCase()}/atoms?label_asym_id=${encodeURIComponent(e.labelAsymId)}&model_nums=1&encoding=bcif`
/** BinaryCIF of one chain, fetched once per page session */
/** the ModelServer occasionally answers 502/503/504 or drops a connection: retry those twice with a short backoff (not 4xx) */
async function fetchRetry(url: string, label: string, tries = 3): Promise<Uint8Array<ArrayBuffer>> {
  for (let k = 1; ; k++) {
    let r: Response | null = null
    try { r = await fetch(url) } catch (err) { if (k >= tries) throw new Error(`${label}: ${String((err as Error)?.message ?? err)}`) }
    if (r?.ok) return new Uint8Array(await r.arrayBuffer())
    if (r && (r.status < 500 || k >= tries)) throw new Error(`${label}: HTTP ${r.status}`)
    await new Promise((res) => setTimeout(res, 800 * k))
  }
}
export function chainBytes(e: SceneEntry): Promise<Uint8Array<ArrayBuffer>> {
  const key = `${e.pdbId.toLowerCase()}|${e.labelAsymId}`
  let p = bytesCache.get(key)
  if (!p) {
    fetches++
    p = fetchRetry(chainUrl(e), e.label)
    bytesCache.set(key, p); p.catch(() => bytesCache.delete(key))
  }
  return p
}
export const sceneFetchCount = () => fetches

/** everything within `radius` Å of the chain (the chain itself, ligands, ions, water, other chains' residues): the source of "Show ligands" */
export const surroundingsUrl = (e: Pick<SceneEntry, 'pdbId' | 'labelAsymId'>, radius: number) =>
  `${MODEL_SERVER}/${e.pdbId.toLowerCase()}/residueSurroundings?label_asym_id=${encodeURIComponent(e.labelAsymId)}&radius=${radius}&model_nums=1&encoding=bcif`
const surroundCache = new Map<string, Promise<Uint8Array<ArrayBuffer>>>()
export function surroundingsBytes(e: SceneEntry, radius: number): Promise<Uint8Array<ArrayBuffer>> {
  const key = `${e.pdbId.toLowerCase()}|${e.labelAsymId}|${radius}`
  let p = surroundCache.get(key)
  if (!p) { p = fetchRetry(surroundingsUrl(e, radius), `${e.label} surroundings`); surroundCache.set(key, p); p.catch(() => surroundCache.delete(key)) }
  return p
}
