import { useEffect, useState } from 'react'
import { loadLigands, loadModificationTracks, loadPocketTracks, loadSasaTracks, loadTracks } from '../api/data'
import type { LigandsDataset, ModificationTracksDataset, PocketTracksDataset, SasaTracksDataset, TracksDataset } from '../types/backend'

export type AnnotationStatus = { state: 'loading' } | { state: 'error'; message: string } | { state: 'ready'; tracks: TracksDataset; ligands: LigandsDataset; sasa: SasaTracksDataset | null; pockets: PocketTracksDataset | null; modifications: ModificationTracksDataset | null }
/** Lazily fetch tracks.json and ligands.json for the Annotations tab (memoised in api/data). */
export function useAnnotations(acc: string): AnnotationStatus {
  const [status, setStatus] = useState<AnnotationStatus>({ state: 'loading' })
  useEffect(() => {
    let alive = true
    setStatus({ state: 'loading' })
    Promise.all([loadTracks(acc), loadLigands(acc), loadSasaTracks(acc), loadPocketTracks(acc), loadModificationTracks(acc)])
      .then(([tracks, ligands, sasa, pockets, modifications]) => { if (alive) setStatus({ state: 'ready', tracks, ligands, sasa, pockets, modifications }) })
      .catch((e: Error) => { if (alive) setStatus({ state: 'error', message: e.message }) })
    return () => { alive = false }
  }, [acc])
  return status
}
