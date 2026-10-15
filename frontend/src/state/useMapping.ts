import { useEffect, useState } from 'react'
import { loadMapping } from '../api/data'
import type { ChainMapping, InstanceRow } from '../types/backend'

export type MappingStatus = { state: 'idle' } | { state: 'loading' } | { state: 'error'; message: string } | { state: 'ready'; mapping: ChainMapping }
/** Lazily fetch the mapping of one chain instance (memoised in api/data). */
export function useMapping(acc: string, instance: InstanceRow | null): MappingStatus {
  const [status, setStatus] = useState<MappingStatus>({ state: 'idle' })
  useEffect(() => {
    if (!instance) { setStatus({ state: 'idle' }); return }
    let alive = true
    setStatus({ state: 'loading' })
    loadMapping(acc, instance.mapping_file)
      .then((mapping) => { if (alive) setStatus({ state: 'ready', mapping }) })
      .catch((e: Error) => { if (alive) setStatus({ state: 'error', message: e.message }) })
    return () => { alive = false }
  }, [acc, instance])
  return status
}
