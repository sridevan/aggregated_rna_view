import { createContext, useContext, useMemo, useState, type ReactNode } from 'react'
import { useFamily } from './family'
import { useMapping } from './useMapping'
import { buildNumbering, type Numbering } from '../lib/numbering'

/** Position numbering (display only), shared by the Family Overview and Annotations tabs: Rfam consensus (null) or one chain's author numbering.
 *  Independent of the active instance; the reference chain's validated mapping is loaded once (memoised loader) and kept across tabs. */
interface NumberingApi { refSeq: string | null; setRef: (seqName: string | null) => void; numbering: Numbering | null; state: 'rfam' | 'loading' | 'ready' | 'error' }
const Ctx = createContext<NumberingApi | null>(null)

export function NumberingProvider({ children }: { children: ReactNode }) {
  const { acc, family, instances } = useFamily()
  const [refSeq, setRef] = useState<string | null>(null)
  const inst = refSeq ? instances.find((i) => i.seq_name === refSeq) ?? null : null
  const ms = useMapping(acc, inst)
  const mapping = ms.state === 'ready' ? ms.mapping : null
  const numbering = useMemo(() => (inst && mapping ? buildNumbering(inst.seq_name, inst.pdb_id, inst.auth_asym_id, mapping, family.clen) : null), [inst, mapping, family.clen])
  const state = !inst ? 'rfam' : numbering ? 'ready' : ms.state === 'error' ? 'error' : 'loading'
  const api = useMemo(() => ({ refSeq, setRef, numbering, state }) as NumberingApi, [refSeq, numbering, state])
  return <Ctx.Provider value={api}>{children}</Ctx.Provider>
}
export function useNumbering(): NumberingApi {
  const v = useContext(Ctx)
  if (!v) throw new Error('useNumbering outside NumberingProvider')
  return v
}

/** The "Position numbering" selector (same control in both tabs). */
export function NumberingSelect() {
  const { instances } = useFamily()
  const { refSeq, setRef, state } = useNumbering()
  return (
    <label className="flex items-center gap-1 text-[11px] text-stone-700" data-testid="numbering-control" data-state={state}>
      Position numbering
      <select data-testid="numbering-select" value={refSeq ?? ''} onChange={(e) => setRef(e.target.value || null)}
        className="max-w-[16rem] rounded border border-stone-300 bg-white px-1 py-px text-[11px]" title="Changes residue labels only; alignment, aggregation and selections use Rfam consensus positions.">
        <option value="">Rfam consensus</option>
        {instances.map((i) => <option key={i.seq_name} value={i.seq_name}>{i.pdb_id.toUpperCase()} {i.auth_asym_id} · author numbering</option>)}
      </select>
      {state === 'loading' && <span className="text-stone-500">loading…</span>}
      {state === 'error' && <span className="text-amber-700">mapping unavailable</span>}
    </label>
  )
}
