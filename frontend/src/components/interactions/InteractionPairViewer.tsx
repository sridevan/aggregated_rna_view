import { lazy, Suspense, useEffect, useMemo, useRef, useState } from 'react'
import { Mat4 } from 'molstar/lib/mol-math/linear-algebra'
import { sceneFetchCount, type SceneApi, type SceneAtom, type SceneEntry } from '../../lib/molstarScene'

const MolstarSuperposition = lazy(() => import('../molstar/MolstarSuperposition').then((m) => ({ default: m.MolstarSuperposition })))
import { baseCorrespondence, kabsch } from '../../lib/interactions'

/*
 * Two RNA chains of the Interactions tab compare mode in the shared Mol* scene (MolstarSuperposition: the Structures tab's chain loading,
 * orientation.json family frame, uniform cartoons; the Family Overview's engine and camera convention).
 *
 * Without a selected interaction both chains are shown in the family frame (the Structures tab's whole-RNA superposition). With one, its
 * nucleotides are drawn as sticks over faded cartoons and, when `fit` is given, chain B is moved by the least-squares fit of the base heavy
 * atoms of ONE nucleotide (the interaction's first nucleotide) onto A's, applied as an extra rigid transform after B's orientation: the
 * partner nucleotide keeps its position relative to the fitted base, so differences in base-pair or stacking geometry stay visible instead
 * of being fitted away. When no valid fit exists (missing residue or fewer than 3 corresponding atoms) the previous valid superposition is
 * kept and `onFit` reports why. Changing the interaction never refetches or re-parses a chain (only B's transform matrix is updated).
 */
type V3 = [number, number, number]
export interface FitResult { state: 'none' | 'local' | 'kept' | 'failed'; message?: string; rmsd?: number; nAtoms?: number }
type RT = { R: number[][]; t: V3 }
const parentOf = (names: Set<string>) => (names.has('N9') ? (names.has('O6') ? 'G' : 'A') : names.has('N4') ? 'C' : 'U')

/** Local fit of B onto A on the base atoms of one corresponding nucleotide (both in the family frame). */
export function localFit(a: SceneAtom[], b: SceneAtom[]): { T: RT; rmsd: number; nAtoms: number } | { error: string } {
  if (!a.length || !b.length) return { error: 'the fitted nucleotide has no coordinates in one structure' }
  const na = new Set(a.map((x) => x.name)), nb = new Set(b.map((x) => x.name))
  const pairs = baseCorrespondence({ names: na, parent: parentOf(na) }, { names: nb, parent: parentOf(nb) })
  if (pairs.length < 3) return { error: `only ${pairs.length} corresponding base atoms` }
  const at = (xs: SceneAtom[], n: string) => { const q = xs.find((x) => x.name === n)!; return [q.x, q.y, q.z] as V3 }
  const f = kabsch(pairs.map(([, y]) => at(b, y)), pairs.map(([x]) => at(a, x)))
  return f ? { T: { R: f.R, t: f.t }, rmsd: f.rmsd, nAtoms: pairs.length } : { error: 'fit failed' }
}
/** column-major 4x4 of x' = R x + t */
const toMat = (T: RT) => { const m = Mat4.identity(); for (let i = 0; i < 3; i++) { for (let j = 0; j < 3; j++) m[j * 4 + i] = T.R[i][j]; m[12 + i] = T.t[i] } return m }
const applyRT = (T: RT | null, a: SceneAtom): V3 => (T ? [0, 1, 2].map((i) => T.R[i][0] * a.x + T.R[i][1] * a.y + T.R[i][2] * a.z + T.t[i]) as V3 : [a.x, a.y, a.z])

export function InteractionPairViewer({ entries, focus, onFit, ligandClasses, height = 460 }: {
  entries: SceneEntry[]
  /** selected interaction: label_seq_ids per structure (sticks) and, for the local superposition, the label_seq_id of the fitted nucleotide per structure */
  focus: { residues: Record<string, number[]>; fit: Record<string, number> | null } | null
  onFit?: (r: FitResult) => void
  /** family ligand classification for Show ligands */
  ligandClasses?: () => Promise<Record<string, string>>
  height?: number
}) {
  const [api, setApi] = useState<SceneApi | null>(null)
  const [fitState, setFitState] = useState<FitResult>({ state: 'none' })
  const [extraT, setExtraT] = useState<RT | null>(null)        // last valid local superposition of B (kept when a new one cannot be made)
  const onFitRef = useRef(onFit); onFitRef.current = onFit
  const key = entries.map((e) => e.seqName).join(',')
  useEffect(() => { setExtraT(null) }, [key])
  const focusKey = JSON.stringify(focus)
  useEffect(() => {
    if (!api || entries.length < 2) return
    const [A, B] = entries
    let res: FitResult = { state: 'none' }, T = extraT
    if (focus?.fit && focus.fit[A.seqName] != null && focus.fit[B.seqName] != null) {
      const f = localFit(api.residueAtoms(A.seqName, focus.fit[A.seqName]), api.residueAtoms(B.seqName, focus.fit[B.seqName]))
      if ('T' in f) { T = f.T; res = { state: 'local', rmsd: f.rmsd, nAtoms: f.nAtoms } }
      else res = { state: T ? 'kept' : 'failed', message: f.error }
    } else if (!focus) T = null
    else if (T) res = { state: 'kept', message: 'no corresponding nucleotide to fit in both structures' }
    setExtraT(T); setFitState(res); onFitRef.current?.(res)
    // geometry of the selected nucleotides as displayed (test hook): C1'-C1' distance per structure, before and after the local fit
    const c1 = (sn: string, s?: number) => (s == null ? undefined : api.residueAtoms(sn, s).find((a) => a.name === "C1'"))
    const geometry = Object.fromEntries(entries.map((e) => {
      const [s1, s2] = focus?.residues[e.seqName] ?? []; const a0 = c1(e.seqName, s1), b0 = c1(e.seqName, s2), M = e.seqName === B.seqName ? T : null
      const d = (p?: V3, q?: V3) => (p && q ? Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]) : null)
      return [e.seqName, { raw: d(a0 && [a0.x, a0.y, a0.z], b0 && [b0.x, b0.y, b0.z]), shown: d(a0 && applyRT(M, a0), b0 && applyRT(M, b0)) }]
    }))
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), interactionViewer: { loaded: entries.map((e) => e.seqName), fetches: sceneFetchCount(), fit: res, geometry,
      sticks: Object.fromEntries(entries.map((e) => [e.seqName, (focus?.residues[e.seqName] ?? []).filter((s) => api.residueAtoms(e.seqName, s).length > 0)])) } }
  }, [api, focusKey, key])   // eslint-disable-line react-hooks/exhaustive-deps
  const extra = useMemo(() => (entries[1] && extraT ? { [entries[1].seqName]: toMat(extraT) } : {}), [entries, extraT])
  return (
    <Suspense fallback={<div className="relative rounded border border-stone-200 bg-white" style={{ height }} data-testid="interaction-pair-viewer" data-status="loading" data-structures={key} />}>
    <MolstarSuperposition entries={entries} extra={extra} focus={focus?.residues ?? null} onApi={setApi} ligandClasses={ligandClasses} testId="interaction-pair-viewer" style={{ height }}
      dataAttrs={{ 'data-fit': fitState.state, 'data-fit-rmsd': fitState.rmsd?.toFixed(3) ?? '', 'data-fit-atoms': String(fitState.nAtoms ?? '') }}>
      <div className="pointer-events-none absolute right-2 top-2 z-10 flex flex-col items-end gap-0.5 text-[11px]">
        {entries.map((e, k) => <span key={e.seqName} className="rounded bg-white/90 px-1.5 py-px shadow-sm" data-testid={`pair-viewer-label-${k === 0 ? 'a' : 'b'}`}><span className="mr-1 inline-block h-2 w-2 rounded-full align-middle" style={{ background: e.color }} />{k === 0 ? 'A' : 'B'} · {e.label}</span>)}
      </div>
    </MolstarSuperposition>
    </Suspense>
  )
}
