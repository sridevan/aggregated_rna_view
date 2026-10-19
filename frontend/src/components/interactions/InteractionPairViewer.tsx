import { useEffect, useRef, useState } from 'react'
import { applyRT, BACKBONE, cartoon, drawableNucleotides, load3Dmol, parseAtomSite, principal, residueLabel, toPdb, tubeScale, type Atom, type ViewerEntry } from '../structures/AlignmentViewer'
import { baseCorrespondence, kabsch } from '../../lib/interactions'

/*
 * Two RNA chains of the Interactions tab compare mode, with the Structures tab's 3D helpers (same 3Dmol build, ModelServer chain fetch,
 * orientation.json transforms into the family frame, principal-axes view, backbone tube). The Structures tab viewer itself is unchanged.
 *
 * Without a selected interaction both chains are shown in the family frame (the whole-RNA superposition of the Structures tab). With one,
 * its nucleotides are drawn as sticks over a faded tube and, when `fit` is given, chain B is moved by the least-squares fit of the base heavy
 * atoms of ONE nucleotide (the interaction's first nucleotide) onto A's: the partner nucleotide keeps its position relative to the fitted base,
 * so differences in base-pair or stacking geometry stay visible instead of being fitted away. When no valid fit exists (missing residue or
 * fewer than 3 corresponding atoms) the previous valid superposition is kept and `onFit` reports why. Changing the interaction never refetches
 * a chain (coordinates are cached per chain); only the scene is rebuilt from the cached atoms.
 */
type V3 = [number, number, number]
type Model = { getID: () => number }
type Viewer3 = { clear: () => void; addModel: (t: string, f: string) => Model; removeModel: (m: Model) => void; setStyle: (s: Record<string, unknown>, st: Record<string, unknown>) => void
  zoomTo: (s?: Record<string, unknown>) => void; zoom: (f: number) => void; render: () => void; resize: () => void; setHoverDuration: (ms: number) => void
  setHoverable: (sel: Record<string, unknown>, on: boolean, hover: (a: { model: number; serial: number }) => void, unhover: () => void) => void }
const MODEL_SERVER = 'https://www.ebi.ac.uk/pdbe/model-server/v1'
const atomCache = new Map<string, Promise<Atom[]>>()   // seq_name -> heavy atoms in the family frame
let fetches = 0
function chainAtoms(e: ViewerEntry): Promise<Atom[]> {
  let p = atomCache.get(e.seqName)
  if (!p) {
    fetches++
    p = fetch(`${MODEL_SERVER}/${e.pdbId.toLowerCase()}/atoms?label_asym_id=${encodeURIComponent(e.labelAsymId)}&model_nums=1&encoding=cif`).then(async (r) => {
      if (!r.ok) throw new Error(`${e.label}: HTTP ${r.status}`)
      const parsed = parseAtomSite(await r.text()).filter((a) => !a.het && a.el !== 'H')
      if (!parsed.length) throw new Error(`${e.label}: no atoms returned`)
      const T = e.transform?.rotation && e.transform.translation ? applyRT(e.transform) : null
      return T ? parsed.map((a) => { const [x, y, z] = T([a.x, a.y, a.z]); return { ...a, x, y, z } }) : parsed
    })
    atomCache.set(e.seqName, p); p.catch(() => atomCache.delete(e.seqName))
  }
  return p
}
const parentOf = (names: Set<string>) => (names.has('N9') ? (names.has('O6') ? 'G' : 'A') : names.has('N4') ? 'C' : 'U')
const residue = (atoms: Atom[], seq: number) => atoms.filter((a) => a.seq === seq)

export interface FitResult { state: 'none' | 'local' | 'kept' | 'failed'; message?: string; rmsd?: number; nAtoms?: number }
type RT = { R: number[][]; t: V3 }
const apply = (T: RT | null) => (p: V3): V3 => (T ? [0, 1, 2].map((i) => T.R[i][0] * p[0] + T.R[i][1] * p[1] + T.R[i][2] * p[2] + T.t[i]) as V3 : p)

/** Local fit of B onto A on the base atoms of one corresponding nucleotide. */
export function localFit(a: Atom[], b: Atom[]): { T: RT; rmsd: number; nAtoms: number } | { error: string } {
  if (!a.length || !b.length) return { error: 'the fitted nucleotide has no coordinates in one structure' }
  const na = new Set(a.map((x) => x.name)), nb = new Set(b.map((x) => x.name))
  const pairs = baseCorrespondence({ names: na, parent: parentOf(na) }, { names: nb, parent: parentOf(nb) })
  if (pairs.length < 3) return { error: `only ${pairs.length} corresponding base atoms` }
  const at = (xs: Atom[], n: string) => { const q = xs.find((x) => x.name === n)!; return [q.x, q.y, q.z] as V3 }
  const f = kabsch(pairs.map(([, y]) => at(b, y)), pairs.map(([x]) => at(a, x)))
  return f ? { T: { R: f.R, t: f.t }, rmsd: f.rmsd, nAtoms: pairs.length } : { error: 'fit failed' }
}

export function InteractionPairViewer({ entries, focus, onFit, height = 460 }: {
  entries: ViewerEntry[]
  /** selected interaction: label_seq_ids per structure (sticks) and, for the local superposition, the label_seq_id of the fitted nucleotide per structure */
  focus: { residues: Record<string, number[]>; fit: Record<string, number> | null } | null
  onFit?: (r: FitResult) => void
  height?: number
}) {
  const hostRef = useRef<HTMLDivElement>(null)
  const viewerRef = useRef<Viewer3 | null>(null)
  const [atoms, setAtoms] = useState<Record<string, Atom[]> | null>(null)
  const [status, setStatus] = useState<'loading' | 'ready' | 'error'>('loading')
  const [err, setErr] = useState<string | null>(null)
  const [fitState, setFitState] = useState<FitResult>({ state: 'none' })
  const fitRef = useRef<RT | null>(null)        // last valid local superposition of B (kept when a new one cannot be made)
  const shownRef = useRef(new Map<number, { e: ViewerEntry; atoms: Atom[] }>())
  const [hover, setHover] = useState<string | null>(null)
  const onFitRef = useRef(onFit); onFitRef.current = onFit
  const key = entries.map((e) => e.seqName).join(',')

  useEffect(() => () => { try { viewerRef.current?.clear() } catch { /* noop */ } viewerRef.current = null }, [])
  useEffect(() => {
    let alive = true
    setStatus('loading'); setAtoms(null); fitRef.current = null
    Promise.all(entries.map(chainAtoms)).then((list) => { if (alive) { setAtoms(Object.fromEntries(entries.map((e, k) => [e.seqName, list[k]]))); setErr(null) } })
      .catch((e: Error) => { if (alive) { setErr(e.message); setStatus('error') } })
    return () => { alive = false }
  }, [key])   // eslint-disable-line react-hooks/exhaustive-deps

  const focusKey = JSON.stringify(focus)
  useEffect(() => {
    if (!atoms) return
    let alive = true
    ;(async () => {
      const $3Dmol = await load3Dmol(); if (!alive || !hostRef.current) return
      viewerRef.current ??= $3Dmol.createViewer(hostRef.current, { backgroundColor: 'white' }) as unknown as Viewer3
      const v = viewerRef.current
      const [A, B] = entries
      // local superposition of B (or keep the previous valid one / the family frame)
      let res: FitResult = { state: 'none' }
      if (focus?.fit && B && focus.fit[A.seqName] != null && focus.fit[B.seqName] != null) {
        const f = localFit(residue(atoms[A.seqName], focus.fit[A.seqName]), residue(atoms[B.seqName], focus.fit[B.seqName]))
        if ('T' in f) { fitRef.current = f.T; res = { state: 'local', rmsd: f.rmsd, nAtoms: f.nAtoms } }
        else res = { state: fitRef.current ? 'kept' : 'failed', message: f.error }
      } else if (!focus) fitRef.current = null
      else if (fitRef.current) res = { state: 'kept', message: 'no corresponding nucleotide to fit in both structures' }
      const T = fitRef.current
      const place = (sn: string) => (B && sn === B.seqName ? atoms[sn].map((a) => { const [x, y, z] = apply(T)([a.x, a.y, a.z]); return { ...a, x, y, z } }) : atoms[sn])
      const placed = Object.fromEntries(entries.map((e) => [e.seqName, place(e.seqName)]))
      const o = principal(atoms[A.seqName].filter((a) => a.name === "C4'").map((a) => [a.x, a.y, a.z]))
      const view = o ? (p: V3): V3 => { const d = [p[0] - o.c[0], p[1] - o.c[1], p[2] - o.c[2]]; return [0, 1, 2].map((i) => o.R[i][0] * d[0] + o.R[i][1] * d[1] + o.R[i][2] * d[2]) as V3 } : undefined
      v.clear(); shownRef.current.clear(); setHover(null)
      const k = tubeScale(Math.max(...entries.map((e) => atoms[e.seqName].filter((a) => a.name === "C4'").length)))
      const sticks: Model[] = []
      for (const e of entries) {
        const bb = placed[e.seqName].filter((a) => BACKBONE.has(a.name))
        const m = v.addModel(toPdb(drawableNucleotides(bb), view), 'pdb'); shownRef.current.set(m.getID(), { e, atoms: bb })
        v.setStyle({ model: m }, cartoon(e.color, focus ? 'faded' : 'normal', k))
        const seqs = focus?.residues[e.seqName] ?? []
        const sel = placed[e.seqName].filter((a) => a.seq != null && seqs.includes(a.seq))
        if (sel.length) {
          const s = v.addModel(toPdb(sel, view), 'pdb'); shownRef.current.set(s.getID(), { e, atoms: sel }); sticks.push(s)
          v.setStyle({ model: s }, { stick: { color: e.color, radius: 0.22 } })
        }
      }
      v.setHoverDuration(80)
      v.setHoverable({}, true, (a) => { const m = shownRef.current.get(a.model), at = m?.atoms[a.serial - 1]; if (m && at) setHover(residueLabel(m.e.pdbId, at)) }, () => setHover(null))
      v.resize()
      if (sticks.length) { v.zoomTo({ model: sticks }); v.zoom(0.7) } else { v.zoomTo(); v.zoom(0.85) }
      v.render()
      // geometry of the selected nucleotides as displayed (test hook): C1'-C1' distance per structure, before and after the local fit
      const c1 = (xs: Atom[], s: number) => xs.find((a) => a.seq === s && a.name === "C1'")
      const geometry = Object.fromEntries(entries.map((e) => {
        const [s1, s2] = focus?.residues[e.seqName] ?? []; const a0 = s1 != null ? c1(atoms[e.seqName], s1) : undefined, b0 = s2 != null ? c1(atoms[e.seqName], s2) : undefined
        const a1 = s1 != null ? c1(placed[e.seqName], s1) : undefined, b1 = s2 != null ? c1(placed[e.seqName], s2) : undefined
        const d = (p?: Atom, q?: Atom) => (p && q ? Math.hypot(p.x - q.x, p.y - q.y, p.z - q.z) : null)
        return [e.seqName, { raw: d(a0, b0), shown: d(a1, b1) }]
      }))
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), interactionViewer: { loaded: entries.map((e) => e.seqName), fetches, fit: res, geometry,
        sticks: Object.fromEntries(entries.map((e) => [e.seqName, [...new Set(placed[e.seqName].filter((a) => a.seq != null && (focus?.residues[e.seqName] ?? []).includes(a.seq)).map((a) => a.seq))]])) } }
      setFitState(res); onFitRef.current?.(res); setStatus('ready')
    })().catch((e: Error) => { if (alive) { setErr(e.message); setStatus('error') } })
    return () => { alive = false }
  }, [atoms, focusKey])   // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <div className="relative rounded border border-stone-200 bg-white" style={{ height }} data-testid="interaction-pair-viewer" data-status={status} data-structures={key}
      data-fit={fitState.state} data-fit-rmsd={fitState.rmsd?.toFixed(3) ?? ''} data-fit-atoms={fitState.nAtoms ?? ''}>
      <div ref={hostRef} className="absolute inset-0" />
      {hover && <div className="pointer-events-none absolute bottom-2 left-2 rounded border border-stone-300 bg-white/95 px-1.5 py-0.5 font-mono text-[11px] text-stone-800 shadow" data-testid="pair-viewer-hover">{hover}</div>}
      <div className="pointer-events-none absolute right-2 top-2 flex flex-col items-end gap-0.5 text-[11px]">
        {entries.map((e, k) => <span key={e.seqName} className="rounded bg-white/90 px-1.5 py-px shadow-sm" data-testid={`pair-viewer-label-${k === 0 ? 'a' : 'b'}`}><span className="mr-1 inline-block h-2 w-2 rounded-full align-middle" style={{ background: e.color }} />{k === 0 ? 'A' : 'B'} · {e.label}</span>)}
      </div>
      {status === 'loading' && <div className="absolute left-2 top-2 rounded bg-white/90 px-2 py-0.5 text-xs text-stone-600">Loading…</div>}
      {err && <div className="absolute inset-x-2 bottom-2 rounded bg-red-50 px-2 py-1 text-xs text-red-800" data-testid="pair-viewer-error">{err}</div>}
    </div>
  )
}
