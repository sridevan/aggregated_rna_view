import { useEffect, useRef, useState } from 'react'
import type { OrientationTransform } from '../../types/backend'

/*
 * 3D alignment view with 3Dmol.js, ported from antibody_antigen_interface/app-react/src/states/SuperpositionViewer.jsx (pinned 3Dmol 2.5.2
 * from the CDN, as in the demo; even oval ribbon per structure; one principal-axes rotation for the whole scene; Fullscreen and Centre
 * view; cached coordinate text). RNA changes: each chain is fetched as mmCIF from the PDBe ModelServer (label_asym_id, model 1) and moved
 * into the family frame with its existing orientation.json transform (Kabsch on C4' onto the family reference, x' = R x + t) before
 * loading. That superposition is for visual comparison; it is not the base-centre fit behind the geometric discrepancy.
 */
const THREEDMOL_URL = 'https://cdn.jsdelivr.net/npm/3dmol@2.5.2/build/3Dmol-min.js'
const MODEL_SERVER = 'https://www.ebi.ac.uk/pdbe/model-server/v1'
type ThreeDmol = { createViewer: (el: HTMLElement, o: Record<string, unknown>) => Viewer }
type HoverAtom = { model: number; serial: number }
type Viewer = { clear: () => void; addModel: (t: string, f: string) => { getID: () => number };
  setHoverable: (sel: Record<string, unknown>, on: boolean, hover: (a: HoverAtom) => void, unhover: (a: HoverAtom) => void) => void; setHoverDuration: (ms: number) => void;
  modelToScreen: (p: { x: number; y: number; z: number }) => { x: number; y: number }; setStyle: (s: Record<string, unknown>, st: Record<string, unknown>) => void; zoomTo: () => void; zoom: (f: number) => void; render: () => void; resize: () => void }
declare global { interface Window { $3Dmol?: ThreeDmol } }

let loader: Promise<ThreeDmol> | null = null
export function load3Dmol(): Promise<ThreeDmol> {
  if (window.$3Dmol) return Promise.resolve(window.$3Dmol)
  loader ??= new Promise<ThreeDmol>((res, rej) => {
    const s = document.createElement('script'); s.src = THREEDMOL_URL; s.async = true
    s.onload = () => (window.$3Dmol ? res(window.$3Dmol) : rej(new Error('3Dmol.js did not initialise')))
    s.onerror = () => { loader = null; rej(new Error('3Dmol.js could not be loaded from the CDN (needs internet)')) }
    document.head.appendChild(s)
  })
  return loader
}

export interface ViewerEntry { seqName: string; label: string; color: string; pdbId: string; labelAsymId: string; transform: OrientationTransform | null }

/** One atom of the chain (from the mmCIF atom_site loop of the PDBe ModelServer response). */
export interface Atom { het: boolean; name: string; resn: string; chain: string; resi: string; icode: string; x: number; y: number; z: number; occ: number; b: number; el: string }

/** Parses the atom_site loop by its column names (first alternate location only). 3Dmol.js 2.5.2 does not read the ModelServer mmCIF, so
 *  the atoms are re-written as PDB records for it (toPdb). Polymer residues are written as ATOM (also nonstandard ones) so the ribbon is continuous. */
export function parseAtomSite(cif: string): Atom[] {
  const lines = cif.split('\n'); const cols: string[] = []; const out: Atom[] = []; let inLoop = false
  for (const raw of lines) {
    const l = raw.trim()
    if (l.startsWith('_atom_site.')) { cols.push(l); inLoop = true; continue }
    if (!inLoop) continue
    if (l.startsWith('ATOM') || l.startsWith('HETATM')) {
      const t = l.split(/\s+/); const g = (k: string) => { const i = cols.indexOf(`_atom_site.${k}`); return i < 0 ? undefined : t[i] }
      const alt = g('label_alt_id'); if (alt && alt !== '.' && alt !== '?' && alt !== 'A') continue
      const unq = (v: string | undefined) => { const x = v ?? ''; return x.length > 1 && (x[0] === '"' || x[0] === "'") && x[x.length - 1] === x[0] ? x.slice(1, -1) : x }   // only enclosing quotes (C4' keeps its prime)
      const seq = g('label_seq_id')
      out.push({ het: seq === '.' || seq === undefined, name: unq(g('auth_atom_id') ?? g('label_atom_id')), resn: unq(g('auth_comp_id') ?? g('label_comp_id')), chain: g('auth_asym_id') ?? g('label_asym_id') ?? 'A',
        resi: g('auth_seq_id') ?? seq ?? '0', icode: (g('pdbx_PDB_ins_code') ?? '?').replace(/[?.]/, ''), x: +g('Cartn_x')!, y: +g('Cartn_y')!, z: +g('Cartn_z')!, occ: +(g('occupancy') ?? 1), b: +(g('B_iso_or_equiv') ?? 0), el: g('type_symbol') ?? '' })
    } else if (cols.length && (l.startsWith('#') || l.startsWith('_') || l.startsWith('loop_'))) inLoop = false
  }
  return out
}
/** Copy of the atoms for drawing only: a modified nucleotide (PSU, 2MG, H2U, 4SU ...) gets its parent base name, inferred from its base atoms
 *  (N9 + O6 -> G, N9 -> A, N4 -> C, else U), because 3Dmol only draws a continuous nucleic-acid cartoon through A/C/G/U residues.
 *  Order and count are unchanged (serials stay aligned); hover labels keep the deposited names from the original atoms. */
export function drawableNucleotides(atoms: Atom[]): Atom[] {
  const names = new Map<string, Set<string>>()
  for (const a of atoms) { const k = `${a.chain}|${a.resi}|${a.icode}`; (names.get(k) ?? names.set(k, new Set()).get(k)!).add(a.name) }
  const parent = (n: Set<string>) => (n.has('N9') ? (n.has('O6') ? 'G' : 'A') : n.has('N4') ? 'C' : 'U')
  return atoms.map((a) => (a.het || ['A', 'C', 'G', 'U'].includes(a.resn) ? a : { ...a, resn: parent(names.get(`${a.chain}|${a.resi}|${a.icode}`)!) }))
}

export function toPdb(atoms: Atom[], f: (p: [number, number, number]) => [number, number, number] = (p) => p): string {
  const pad = (s: string | number, n: number, left = false) => { const v = String(s); return v.length >= n ? v.slice(0, n) : left ? v.padEnd(n) : v.padStart(n) }
  return atoms.map((a, i) => {
    const [x, y, z] = f([a.x, a.y, a.z])
    const name = a.name.length < 4 && a.el.length === 1 ? ` ${a.name}` : a.name
    return `${a.het ? 'HETATM' : 'ATOM  '}${pad(i + 1, 5)} ${pad(name, 4, true)} ${pad(a.resn, 3)} ${pad(a.chain.slice(0, 1), 1)}${pad(a.resi, 4)}${pad(a.icode || ' ', 1)}   ${pad(x.toFixed(3), 8)}${pad(y.toFixed(3), 8)}${pad(z.toFixed(3), 8)}${pad(a.occ.toFixed(2), 6)}${pad(a.b.toFixed(2), 6)}          ${pad(a.el, 2)}`
  }).join('\n') + '\nEND\n'
}
const applyRT = (T: OrientationTransform) => (p: [number, number, number]): [number, number, number] => {
  const R = T.rotation!, v = T.translation!
  return [0, 1, 2].map((i) => R[i][0] * p[0] + R[i][1] * p[1] + R[i][2] * p[2] + v[i]) as [number, number, number]
}
// Principal axes (cyclic Jacobi), from the reference demo: longest axis up the screen, view down the shortest; deterministic signs.
export function principal(coords: number[][]): { R: number[][]; c: number[] } | null {
  const n = coords.length; if (n < 12) return null
  const c = [0, 1, 2].map((k) => coords.reduce((s, p) => s + p[k], 0) / n)
  const a = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
  for (const p of coords) for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) a[i][j] += ((p[i] - c[i]) * (p[j] - c[j])) / n
  const V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
  for (let sweep = 0; sweep < 24; sweep++) {
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 3; q++) {
      if (Math.abs(a[p][q]) < 1e-15) continue
      const th = (a[q][q] - a[p][p]) / (2 * a[p][q]); const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)); const co = 1 / Math.sqrt(t * t + 1), s = t * co
      for (let k = 0; k < 3; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = co * akp - s * akq; a[k][q] = s * akp + co * akq }
      for (let k = 0; k < 3; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = co * apk - s * aqk; a[q][k] = s * apk + co * aqk }
      for (let k = 0; k < 3; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = co * vkp - s * vkq; V[k][q] = s * vkp + co * vkq }
    }
  }
  const e = [0, 1, 2].map((i) => ({ val: a[i][i], vec: [V[0][i], V[1][i], V[2][i]] })).sort((x, y) => y.val - x.val)
  const fix = (v: number[]) => { let k = 0; for (let i = 1; i < 3; i++) if (Math.abs(v[i]) > Math.abs(v[k])) k = i; return v[k] < 0 ? v.map((x) => -x) : v }
  const y = fix(e[0].vec), x = fix(e[1].vec), z = [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]]
  return { R: [x, y, z], c }
}
const FIT = 0.85
const fit = (v: Viewer) => { v.zoomTo(); v.zoom(FIT) }

/** Without `height` the viewer flexes into its parent (min 320 px); the ResizeObserver below refits the scene whenever the box changes. */
/** overall chain geometry only: a smooth tube along the backbone (as in the protein-complex views). 3Dmol draws base ladders from the base
 *  atoms, so only the sugar-phosphate backbone atoms are passed to the scene (BACKBONE). With an emphasised structure among several, it is
 *  thicker and opaque while the others are faded */
export const cartoon = (color: string, mode: 'normal' | 'emphasis' | 'faded', k = 1) => ({ cartoon: { color, style: 'oval', arrows: false,
  thickness: (mode === 'emphasis' ? 0.9 : 0.55) * k, width: (mode === 'emphasis' ? 0.9 : 0.55) * k, opacity: mode === 'faded' ? 0.35 : 1 } })
/** tube scale for large chains (e.g. 16S rRNA, ~1,500 nt), so the tube keeps a similar on-screen weight as for riboswitch-sized chains */
export const tubeScale = (nResidues: number) => (nResidues > 800 ? 1.9 : nResidues > 300 ? 1.4 : 1)
export const BACKBONE = new Set(["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"])

/** Hover label of a residue: PDB id, author chain, residue name and author number (with insertion code), e.g. "5J7L.AA G230 (author)". */
export const residueLabel = (pdbId: string, a: Pick<Atom, 'chain' | 'resn' | 'resi' | 'icode'>) => `${pdbId.toUpperCase()}.${a.chain} ${a.resn}${a.resi}${a.icode} (author)`

export function AlignmentViewer({ entries, hidden, height, emphasis = null }: { entries: ViewerEntry[]; hidden: string[]; height?: number
  /** seq_name of the structure to emphasise (the selected instance) when it is displayed with others */
  emphasis?: string | null }) {
  const hostRef = useRef<HTMLDivElement>(null), frameRef = useRef<HTMLDivElement>(null)
  const viewerRef = useRef<Viewer | null>(null)
  const modelsRef = useRef(new Map<string, unknown>())   // seq_name -> 3Dmol model of the current scene
  const nResRef = useRef(0)                                // largest chain length (C4' atoms) in the scene
  const cacheRef = useRef(new Map<string, Atom[]>())   // seq_name -> atoms already in the family frame
  const orientRef = useRef<{ R: number[][]; c: number[] } | null>(null)
  const [status, setStatus] = useState<'loading' | 'ready' | 'error'>('loading')
  const [err, setErr] = useState<string | null>(null)
  const [isFull, setIsFull] = useState(false)
  const shownAtomsRef = useRef(new Map<number, { entry: ViewerEntry; atoms: Atom[] }>())   // 3Dmol model id -> its entry and the atoms written (serial = index + 1)
  const pointer = useRef({ x: 0, y: 0, flip: false })   // pointer in the viewer; flip = label to the left (near the right edge)
  const [hover, setHover] = useState<{ text: string; x: number; y: number; flip: boolean } | null>(null)
  const shown = entries.filter((e) => !hidden.includes(e.seqName))
  const key = shown.map((e) => `${e.seqName}:${e.color}`).join(',')
  useEffect(() => {
    const onChange = () => { setIsFull(document.fullscreenElement === frameRef.current); const v = viewerRef.current; if (v) setTimeout(() => { try { v.resize(); v.render() } catch { /* torn down */ } }, 80) }
    document.addEventListener('fullscreenchange', onChange); return () => document.removeEventListener('fullscreenchange', onChange)
  }, [])
  useEffect(() => {
    const el = hostRef.current
    if (!el || typeof ResizeObserver === 'undefined') return
    let frame = 0
    const ro = new ResizeObserver(() => { cancelAnimationFrame(frame); frame = requestAnimationFrame(() => { const v = viewerRef.current; if (v) try { v.resize(); fit(v); v.render() } catch { /* noop */ } }) })
    ro.observe(el); return () => { cancelAnimationFrame(frame); ro.disconnect() }
  }, [])
  useEffect(() => () => { try { viewerRef.current?.clear() } catch { /* noop */ } viewerRef.current = null; if (hostRef.current) hostRef.current.innerHTML = '' }, [])   // release the WebGL scene on unmount
  useEffect(() => {
    let cancelled = false
    ;(async () => {
      setStatus('loading')
      try {
        const $3Dmol = await load3Dmol()
        if (cancelled || !hostRef.current) return
        viewerRef.current ??= $3Dmol.createViewer(hostRef.current, { backgroundColor: 'white' })
        const v = viewerRef.current
        v.clear(); modelsRef.current.clear(); shownAtomsRef.current.clear(); nResRef.current = 0; setHover(null)
        for (const e of shown) {
          let atoms = cacheRef.current.get(e.seqName)
          if (atoms == null) {
            const r = await fetch(`${MODEL_SERVER}/${e.pdbId.toLowerCase()}/atoms?label_asym_id=${encodeURIComponent(e.labelAsymId)}&model_nums=1&encoding=cif`)
            if (!r.ok) throw new Error(`${e.label}: HTTP ${r.status}`)
            const parsed = parseAtomSite(await r.text())
            if (!parsed.length) throw new Error(`${e.label}: no atoms returned`)
            const T = e.transform?.rotation && e.transform.translation ? applyRT(e.transform) : null
            atoms = T ? parsed.map((a) => { const [x, y, z] = T([a.x, a.y, a.z]); return { ...a, x, y, z } }) : parsed
            cacheRef.current.set(e.seqName, atoms)
          }
          if (cancelled) return
          orientRef.current ??= principal(atoms.filter((a) => a.name === "C4'").map((a) => [a.x, a.y, a.z]))
          const o = orientRef.current
          const drawn = atoms.filter((a) => BACKBONE.has(a.name))
          const scene = toPdb(drawableNucleotides(atoms).filter((a) => BACKBONE.has(a.name)), o ? (p) => { const d = [p[0] - o.c[0], p[1] - o.c[1], p[2] - o.c[2]]; return [0, 1, 2].map((i) => o.R[i][0] * d[0] + o.R[i][1] * d[1] + o.R[i][2] * d[2]) as [number, number, number] } : undefined)
          const model = v.addModel(scene, 'pdb')
          modelsRef.current.set(e.seqName, model)
          shownAtomsRef.current.set(model.getID(), { entry: e, atoms: drawn })
          nResRef.current = Math.max(nResRef.current, atoms.filter((a) => a.name === "C4'").length)
        }
        if (cancelled) return
        applyStyles(v)
        // residue hover label from our own parsed atoms (the PDB records written for 3Dmol truncate multi-letter author chain ids)
        v.setHoverDuration(80)
        v.setHoverable({}, true, (a) => {
          const m = shownAtomsRef.current.get(a.model), at = m?.atoms[a.serial - 1]
          if (m && at) setHover({ text: residueLabel(m.entry.pdbId, at), ...pointer.current })
        }, () => setHover(null))
        v.resize(); fit(v); v.render()
        setErr(null); setStatus('ready')
        window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), alignment: { loaded: shown.map((e) => e.seqName), transformed: shown.map((e) => !!e.transform?.rotation) },
          // test hook: screen position (page coordinates) and expected label of the i-th P atom (on the tube's trace) of a displayed structure
          alignmentProbe: (seqName: string, i: number) => {
            const id = modelsRef.current.get(seqName) as { getID: () => number } | undefined, m = id && shownAtomsRef.current.get(id.getID())
            const k = m ? m.atoms.map((a, j) => [a, j] as const).filter(([a]) => a.name === 'P')[i]?.[1] : undefined
            if (!m || k == null) return null
            const sel = (v as unknown as { getModel: (id: number) => { selectedAtoms: (s: Record<string, unknown>) => { x: number; y: number; z: number }[] } }).getModel(id.getID()).selectedAtoms({ serial: k + 1 })[0]
            return { ...v.modelToScreen(sel), label: residueLabel(m.entry.pdbId, m.atoms[k]) }
          } }
      } catch (e) {
        if (!cancelled) { setErr(String((e as Error).message ?? e)); setStatus('error') }
      }
    })()
    return () => { cancelled = true }
  }, [key])   // eslint-disable-line react-hooks/exhaustive-deps
  /** styles of the loaded models (no reload): emphasis only when the emphasised structure is shown together with others */
  const emphasised = emphasis && shown.length > 1 && shown.some((e) => e.seqName === emphasis) ? emphasis : null
  function applyStyles(v: Viewer) {
    const k = tubeScale(nResRef.current)
    for (const e of shown) { const model = modelsRef.current.get(e.seqName); if (model != null) v.setStyle({ model }, cartoon(e.color, !emphasised ? 'normal' : e.seqName === emphasised ? 'emphasis' : 'faded', k)) }
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), alignmentEmphasis: emphasised }
  }
  useEffect(() => { const v = viewerRef.current; if (v && status === 'ready') { applyStyles(v); v.render() } }, [emphasised, status])   // eslint-disable-line react-hooks/exhaustive-deps
  const toggleFull = () => { const el = frameRef.current; if (!el) return; if (document.fullscreenElement) void document.exitFullscreen?.(); else void el.requestFullscreen?.() }
  const recenter = () => { const v = viewerRef.current; if (v) { v.resize(); fit(v); v.render() } }
  return (
    <div ref={frameRef} className={`relative rounded border border-stone-200 bg-white ${height ? '' : 'min-h-[320px] flex-1 lg:min-h-[340px]'}`} style={height ? { height } : undefined} data-testid="alignment-viewer" data-status={status} data-structures={shown.map((e) => e.seqName).join(',')} data-colors={shown.map((e) => e.color).join(',')} data-emphasis={emphasised ?? ''}>
      <div ref={hostRef} className="absolute inset-0" onMouseMove={(e) => { const r = e.currentTarget.getBoundingClientRect(); pointer.current = { x: e.clientX - r.left, y: e.clientY - r.top, flip: e.clientX - r.left > r.width - 200 } }} onMouseLeave={() => setHover(null)} />
      {hover && <div role="tooltip" className="pointer-events-none absolute z-10 whitespace-nowrap rounded border border-stone-300 bg-white/95 px-1.5 py-0.5 font-mono text-[11px] text-stone-800 shadow" style={{ left: hover.flip ? undefined : hover.x + 12, right: hover.flip ? `calc(100% - ${hover.x - 12}px)` : undefined, top: hover.y + 12 }} data-testid="viewer-residue-hover">{hover.text}</div>}
      <div className="absolute right-2 top-2 flex gap-1">
        <button type="button" className="rounded border border-stone-300 bg-white px-2 py-0.5 text-xs shadow-sm hover:bg-stone-50" onClick={toggleFull} data-testid="viewer-fullscreen">{isFull ? 'Exit fullscreen' : 'Fullscreen'}</button>
        <button type="button" className="rounded border border-stone-300 bg-white px-2 py-0.5 text-xs shadow-sm hover:bg-stone-50" onClick={recenter} data-testid="viewer-centre">Centre view</button>
      </div>
      {status === 'loading' && <div className="absolute left-2 top-2 rounded bg-white/90 px-2 py-0.5 text-xs text-stone-600">Loading…</div>}
      {err && <div className="absolute inset-x-2 bottom-2 rounded bg-red-50 px-2 py-1 text-xs text-red-800" data-testid="viewer-error">{err}</div>}
      {!shown.length && <div className="absolute inset-0 flex items-center justify-center text-sm text-stone-500">No structure displayed.</div>}
    </div>
  )
}
