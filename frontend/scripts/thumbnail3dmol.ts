/*
 * 3Dmol.js helpers for the static home-page thumbnails only (scripts/make-thumbnails.mjs loads this file through the Vite dev server).
 * The app's 3D viewers are Mol* (components/molstar/MolstarSuperposition.tsx); these helpers are the former Structures-tab 3Dmol code
 * (AlignmentViewer.tsx before the Mol* unification), kept so all thumbnails are rendered the same way: 3Dmol 2.5.2 from the CDN, chain
 * mmCIF from the PDBe ModelServer re-written as PDB records, principal-axes orientation. Never imported by the app.
 */
const THREEDMOL_URL = 'https://cdn.jsdelivr.net/npm/3dmol@2.5.2/build/3Dmol-min.js'

type ThreeDmol = { createViewer: (el: HTMLElement, o: Record<string, unknown>) => Viewer }
type Viewer = { addModel: (t: string, f: string) => unknown; setStyle: (s: Record<string, unknown>, st: Record<string, unknown>) => void; zoomTo: () => void; zoom: (f: number) => void; render: () => void }
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
/** One atom of the chain (from the mmCIF atom_site loop of the PDBe ModelServer response). */
export interface Atom { het: boolean; name: string; resn: string; chain: string; resi: string; icode: string; x: number; y: number; z: number; occ: number; b: number; el: string
  /** label_seq_id (polymer residues; the Interactions tab addresses residues by it) */
  seq?: number }

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
        resi: g('auth_seq_id') ?? seq ?? '0', icode: (g('pdbx_PDB_ins_code') ?? '?').replace(/[?.]/, ''), x: +g('Cartn_x')!, y: +g('Cartn_y')!, z: +g('Cartn_z')!, occ: +(g('occupancy') ?? 1), b: +(g('B_iso_or_equiv') ?? 0), el: g('type_symbol') ?? '', seq: seq && seq !== '.' && seq !== '?' ? +seq : undefined })
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
/** overall chain geometry only: a smooth tube along the backbone (as in the protein-complex views). 3Dmol draws base ladders from the base
 *  atoms, so only the sugar-phosphate backbone atoms are passed to the scene (BACKBONE). With an emphasised structure among several, it is
 *  thicker and opaque while the others are faded */
export const cartoon = (color: string, mode: 'normal' | 'emphasis' | 'faded', k = 1) => ({ cartoon: { color, style: 'oval', arrows: false,
  thickness: (mode === 'emphasis' ? 0.9 : 0.55) * k, width: (mode === 'emphasis' ? 0.9 : 0.55) * k, opacity: mode === 'faded' ? 0.35 : 1 } })
/** tube scale for large chains (e.g. 16S rRNA, ~1,500 nt), so the tube keeps a similar on-screen weight as for riboswitch-sized chains */
export const tubeScale = (nResidues: number) => (nResidues > 800 ? 1.9 : nResidues > 300 ? 1.4 : 1)
