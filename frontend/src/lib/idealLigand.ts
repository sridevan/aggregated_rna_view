/*
 * Ideal (CCD) ligand coordinates for the Ligands tab 3D viewer: source URLs, a minimal molfile check and the principal axes used to orient
 * the view. Pure functions (no Mol*), so they are unit-tested without a browser.
 */
export const idealSdfUrls = (ccd: string) => [
  `https://files.rcsb.org/ligands/download/${encodeURIComponent(ccd)}_ideal.sdf`,
  `https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/${encodeURIComponent(ccd)}_ideal.sdf`,
]
export const IDEAL_INFO = 'This viewer displays idealized coordinates from the PDB Chemical Component Dictionary. These coordinates represent a reference molecular geometry and may differ from the ligand conformation observed in an experimental RNA complex.'

export type SdfCheck = { ok: true; nAtoms: number; xyz: number[][] } | { ok: false; reason: 'invalid' | 'no_ideal_coordinates' }
/** Minimal V2000/V3000 molfile check: an atom block that can be read, and coordinates that are not all at one point (missing ideal coordinates). */
export function checkIdealSdf(text: string): SdfCheck {
  const lines = text.split(/\r?\n/)
  if (lines.length < 4 || !lines.some((l) => l.startsWith('M  END'))) return { ok: false, reason: 'invalid' }
  const xyz: number[][] = []
  if (/V2000/.test(lines[3])) {
    const n = parseInt(lines[3].slice(0, 3), 10)
    if (!(n > 0) || lines.length < 4 + n) return { ok: false, reason: 'invalid' }
    for (const l of lines.slice(4, 4 + n)) xyz.push([l.slice(0, 10), l.slice(10, 20), l.slice(20, 30)].map(Number))
  } else if (/V3000/.test(lines[3])) {
    const a = lines.findIndex((l) => /BEGIN ATOM/.test(l)), b = lines.findIndex((l) => /END ATOM/.test(l))
    if (a < 0 || b <= a + 1) return { ok: false, reason: 'invalid' }
    for (const l of lines.slice(a + 1, b)) { const f = l.trim().split(/\s+/); xyz.push([f[4], f[5], f[6]].map(Number)) }
  } else return { ok: false, reason: 'invalid' }
  if (xyz.some((p) => p.some((v) => !Number.isFinite(v)))) return { ok: false, reason: 'invalid' }
  const spread = Math.max(...[0, 1, 2].map((k) => Math.max(...xyz.map((p) => p[k])) - Math.min(...xyz.map((p) => p[k]))))
  return xyz.length > 1 && spread < 1e-3 ? { ok: false, reason: 'no_ideal_coordinates' } : { ok: true, nAtoms: xyz.length, xyz }
}

/** principal axes of a point set (largest variance first), by power iteration with deflation on the 3x3 covariance */
export function principalAxes(xyz: number[][]): { centre: number[]; axes: number[][] } {
  const n = xyz.length, centre = [0, 1, 2].map((k) => xyz.reduce((a, p) => a + p[k], 0) / n)
  let C = [0, 1, 2].map((i) => [0, 1, 2].map((j) => xyz.reduce((a, p) => a + (p[i] - centre[i]) * (p[j] - centre[j]), 0) / n))
  const norm = (v: number[]) => { const l = Math.hypot(...v) || 1; return v.map((x) => x / l) }
  const axes: number[][] = []
  for (let k = 0; k < 2; k++) {
    let v = norm(k === 0 ? [1, 0.7, 0.3] : [0.3, -0.5, 1].map((x, i) => x - axes[0][i] * (0.3 * axes[0][0] - 0.5 * axes[0][1] + axes[0][2])))
    for (let it = 0; it < 100; it++) v = norm(C.map((row) => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]).map((x, i) => (Math.abs(x) < 1e-12 && it === 0 ? v[i] : x)))
    const lam = C.reduce((a, row, i) => a + v[i] * (row[0] * v[0] + row[1] * v[1] + row[2] * v[2]), 0)
    axes.push(v); C = C.map((row, i) => row.map((x, j) => x - lam * v[i] * v[j]))
  }
  const [a, b] = axes
  axes.push([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])
  return { centre, axes }
}
