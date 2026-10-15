import { Mat4, Vec3 } from 'molstar/lib/mol-math/linear-algebra'
import type { OrientationTransform } from '../types/backend'

/** Build the Mol* 4x4 matrix for x' = R·x + t. Mol* stores Mat4 column-major; Mat4.ofRows takes row-major rows, so the rows of the
 *  affine matrix are written directly. Applied with TransformStructureConformation (transpose: false). */
export function toMat4(t: Pick<OrientationTransform, 'rotation' | 'translation'>): Mat4 | null {
  const R = t.rotation, v = t.translation
  if (!R || !v) return null
  return Mat4.ofRows([[R[0][0], R[0][1], R[0][2], v[0]], [R[1][0], R[1][1], R[1][2], v[1]], [R[2][0], R[2][1], R[2][2], v[2]], [0, 0, 0, 1]])
}
export function applyMat4(m: Mat4, p: number[]): number[] {
  const out = Vec3.transformMat4(Vec3(), Vec3.create(p[0], p[1], p[2]), m)
  return [out[0], out[1], out[2]]
}
/** A transform is applied only when the fit is usable; 'insufficient' (too few shared coordinates) leaves the structure as deposited. */
export const isApplicable = (t: OrientationTransform | undefined | null): t is OrientationTransform => !!t && t.rotation != null && t.translation != null && t.fit_quality !== 'insufficient'
