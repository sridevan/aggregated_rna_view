import type { StructureClusteringDataset, StructureGroupRecord, StructureMetric, StructureMetricId } from '../types/backend'

/*
 * Structural groups for the Structures tab. Adapted from antibody_antigen_interface/app-react/src/groups/clustering.js: the tree, the
 * cut and the medoids are NOT computed here (they are precomputed by rna_explorer/structure_clustering.py with the same algorithm);
 * this module rebuilds the drawing tree from the stored linkage and carries the demo's colour and display rules.
 */

/** A tree node; `leaves` are matrix indices in drawing order (left subtree first), as stored by the backend. */
export interface TreeNode { id: number; left: TreeNode | null; right: TreeNode | null; height: number; leaves: number[] }
export interface Group extends StructureGroupRecord { color: string; nodeRef: TreeNode; from: number; to: number }
export interface Clustering { metricId: StructureMetricId; metric: StructureMetric; ids: string[]; root: TreeNode; order: string[]; groups: Group[]; groupOf: Record<string, number>; cutHeight: number; maxHeight: number; floor: number | null; dist: number[][]; value: number[][] }

/** Group colours from the reference demo: one named colour per multi-member group (by size), grey for single instances and beyond the
 *  fourth group. They are identities, kept away from the viridis ramp used for pairwise values. */
export const GROUP_COLORS = ['#E8685A', '#B0327A', '#8C4A1E', '#3FA7D6']
export const UNPLACED = '#6b7480'
/** Per-structure colours of the 3D overlay (Okabe-Ito first; reference demo). The representative always takes the first. */
export const SERIES = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#56B4E9', '#E69F00', '#9F3960', '#815C0A']
/** Maximum number of structures of one structural group in the 3D overlay (representative included). Same-group structures are told
 *  apart by shades of the group colour, which stay distinguishable for about five. */
export const OVERLAY_MAX = 5
/** Number of displayed structures that belong to the same structural group as `sn` (non-clustered instances count together). */
export const groupCount = (displayed: string[], sn: string, groupOf: Record<string, number | undefined>) => displayed.filter((x) => groupOf[x] === groupOf[sn]).length

export function buildTree(linkage: [number, number, number, number][], n: number): TreeNode {
  const nodes = new Map<number, TreeNode>()
  for (let i = 0; i < n; i++) nodes.set(i, { id: i, left: null, right: null, height: 0, leaves: [i] })
  linkage.forEach(([a, b, h], k) => {
    const L = nodes.get(a)!, R = nodes.get(b)!
    nodes.set(n + k, { id: n + k, left: L, right: R, height: h, leaves: [...L.leaves, ...R.leaves] })
  })
  return nodes.get(n + linkage.length - 1) ?? nodes.get(0)!
}

export function findNode(root: TreeNode, id: number): TreeNode | null {
  if (root.id === id) return root
  if (!root.left || !root.right) return null
  return findNode(root.left, id) ?? findNode(root.right, id)
}

/** Full square matrix from the stored upper triangle. */
export function squareMatrix(upper: number[], n: number, diag: number): number[][] {
  const M = Array.from({ length: n }, () => new Array<number>(n).fill(diag))
  let k = 0
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { M[i][j] = M[j][i] = upper[k++] }
  return M
}

export function clusteringFor(data: StructureClusteringDataset, metricId: StructureMetricId): Clustering {
  const metric = data.metrics[metricId]
  const ids = data.dataset.instances
  const n = ids.length
  const c = metric.clustering
  const root = buildTree(c.linkage, n)
  if (root.leaves.join() !== c.leaf_order.join()) throw new Error('stored leaf order does not match the stored linkage')
  const value = squareMatrix(metric.values_upper, n, metric.kind === 'similarity' ? 1 : 0)
  const dist = metric.kind === 'similarity' ? value.map((r) => r.map((v) => 1 - v)) : value
  let coloured = 0
  const groups: Group[] = c.groups.map((g) => ({ ...g, color: g.n_instances > 1 ? (GROUP_COLORS[coloured++] ?? UNPLACED) : UNPLACED, nodeRef: findNode(root, g.node)!, from: g.leaf_range[0], to: g.leaf_range[1] }))
  return { metricId, metric, ids, root, order: c.leaf_order.map((i) => ids[i]), groups, groupOf: c.assignments, cutHeight: c.cut.height, maxHeight: c.max_height, floor: c.floor.value, dist, value }
}

/** Largest group (groups are numbered by size, so group 1). */
export const largestGroup = (c: Clustering) => c.groups[0]

/** Default 3D content of a group: every member when the group has 2..OVERLAY_MAX members (representative first), else the representative alone. */
export function defaultDisplayed(g: Pick<StructureGroupRecord, 'members' | 'representative'>): string[] {
  if (g.members.length > 1 && g.members.length <= OVERLAY_MAX) return [g.representative, ...g.members.filter((m) => m !== g.representative)]
  return [g.representative]
}

/** Viridis (matplotlib control points, linear in sRGB), as in the reference demo's heatmap. */
const VIRIDIS = ['#440154', '#482878', '#3e4a89', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6dcd59', '#b4de2c', '#fde725'].map((h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16)))
export function viridis(t: number): string {
  const x = Math.max(0, Math.min(1, t)) * (VIRIDIS.length - 1)
  const i = Math.min(VIRIDIS.length - 2, Math.floor(x)), f = x - i
  return `rgb(${VIRIDIS[i].map((c, k) => Math.round(c + (VIRIDIS[i + 1][k] - c) * f)).join(',')})`
}
export const DIAGONAL = '#d7dae0'

/** "4fe5_1_B_A" -> "4FE5 B" (PDB entry + author chain, for display; the full seq_name stays the identity). */
export const instanceLabel = (seqName: string) => { const [pdb, , auth] = seqName.split('_'); return `${pdb.toUpperCase()} ${auth}` }
export const fmt = (v: number) => (Math.abs(v) >= 1 ? v.toFixed(2) : v.toFixed(3))

/** Node ids on the route through the tree between two nodes: up from each to their lowest common ancestor, that merge included
 *  (ported from antibody_antigen_interface/app-react/src/groups/clustering.js pathBetween). Null when either node is not in this tree. */
export function pathBetween(root: TreeNode, idA: number, idB: number): Set<number> | null {
  const chain = (target: number): number[] | null => {
    const walk = (node: TreeNode, trail: number[]): number[] | null => {
      const here = [...trail, node.id]
      if (node.id === target) return here
      if (!node.left || !node.right) return null
      return walk(node.left, here) ?? walk(node.right, here)
    }
    return walk(root, [])
  }
  const a = chain(idA), b = chain(idB)
  if (!a || !b) return null
  let k = 0
  while (k < a.length && k < b.length && a[k] === b[k]) k++
  return new Set([...a.slice(k - 1), ...b.slice(k)])
}

/** Shades of a group colour for the non-representative members shown in 3D (the representative keeps the exact group colour).
 *  Alternating lighter / darker steps so neighbours in the legend differ in lightness; mixing in sRGB with white or near-black. */
const SHADE_STEPS = [0.45, -0.45, 0.72, -0.68, 0.22, -0.22, 0.86]
export function groupShade(hex: string, k: number): string {
  const t = SHADE_STEPS[k % SHADE_STEPS.length]
  const rgb = [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16))
  const target = t > 0 ? 255 : 20
  const out = rgb.map((c) => Math.round(c + (target - c) * Math.abs(t)))
  return `#${out.map((c) => c.toString(16).padStart(2, '0')).join('')}`
}
/** 3D colour of each displayed structure: a group's representative takes its group colour; other members take shades of their own
 *  group's colour, numbered per group in display order. */
export function structureColours(displayed: string[], c: Pick<Clustering, 'groups' | 'groupOf'>): Record<string, string> {
  const seen: Record<number, number> = {}, out: Record<string, string> = {}
  for (const sn of displayed) {
    const g = c.groups.find((x) => x.group_id === c.groupOf[sn])
    if (!g) { out[sn] = UNPLACED; continue }
    if (g.representative === sn) { out[sn] = g.color; continue }
    const k = seen[g.group_id] ?? 0; seen[g.group_id] = k + 1
    out[sn] = groupShade(g.color, k)
  }
  return out
}
