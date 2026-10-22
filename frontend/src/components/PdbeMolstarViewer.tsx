import { useEffect, useRef, useState } from 'react'
import { PDBeMolstarPlugin } from 'pdbe-molstar/lib/viewer'
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder'
import { StateTransforms } from 'molstar/lib/mol-plugin-state/transforms'
import { createStructureRepresentationParams } from 'molstar/lib/mol-plugin-state/helpers/structure-representation-params'
import { PluginCommands } from 'molstar/lib/mol-plugin/commands'
import { setSubtreeVisibility } from 'molstar/lib/mol-plugin/behavior/static/state'
import { StructureFocusRepresentation } from 'molstar/lib/mol-plugin/behavior/dynamic/selection/structure-focus-representation'
import { Script } from 'molstar/lib/mol-script/script'
import { StructureSelection, StructureElement, StructureProperties } from 'molstar/lib/mol-model/structure'
import { Vec3 } from 'molstar/lib/mol-math/linear-algebra'
import { Color } from 'molstar/lib/mol-util/color'
import { OrderedSet } from 'molstar/lib/mol-data/int'
import { StructureSelectionQueries } from 'molstar/lib/mol-plugin-state/helpers/structure-selection-query'
import type { Expression } from 'molstar/lib/mol-script/language/expression'
import type { Structure } from 'molstar/lib/mol-model/structure'
import type { OrientationTransform } from '../types/backend'
import { isApplicable, toMat4 } from '../lib/orientation'
import 'pdbe-molstar/build/pdbe-molstar-light.css'

export interface MolstarSelectionItem { label_seq_id: number }
export type ViewMode = 'region' | 'complete'
interface Props {
  pdbId: string
  labelAsymId: string
  authAsymId: string
  region?: { start: number; end: number }
  selection: MolstarSelectionItem[]
  /** centre the camera on the selection whenever it changes (single-residue / group exploration); the view direction is kept */
  focusOnSelection?: boolean
  /** precomputed rigid-body transform into the family reference frame (orientation.json); null/insufficient = shown as deposited */
  orientation?: OrientationTransform | null
  onResidueClick?: (labelSeqId: number) => void
  /** changes on every explicit deselection: clear highlight and focus, restore the default whole-structure (region) camera */
  resetKey?: number
  /** one fixed colour (hex) for the polymer / Rfam-region representation of every structure; undefined = PDBe Mol* default colouring */
  structureColor?: number
  /** "Clear selection" in the viewer controls: clears the shared selection (the parent then restores the default 2D and 3D views) */
  onClearSelection?: () => void
  canClear?: boolean
  className?: string
  initialMode?: ViewMode
  /** CCD id -> ligand class of the family (ligands.json by_ccd; loaded lazily when "Show ligands" is first switched on) */
  ligandClasses?: () => Promise<Record<string, string>>
}
interface MolstarClickEvent extends Event { eventData?: { seq_id?: number; auth_seq_id?: number; auth_asym_id?: string; entity_id?: string; ins_code?: string } }
/** selection (nucleotide or whole loop): violet #7C3AED, distinct from the slate structure colour; identical for every instance */
const SELECT_HEX = 0x7c3aed
const SELECT_COLOR = { r: 0x7c, g: 0x3a, b: 0xed }
const WHITE = { r: 255, g: 255, b: 255 }
const REGION_LABEL = 'Rfam-mapped region'
/**
 * Complete structure: the selected RNA chain stays prominent (its own opaque cartoon in the structure colour) while the rest of the
 * deposited model gives context: other nucleic-acid chains and proteins as subdued, semi-transparent grey cartoons (proteins lighter);
 * ligands, non-standard residues and ions keep Mol*'s atom-coloured representations; water stays hidden. Same rules for every family.
 */
const CHAIN_LABEL = 'Selected RNA chain', CONTEXT_RNA_LABEL = 'Other nucleic-acid chains', CONTEXT_PROTEIN_LABEL = 'Protein chains'
const CONTEXT_RNA = { color: 0xa8a29e, alpha: 0.25 }, CONTEXT_PROTEIN = { color: 0xd6d3d1, alpha: 0.2 }
/**
 * "Show ligands" (off by default, independent of the view mode). PDBe binding-site annotations are entity-level (bound-molecule ids, no
 * chain-specific ligand instance), so association is spatial and per instance: a small-molecule residue (non-polymer, not water, not an
 * ion, saccharide or lipid entity) with any atom within LIGAND_RADIUS of the selected RNA chain in the loaded coordinates. CCDs the family
 * classifies as crystallisation additives or metal ions (ligands.json by_ccd) are not part of the toggle; ions keep their own component.
 */
const LIGAND_RADIUS = 4.0
const LIGAND_CARBON_DEFAULT = 0x64748b   // ligand carbons follow the chain colour (structureColor); slate when none is given
const LIGAND_EXCLUDED_CLASSES = ['additive', 'metal_ion']
const ASSOC_LIGANDS_LABEL = 'Ligands near the selected chain', OTHER_LIGANDS_LABEL = 'Other ligands'
const ORIENT_TAG = 'rna-explorer-orientation'
/** simplified display: the focus representation shows only the selected nucleotides (ball-and-stick); no 5 A surroundings, no contact lines */
const FOCUS_COMPONENTS = ['target'] as const
/** water is hidden (visibility only: the component and the coordinates stay in the structure) */
const isWater = (c: { key?: string; cell: { obj?: { label?: string } } }) => c.key === 'water' || c.cell.obj?.label === 'Water'
/** PDBe Mol*'s default ligand component: replaced by the "Show ligands" layers (hidden in both modes) */
const isLigandComponent = (c: { key?: string; cell: { obj?: { label?: string } } }) => c.key === 'ligand' || c.cell.obj?.label === 'Ligand'
/** PDBe Mol*'s default polymer component (all chains); replaced by the three Complete-structure layers once they exist */
const isPolymer = (c: { key?: string; cell: { obj?: { label?: string } } }) => c.key === 'polymer' || c.cell.obj?.label === 'Polymer'
const VIEW_UP = Vec3.create(0, 1, 0), VIEW_DIR = Vec3.create(0, 0, -1)   // fixed camera orientation in the reference frame

/**
 * Reusable PDBe Mol* wrapper. One plugin per mount; structures are switched with `visual.update(..., fullLoad=false)` (load the new
 * entry, then remove the old one) so the canvas and its white background are never re-initialised. A precomputed rigid-body transform
 * (TransformStructureConformation, inserted below the structure node) moves each structure into the family reference frame; the
 * camera is then set from a fixed view direction, so all instances appear in the same orientation. Selection colours residues; with
 * `focusOnSelection` the camera is also centred on them (direction kept, context margin), and clearing never moves the camera.
 */
export function PdbeMolstarViewer({ pdbId, labelAsymId, authAsymId, region, selection, focusOnSelection = false, orientation, onResidueClick, resetKey = 0, structureColor, onClearSelection, canClear = false, className = '', initialMode = 'region', ligandClasses }: Props) {
  const host = useRef<HTMLDivElement>(null)
  const outer = useRef<HTMLDivElement>(null)
  const plugin = useRef<PDBeMolstarPlugin | null>(null)
  const queue = useRef<Promise<void>>(Promise.resolve())
  const regionRef = useRef<string | null>(null)
  /** state refs of the Complete-structure layers of the loaded structure (selected chain, other nucleic acids, proteins) */
  const layersRef = useRef<string[] | null>(null)
  /** "Show ligands": [near the selected chain, other ligands] layer refs of the loaded structure */
  const ligandLayersRef = useRef<string[] | null>(null)
  const [showLigands, setShowLigands] = useState(false)
  const showLigandsRef = useRef(false); showLigandsRef.current = showLigands
  const modeRef = useRef<ViewMode>('region')
  /** spatial candidates near the selected chain (before the class filter) and, once the layers exist, the ligand instances shown */
  const [ligandInfo, setLigandInfo] = useState<{ candidates: string[]; shown: string[] | null; excludedNearby: string[] } | null>(null)
  const classesRef = useRef(ligandClasses); classesRef.current = ligandClasses
  const [status, setStatus] = useState<'loading' | 'ready' | 'error'>('loading')
  const [loadedPdb, setLoadedPdb] = useState<string | null>(null)
  const [pluginGen, setPluginGen] = useState(0)        // bumped when a plugin instance exists: drives the load effect deterministically
  const [retry, setRetry] = useState(0)
  const [errorMsg, setErrorMsg] = useState<string | null>(null)
  const [mode, setMode] = useState<ViewMode>(region ? initialMode : 'complete')
  const clickRef = useRef(onResidueClick); clickRef.current = onResidueClick
  const authRef = useRef(authAsymId); authRef.current = authAsymId
  const selRef = useRef<{ label_asym_id: string; label_seq_id: number; color: typeof SELECT_COLOR }[]>([])
  const applied = isApplicable(orientation) ? orientation : null

  const structureOf = (p: PDBeMolstarPlugin) => {
    const all = p.plugin.managers.structure.hierarchy.current.structures
    return [...all].sort((a, b) => b.components.length - a.components.length)[0]
  }
  const visibleStructureData = (p: PDBeMolstarPlugin): Structure | undefined => {
    const s = structureOf(p)
    return (s?.transform?.cell.obj?.data as Structure | undefined) ?? (s?.cell.obj?.data as Structure | undefined)
  }

  /** Insert (or update) the conformation transform below the structure node, mirroring Mol*'s own superposition helper. */
  const applyOrientation = async (p: PDBeMolstarPlugin, t: OrientationTransform | null) => {
    const ctx = p.plugin
    const s = structureOf(p)
    if (!s) return
    const existing = ctx.state.data.selectQ((q) => q.byRef(s.cell.transform.ref).subtree().withTransformer(StateTransforms.Model.TransformStructureConformation))[0]
    const m = t ? toMat4(t) : null
    if (!m) { if (existing) await PluginCommands.State.RemoveObject(ctx, { state: ctx.state.data, ref: existing.transform.ref }); return }
    const params = { transform: { name: 'matrix' as const, params: { data: m, transpose: false } } }
    const b = existing ? ctx.state.data.build().to(existing).update(params) : ctx.state.data.build().to(s.cell.transform.ref).insert(StateTransforms.Model.TransformStructureConformation, params, { tags: ORIENT_TAG })
    await ctx.runTask(ctx.state.data.updateTree(b))
  }

  /** Complete-structure layers (built once per loaded structure, then only shown or hidden): the selected chain as an opaque cartoon in the
   *  structure colour, other nucleic-acid chains and proteins as subdued semi-transparent grey cartoons. */
  const ensureLayers = async (p: PDBeMolstarPlugin, s: NonNullable<ReturnType<typeof structureOf>>) => {
    const ctx = p.plugin
    if (layersRef.current?.every((r) => ctx.state.data.cells.has(r))) return
    let structureData = visibleStructureData(p)
    for (let i = 0; i < 60 && !structureData; i++) { await new Promise((r) => setTimeout(r, 50)); structureData = visibleStructureData(p) }
    if (!structureData) { layersRef.current = null; return }
    const chain = MS.struct.generator.atomGroups({ 'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), labelAsymId]) })
    const others = (q: typeof StructureSelectionQueries.protein) => MS.struct.modifier.exceptBy({ 0: q.expression, by: chain })
    const parentRef = s.transform?.cell.transform.ref ?? s.cell.transform.ref
    const update = ctx.build()
    const refs: string[] = []
    const layer = (expression: Expression, label: string, color: number | null, alpha: number) => {
      const sel = update.to(parentRef).apply(StateTransforms.Model.StructureSelectionFromExpression, { expression, label })
      sel.apply(StateTransforms.Representation.StructureRepresentation3D, createStructureRepresentationParams(ctx, structureData, {
        type: 'cartoon', typeParams: { alpha }, ...(color != null ? { color: 'uniform' as const, colorParams: { value: Color(color) } } : {}) }))
      refs.push(sel.ref)
    }
    layer(chain, CHAIN_LABEL, structureColor ?? null, 1)
    layer(others(StructureSelectionQueries.nucleic), CONTEXT_RNA_LABEL, CONTEXT_RNA.color, CONTEXT_RNA.alpha)
    layer(others(StructureSelectionQueries.protein), CONTEXT_PROTEIN_LABEL, CONTEXT_PROTEIN.color, CONTEXT_PROTEIN.alpha)
    await update.commit()
    layersRef.current = refs
  }

  const chainExpr = () => MS.struct.generator.atomGroups({ 'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), labelAsymId]) })
  /** small-molecule residues (one atom group per residue), optionally without the given CCD ids */
  const smallMolecules = (exclude: string[]) => MS.struct.generator.atomGroups({
    'entity-test': MS.core.logic.and([MS.core.rel.eq([MS.ammp('entityType'), 'non-polymer']), MS.core.logic.not([MS.core.str.match([MS.re('(ion|oligosaccharide|lipid)', 'i'), MS.ammp('entitySubtype')])])]),
    ...(exclude.length ? { 'residue-test': MS.core.logic.not([MS.core.set.has([MS.set(...exclude), MS.ammp('label_comp_id')])]) } : {}),
    'group-by': MS.ammp('residueKey'),
  })
  const nearChain = (exclude: string[]) => MS.struct.filter.within({ 0: smallMolecules(exclude), target: chainExpr(), 'max-radius': LIGAND_RADIUS })
  /** ligand instances of an expression in the loaded structure: "CCD auth_chain auth_seq" (instance identifiers kept) */
  const instancesOf = (p: PDBeMolstarPlugin, expr: Expression): string[] => {
    const data = visibleStructureData(p); if (!data) return []
    const loci = StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(expr, data))
    const out = new Set<string>()
    StructureElement.Loci.forEachLocation(loci, (l) => { out.add(`${StructureProperties.atom.label_comp_id(l)} ${StructureProperties.chain.auth_asym_id(l)}${StructureProperties.residue.auth_seq_id(l)}${StructureProperties.residue.pdbx_PDB_ins_code(l) || ''}`) })
    return [...out].sort()
  }
  /** build the two ligand layers once per loaded structure (class filter from the family classification), then show / hide them */
  const applyLigands = async (p: PDBeMolstarPlugin, m: ViewMode) => {
    const ctx = p.plugin
    if (ligandLayersRef.current && !ligandLayersRef.current.every((r) => ctx.state.data.cells.has(r))) ligandLayersRef.current = null
    const show = showLigandsRef.current
    if (show && !ligandLayersRef.current) {
      const s = structureOf(p), data = visibleStructureData(p)
      if (!s || !data) return
      const classes = classesRef.current ? await classesRef.current().catch(() => ({} as Record<string, string>)) : {}
      const exclude = Object.entries(classes).filter(([, c]) => LIGAND_EXCLUDED_CLASSES.includes(c)).map(([k]) => k)
      const near = nearChain(exclude), other = MS.struct.modifier.exceptBy({ 0: smallMolecules(exclude), by: near })
      const parentRef = s.transform?.cell.transform.ref ?? s.cell.transform.ref
      const update = ctx.build(); const refs: string[] = []
      for (const [expression, label] of [[near, ASSOC_LIGANDS_LABEL], [other, OTHER_LIGANDS_LABEL]] as const) {
        const sel = update.to(parentRef).apply(StateTransforms.Model.StructureSelectionFromExpression, { expression, label })
        sel.apply(StateTransforms.Representation.StructureRepresentation3D, createStructureRepresentationParams(ctx, data, {
          // sticks (equal sphere and bond radii); carbons in the RNA chain's colour, other elements CPK (N blue, O red, P orange, S yellow)
          type: 'ball-and-stick', typeParams: { sizeFactor: 0.22, sizeAspectRatio: 1 }, color: 'element-symbol',
          colorParams: { carbonColor: { name: 'uniform', params: { value: Color(structureColor ?? LIGAND_CARBON_DEFAULT), saturation: 0, lightness: 0 } }, lightness: 0 } }))
        refs.push(sel.ref)
      }
      await update.commit()
      ligandLayersRef.current = refs
      const shown = instancesOf(p, near), all = instancesOf(p, nearChain([]))
      setLigandInfo((cur) => ({ candidates: cur?.candidates ?? all, shown, excludedNearby: all.filter((x) => !shown.includes(x)) }))
    }
    const [near, other] = ligandLayersRef.current ?? []
    if (near) setSubtreeVisibility(ctx.state.data, near, !show)
    if (other) setSubtreeVisibility(ctx.state.data, other, !(show && m === 'complete'))
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), molstarLigands: { show, mode: m, built: !!ligandLayersRef.current, nearVisible: !!near && show, otherVisible: !!other && show && m === 'complete',
      near: near ? instancesOf(p, (ctx.state.data.cells.get(near)?.transform.params as { expression: Expression }).expression) : [],
      style: (() => {   // representation and colouring of the ligand sticks (test hook)
        const repr = near ? ctx.state.data.tree.children.get(near).first() : undefined
        const q = repr ? ctx.state.data.cells.get(repr)?.transform.params as { type?: { name: string; params: { sizeAspectRatio?: number } }; colorTheme?: { name: string; params: { carbonColor?: { name: string; params: { value?: number } } } } } : undefined
        return q ? { type: q.type?.name, aspect: q.type?.params.sizeAspectRatio, color: q.colorTheme?.name, carbon: q.colorTheme?.params.carbonColor?.name, carbonValue: q.colorTheme?.params.carbonColor?.params.value } : null
      })() } }
  }

  const applyMode = async (p: PDBeMolstarPlugin, m: ViewMode, fit = true) => {
    modeRef.current = m
    const ctx = p.plugin
    for (let i = 0; i < 60 && ctx.managers.structure.hierarchy.current.structures.length === 0; i++) await new Promise((r) => setTimeout(r, 50))
    const s = structureOf(p)
    if (!s) return
    const comps = () => ctx.managers.structure.hierarchy.current.structures.flatMap((x) => x.components)
    if (regionRef.current && !ctx.state.data.cells.has(regionRef.current)) regionRef.current = null   // region of a structure that has been removed
    if (layersRef.current && !layersRef.current.every((r) => ctx.state.data.cells.has(r))) layersRef.current = null
    if (m === 'region' && region) {
      if (!regionRef.current) {
        const expr = MS.struct.generator.atomGroups({
          'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), labelAsymId]),
          'residue-test': MS.core.rel.inRange([MS.struct.atomProperty.macromolecular.label_seq_id(), region.start, region.end]),
        })
        const parentRef = s.transform?.cell.transform.ref ?? s.cell.transform.ref   // below the orientation transform so the region is transformed too
        // the structure (or its transformed copy) may not have resolved yet right after a switch: wait for it rather than skipping the region
        let structureData = visibleStructureData(p)
        for (let i = 0; i < 60 && !structureData; i++) { await new Promise((r) => setTimeout(r, 50)); structureData = visibleStructureData(p) }
        if (structureData) {
          const update = ctx.build()
          const sel = update.to(parentRef).apply(StateTransforms.Model.StructureSelectionFromExpression, { expression: expr, label: REGION_LABEL })
          // cartoon only (Mol* draws the bases as rings/ladders): the selected nucleotide then stands out in its atom-level focus representation
          sel.apply(StateTransforms.Representation.StructureRepresentation3D, createStructureRepresentationParams(ctx, structureData, { type: 'cartoon', ...uniformColor() }))
          await update.commit()
          regionRef.current = sel.ref
        } else {
          console.warn('[PdbeMolstarViewer] region component could not be built; showing the complete structure')
        }
      }
      // never leave an empty scene: without a region component the default components stay visible (Mol* drops camera updates when nothing is visible)
      for (const c of comps()) setSubtreeVisibility(ctx.state.data, c.cell.transform.ref, regionRef.current ? c.cell.transform.ref !== regionRef.current : isWater(c))
      // Complete-structure layers (also empty ones, which Mol* does not list as components) stay hidden in Rfam region mode
      if (regionRef.current) for (const r of layersRef.current ?? []) setSubtreeVisibility(ctx.state.data, r, true)
    } else {
      await ensureLayers(p, s)
      const layers = layersRef.current ?? []
      // the default polymer component is replaced by the layers (when they could be built); the region and water stay hidden
      for (const c of comps()) setSubtreeVisibility(ctx.state.data, c.cell.transform.ref, c.cell.transform.ref === regionRef.current || isWater(c) || isLigandComponent(c) || (layers.length > 0 && isPolymer(c) && !layers.includes(c.cell.transform.ref)))
      for (const r of layers) setSubtreeVisibility(ctx.state.data, r, false)
      // a selection made before the layers existed: re-apply its highlight so the selected chain layer shows it too
      if (layers.length && selRef.current.length) await p.visual.select({ data: selRef.current })
    }
    // the same fixed colour on the default polymer component (shown only when the Complete-structure layers could not be built)
    if (structureColor != null) {
      const polymers = comps().filter(isPolymer)
      if (polymers.length) await ctx.managers.structure.component.updateRepresentationsTheme(polymers, { color: 'uniform', colorParams: { value: Color(structureColor) } })
    }
    const regionCell = regionRef.current ? ctx.state.data.cells.get(regionRef.current) : undefined
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), molstarColours: representationColours(p) }
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), molstarComponents: [...comps().map((c) => ({ label: c.cell.obj?.label, hidden: !!c.cell.state.isHidden })), { label: REGION_LABEL, hidden: regionCell?.state.isHidden ?? null, present: !!regionCell }],
      molstarLayers: (layersRef.current ?? []).map((r) => { const cell = ctx.state.data.cells.get(r); const rep = cell ? ctx.state.data.tree.children.get(r).first() : undefined
        const params = rep ? ctx.state.data.cells.get(rep)?.transform.params as { type?: { params?: { alpha?: number } }; colorTheme?: { name: string; params?: { value?: number } } } | undefined : undefined
        return { label: cell?.obj?.label ?? null, hidden: !!cell?.state.isHidden, alpha: params?.type?.params?.alpha ?? null, theme: params?.colorTheme?.name ?? null, value: params?.colorTheme?.params?.value ?? null, elements: (cell?.obj?.data as Structure | undefined)?.elementCount ?? null } }) }
    await applyLigands(p, m)
    // a mode switch keeps the current camera (no zoom-out to a whole assembly); a new structure gets the default RNA-focused view
    if (fit) await resetView(p, m)
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), molstarModeApplied: { mode: m, fit, n: Number((window.__rnaExplorer?.molstarModeApplied as { n?: number } | undefined)?.n ?? 0) + 1 } }
  }

  /**
   * Deterministic, tight camera: the visible atoms (region or whole chain) are projected onto the fixed view axes (dir, up) of the
   * reference frame; the camera is centred on the projected bounding box and placed at the distance at which that box fills the
   * canvas (given the camera fov and the actual canvas aspect ratio) with a small margin. No arbitrary zoom factor; hidden
   * components never enter the calculation; the view direction is identical for every instance.
   */
  const resetView = async (p: PDBeMolstarPlugin, m: ViewMode) => {
    const ctx = p.plugin
    const data = visibleStructureData(p)
    if (!data || !ctx.canvas3d) return
    const expr = m === 'region' && region
      ? MS.struct.generator.atomGroups({ 'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), labelAsymId]), 'residue-test': MS.core.rel.inRange([MS.struct.atomProperty.macromolecular.label_seq_id(), region.start, region.end]) })
      : MS.struct.generator.atomGroups({ 'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), labelAsymId]) })   // Complete structure: the selected RNA chain, not the whole deposited model
    const loci = StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(expr, data))
    const sphere = StructureElement.Loci.getBoundary(loci).sphere
    // view axes
    const dir = Vec3.normalize(Vec3(), VIEW_DIR)
    const right = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), dir, VIEW_UP))
    const up = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), right, dir))
    // projected coordinates of every visible atom in the view frame (x: right, y: up, z: towards the camera = -dir)
    const xs: number[] = [], ys: number[] = [], zs: number[] = []
    const v = Vec3(), d = Vec3()
    for (const e of loci.elements) {
      const unit = e.unit
      OrderedSet.forEach(e.indices, (i) => {
        unit.conformation.position(unit.elements[i], v); Vec3.sub(d, v, sphere.center)
        xs.push(Vec3.dot(d, right)); ys.push(Vec3.dot(d, up)); zs.push(-Vec3.dot(d, dir))
      })
    }
    const n = xs.length
    if (!n) return
    const cam = ctx.canvas3d.camera
    const fov = cam.state.fov
    const vp = cam.viewport
    const aspect = vp.width > 0 && vp.height > 0 ? vp.width / vp.height : 1
    const tanHalf = Math.tan(fov / 2)
    const MARGIN = 1.06                                   // 6 % breathing space so sticks at the silhouette are never clipped
    const OVERLAY_PX = 40                                 // top controls / bottom badge overlay the canvas: keep the molecule out of those bands
    const usableH = Math.max(vp.height - 2 * OVERLAY_PX, vp.height * 0.5)
    const tanHalfY = tanHalf * usableH / vp.height        // effective vertical half-angle of the usable band
    const xmin = Math.min(...xs), xmax = Math.max(...xs), ymin = Math.min(...ys), ymax = Math.max(...ys)
    const cx = (xmin + xmax) / 2, cy = (ymin + ymax) / 2
    // exact perspective fit: an atom at (x, y) displaced z towards the camera is inside the frustum when |y| <= (d - z)·tan(fov/2)
    // and |x| <= (d - z)·tan(fov/2)·aspect, so the camera distance is the maximum over atoms of the required d
    let distance = 0
    for (let i = 0; i < n; i++) {
      const need = Math.max(Math.abs(ys[i] - cy) * MARGIN / tanHalfY, Math.abs(xs[i] - cx) * MARGIN / (tanHalf * aspect)) + zs[i]
      if (need > distance) distance = need
    }
    distance = Math.max(distance, 5)
    const target = Vec3.scaleAndAdd(Vec3(), Vec3.scaleAndAdd(Vec3(), Vec3.clone(sphere.center), right, cx), up, cy)
    const position = Vec3.sub(Vec3(), target, Vec3.scale(Vec3(), dir, distance))
    // Camera.Reset re-derives the distance from the bounding radius (Camera.targetDistance), discarding the position; SetSnapshot applies ours as is
    const snapshot = { target, position, up: Vec3.clone(VIEW_UP), radius: Math.max(sphere.radius, 5), radiusMax: Math.max(sphere.radius, 5) }
    const halfW = (xmax - xmin) / 2, halfH = (ymax - ymin) / 2
    // camera resets requested by the loader are resolved on the next frames; let them pass, apply ours, then verify once
    const frames = (n: number) => new Promise<void>((res) => { const step = (k: number) => (k <= 0 ? res() : requestAnimationFrame(() => step(k - 1))); step(n) })
    await frames(2)
    await PluginCommands.Camera.SetSnapshot(ctx, { snapshot, durationMs: 0 })
    await frames(2)
    const st = ctx.canvas3d.camera.state
    if (Math.abs(Vec3.distance(st.position, st.target) - distance) > 0.01 * distance) await PluginCommands.Camera.SetSnapshot(ctx, { snapshot, durationMs: 0 })
    const occupancy = { width: halfW / (distance * tanHalf * aspect), height: halfH / (distance * tanHalf) }   // projected extent at the box centre depth / visible extent
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), lastView: { center: [...sphere.center], radius: sphere.radius, up: [...VIEW_UP], dir: [...VIEW_DIR], atoms: n, distance, aspect, fov, occupancy } }
  }

  useEffect(() => {
    const el = host.current, wrap = outer.current
    if (!el || !wrap) return
    let alive = true
    let p: PDBeMolstarPlugin | null = null
    const handler = (e: Event) => {
      const d = (e as MolstarClickEvent).eventData
      if (!d || d.seq_id === undefined || d.auth_asym_id !== authRef.current) return
      clickRef.current?.(d.seq_id)
    }
    wrap.addEventListener('PDB.molstar.click', handler)
    // Mol* only follows window resizes: follow the canvas area too (toolbar wrapping, panels that change height)
    const ro = new ResizeObserver(() => { try { p?.plugin?.canvas3d?.handleResize() } catch { /* not rendered or torn down */ } })
    ro.observe(el)
    const timer = window.setTimeout(() => {
      if (!alive) return
      p = new PDBeMolstarPlugin()
      plugin.current = p
      const inst = p
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), molstar: inst, getCamera: () => { const c = inst.plugin.canvas3d?.camera.state; return c ? { position: [...c.position], target: [...c.target], up: [...c.up], radius: c.radius } : null } }
      setPluginGen((g) => g + 1)
    }, 0)
    return () => {
      alive = false
      window.clearTimeout(timer)
      wrap.removeEventListener('PDB.molstar.click', handler); ro.disconnect()
      const inst = p
      plugin.current = null
      queue.current = queue.current.catch(() => undefined).then(() => inst?.clear()).catch(() => undefined).then(() => { el.innerHTML = '' })
    }
  }, [])

  // load / switch structure without re-initialising the canvas, then orientation + mode
  useEffect(() => {
    const p = plugin.current, el = host.current
    if (!p || !el) return
    let alive = true
    setStatus('loading'); setErrorMsg(null)
    // watchdog: never leave an indefinite blank viewer
    const watchdog = window.setTimeout(() => { if (alive) { setErrorMsg(`Loading ${pdbId.toUpperCase()} timed out after 60 s`); setStatus('error') } }, 60_000)
    queue.current = queue.current.catch(() => undefined).then(async () => {
      if (!alive || plugin.current !== p) return
      regionRef.current = null; layersRef.current = null; ligandLayersRef.current = null   // reset inside the queue: a still-running previous load must not leave its refs behind
      if (loadedPdb === null && !el.querySelector('canvas')) {
        await p.render(el, { moleculeId: pdbId, subscribeEvents: true, hideControls: true, sequencePanel: false, bgColor: WHITE, pdbeLink: false, expanded: false, loadingOverlay: false,
          hideCanvasControls: ['expand', 'controlToggle', 'controlInfo', 'animation', 'trajectory'] })
        await p.canvas.setBgColor(WHITE)
        p.plugin.canvas3d?.setProps({ camera: { manualReset: true } })   // no automatic camera resets on state changes: the view is set explicitly below
        await p.canvas.setBgColor(WHITE)                                   // setProps may re-derive renderer defaults: re-assert the white background
        await setFocusComponents(p)
      } else {
        const ctx = p.plugin
        const old = ctx.managers.structure.hierarchy.current.trajectories.map((t) => t.cell.transform.ref)
        p.plugin.managers.structure.focus.clear()                     // no focus representation of the previous structure survives the switch
        await p.visual.update({ moleculeId: pdbId }, false)          // fullLoad=false: no plugin.clear(), canvas settings untouched
        const oldRefs = ctx.managers.structure.hierarchy.current.trajectories.filter((t) => old.includes(t.cell.transform.ref))
        if (oldRefs.length) await ctx.managers.structure.hierarchy.remove(oldRefs, false)
        await p.canvas.setBgColor(WHITE)
      }
      await applyOrientation(p, applied)
      await applyMode(p, region ? mode : 'complete')
      // "Show ligands" availability: small-molecule instances within the radius of the selected chain (before the class filter)
      const candidates = instancesOf(p, nearChain([]))
      if (alive && plugin.current === p) setLigandInfo({ candidates, shown: null, excludedNearby: [] })
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), structuresLoaded: p.plugin.managers.structure.hierarchy.current.structures.length }
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), orientationApplied: applied ? { pdbId, rmsd: applied.rmsd, n_pairs: applied.n_pairs, fit_quality: applied.fit_quality } : null, bgColor: p.plugin.canvas3d?.props.renderer.backgroundColor }
      if (alive && plugin.current === p) { window.clearTimeout(watchdog); setLoadedPdb(pdbId); setStatus('ready') }
    }).catch((err: unknown) => { if (alive) { console.error('[PdbeMolstarViewer] load', err); window.clearTimeout(watchdog); setErrorMsg(`Could not load ${pdbId.toUpperCase()}: ${err instanceof Error ? err.message : String(err)}`); setStatus('error') } })
    return () => { alive = false; window.clearTimeout(watchdog) }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pdbId, labelAsymId, region?.start, region?.end, applied, pluginGen, retry])

  // another RNA instance: the ligand toggle starts off again and its availability is re-evaluated after loading
  useEffect(() => { setShowLigands(false); setLigandInfo(null) }, [pdbId, labelAsymId])
  useEffect(() => {
    const p = plugin.current
    if (!p || status !== 'ready' || loadedPdb !== pdbId) return
    queue.current = queue.current.catch(() => undefined).then(() => applyLigands(p, modeRef.current)).catch((err: unknown) => console.error('[PdbeMolstarViewer] ligands', err))
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [showLigands])
  useEffect(() => {
    const p = plugin.current
    if (!p || status !== 'ready' || loadedPdb !== pdbId) return
    queue.current = queue.current.catch(() => undefined).then(() => applyMode(p, mode, false)).catch((err: unknown) => console.error('[PdbeMolstarViewer] mode', err))
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [mode])

  const key = selection.map((s) => s.label_seq_id).sort((a, b) => a - b).join(',')
  useEffect(() => {
    const p = plugin.current
    if (!p || status !== 'ready' || loadedPdb !== pdbId) return
    let alive = true
    queue.current = queue.current.catch(() => undefined).then(async () => {
      if (!alive || plugin.current !== p) return
      await p.visual.clearSelection()
      const data = key ? key.split(',').map((n) => ({ label_asym_id: labelAsymId, label_seq_id: +n, color: SELECT_COLOR })) : []
      selRef.current = data
      // atom-level close-up: Mol*'s structure focus restricted to the selected nucleotides (ball-and-stick, no surroundings or contact lines; see
      // FOCUS_COMPONENTS) replaces any previous focus, so repeated selections never accumulate representations; cleared when nothing is selected
      const loci = data.length ? lociFor(p, data) : null
      if (loci) p.plugin.managers.structure.focus.setFromLoci(loci); else p.plugin.managers.structure.focus.clear()
      if (data.length) await p.visual.select({ data })
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), lastMolstarSelection: { pdbId, labelAsymId, labelSeqIds: data.map((d) => d.label_seq_id), color: SELECT_HEX }, molstarMode: mode, focusRepresentation: focusedSeqIds(p) }
      if (loci && focusOnSelection) await focusOn(p, data)   // clearing never touches the camera
    }).catch((err: unknown) => console.error('[PdbeMolstarViewer] selection', err))
    return () => { alive = false }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [key, labelAsymId, pdbId, status, loadedPdb])

  // explicit deselection (same nucleotide/loop clicked again, Clear selection): no highlight, no focus, default camera of the current mode
  const lastReset = useRef(resetKey)
  useEffect(() => {
    if (resetKey === lastReset.current) return
    lastReset.current = resetKey
    const p = plugin.current
    if (!p || status !== 'ready' || loadedPdb !== pdbId) return
    queue.current = queue.current.catch(() => undefined).then(async () => {
      if (plugin.current !== p) return
      p.plugin.managers.structure.focus.clear()
      await p.visual.clearSelection()
      selRef.current = []
      await resetView(p, region ? mode : 'complete')
      window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), focusRepresentation: focusedSeqIds(p), molstarReset: Number(window.__rnaExplorer?.molstarReset ?? 0) + 1 }
    }).catch((err: unknown) => console.error('[PdbeMolstarViewer] reset', err))
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [resetKey])

  const uniformColor = () => (structureColor != null ? { color: 'uniform' as const, colorParams: { value: Color(structureColor) } } : {})
  /** colour theme of every representation in the scene (test hook): parent component label, theme name, uniform value, visibility */
  const representationColours = (p: PDBeMolstarPlugin) => p.plugin.managers.structure.hierarchy.current.structures.flatMap((st) => st.components.flatMap((c) => c.representations.map((r) => {
    const ct = (r.cell.transform.params as { colorTheme?: { name: string; params?: { value?: number } } } | undefined)?.colorTheme
    return { component: c.cell.obj?.label ?? '', theme: ct?.name ?? null, value: ct?.params?.value ?? null, hidden: !!c.cell.state.isHidden }
  })))

  /** Restrict Mol*'s focus representation to the target (selected nucleotides, heavy atoms): no surroundings, no non-covalent contact lines. */
  const setFocusComponents = async (p: PDBeMolstarPlugin) => {
    const ctx = p.plugin
    const cell = ctx.state.behaviors.selectQ((q) => q.ofTransformer(StructureFocusRepresentation))[0]
    if (cell) await ctx.state.updateBehavior(StructureFocusRepresentation, (old) => {
      old.components = [...FOCUS_COMPONENTS]; old.ignoreHydrogens = true            // heavy atoms only: less clutter on large loops
      // contrasting target: uniform dark selection colour, slightly thicker ball-and-stick than Mol*'s default
      old.targetParams = createStructureRepresentationParams(ctx, undefined, { type: 'ball-and-stick', typeParams: { sizeFactor: 0.32 }, color: 'uniform', colorParams: { value: Color(SELECT_HEX) } })
    })
    const now = ctx.state.behaviors.selectQ((q) => q.ofTransformer(StructureFocusRepresentation))[0]
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), focusComponents: (now?.transform.params as { components?: string[] } | undefined)?.components ?? null, focusIgnoreHydrogens: (now?.transform.params as { ignoreHydrogens?: boolean } | undefined)?.ignoreHydrogens ?? null }
  }

  /** Loci of the given residues in the displayed (transformed) structure; null when none of them is present (no coordinates). */
  const lociFor = (p: PDBeMolstarPlugin, items: { label_asym_id: string; label_seq_id: number }[]) => {
    const data = visibleStructureData(p)
    if (!data || !items.length) return null
    const expr = MS.struct.generator.atomGroups({
      'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), items[0].label_asym_id]),
      'residue-test': MS.core.set.has([MS.set(...items.map((i) => i.label_seq_id)), MS.struct.atomProperty.macromolecular.label_seq_id()]),
    })
    const loci = StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(expr, data))
    return StructureElement.Loci.isEmpty(loci) ? null : loci
  }
  /** label_seq_ids of the residues in Mol*'s current focus representation (test hook). */
  const focusedSeqIds = (p: PDBeMolstarPlugin): number[] => {
    const cur = p.plugin.managers.structure.focus.current
    if (!cur) return []
    const ids = new Set<number>()
    StructureElement.Loci.forEachLocation(cur.loci, (l) => { ids.add(StructureProperties.residue.label_seq_id(l)) })
    return [...ids].sort((a, b) => a - b)
  }
  /** Close-up of the selected residues, keeping the current view direction: the nucleotide fills the view with its immediate neighbours as context. */
  const focusOn = async (p: PDBeMolstarPlugin, items: { label_asym_id: string; label_seq_id: number }[]) => {
    const ctx = p.plugin
    const loci = lociFor(p, items)
    if (!loci || !ctx.canvas3d) return
    const sphere = StructureElement.Loci.getBoundary(loci).sphere
    // focusLoci keeps the camera direction/up and moves target + distance (smooth transition); the margin keeps the 5 A surroundings in view
    ctx.managers.camera.focusLoci(loci, { extraRadius: 7, minRadius: 12, durationMs: 250 })
    window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), lastFocus: items.map((d) => d.label_seq_id), lastFocusCenter: [...sphere.center], lastFocusRadius: sphere.radius }
  }
  const retryLoad = () => {
    const p = plugin.current
    if (!p) return
    queue.current = queue.current.catch(() => undefined).then(() => p.clear()).catch(() => undefined)
    setLoadedPdb(null); setRetry((r) => r + 1)
  }
  const resetOrientation = () => { const p = plugin.current; if (p) queue.current = queue.current.catch(() => undefined).then(() => resetView(p, region ? mode : 'complete')) }

  return (
    <div ref={outer} className={`relative flex flex-col overflow-hidden rounded border border-stone-200 bg-white ${className}`} data-testid="molstar-viewer" data-pdb={status === 'ready' ? loadedPdb ?? pdbId : pdbId} data-status={status} data-mode={region ? mode : 'complete'} data-oriented={applied ? 'true' : 'false'}>
      {/* controls in a toolbar above the canvas (not over it), so they never cover the structure or Mol*'s own canvas buttons */}
      <div className="flex shrink-0 flex-wrap items-center gap-1.5 border-b border-stone-200 bg-stone-50 px-2 py-1.5 text-xs" data-testid="molstar-toolbar">
        {region && (
          <div className="flex overflow-hidden rounded border border-stone-300 bg-white" role="group" aria-label="View mode">
            <button data-testid="mode-region" onClick={() => setMode('region')} aria-pressed={mode === 'region'} className={`px-2 py-0.5 ${mode === 'region' ? 'bg-blue-600 text-white' : 'text-stone-700'}`}>Rfam region</button>
            <button data-testid="mode-complete" onClick={() => setMode('complete')} aria-pressed={mode === 'complete'} className={`px-2 py-0.5 ${mode === 'complete' ? 'bg-blue-600 text-white' : 'text-stone-700'}`}>Complete structure</button>
          </div>
        )}
        {(() => {
          const none = ligandInfo !== null && ligandInfo.candidates.length === 0
          const shown = ligandInfo?.shown
          const title = none ? `No small-molecule ligand within ${LIGAND_RADIUS} Å of this RNA chain in the loaded coordinates`
            : `Show ligand instances within ${LIGAND_RADIUS} Å of the selected RNA as atom-coloured sticks. Complete structure also shows other model ligands. Excludes water, ions and crystallisation additives.`
          // the shown ligands are named in the tooltip only (not next to the checkbox), so the toolbar stays compact
          const codes = shown ? [...new Set(shown.map((x) => x.split(' ')[0]))] : []
          // status notes are in the tooltip too (not beside the checkbox): the toolbar stays on one line in every family
          const onlyExcluded = showLigands && !!shown && shown.length === 0
          const tip = title + (showLigands && codes.length ? `\nShown: ${codes.join(', ')}` : '') + (onlyExcluded ? `\nOnly ions or crystallisation additives are within ${LIGAND_RADIUS} Å.` : '')
          return <label className={`flex items-center gap-1 rounded border border-stone-300 bg-white px-2 py-0.5 ${none ? 'text-stone-400' : 'text-stone-700'}`} title={tip} data-testid="show-ligands-control" data-only-excluded={String(onlyExcluded)}
            data-candidates={ligandInfo?.candidates.length ?? ''} data-shown={shown ? shown.join(';') : ''}>
            <input type="checkbox" checked={showLigands} disabled={none || ligandInfo === null || status !== 'ready'} onChange={(e) => setShowLigands(e.target.checked)} data-testid="show-ligands" />
            Show ligands
          </label>
        })()}
        <button data-testid="reset-orientation" onClick={resetOrientation} className="rounded border border-stone-300 bg-white px-2 py-0.5 text-stone-700" title="Reset orientation and fit the view">Reset view</button>
        {onClearSelection && <button data-testid="molstar-clear-selection" onClick={onClearSelection} disabled={!canClear} className="rounded border border-stone-300 bg-white px-2 py-0.5 text-stone-700 disabled:opacity-40">Clear selection</button>}
        {status === 'error' && <span className="flex items-center gap-2 rounded bg-red-50 px-2 py-0.5 text-red-700" data-testid="molstar-error">{errorMsg ?? `Could not load ${pdbId.toUpperCase()}`}<button data-testid="molstar-retry" onClick={retryLoad} className="rounded border border-red-300 bg-white px-1.5 py-0.5 text-red-700">Retry</button></span>}
      </div>
      <div className="relative min-h-0 flex-1" data-testid="molstar-canvas-area">
        <div ref={host} className="absolute inset-0 bg-white" />
        {status !== 'ready' && <div className="pointer-events-none absolute inset-0 z-[5] bg-white" data-testid="molstar-loading-overlay" />}
        {/* on the canvas, not in the toolbar: the toolbar height (and so the canvas size the camera is fitted to) never changes while loading */}
        {status === 'loading' && <span className="absolute left-2 top-2 z-10 rounded bg-white/90 px-2 py-0.5 text-xs text-stone-600">Loading {pdbId.toUpperCase()}…</span>}
      </div>
    </div>
  )
}
