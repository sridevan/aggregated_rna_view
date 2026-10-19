import { lazy, Suspense } from 'react'
import { residueLabel, type SceneEntry } from '../../lib/molstarScene'

/** the Mol* scene is loaded on demand (keeps Mol* out of the initial bundle and of modules that only need the helpers) */
const MolstarSuperposition = lazy(() => import('../molstar/MolstarSuperposition').then((m) => ({ default: m.MolstarSuperposition })))

/*
 * Structures tab 3D view: one structure, or two superposed, in the family frame. Mol* (shared MolstarSuperposition scene, the same
 * engine and camera convention as the Family Overview viewer): each chain is fetched as BinaryCIF from the PDBe ModelServer
 * (label_asym_id, model 1) and placed with its orientation.json transform (Kabsch on C4' onto the family reference, x' = R x + t). That
 * superposition is for visual comparison; it is not the base-centre fit behind the geometric discrepancy. Replaces the former 3Dmol.js view.
 */
export type ViewerEntry = SceneEntry
export { residueLabel }

export function AlignmentViewer({ entries, hidden, height, emphasis = null, ligandClasses }: { entries: ViewerEntry[]; hidden: string[]; height?: number
  /** family ligand classification for Show ligands (ions and additives left out) */
  ligandClasses?: () => Promise<Record<string, string>>
  /** seq_name of the structure to emphasise (the selected instance) when it is displayed with others */
  emphasis?: string | null }) {
  const shown = entries.filter((e) => !hidden.includes(e.seqName))
  const cls = height ? '' : 'min-h-[320px] flex-1 lg:min-h-[340px]'
  return (
    <Suspense fallback={<div className={`relative rounded border border-stone-200 bg-white ${cls}`} style={height ? { height } : undefined} data-testid="alignment-viewer" data-status="loading" data-structures={shown.map((e) => e.seqName).join(',')} />}>
    <MolstarSuperposition entries={shown} emphasis={emphasis} ligandClasses={ligandClasses} testId="alignment-viewer"
      className={cls} style={height ? { height } : undefined}
      onHooks={(h) => { window.__rnaExplorer = { ...(window.__rnaExplorer ?? {}), alignment: { loaded: h.loaded, transformed: h.transformed }, alignmentProbe: h.probe } }} />
    </Suspense>
  )
}
