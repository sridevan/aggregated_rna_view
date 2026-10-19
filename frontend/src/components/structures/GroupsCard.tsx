import { useState } from 'react'
import type { Clustering } from '../../lib/structureGroups'
import { instanceLabel } from '../../lib/structureGroups'
import type { StructureClusteringDataset } from '../../types/backend'
import { InfoTip } from '../InfoTip'

/*
 * Structural groups panel, ported from antibody_antigen_interface/app-react/src/groups/GroupPanels.jsx (GroupsCard): group list (largest
 * three, the rest folded; the selected group is never folded out of sight), "largest"/"small"/"single" tags, then the selected group's size,
 * representative (medoid), within-group and across-group dissimilarity with the ratio of the figures as shown, and the clustering details.
 * All figures come from structure_clustering.json (metric-specific); nothing is recomputed in the browser.
 */
export const GROUP_NOTE = 'Structural groups are derived computationally from pairwise structural comparison. They are descriptive and do not necessarily correspond to biological conformational states; differences may reflect sequence, ligands, construct design, crystal environment or modelling.'
const MAX_LISTED = 3
const plural = (n: number, one: string, many: string) => `${n} ${n === 1 ? one : many}`

export function GroupsCard({ c, data, selectedId, onSelectGroup }: { c: Clustering; data: StructureClusteringDataset; selectedId: number; onSelectGroup: (id: number) => void }) {
  const [showAll, setShowAll] = useState(false)
  const groups = c.groups, small = data.visualisation.small_group
  const folds = groups.length > MAX_LISTED + 1
  const hidden = folds ? groups.slice(MAX_LISTED) : []
  const selHidden = hidden.some((g) => g.group_id === selectedId)
  const listed = folds && !showAll ? [...groups.slice(0, MAX_LISTED), ...(selHidden ? [groups.find((g) => g.group_id === selectedId)!] : [])] : groups
  const rest = hidden.filter((g) => g.group_id !== selectedId)
  const sel = groups.find((g) => g.group_id === selectedId)!
  const sim = c.metric.kind === 'similarity'
  // display only: dissimilarities to two decimals; the ratio from the stored (unrounded) means, to one decimal
  const fmt2 = (v: number) => v.toFixed(2)
  const ratio = sel.within.mean != null && sel.across.mean != null && sel.within.mean > 0 ? sel.across.mean / sel.within.mean : null
  const mName = sim ? '1 − TM-score' : 'geometric discrepancy'
  const cl = c.metric.clustering
  return (
    <section className="flex h-full flex-col rounded border border-stone-200 bg-white p-3" data-testid="groups-card">
      <h2 className="flex items-center text-base font-semibold text-stone-900">Structural groups<InfoTip label="About structural groups" text={`Groups of RNA instances with similar 3D structure under the selected metric. They are descriptive and do not necessarily correspond to biological conformational states.`} /></h2>
      <div className="mt-2 space-y-1.5" role="group" aria-label="Structural groups">
        {listed.map((g) => {
          const on = g.group_id === selectedId
          return (
            <button key={g.group_id} type="button" aria-pressed={on} onClick={() => onSelectGroup(g.group_id)} data-testid={`group-chip-${g.group_id}`}
              className={`flex w-full items-center gap-2 rounded border px-2 py-1 text-left text-xs ${on ? 'border-2' : 'border-stone-200 hover:bg-stone-50'}`}
              style={on ? { borderColor: g.color, background: `color-mix(in srgb, ${g.color} 12%, white)` } : undefined}>
              <span className="inline-block h-2.5 w-2.5 shrink-0 rounded-full" style={{ background: g.color }} />
              <b className="text-stone-900">{g.label}</b>
              <span className="text-stone-600">{plural(g.n_instances, 'RNA instance', 'RNA instances')}</span>
              <span className="ml-auto flex gap-1">
                {g.group_id === groups[0].group_id && <span className="rounded-full border border-stone-400 px-1.5 text-[10px] text-stone-700" title="Largest refers only to the number of RNA instances in the group; it does not imply biological importance.">largest</span>}
                {g.n_instances < small && <span className="rounded-full border border-amber-400 bg-amber-50 px-1.5 text-[10px] text-amber-800" title={`Fewer than ${small} RNA instances: possibly an outlier rather than a population; its statistics rest on few pairs.`}>{g.n_instances === 1 ? 'single' : 'small'}</span>}
              </span>
            </button>
          )
        })}
      </div>
      {folds && <button type="button" className="mt-1 self-start text-xs text-blue-700 underline" aria-expanded={showAll} onClick={() => setShowAll(!showAll)}>{showAll ? 'Show fewer groups' : `${rest.length} more groups`}</button>}
      <p className="mt-1 text-[11px] text-stone-500">Largest group selected by default.</p>

      <div className="mt-3 border-t border-stone-200 pt-2 text-xs" data-testid="selected-group">
        <p className="text-[10px] font-semibold uppercase tracking-wide text-stone-500">Selected structural group</p>
        <p className="mt-0.5 flex items-center gap-1.5 text-base font-semibold text-stone-900"><span className="inline-block h-3 w-3 rounded-full" style={{ background: sel.color }} />{sel.label}</p>
        <p className="text-stone-600" data-testid="selected-group-counts">{plural(sel.n_instances, 'RNA instance', 'RNA instances')} · {plural(sel.n_entries, 'PDB entry', 'PDB entries')}</p>
        {sel.n_instances < small && <p className="mt-1 rounded bg-amber-50 px-2 py-1 text-[11px] text-amber-900" role="note">{sel.n_instances === 1 ? 'Single RNA instance: no within-group statistics.' : `Small group of ${sel.n_instances}: its statistics rest on ${sel.within.n_pairs} pairs.`}</p>}
        <p className="mt-2 font-semibold text-stone-800">Representative<InfoTip label="About the representative" text="Group member with the lowest mean dissimilarity to the others, selected among chains without backbone breaks." /></p>
        <p><span className="font-mono" data-testid="group-representative">{instanceLabel(sel.representative)}</span> <span className="rounded bg-stone-100 px-1 text-[10px] text-stone-600">Medoid</span></p>
        {sel.medoid && sel.medoid !== sel.representative && <p className="text-[11px] text-stone-500" data-testid="group-medoid-note">The overall medoid, {instanceLabel(sel.medoid)}, has unmodelled nucleotides inside the chain.</p>}
        {!!sel.representative_internal_gaps?.length && <p className="text-[11px] text-amber-800" data-testid="group-representative-gaps">Every member has unmodelled nucleotides inside the chain; this one lacks {sel.representative_internal_gaps.length} ({sel.representative_internal_gaps.slice(0, 6).join(', ')}{sel.representative_internal_gaps.length > 6 ? ', …' : ''}).</p>}
        <p className="mt-2 font-semibold text-stone-800" ><span data-testid="within-heading">Within-group {mName}</span><InfoTip label="About within-group dissimilarity" text={`Mean pairwise ${mName} within this group; lower indicates greater similarity.`} /></p>
        {sel.within.mean == null ? <p className="text-stone-500">Single RNA instance</p> : (
          <div data-testid="within-stats"><p className="flex justify-between"><span>Mean</span><b className="tabular-nums" data-testid="within-mean">{fmt2(sel.within.mean)}</b></p><p className="flex justify-between"><span>Maximum</span><b className="tabular-nums" data-testid="within-max">{fmt2(sel.within.max!)}</b></p></div>)}
        {sel.across.mean != null && (
          <>
            <p className="mt-2 font-semibold text-stone-800" ><span data-testid="across-heading">Across-group {mName}</span><InfoTip label="About across-group dissimilarity" text={`Mean pairwise ${mName} between this group and all other clustered instances.`} /></p>
            <p className="flex justify-between" data-testid="across-stats"><span>Mean</span><b className="tabular-nums" data-testid="across-mean">{fmt2(sel.across.mean)}</b></p>
            {ratio != null && <p className="text-stone-500" data-testid="across-ratio">{ratio.toFixed(1)}× the within-group mean</p>}
          </>
        )}
        <details className="mt-2" data-testid="clustering-details">
          <summary className="cursor-pointer text-blue-700">Clustering details</summary>
          <div className="mt-1 space-y-0.5 text-[11px] text-stone-700">
            <p>Metric: {c.metric.definition}{sim ? '' : ' Dimensionless; not an RMSD.'}</p>
            <p>Method: {cl.method}, {cl.linkage_method} linkage on {sim ? '1 − TM-score' : 'geometric discrepancy'}</p>
            <p>Fixed cut: {cl.cut.k} groups at {cl.cut.height.toFixed(3)} ({cl.cut.rule})</p>
            <p>Status: <b>{cl.validation}</b>. {cl.rationale}</p>
            <p>Cophenetic correlation {cl.cophenetic_correlation}. Duplicate floor {cl.floor.value?.toFixed(3) ?? 'n/a'} (dotted orange line): {cl.floor.definition}.</p>
            <p>Statistics: within-group {sel.within.n_pairs} pairs, across-group {sel.across.n_pairs} pairs; the ratio uses the unrounded means and is not a statistical test.</p>
            <p>The instance-level matrix and dendrogram are drawn for at most {data.visualisation.detail_max} RNA instances; larger sets are shown as collapsed groups.</p>
            <p data-testid="details-excluded">Not clustered: {data.dataset.excluded.map((e) => `${instanceLabel(e.seq_name)} (${e.reason_text})`).join('; ')}.</p>
            <p>3D alignment view: each RNA chain is placed in the family reference frame of the Family Overview by Kabsch superposition of C4′ atoms at shared Rfam consensus positions onto the family reference structure. It is the same for both metrics and is not the base-centre superposition used for geometric discrepancy; fits with C4′ RMSD above 4 Å are flagged as poor.</p>
            <p>{GROUP_NOTE}</p>
            <p>{data.dataset.n_clustered_instances} of {data.dataset.n_mapped_instances} mapped RNA instances clustered ({data.dataset.n_clustered_entries} of {data.dataset.n_mapped_entries} PDB entries); {data.dataset.n_excluded_instances} not clustered. {data.dataset.inclusion_rule}.</p>
          </div>
        </details>
      </div>
    </section>
  )
}
