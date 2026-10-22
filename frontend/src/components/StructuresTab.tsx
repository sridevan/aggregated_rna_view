import { useEffect, useMemo, useRef, useState } from 'react'
import { useFamily } from '../state/family'
import { useSelection } from '../state/selection'
import { useStructures } from '../state/structures'
import { displayOrganism } from '../lib/format'
import { familyFileUrl, loadLigandClasses, loadLigands, loadModifications, loadStructureClustering } from '../api/data'
import { clusteringFor, dissimilarityMatrix, fmt, instanceLabel, viewerColours, type Clustering } from '../lib/structureGroups'
import type { InstanceRow, LigandsDataset, ModificationsDataset, StructureClusteringDataset, StructureMetricId } from '../types/backend'
import { InstanceTable, type TableGroups } from './InstanceTable'
import { InfoTip } from './InfoTip'
import { GroupsCard } from './structures/GroupsCard'
import { GroupOverview } from './structures/GroupOverview'
import { PairwiseHeatmap } from './structures/PairwiseHeatmap'
import { AlignmentViewer, type ViewerEntry } from './structures/AlignmentViewer'
import { AssociatedFeatures } from './structures/AssociatedFeatures'
import { DissimilarityHeatmap } from './structures/DissimilarityHeatmap'

/*
 * Structures tab (layout of the PDBe-KB Complexes "Structural groups" prototype, antibody_antigen_interface/app-react/src/groups/
 * StructuralGroupsApp.jsx): A = RNA instance table + Structural groups panel (compact) or full-width table (expanded); B = pairwise
 * structural comparison (metric selector; collapsed group dendrogram above detail_max instances, with drill-down into groups of
 * <= detail_max) + 3D alignment view; C = associated features. Groups, cut, medoids and statistics come precomputed from
 * structure_clustering.json; switching metric switches to that metric's own groups.
 * Two modes (segmented control): Structural Dissimilarity (default; no groups panel, full-width table without the group column,
 * the geometric discrepancy matrix in the stored TPPL seriation order instead of the dendrogram/heatmap, the same 3D alignment view; a
 * click below the diagonal displays that pair through the existing displayed-structures state) and Clustering (all of the above).
 * Without a stored seriation the tab falls back to Clustering.
 */
const METRIC_IDS: StructureMetricId[] = ['discrepancy', 'tm']

function useClusteringData(acc: string) {
  const [s, set] = useState<{ data: StructureClusteringDataset | null; ligands: LigandsDataset | null; mods: ModificationsDataset | null; ready: boolean }>({ data: null, ligands: null, mods: null, ready: false })
  useEffect(() => {
    let alive = true
    Promise.all([loadStructureClustering(acc), loadLigands(acc).catch(() => null), loadModifications(acc)]).then(([data, ligands, mods]) => { if (alive) set({ data, ligands, mods, ready: true }) })
    return () => { alive = false }
  }, [acc])
  return s
}

export function StructuresTab() {
  const { acc, family, instances, orientation } = useFamily()
  const sel = useSelection()
  const st = useStructures()
  const { data, ligands, mods, ready } = useClusteringData(acc)
  const current = instances.find((i) => i.seq_name === sel.instance) ?? family.display_representative ?? family.representative
  const clusterings = useMemo(() => (data ? Object.fromEntries(METRIC_IDS.map((m) => [m, clusteringFor(data, m)])) as Record<StructureMetricId, Clustering> : null), [data])
  const c = clusterings?.[st.metric] ?? null
  // built once per loaded dataset (not on mode switches or renders); counted for the tests
  const dm = useMemo(() => {
    if (!data) return null
    if (typeof window !== 'undefined') { const w = window.__rnaExplorer ?? (window.__rnaExplorer = {}); w.dissimilarityBuilds = ((w.dissimilarityBuilds as number) ?? 0) + 1 }
    return dissimilarityMatrix(data)
  }, [data])
  const dissim = st.view === 'dissimilarity' && !!dm
  // Structural Dissimilarity mode has its own 3D content (st.dissimDisplayed), never the Clustering one: by default the first instance
  // of the stored TPPL order (also the table's first row and the first diagonal cell), one instance chosen in the table or on the
  // diagonal, or exactly one pair chosen below the diagonal. Selections are seq_names, never row or matrix positions.
  const rep = family.display_representative?.seq_name ?? null                      // the family's designated representative (label tag only)
  const first = dm?.ids[0] ?? null
  const dissimShown = st.dissimDisplayed ?? (first ? [first] : [])
  const pair: [string, string] | null = dissimShown.length === 2 ? [dissimShown[0], dissimShown[1]] : null
  const [holdPage, setHoldPage] = useState(false)                                  // Reset keeps the table page
  const selectPair = (a: string, b: string) => { setHoldPage(false); st.setDissimDisplayed([a, b]) }
  const showOne = (sn: string) => { setHoldPage(false); sel.setInstance(sn); st.setDissimDisplayed([sn]) }
  // the "Selected:" line names what is displayed in Structural Dissimilarity mode (one instance or the pair; no structural groups there)
  const hdr = dissim ? (pair ? null : instances.find((r) => r.seq_name === dissimShown[0]) ?? null) : current
  const ho = hdr ? displayOrganism(hdr) : null
  const resetDissim = () => { setHoldPage(true); st.setDissimDisplayed(null); if (first) sel.setInstance(first) }
  const groupId = c ? (st.groupByMetric[st.metric] ?? c.groups[0].group_id) : null
  const group = c?.groups.find((g) => g.group_id === groupId) ?? c?.groups[0] ?? null
  const excluded = useMemo(() => new Map((data?.dataset.excluded ?? []).map((e) => [e.seq_name, e.reason_text])), [data])
  const tableGroups: TableGroups | undefined = c ? {
    of: (sn) => { const g = c.groupOf[sn]; if (g != null) { const gr = c.groups.find((x) => x.group_id === g)!; return { kind: 'group', id: g, label: gr.label, color: gr.color, representative: gr.representative === sn, order: c.order.indexOf(sn) } } return excluded.has(sn) ? { kind: 'excluded', reason: excluded.get(sn)! } : null },
    selectedId: groupId, onSelectGroup: (id) => st.selectGroup(st.metric, id),
  } : undefined
  const expanded = st.tableMode === 'expanded' || dissim   // Structural Dissimilarity: the table takes the groups panel's space
  const [groupOnly, setGroupOnly] = useState(false)
  const restrict = useMemo(() => (groupOnly && c && group ? (r: { seq_name: string }) => c.groupOf[r.seq_name] === group.group_id : null), [groupOnly, c, group])
  const curGroup = current && c ? c.groupOf[current.seq_name] : undefined
  /**
   * Clustering mode: selecting an instance (table row, heatmap or dendrogram leaf) shows that RNA alone in the 3D view. Its structural
   * group becomes the selected group (membership and the displayed instance stay separate states); any representative comparison ends.
   * A group's own representative is shown as "the representative" (displayed = null).
   */
  const followInstance = (sn: string) => {
    sel.setInstance(sn)
    if (!c || !group) return
    const gid = c.groupOf[sn]
    if (gid != null && gid !== group.group_id) st.selectGroup(st.metric, gid)   // also ends a comparison
    else if (st.compareGroup != null) st.setCompareGroup(null)
    st.setDisplayed(gid != null && c.groups.find((g) => g.group_id === gid)!.representative === sn ? null : [sn])
  }
  // what the one 3D viewer shows, in either mode: one structure or exactly two, a caption, the stored pair value, and what Reset restores
  const viewerModel = ((): ViewerModel | null => {
    if (!c || !group || !data) return null
    if (dissim) {
      const pairLine = pair ? comparisonSummary(clusterings!.discrepancy, pair[0], pair[1], data.dataset.n_shared_positions) : null
      return { displayed: dissimShown, caption: pair ? `Comparing ${instanceLabel(pair[0])} and ${instanceLabel(pair[1])}` : 'Select a pair from the heatmap to compare structures.',
        metric: pairLine && { ...pairLine, testId: 'pair-summary' }, isRep: (sn) => sn === rep, onReset: resetDissim, notClustered: null }
    }
    const compare = st.compareGroup != null ? c.groups.find((g) => g.group_id === st.compareGroup) ?? null : null
    if (compare) {
      const line = comparisonSummary(c, group.representative, compare.representative, data.dataset.n_shared_positions)
      return { displayed: [group.representative, compare.representative], caption: `Comparing representatives of structural groups ${group.group_id} and ${compare.group_id}`,
        metric: { ...line, testId: 'compare-summary' }, isRep: () => false, onReset: () => st.setCompareGroup(null), notClustered: null }
    }
    const one = st.displayed?.[0] ?? group.representative
    const gid = c.groupOf[one], og = gid != null ? c.groups.find((g) => g.group_id === gid)! : null
    const caption = !og ? `${instanceLabel(one)} is not clustered; shown on its own` : og.representative === one ? `${og.label} · Showing representative` : `${og.label} · Showing member ${instanceLabel(one)}`
    return { displayed: [one], caption, metric: null, isRep: (sn) => c.groups.some((g) => g.representative === sn), onReset: () => { st.setCompareGroup(null); st.setDisplayed(null) },
      notClustered: !og && excluded.has(one) ? `${instanceLabel(one)} is not clustered (${excluded.get(one)}).` : null }
  })()
  const table = (
    <InstanceTable rows={instances} selected={dissim ? (pair ? null : dissimShown[0] ?? null) : current?.seq_name ?? null} onSelect={dissim ? showOne : followInstance}
      seriation={dissim ? dm!.ids : null} marked={dissim && pair ? pair : []} jumpOnSelect={!(dissim && holdPage)} clen={family.clen} groups={dissim ? undefined : tableGroups} compact={!!c && !expanded} restrict={restrict}
      toolbar={c && <>{dissim && pair && <span className="ml-auto text-xs text-stone-600" data-testid="table-pair-note">Comparing <b className="font-mono">{instanceLabel(pair[0])}</b> ↔ <b className="font-mono">{instanceLabel(pair[1])}</b> (rows marked)</span>}{!dissim && c.ids.length > (data?.visualisation.detail_max ?? 50) && <label className="ml-auto flex items-center gap-1 text-xs text-stone-700"><input type="checkbox" checked={groupOnly} onChange={(e) => setGroupOnly(e.target.checked)} data-testid="group-only" />selected group only</label>}{!dissim && <button type="button" className="rounded border border-stone-300 px-2 py-1 text-xs hover:bg-stone-50 first:ml-auto" onClick={() => st.setTableMode(expanded ? 'compact' : 'expanded')} data-testid="table-mode-toggle">{expanded ? 'Compact view' : 'Expand table'}</button>}</>} />
  )
  return (
    <div className="space-y-4" data-testid="structures-tab" data-view={st.view}>
      {dm && (
        <div className="inline-flex overflow-hidden rounded border border-stone-300 text-sm" role="radiogroup" aria-label="Structures view" data-testid="view-selector">
          {([['dissimilarity', 'Structural Dissimilarity'], ['clustering', 'Clustering']] as const).map(([v, label]) => (
            <button key={v} type="button" role="radio" aria-checked={st.view === v} onClick={() => st.setView(v)} data-testid={`view-${v}`}
              className={`px-3 py-1 ${st.view === v ? 'bg-stone-800 text-white' : 'bg-white text-stone-700 hover:bg-stone-50'}`}>{label}</button>))}
        </div>
      )}
      <div className="flex flex-wrap items-baseline justify-between gap-2">
        {dissim && pair ? <p className="text-sm" data-testid="structures-current">Selected: <span className="font-mono">{instanceLabel(pair[0])}</span> ↔ <span className="font-mono">{instanceLabel(pair[1])}</span></p>
        : <p className="text-sm" data-testid="structures-current">Selected: <span className="font-mono">{hdr?.pdb_id.toUpperCase()}</span> chain {hdr?.auth_asym_id} · <span className={ho?.isSpecies ? 'italic' : ''}>{ho?.name}</span> · {hdr?.experimental_method}{hdr?.resolution != null ? ` ${hdr.resolution} Å` : ''}
          {!dissim && c && current && <> · {curGroup != null ? <span data-testid="current-group">{c.groups.find((g) => g.group_id === curGroup)?.label}</span> : <span className="text-stone-500" data-testid="current-group">not clustered</span>}</>}</p>}
        {data && <p className="text-xs text-stone-600" data-testid="structures-summary">{data.dataset.n_mapped_instances} RNA instances from {data.dataset.n_mapped_entries} PDB entries · {data.dataset.n_clustered_instances} clustered ({data.dataset.n_clustered_entries} PDB entries) · {data.dataset.n_excluded_instances} not clustered<InfoTip label="About clustering coverage" text={`${data.dataset.n_clustered_instances} of ${data.dataset.n_mapped_instances} RNA instances were clustered. Inclusion requires at least ${data.dataset.coverage_threshold} core positions with complete base, C4′ and C3′ coordinates. Comparisons use ${data.dataset.n_shared_positions} positions shared by all included instances.`} /></p>}
      </div>

      {/* A: the table keeps its place in the tree in both modes, so its search, sort and page survive switching. */}
      <div className={`grid gap-4 ${c && !expanded ? 'lg:grid-cols-[minmax(0,1fr)_340px]' : ''}`} data-testid="structures-section-a" data-mode={expanded ? 'expanded' : 'compact'}>
        <div className="min-w-0">{table}</div>
        {c && group && !expanded && <GroupsCard c={c} data={data!} selectedId={group.group_id} onSelectGroup={(id) => st.selectGroup(st.metric, id)} />}
      </div>
      {ready && !data && <p className="rounded border border-stone-200 bg-stone-50 p-3 text-sm text-stone-600" data-testid="no-clustering">No structural clustering is available for this family yet.</p>}

      {c && group && data && (
        <>
          <div className="grid items-stretch gap-4 lg:grid-cols-2" data-testid="structures-section-b">
            {dissim ? <DissimilarityPanel m={dm!} nShared={data.dataset.n_shared_positions} rowOf={(sn) => instances.find((r) => r.seq_name === sn)} shown={dissimShown} onSelectPair={selectPair} onSelectOne={showOne} />
              : <PairwisePanel c={c} data={data} groupId={group.group_id} onSelectGroup={(id) => st.selectGroup(st.metric, id)} metric={st.metric} setMetric={st.setMetric} selectedInstance={current?.seq_name ?? null} onSelectInstance={followInstance} />}
            <AlignmentPanel model={viewerModel!} compareControl={!dissim} c={c} group={group} acc={acc} frameReference={orientation?.reference ? instanceLabel(orientation.reference.seq_name) : null} instanceRow={(sn) => instances.find((r) => r.seq_name === sn)} transforms={orientation?.transforms ?? {}} />
          </div>
          {!dissim && <AssociatedFeatures c={c} group={group} compared={st.compareGroup != null ? c.groups.find((g) => g.group_id === st.compareGroup) ?? null : null} instances={instances} ligands={ligands} mods={mods} clen={family.clen} />}
        </>
      )}
    </div>
  )
}

/** Structural Dissimilarity mode: the discrepancy matrix in TPPL seriation order (no groups, no dendrogram). */
function DissimilarityPanel({ m, nShared, rowOf, shown, onSelectPair, onSelectOne }: { m: NonNullable<ReturnType<typeof dissimilarityMatrix>>; nShared: number; rowOf: (sn: string) => InstanceRow | undefined; shown: string[]; onSelectPair: (a: string, b: string) => void; onSelectOne: (sn: string) => void }) {
  const info = `Pairwise geometric discrepancy for ${m.ids.length} RNA instances at ${nShared} shared Rfam positions; lower means more similar. Rows and columns are ordered by tree-penalized path length seriation (TPPL), independently of structural groups.`
  return (
    <section className="flex min-w-0 flex-col rounded border border-stone-200 bg-white p-3" data-testid="dissimilarity-panel">
      <h2 className="flex items-center text-base font-semibold text-stone-900">Structural dissimilarity<InfoTip label="About the structural dissimilarity matrix" text={info} /></h2>
      <p className="text-xs text-stone-600" data-testid="dissimilarity-caption">{m.ids.length} RNA instances · geometric discrepancy · ordered by similarity (TPPL seriation)</p>
      <div className="mt-2 min-w-0"><DissimilarityHeatmap m={m} nShared={nShared} rowOf={rowOf} shown={shown} onSelectPair={onSelectPair} onSelectOne={onSelectOne} /></div>
      {/* methods kept out of the info popover (short meaning only) */}
      <details className="mt-2 text-xs text-stone-600" data-testid="dissimilarity-methods">
        <summary className="cursor-pointer text-blue-700">Ordering method</summary>
        <p className="mt-1">Tree-penalized path length seriation (TPPL; Aliyev &amp; Zirbel 2023), {m.seriation.repetitions} repetitions, seed {m.seriation.seed}. Similar structures are placed next to each other; the order uses the dissimilarity matrix only, not the structural groups.</p>
      </details>
      <p className="mt-2 text-xs text-stone-500" data-testid="dissimilarity-hint">Click a diagonal cell to show that structure, or a cell below the diagonal to superpose the pair in 3D.</p>
    </section>
  )
}

function PairwisePanel({ c, data, groupId, onSelectGroup, metric, setMetric, selectedInstance, onSelectInstance }: {
  c: Clustering; data: StructureClusteringDataset; groupId: number; onSelectGroup: (id: number) => void; metric: StructureMetricId; setMetric: (m: StructureMetricId) => void
  selectedInstance: string | null; onSelectInstance: (sn: string) => void
}) {
  const [drill, setDrill] = useState(false)
  useEffect(() => setDrill(false), [groupId, metric])
  const { compareGroup } = useStructures()
  // representative comparison (temporary view): the two groups' route is highlighted in the tree, their representatives marked in their colours
  const comparePair: [number, number] | null = compareGroup != null && compareGroup !== groupId && c.groups.some((g) => g.group_id === compareGroup) ? [groupId, compareGroup] : null
  const compareLeaves = (comparePair ?? []).map((id) => { const g = c.groups.find((x) => x.group_id === id)!; return { leaf: c.ids.indexOf(g.representative), color: g.color, label: `${g.label} · ${instanceLabel(g.representative)}` } })
  const figRef = useRef<HTMLDivElement>(null)
  const detailMax = data.visualisation.detail_max
  const n = c.ids.length
  const full = n <= detailMax
  const group = c.groups.find((g) => g.group_id === groupId)!
  const sim = c.metric.kind === 'similarity'
  const q = sim ? '1 − TM-score' : 'geometric discrepancy'
  const cl = c.metric.clustering
  // short meaning only; methods, cut rationale, thresholds and caveats are in Clustering details (Structural groups panel)
  const info = sim
    ? 'TM-score measures 3D similarity at shared Rfam positions (1 = identical). Clustering uses 1 − TM-score, where lower values indicate greater similarity.'
    : 'Geometric discrepancy compares nucleotide positions and orientations at shared Rfam positions. Lower values indicate greater similarity (0 = identical).'
  const dendroGroups = c.groups.map((g) => ({ id: g.group_id, label: g.label, color: g.color, node: g.nodeRef, from: g.from, to: g.to, size: g.n_instances }))
  const download = (name: string, text: string, type: string) => { const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000) }
  const pairsCsv = () => { const rows = ['instance_a,instance_b,' + (sim ? 'tm_score,one_minus_tm' : 'geometric_discrepancy')]; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) rows.push(`${c.ids[i]},${c.ids[j]},` + (sim ? `${c.value[i][j]},${c.dist[i][j].toFixed(6)}` : `${c.dist[i][j]}`)); return rows.join('\n') }
  const groupsCsv = () => ['seq_name,pdb_id,structural_group,representative,dendrogram_order,status', ...data.dataset.instances.map((sn) => { const g = c.groupOf[sn]; const gr = c.groups.find((x) => x.group_id === g)!; return `${sn},${sn.split('_')[0]},${g},${gr.representative === sn},${c.order.indexOf(sn) + 1},clustered` }),
    ...data.dataset.excluded.map((e) => `${e.seq_name},${e.pdb_id},,,,"not clustered: ${e.reason_text}"`)].join('\n')
  const svgOut = () => { const svgs = figRef.current?.querySelectorAll('svg'); if (!svgs?.length) return; download(`${data.rfam_acc}_${metric}_${full || drill ? 'pairwise' : 'groups'}.svg`, new XMLSerializer().serializeToString(svgs[svgs.length > 1 ? 1 : 0]), 'image/svg+xml') }
  return (
    <section className="flex min-w-0 flex-col rounded border border-stone-200 bg-white p-3" data-testid="pairwise-panel" data-metric={metric}>
      <div className="flex flex-wrap items-center gap-2">
        <h2 className="flex items-center text-base font-semibold text-stone-900">Pairwise structural comparison<InfoTip label="About the comparison metric" text={info} /></h2>
        <div className="ml-auto flex items-center gap-2">
          <div className="inline-flex overflow-hidden rounded border border-stone-300 text-xs" role="radiogroup" aria-label="Comparison metric" data-testid="metric-selector">
            {METRIC_IDS.map((m) => <button key={m} type="button" role="radio" aria-checked={metric === m} onClick={() => setMetric(m)} data-testid={`metric-${m}`}
              className={`px-2 py-1 ${metric === m ? 'bg-stone-800 text-white' : 'bg-white text-stone-700 hover:bg-stone-50'}`}>{data.metrics[m].label}</button>)}
          </div>
          <details className="relative" data-testid="structures-download">
            <summary className="cursor-pointer list-none rounded border border-stone-300 px-2 py-1 text-xs hover:bg-stone-50">Download ▾</summary>
            <div className="absolute right-0 z-30 mt-1 w-64 rounded border border-stone-300 bg-white p-1 text-xs shadow">
              <a className="block rounded px-2 py-1 hover:bg-stone-100" href={familyFileUrl(data.rfam_acc, 'structure_clustering.json')} download={`${data.rfam_acc}_structure_clustering.json`} data-testid="dl-json">Clustering data, both metrics (JSON)</a>
              <button type="button" className="block w-full rounded px-2 py-1 text-left hover:bg-stone-100" onClick={() => download(`${data.rfam_acc}_${metric}_pairs.csv`, pairsCsv(), 'text/csv')} data-testid="dl-pairs">Pairwise values, {c.metric.label} (CSV)</button>
              <button type="button" className="block w-full rounded px-2 py-1 text-left hover:bg-stone-100" onClick={() => download(`${data.rfam_acc}_${metric}_groups.csv`, groupsCsv(), 'text/csv')} data-testid="dl-groups">Structural groups and not-clustered instances (CSV)</button>
              <button type="button" className="block w-full rounded px-2 py-1 text-left hover:bg-stone-100" onClick={svgOut} data-testid="dl-svg">Current figure (SVG)</button>
            </div>
          </details>
        </div>
      </div>
      <p className="mt-1 text-xs text-stone-600" data-testid="pairwise-caption">
        {n} RNA instances · {c.groups.length} structural groups · Fixed cut: {cl.cut.height.toFixed(3)}{cl.validation === 'provisional' && <span className="ml-1.5 rounded border border-amber-300 bg-amber-50 px-1 text-[10px] text-amber-800" title={cl.rationale}>provisional</span>}
      </p>
      <div ref={figRef} className="mt-2 min-w-0">
        {full ? (
          <PairwiseHeatmap c={c} root={c.root} groups={dendroGroups} selectedId={groupId} onSelectGroup={onSelectGroup} selectedInstance={selectedInstance} onSelectInstance={onSelectInstance} nPositions={data.dataset.n_shared_positions} compareLeaves={compareLeaves} />
        ) : (
          <>
            <GroupOverview c={c} selectedId={groupId} onSelectGroup={onSelectGroup} quantity={q} selectedInstance={selectedInstance} compare={comparePair} />
            {drill && group.instance_level_view && (
              <div className="mt-2 border-t border-stone-200 pt-2" data-testid="drilldown-view">
                <PairwiseHeatmap c={c} root={group.nodeRef} groups={[{ id: group.group_id, label: group.label, color: group.color, node: group.nodeRef, from: 0, to: group.n_instances - 1, size: group.n_instances }]} selectedId={groupId}
                  onSelectGroup={onSelectGroup} selectedInstance={selectedInstance} onSelectInstance={onSelectInstance} nPositions={data.dataset.n_shared_positions} compareLeaves={compareLeaves} />
              </div>
            )}
          </>
        )}
      </div>
      {/* one short contextual sentence directly below the figure (plus the drill-down control when the selected group is eligible) */}
      <div className="mt-2 flex flex-wrap items-center gap-x-3 gap-y-1 text-xs text-stone-500" data-testid="drilldown">
        <span data-testid="pairwise-note">{full ? `Pairwise ${q} between all clustered RNA instances, in dendrogram order.`
          : drill && group.instance_level_view ? `Pairwise ${q} within ${group.label}, in dendrogram order.`
          : `Groups with more than ${detailMax} instances are shown in collapsed form.`}</span>
        {!full && group.instance_level_view && (
          <button type="button" className="rounded border border-stone-300 px-2 py-0.5 text-stone-700 hover:bg-stone-50" onClick={() => setDrill(!drill)} data-testid="drilldown-toggle" aria-expanded={drill}>
            {drill ? 'Hide the matrix' : `Show the pairwise matrix of ${group.label}`}</button>
        )}
      </div>
    </section>
  )
}

/**
 * The one-line representative comparison: the active metric's stored pairwise value (structure_clustering.json, never recomputed) and the
 * shared positions it was computed on. The chains and groups are already named by the comparison controls and the structure chips.
 */
export function comparisonSummary(c: Pick<Clustering, 'ids' | 'value' | 'metric'>, a: string, b: string, nShared: number) {
  const i = c.ids.indexOf(a), j = c.ids.indexOf(b)
  const v = i >= 0 && j >= 0 ? c.value[i][j] : null
  const metric = v == null ? 'pairwise value not available' : c.metric.kind === 'similarity' ? `TM-score ${fmt(v)} (1 − TM-score ${fmt(1 - v)})` : `${c.metric.label} ${fmt(v)}`
  return { metric: `${metric} · ${nShared} shared Rfam positions`, value: v }
}

/** The 3D alignment ⓘ: the visual superposition onto the family frame with each displayed chain's fit (orientation.json), kept apart from the metric's alignment */
export function alignmentInfo(displayed: string[], fit: (sn: string) => { quality: string; rmsd: number | null } | null, frameReference: string | null, metricKind: string) {
  const fits = displayed.map((sn) => { const f = fit(sn); return `• ${instanceLabel(sn)}: ${f ? `${f.quality} fit${f.rmsd != null ? `, RMSD ${f.rmsd.toFixed(2)} Å` : ''}` : 'no alignment'}` })
  return [`Superposition: each chain is fitted onto ${frameReference ?? 'the family reference'} using C4′ atoms at shared Rfam positions.`, '', 'Fit to the reference (C4′ RMSD):', ...fits, '',
    `This visual superposition may differ from the alignment used to calculate the ${metricKind === 'similarity' ? 'TM-score' : 'geometric discrepancy'}. Poor fits are flagged on the structure chips.`].join('\n')
}

/** What the 3D viewer shows (Clustering or Structural Dissimilarity): never more than two structures. */
interface ViewerModel {
  displayed: string[]
  caption: string
  /** the stored pairwise value of the two displayed structures (null for one structure) */
  metric: { metric: string; value: number | null; testId: 'compare-summary' | 'pair-summary' } | null
  isRep: (seqName: string) => boolean
  onReset: () => void
  notClustered: string | null
}

/**
 * The one 3D panel of both modes: title (3D structure / 3D superposition), caption and stored pair value, the Compare representatives
 * control (Clustering only), the shared viewer and superposition, and read-only colour-coded labels with Reset. No member, visibility or
 * remove controls: the content is always one structure or exactly two.
 */
function AlignmentPanel({ model, compareControl, c, group, acc, frameReference, instanceRow, transforms }: {
  model: ViewerModel; compareControl: boolean
  c: Clustering; group: Clustering['groups'][number]; acc: string; frameReference: string | null; instanceRow: (sn: string) => { pdb_id: string; label_asym_id: string; auth_asym_id: string } | undefined
  transforms: Record<string, import('../types/backend').OrientationTransform>
}) {
  const st = useStructures()
  const [pick, setPick] = useState<number | ''>('')
  const displayed = model.displayed.slice(0, 2)
  const pairShown = displayed.length === 2
  // a group's representative always takes its group colour (as in the dendrogram, groups panel and table); other members take shades of their group colour
  const colours = viewerColours(displayed, c)
  const entries: ViewerEntry[] = displayed.map((sn) => { const r = instanceRow(sn)!; return { seqName: sn, label: instanceLabel(sn), color: colours[sn], pdbId: r.pdb_id, labelAsymId: r.label_asym_id, transform: transforms[sn] ?? null } })
  const quality = (sn: string): string => transforms[sn]?.fit_quality ?? 'missing'
  const others = c.groups.filter((g) => g.group_id !== group.group_id)
  const comparing = compareControl && st.compareGroup != null
  return (
    <section className="flex min-w-0 flex-col rounded border border-stone-200 bg-white p-3" data-testid="alignment-panel" data-state={pairShown ? 'pair' : 'single'}>
      <h2 className="flex items-center text-base font-semibold text-stone-900"><span data-testid="alignment-title">{pairShown ? '3D superposition' : '3D structure'}</span><InfoTip label="About the 3D alignment" title="3D alignment"
        text={alignmentInfo(displayed, (sn) => (transforms[sn] ? { quality: transforms[sn].fit_quality, rmsd: transforms[sn].rmsd } : null), frameReference, c.metric.kind)} /></h2>
      <p className="text-xs text-stone-600" data-testid="alignment-caption">{model.caption}</p>
      {model.metric && <p className="text-xs font-medium text-stone-800" data-testid={model.metric.testId} data-value={model.metric.value ?? ''}><span data-testid={`${model.metric.testId}-metric`}>{model.metric.metric}</span></p>}
      {compareControl && others.length > 0 && (
        <div className="mt-1.5 flex flex-wrap items-center gap-x-2 gap-y-1 rounded border border-stone-200 bg-stone-50 px-2 py-1 text-xs" data-testid="compare-box">
          <span className="text-[10px] font-semibold uppercase tracking-wide text-stone-500">Compare representatives</span>
          <span className="inline-flex items-center gap-1"><span className="inline-block h-2 w-2 rounded-full" style={{ background: group.color }} />{group.label.replace('Structural group', 'Group')} · {instanceLabel(group.representative)}</span>
          <select className="rounded border border-stone-300 bg-white px-1 py-0.5" value={pick === '' ? others[0].group_id : pick} onChange={(e) => setPick(+e.target.value)} data-testid="compare-select" aria-label="Compare with">
            {others.map((g) => <option key={g.group_id} value={g.group_id}>{g.label.replace('Structural group', 'Group')} · {instanceLabel(g.representative)}</option>)}
          </select>
          {comparing ? <button type="button" className="rounded bg-stone-700 px-2 py-0.5 text-white" onClick={() => st.setCompareGroup(null)} data-testid="compare-off">Back to {group.label.replace('Structural group', 'group')}</button>
            : <button type="button" className="rounded bg-blue-700 px-2 py-0.5 text-white" onClick={() => st.setCompareGroup(pick === '' ? others[0].group_id : pick)} data-testid="compare-on" title="Superpose the two group representatives">Compare</button>}
        </div>
      )}
      {/* the viewer takes the vertical space the shared grid row gives the panel (equal outer heights with the pairwise panel) */}
      <div className="mt-1.5 flex min-h-[320px] flex-1 flex-col lg:min-h-[340px]"><AlignmentViewer key={acc} entries={entries} hidden={[]} ligandClasses={() => loadLigandClasses(acc)} /></div>
      {/* read-only labels of what is displayed (colours as in the viewer) and Reset */}
      <div className="mt-1.5 flex flex-wrap items-center gap-1.5 text-[11px]" data-testid="viewer-labels" data-state={pairShown ? 'pair' : 'single'}>
        <span className="font-semibold text-stone-600">{pairShown ? 'Comparing:' : 'Displayed:'}</span>
        {entries.map((e, k) => {
          const r = instanceRow(e.seqName), qy = quality(e.seqName)
          return (
            <span key={e.seqName} className="inline-flex items-center gap-1">
              {k === 1 && <span className="px-0.5 text-stone-500" aria-hidden>↔</span>}
              <span className="inline-flex items-center gap-1 rounded border border-stone-200 bg-white px-1 py-px" data-testid={`viewer-label-${e.seqName}`} title={e.seqName}>
                <span className="inline-block h-2.5 w-2.5 rounded-sm ring-1 ring-inset ring-black/10" style={{ background: e.color }} data-testid={`swatch-${e.seqName}`} data-color={e.color} />
                <span className="font-mono">{r ? `${r.pdb_id.toUpperCase()} ${r.auth_asym_id}` : e.label}</span>
                {!pairShown && model.isRep(e.seqName) && <span className="text-stone-500" data-testid="viewer-rep-label">· Representative</span>}
                {(qy === 'poor' || qy === 'insufficient' || qy === 'missing') && <span className="rounded border border-amber-400 bg-amber-50 px-0.5 text-[10px] text-amber-800" title="The superposition onto the family frame is poor for this chain (C4′ RMSD > 4 Å) or unavailable; its placement may look offset." data-testid={`fit-warning-${e.seqName}`}>{qy === 'missing' ? 'no alignment' : `${qy} fit`}</span>}
              </span>
            </span>
          )
        })}
        <button type="button" className="ml-auto rounded border border-stone-300 px-1.5 py-px text-stone-700 hover:bg-stone-50" onClick={model.onReset} data-testid="viewer-reset" title="Show the representative alone">Reset</button>
        {model.notClustered && <span className="basis-full text-stone-500" data-testid="selected-not-clustered">{model.notClustered}</span>}
      </div>
    </section>
  )
}
