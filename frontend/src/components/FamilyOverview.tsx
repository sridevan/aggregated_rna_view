import { useFamily } from '../state/family'
import { useSelection } from '../state/selection'
import { displayOrganism, pct } from '../lib/format'
import { StructureExplorer } from './StructureExplorer'
import { Stat } from './Stat'
import { SpeciesChart } from './SpeciesChart'
import { SectionCard } from './SectionCard'
import { InfoTip } from './InfoTip'
import type { InstanceRow } from '../types/backend'

function RepCard({ rep, rule, testId, className = '' }: { rep: InstanceRow | null | undefined; rule: string; testId: string; className?: string }) {
  const o = rep ? displayOrganism(rep) : null
  const Row = ({ k, v, testId: t }: { k: string; v: React.ReactNode; testId?: string }) => <><dt className="text-xs text-stone-500">{k}</dt><dd className="text-right text-xs text-stone-900" data-testid={t}>{v}</dd></>
  return (
    <SectionCard eyebrow="Structure" title="Representative structure" testId={testId} className={className}
      actions={<InfoTip label="How the representative is chosen" text={`Selection rule: ${rule}. The coverage-based representative (most matched residues with coordinates) is kept in the data as a fallback.`} />}
      description={rep ? undefined : 'None eligible.'}>
      {rep && o && (
        <>
          <p className="text-base font-semibold leading-tight" data-testid="rep-id"><span className="font-mono">{rep.pdb_id.toUpperCase()}</span> chain {rep.auth_asym_id} · <span className={o.isSpecies ? 'italic font-normal' : 'font-normal'}>{o.name}</span></p>
          <dl className="mt-1.5 grid grid-cols-[auto_1fr] gap-x-4 gap-y-0.5">
            <Row k="Method" v={rep.experimental_method ?? '–'} />
            <Row k="Resolution" v={rep.resolution != null ? `${rep.resolution} Å` : '–'} />
            <Row k="Rfam model coverage" v={`${pct(rep.fraction_clen_matched, 0)} (${rep.n_match} of ${rep.n_match + rep.n_deleted_in_span + rep.n_uncovered} positions)`} />
            <Row k="Matched residues with coordinates" v={`${rep.n_match_with_coordinates} / ${rep.n_match}`} />
            <Row k="PDBe 2D topology" v={rep.pdbe_2d_available ? 'available' : 'not available'} />
          </dl>
        </>
      )}
    </SectionCard>
  )
}

export function FamilyOverview() {
  const { family, stats, instances } = useFamily()
  const sel = useSelection()
  const display = family.display_representative ?? family.representative
  const current = instances.find((i) => i.seq_name === sel.instance) ?? display
  const methods = Object.entries(family.experimental_methods).sort((a, b) => b[1] - a[1])
  return (
    <div className="space-y-4">
      <SectionCard title="About this family" testId="family-description"
        description={<>Source: Rfam family entry (<code>comment</code> field). <a className="text-blue-700 underline" href={family.rfam_url ?? `https://rfam.org/family/${family.rfam_acc}`} target="_blank" rel="noreferrer" data-testid="rfam-link">View {family.rfam_acc} at Rfam</a></>}>
        {family.comment ? <p className="text-sm leading-5 text-stone-700">{family.comment}</p> : <p className="text-sm text-stone-500">{family.description}</p>}
      </SectionCard>
      <section className="grid grid-cols-2 gap-3 md:grid-cols-3 lg:grid-cols-6 lg:items-stretch" data-testid="stat-row">
        <Stat label="PDB entries" value={family.counts.pdb_entries} />
        <Stat label="RNA instances" value={family.counts.instances} hint={`${family.counts.instances_validation_ok} validated`} />
        <Stat label="Unique sequences" value={family.counts.unique_mapped_sequences} hint="identical mapped regions deduplicated" />
        <Stat label="Model length" value={family.clen} hint="consensus match positions" />
        <Stat label="Mean sequence occupancy" value={pct(family.coverage_summary.sequence_coverage_mean)} hint={`mean over ${family.clen} positions of the fraction of ${stats.sequence_level.n_unique_sequences} unique sequences with a residue`} testId="stat-sequence-occupancy" />
        <Stat label="Mean structural occupancy" value={pct(family.coverage_summary.structural_coverage_mean)} hint={`mean over ${family.clen} positions of the fraction of ${stats.structure_level.n_instances} instances with coordinates`} testId="stat-structural-occupancy" />
      </section>
      {/* Row 3: structure explorer (shared instance selector, coverage tracks, 2D + 3D viewers), full width */}
      <SectionCard eyebrow="Molecular view" testId="viewer-card"
        title={current ? <>{current.seq_name === display?.seq_name ? 'Representative RNA' : 'Selected RNA'} — <span className="font-mono">{current.pdb_id.toUpperCase()}</span> chain {current.auth_asym_id}</> : 'RNA structure'}
        description={current ? <span data-testid="overview-current">{current.pdb_id.toUpperCase()} chain {current.auth_asym_id}: explore RNA sequence and structure using synchronised consensus-position selections across the coverage tracks, 2D viewer and Mol*.</span> : undefined}>
        <StructureExplorer />
      </SectionCard>
      {/* Row 4: metadata (40 %) and species diversity (60 %), equal overall height */}
      <div className="grid gap-4 lg:grid-cols-5 lg:items-start" data-testid="metadata-row">
        <div className="flex flex-col gap-4 lg:col-span-2" data-testid="metadata-panel">
          <RepCard rep={family.display_representative ?? family.representative} testId="display-representative"
            rule={family.display_representative?.selection_rule ?? 'coverage representative: clan-best instance with most consensus-matched residues with coordinates, then best resolution'} />
          <SectionCard eyebrow="Statistics" title="Experimental methods" testId="methods-card"
            actions={<InfoTip label="How experimental methods are counted" text={`Counts RNA chain instances (${family.counts.instances} mapped chains), not PDB entries: an entry with several mapped chains contributes each chain. Percentages are of all chain instances. Sequence-level statistics elsewhere count each unique mapped sequence once (${stats.sequence_level.n_unique_sequences}).`} />}>
            <ul className="space-y-1" data-testid="methods-list">
              {methods.map(([m, n]) => (
                <li key={m} className="grid grid-cols-[auto_1fr_auto] items-center gap-x-3">
                  <span className="whitespace-nowrap text-sm">{m}</span>
                  <span className="relative h-3.5 w-full rounded-sm bg-stone-100"><span className="absolute inset-y-0 left-0 rounded-sm bg-blue-600" style={{ width: `${(n / family.counts.instances) * 100}%` }} /></span>
                  <span className="whitespace-nowrap text-sm tabular-nums text-stone-700"><b className="text-stone-900">{n}</b> <span className="text-xs text-stone-500">({((n / family.counts.instances) * 100).toFixed(0)} %)</span></span>
                </li>))}
            </ul>
            <p className="mt-1.5 text-[11px] text-stone-500">Counts RNA chain instances ({family.counts.instances}).</p>
          </SectionCard>
        </div>
        <div className="lg:col-span-3"><SpeciesChart instances={instances} /></div>
      </div>
    </div>
  )
}
