import { useFamily } from '../state/family'
import { useSelection } from '../state/selection'
import { displayOrganism, pct } from '../lib/format'
import { StructureExplorer } from './StructureExplorer'
import { useLayoutEffect, useRef, useState } from 'react'
import { SpeciesChart } from './SpeciesChart'
import { SectionCard } from './SectionCard'
import { InfoTip } from './InfoTip'
import type { InstanceRow } from '../types/backend'

function RepCard({ rep, rule, testId, className = '' }: { rep: InstanceRow | null | undefined; rule: string; testId: string; className?: string }) {
  const o = rep ? displayOrganism(rep) : null
  const Row = ({ k, v, testId: t }: { k: string; v: React.ReactNode; testId?: string }) => <><dt className="text-xs text-stone-500">{k}</dt><dd className="text-right text-xs text-stone-900" data-testid={t}>{v}</dd></>
  return (
    <SectionCard eyebrow="Structure" title="Representative structure" testId={testId} className={className}
      actions={<InfoTip label="How the representative is chosen" text={`Selected using ${rule}. A coverage-based representative is retained as fallback.`} />}
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

/**
 * Compact family summary (one card): the Rfam description clamped to three lines with Show more / Show less (only when it overflows),
 * the three primary counts, and one lightweight row with the other summary statistics (definitions behind an info icon). Values and their
 * observation definitions are unchanged from family.json / consensus_stats.json.
 */
function FamilySummary() {
  const { family, stats } = useFamily()
  const textRef = useRef<HTMLParagraphElement>(null)
  const [expanded, setExpanded] = useState(false)
  const [overflows, setOverflows] = useState(false)
  useLayoutEffect(() => {
    const el = textRef.current; if (!el) return
    const measure = () => { if (!expanded) setOverflows(el.scrollHeight > el.clientHeight + 1) }
    measure()
    if (typeof ResizeObserver === 'undefined') return
    const ro = new ResizeObserver(measure); ro.observe(el); return () => ro.disconnect()
  }, [expanded, family.comment])
  const Primary = ({ value, label, testId, title, note }: { value: string | number; label: string; testId: string; title?: string; note?: string }) => (
    <div className="flex items-baseline gap-1.5" data-testid={`${testId}-card`} title={title}>
      <span className="text-lg font-semibold leading-6 tabular-nums text-stone-900" data-testid={testId}>{value}</span>
      <span className="text-xs text-stone-600">{label}{note && <span className="text-stone-400"> · {note}</span>}</span>
    </div>)
  const Secondary = ({ value, label, testId }: { value: string | number; label: string; testId: string }) => (
    <span className="whitespace-nowrap" data-testid={`${testId}-card`}><span className="font-medium tabular-nums text-stone-800" data-testid={testId}>{value}</span> {label}</span>)
  const definitions = [
    `RNA instances: ${family.counts.instances} mapped chains (${family.counts.instances_validation_ok} validated); chains from the same PDB entry count separately.`,
    `Unique sequences: ${family.counts.unique_mapped_sequences} distinct mapped-region sequences.`,
    `Mean sequence occupancy: average fraction of ${stats.sequence_level.n_unique_sequences} unique sequences with a residue across ${family.clen} Rfam positions.`,
    `Mean structural occupancy: average fraction of ${stats.structure_level.n_instances} RNA instances with coordinates across ${family.clen} Rfam positions.`,
  ].join('\n\n')
  return (
    <SectionCard title="About this family" testId="family-description"
      description={<>Source: Rfam family entry (<code>comment</code> field). <a className="text-blue-700 underline" href={family.rfam_url ?? `https://rfam.org/family/${family.rfam_acc}`} target="_blank" rel="noreferrer" data-testid="rfam-link">View {family.rfam_acc} at Rfam</a></>}>
      <div data-testid="family-about">
        {family.comment
          ? <p ref={textRef} id="family-about-text" className={`text-sm leading-5 text-stone-700 ${expanded ? '' : 'line-clamp-3'}`} data-testid="family-about-text" data-expanded={String(expanded)}>{family.comment}</p>
          : <p ref={textRef} className="text-sm text-stone-500" data-testid="family-about-text">{family.description}</p>}
        {family.comment && (overflows || expanded) && (
          <button type="button" className="mt-0.5 text-xs text-blue-700 hover:underline" aria-expanded={expanded} aria-controls="family-about-text" onClick={() => setExpanded((e) => !e)} data-testid="family-about-toggle">{expanded ? 'Show less' : 'Show more'}</button>)}
      </div>
      {/* primary counts, then one lightweight row with the remaining statistics (same values and definitions as before) */}
      <div className="mt-2 flex flex-wrap items-baseline gap-x-6 gap-y-1 border-t border-stone-100 pt-2" data-testid="stat-row">
        <Primary value={family.counts.pdb_entries} label="PDB entries" testId="stat-pdb-entries" />
        <Primary value={family.counts.instances} label="RNA chain instances" testId="stat-rna-instances" note={`${family.counts.instances_validation_ok} validated`} />
        <Primary value={family.clen} label="Rfam consensus positions" testId="stat-model-length" title="Match-state positions in the Rfam covariance model." />
        <div className="flex flex-wrap items-center gap-x-4 gap-y-0.5 text-xs text-stone-600 lg:ml-auto" data-testid="stat-secondary">
          <Secondary value={family.counts.unique_mapped_sequences} label="unique sequences" testId="stat-unique-sequences" />
          <Secondary value={pct(family.coverage_summary.sequence_coverage_mean)} label="mean sequence occupancy" testId="stat-sequence-occupancy" />
          <Secondary value={pct(family.coverage_summary.structural_coverage_mean)} label="mean structural occupancy" testId="stat-structural-occupancy" />
          <InfoTip label="How these statistics are defined" text={definitions} />
        </div>
      </div>
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
      <FamilySummary />
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
            actions={<InfoTip label="How experimental methods are counted" text={`Counts ${family.counts.instances} mapped RNA chains, including multiple chains from one PDB entry. Percentages use RNA instances; sequence statistics use ${stats.sequence_level.n_unique_sequences} distinct mapped sequences.`} />}>
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
