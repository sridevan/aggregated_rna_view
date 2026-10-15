import { useMemo, useState } from 'react'
import { summariseSpecies, type CountUnit, type SpeciesRow } from '../lib/species'
import type { InstanceRow } from '../types/backend'
import { SectionCard } from './SectionCard'
import { InfoTip } from './InfoTip'

const COLOR: Record<string, string> = { species: '#2a78d6', other: '#86b6ef', synthetic: '#a8a69f', unresolved: '#c98500', unknown: '#9b6b6b' }

/** Horizontal bar chart: top species by distinct PDB entries (default) or RNA chain instances; synthetic/unresolved/unknown kept apart. */
export function SpeciesChart({ instances, className = '' }: { instances: InstanceRow[]; className?: string }) {
  const [unit, setUnit] = useState<CountUnit>('entries')
  const s = useMemo(() => summariseSpecies(instances, unit, 10), [instances, unit])
  const total = unit === 'entries' ? s.totals.entries : s.totals.instances
  const value = (r: SpeciesRow) => (unit === 'entries' ? r.entries : r.instances)
  const max = Math.max(1, ...s.rows.map(value))
  const unitLabel = unit === 'entries' ? 'PDB entries' : 'RNA chain instances'
  const toggle = (
    <div className="flex overflow-hidden rounded border border-stone-300 text-xs" role="group" aria-label="Counting unit">
      <button data-testid="unit-entries" aria-pressed={unit === 'entries'} onClick={() => setUnit('entries')} className={`px-2 py-0.5 ${unit === 'entries' ? 'bg-blue-600 text-white' : 'text-stone-700'}`}>PDB entries</button>
      <button data-testid="unit-instances" aria-pressed={unit === 'instances'} onClick={() => setUnit('instances')} className={`px-2 py-0.5 ${unit === 'instances' ? 'bg-blue-600 text-white' : 'text-stone-700'}`}>RNA chain instances</button>
    </div>
  )
  const description = (
    <span data-testid="species-summary" className="flex flex-wrap items-center gap-x-2">
      <span><b className="text-base text-stone-900">{s.identifiedSpecies}</b> distinct identified species</span>
      <span>· {s.identifiedEntries} of {s.totals.entries} PDB entries have species-level assignments ({s.identifiedInstances} of {s.totals.instances} chain instances)</span>
      <InfoTip label="About the species statistics" text={`Species come from the backend's NCBI Taxonomy resolution of each chain's source organism (strain and subspecies ids collapsed to the species rank); nothing is inferred from sequence. Synthetic constructs, tax ids without a species-rank name (unresolved) and chains without taxonomy (unknown) are shown separately. Counting units: PDB entries counts an entry once per species even when several of its chains share that species; RNA chain instances counts every mapped chain. Percentages are of all ${unit === 'entries' ? s.totals.entries + ' entries' : s.totals.instances + ' instances'}; the top 10 identified species are listed and the rest are grouped as "Other identified species".`} />
    </span>
  )
  return (
    <SectionCard eyebrow="Diversity" title="Species distribution" description={description} actions={toggle} className={className} testId="species-chart">
      <ul className="space-y-0.5" role="list" aria-label={`Species by ${unitLabel}`}>
        {s.rows.map((r) => {
          const v = value(r), pctv = total ? (v / total) * 100 : 0
          const fill = r.key === 'other' ? COLOR.other : COLOR[r.category]
          const tip = `${r.label}: ${v} ${unitLabel} = ${pctv.toFixed(1)} % of ${total}` + (unit === 'entries' ? ` (${r.instances} chain instances)` : ` (${r.entries} PDB entries)`) + (r.key !== 'other' && r.sourceNames.length && r.sourceNames.join(',') !== r.label ? `. Source names: ${r.sourceNames.join('; ')}` : '') + (r.key === 'other' ? `. Species: ${r.sourceNames.join('; ')}` : '')
          return (
            <li key={r.key} data-testid={`species-row-${r.key}`} title={tip} className="grid grid-cols-1 items-center gap-x-3 gap-y-0.5 rounded px-1 py-0.5 hover:bg-stone-100 sm:grid-cols-[minmax(10rem,14rem)_1fr_3.5rem]">
              <span data-testid="species-label" data-italic={r.italic ? 'true' : 'false'} className={`truncate text-sm leading-5 text-stone-900 sm:text-right ${r.italic ? 'italic' : ''}`} title={r.label}>{r.label}</span>
              <span className="relative h-3.5 w-full rounded-sm bg-stone-100"><span className="absolute inset-y-0 left-0 rounded-sm" style={{ width: `${(v / max) * 100}%`, backgroundColor: fill }} data-testid="species-bar" data-value={v} /></span>
              <span className="text-sm tabular-nums text-stone-700" data-testid="species-count">{v}</span>
            </li>
          )
        })}
      </ul>
      <p className="mt-1.5 flex flex-wrap gap-x-3 text-[11px] text-stone-500">
        <span><i className="mr-1 inline-block h-2.5 w-2.5 rounded-sm align-middle" style={{ background: COLOR.species }} />identified species</span>
        {s.rows.some((r) => r.key === 'other') && <span><i className="mr-1 inline-block h-2.5 w-2.5 rounded-sm align-middle" style={{ background: COLOR.other }} />other identified species</span>}
        {s.rows.some((r) => r.category === 'synthetic') && <span><i className="mr-1 inline-block h-2.5 w-2.5 rounded-sm align-middle" style={{ background: COLOR.synthetic }} />synthetic construct</span>}
        {s.rows.some((r) => r.category === 'unresolved') && <span><i className="mr-1 inline-block h-2.5 w-2.5 rounded-sm align-middle" style={{ background: COLOR.unresolved }} />unresolved taxonomy</span>}
        {s.rows.some((r) => r.category === 'unknown') && <span><i className="mr-1 inline-block h-2.5 w-2.5 rounded-sm align-middle" style={{ background: COLOR.unknown }} />unknown taxonomy</span>}
        <span className="ml-auto">bar length ∝ {unitLabel}; hover a row for percentages</span>
      </p>
    </SectionCard>
  )
}
