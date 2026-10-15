import type { InstanceRow } from '../types/backend'

export type CountUnit = 'entries' | 'instances'
export type SpeciesCategory = 'species' | 'synthetic' | 'unresolved' | 'unknown'
export interface SpeciesRow {
  key: string; label: string; category: SpeciesCategory; italic: boolean
  entries: number; instances: number; speciesTaxId: string | null
  /** for the "Other identified species" row: how many species were folded in */
  groupedSpecies?: number
  /** original organism names seen (for tooltips) */
  sourceNames: string[]
}
export interface SpeciesSummary {
  rows: SpeciesRow[]; unit: CountUnit
  totals: { entries: number; instances: number }
  identifiedSpecies: number; identifiedEntries: number; identifiedInstances: number
}

type Bucket = { key: string; label: string; category: SpeciesCategory; speciesTaxId: string | null; pdbIds: Set<string>; instances: number; sourceNames: Set<string> }

/** Classify one instance. Identification is taken only from the backend's NCBI-resolved species (never inferred): a species-rank tax id
 *  gives 'species'; synthetic constructs are kept apart; a tax id without a species-rank ancestor is 'unresolved'; no tax id = 'unknown'. */
export function classify(r: Pick<InstanceRow, 'organism' | 'organism_species' | 'species_tax_id' | 'tax_id'>): { category: SpeciesCategory; key: string; label: string; speciesTaxId: string | null } {
  const name = (r.organism_species ?? r.organism ?? '').trim()
  if (name.toLowerCase() === 'synthetic construct' || String(r.tax_id) === '32630') return { category: 'synthetic', key: 'synthetic', label: 'Synthetic construct', speciesTaxId: null }
  if (r.species_tax_id && r.organism_species) return { category: 'species', key: `tax:${r.species_tax_id}`, label: r.organism_species, speciesTaxId: String(r.species_tax_id) }
  if (r.tax_id != null && r.tax_id !== '') return { category: 'unresolved', key: 'unresolved', label: 'Unresolved taxonomy (no species-rank name)', speciesTaxId: null }
  return { category: 'unknown', key: 'unknown', label: 'Unknown taxonomy', speciesTaxId: null }
}

/** Aggregate instances per species. PDB entries are counted once per species even when several chains of that entry share the species. */
export function summariseSpecies(rows: InstanceRow[], unit: CountUnit = 'entries', top = 10): SpeciesSummary {
  const buckets = new Map<string, Bucket>()
  const allPdb = new Set<string>()
  for (const r of rows) {
    allPdb.add(r.pdb_id)
    const c = classify(r)
    let b = buckets.get(c.key)
    if (!b) { b = { key: c.key, label: c.label, category: c.category, speciesTaxId: c.speciesTaxId, pdbIds: new Set(), instances: 0, sourceNames: new Set() }; buckets.set(c.key, b) }
    b.pdbIds.add(r.pdb_id); b.instances += 1
    if (r.organism) b.sourceNames.add(r.organism)
  }
  const toRow = (b: Bucket): SpeciesRow => ({ key: b.key, label: b.label, category: b.category, italic: b.category === 'species', entries: b.pdbIds.size, instances: b.instances, speciesTaxId: b.speciesTaxId, sourceNames: [...b.sourceNames].sort() })
  const value = (x: SpeciesRow) => (unit === 'entries' ? x.entries : x.instances)
  const species = [...buckets.values()].filter((b) => b.category === 'species').map(toRow).sort((a, b) => value(b) - value(a) || a.label.localeCompare(b.label))
  const shown = species.slice(0, top)
  const rest = species.slice(top)
  const out: SpeciesRow[] = [...shown]
  if (rest.length) {
    const pdb = new Set<string>()
    for (const b of [...buckets.values()].filter((b) => b.category === 'species' && rest.some((r) => r.key === b.key))) for (const id of b.pdbIds) pdb.add(id)
    out.push({ key: 'other', label: `Other identified species (${rest.length})`, category: 'species', italic: false, entries: pdb.size, instances: rest.reduce((s, r) => s + r.instances, 0), speciesTaxId: null, groupedSpecies: rest.length, sourceNames: rest.map((r) => r.label) })
  }
  for (const k of ['synthetic', 'unresolved', 'unknown']) { const b = buckets.get(k); if (b) out.push(toRow(b)) }
  return {
    rows: out, unit, totals: { entries: allPdb.size, instances: rows.length },
    identifiedSpecies: species.length,
    identifiedEntries: new Set([...buckets.values()].filter((b) => b.category === 'species').flatMap((b) => [...b.pdbIds])).size,
    identifiedInstances: species.reduce((s, r) => s + r.instances, 0),
  }
}
