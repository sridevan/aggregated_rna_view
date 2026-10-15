import type { InstanceRow } from '../types/backend'

/** Species-level name resolved by the backend from NCBI taxonomy (strain/subspecies ids collapsed to the species); falls back to the source name. */
export function displayOrganism(r: Pick<InstanceRow, 'organism' | 'organism_species'>): { name: string | null; isSpecies: boolean } {
  if (r.organism_species) return { name: r.organism_species, isSpecies: r.organism_species !== 'synthetic construct' }
  return { name: r.organism, isSpecies: false }
}
export const pct = (x: number | null | undefined, digits = 1) => (x == null ? '–' : `${(x * 100).toFixed(digits)} %`)
