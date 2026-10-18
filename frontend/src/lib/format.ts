import type { InstanceRow } from '../types/backend'

/** Species-level name resolved by the backend from NCBI taxonomy (strain/subspecies ids collapsed to the species); falls back to the source name. */
export function displayOrganism(r: Pick<InstanceRow, 'organism' | 'organism_species'>): { name: string | null; isSpecies: boolean } {
  if (r.organism_species) return { name: r.organism_species, isSpecies: r.organism_species !== 'synthetic construct' }
  return { name: r.organism, isSpecies: false }
}
export const pct = (x: number | null | undefined, digits = 1) => (x == null ? '–' : `${(x * 100).toFixed(digits)} %`)

/** Display-only abbreviation of experimental methods (the data keep the full PDBe value; unknown methods are shown unchanged). */
const METHOD_SHORT: Record<string, string> = { 'x-ray diffraction': 'X-ray', 'electron microscopy': 'EM', 'solution nmr': 'NMR' }
export const shortMethod = (m: string | null | undefined): string | null => (m == null ? null : METHOD_SHORT[m.trim().toLowerCase()] ?? m)
