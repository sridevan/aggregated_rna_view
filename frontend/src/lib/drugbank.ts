/*
 * DrugBank cross-reference of a chemical component (CCD), fetched on demand for the gallery cards on screen.
 * Source order (user rule: PDBe first, RCSB only when PDBe cannot provide the data):
 *   1. PDBe CCD file pdbechem_v2/{CCD}.cif, category _pdbe_chem_comp_drugbank_details (drugbank_id, name);
 *   2. RCSB chem comp record data.rcsb.org/rest/v1/core/chemcomp/{CCD}, rcsb_chem_comp_related with resource_name DrugBank,
 *      only when PDBe has no DrugBank record or cannot be reached.
 * Checked 2026-10-22 for the 72 gallery CCDs: PDBe maps 27, RCSB agrees on all 27 and adds SAM, TAC, MIY and TPS. Both allow cross-origin
 * requests. The mapping is chemical identity only: DrugBank targets (proteins, UniProt) are deliberately not shown.
 */
export interface DrugbankXref { id: string; name: string | null; source: 'PDBe' | 'RCSB'; sourceUrl: string }
export type DrugbankResult = DrugbankXref | { id: null; checked: ('PDBe' | 'RCSB')[] }

export const pdbeCcdUrl = (ccd: string) => `https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/${encodeURIComponent(ccd)}.cif`
export const rcsbChemCompUrl = (ccd: string) => `https://data.rcsb.org/rest/v1/core/chemcomp/${encodeURIComponent(ccd)}`
export const drugbankUrl = (id: string) => `https://go.drugbank.com/drugs/${encodeURIComponent(id)}`

/** value of an mmCIF item written as a key-value pair: bare token, 'quoted', "quoted" or a ;-delimited text field */
function kvValue(cif: string, item: string): string | null {
  const lines = cif.split(/\r?\n/)
  const i = lines.findIndex((l) => l.startsWith(item + ' ') || l.trim() === item)
  if (i < 0) return null
  let rest = lines[i].slice(item.length).trim()
  if (!rest) {   // value on the next line(s)
    const next = lines[i + 1] ?? ''
    if (next.startsWith(';')) { const out = [next.slice(1)]; for (let k = i + 2; k < lines.length && !lines[k].startsWith(';'); k++) out.push(lines[k]); return out.join(' ').trim() || null }
    rest = next.trim()
  }
  const m = rest.match(/^'(.*)'$|^"(.*)"$|^(\S+)/)
  const v = m ? (m[1] ?? m[2] ?? m[3]) : null
  return v && v !== '?' && v !== '.' ? v : null
}

/** DrugBank id and name from a PDBe CCD file (key-value form, as in every gallery CCD; a looped category falls back to the first DB accession in it) */
export function parsePdbeDrugbank(cif: string): { id: string; name: string | null } | null {
  const id = kvValue(cif, '_pdbe_chem_comp_drugbank_details.drugbank_id')
  if (id && /^DB\d{5}$/.test(id)) return { id, name: kvValue(cif, '_pdbe_chem_comp_drugbank_details.name') }
  const block = cif.split(/\n(?=loop_|_[a-z])/).find((b) => b.includes('_pdbe_chem_comp_drugbank_details.drugbank_id'))
  const loopId = block?.match(/\bDB\d{5}\b/)?.[0]
  return loopId ? { id: loopId, name: null } : null
}

/** DrugBank accession from an RCSB chem comp record */
export function parseRcsbDrugbank(rec: unknown): { id: string; name: string | null } | null {
  const r = rec as { rcsb_chem_comp_related?: { resource_name?: string; resource_accession_code?: string }[]; drugbank?: { drugbank_info?: { name?: string } } } | null
  const id = r?.rcsb_chem_comp_related?.find((x) => x.resource_name === 'DrugBank' && /^DB\d{5}$/.test(x.resource_accession_code ?? ''))?.resource_accession_code
  return id ? { id, name: r?.drugbank?.drugbank_info?.name ?? null } : null
}

const memo = new Map<string, Promise<DrugbankResult>>()
/** DrugBank cross-reference for one CCD (memoised for the session). Never throws: unreachable sources count as not checked. */
export function loadDrugbank(ccd: string, fetchImpl: typeof fetch = fetch): Promise<DrugbankResult> {
  const hit = memo.get(ccd); if (hit) return hit
  const p = (async (): Promise<DrugbankResult> => {
    const checked: ('PDBe' | 'RCSB')[] = []
    try {
      const r = await fetchImpl(pdbeCcdUrl(ccd))
      if (r.ok) { checked.push('PDBe'); const x = parsePdbeDrugbank(await r.text()); if (x) return { ...x, source: 'PDBe', sourceUrl: pdbeCcdUrl(ccd) } }
    } catch { /* PDBe unreachable: RCSB below */ }
    try {
      const r = await fetchImpl(rcsbChemCompUrl(ccd))
      if (r.ok) { checked.push('RCSB'); const x = parseRcsbDrugbank(await r.json()); if (x) return { ...x, source: 'RCSB', sourceUrl: rcsbChemCompUrl(ccd) } }
    } catch { /* RCSB unreachable */ }
    return { id: null, checked }
  })()
  memo.set(ccd, p)
  p.then((v) => { if (v.id === null && v.checked.length === 0) memo.delete(ccd) })   // retry later when neither source was reachable
  return p
}

/** ligand name shortened to `max` characters (with an ellipsis) when longer */
export const NAME_MAX = 20
export const shortName = (name: string, max = NAME_MAX) => (name.length > max ? `${name.slice(0, max).trimEnd()}…` : name)
