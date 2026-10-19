/** The Rfam family shown is taken from the URL path: /{rfam_acc} (e.g. /RF00167). */
export const DEFAULT_ACC = 'RF00167'
/** Families with generated data in this build (data/families/<ACC>). */
export const AVAILABLE_ACCS = ['RF00167', 'RF00162', 'RF00059', 'RF00005', 'RF00177', 'RF00011']
/** Short names for the family selector (the page heading uses the Rfam description from family.json). */
export const FAMILY_NAMES: Record<string, string> = { RF00167: 'Purine riboswitch', RF00162: 'SAM-I riboswitch', RF00059: 'TPP riboswitch', RF00005: 'tRNA', RF00177: 'Bacterial 16S rRNA', RF00011: 'RNase P RNA (class B)' }

/** Rfam accession from a URL path ("/RF00167", "/rf00167/" -> "RF00167"); null when the path is not a single accession segment. */
export function accFromPath(path: string): string | null {
  const m = path.match(/^\/(RF\d{5})\/?$/i)
  return m ? m[1].toUpperCase() : null
}

/** Canonical path for a location: "/" (home page) for the root and any path that is not a single accession; otherwise the accession in upper case. */
export function canonicalPath(path: string): string {
  const acc = accFromPath(path)
  return acc ? `/${acc}` : '/'
}

/** One-sentence description of each family for the home page (biology only; numbers on the page come from the generated family data). */
export const FAMILY_BLURBS: Record<string, string> = {
  RF00167: 'Bacterial riboswitch that binds adenine or guanine to regulate genes involved in purine metabolism and transport.',
  RF00162: 'Bacterial riboswitch that recognises S-adenosylmethionine (SAM) and regulates genes involved in methionine and sulfur metabolism.',
  RF00059: 'Riboswitch that binds thiamine pyrophosphate (TPP) and regulates genes involved in thiamine metabolism and transport.',
  RF00005: 'Transfer RNA delivers amino acids to the ribosome during protein synthesis and participates in diverse RNA–protein complexes.',
  RF00177: 'RNA component of the bacterial small ribosomal subunit, involved in decoding and translation.',
  RF00011: 'Catalytic RNA component of bacterial RNase P, which processes the 5′ ends of precursor tRNAs.',
}
