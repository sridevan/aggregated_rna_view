/** The Rfam family shown is taken from the URL path: /{rfam_acc} (e.g. /RF00167). */
export const DEFAULT_ACC = 'RF00167'
/** Families with generated data in this build (data/families/<ACC>). */
export const AVAILABLE_ACCS = ['RF00167', 'RF00162', 'RF00059', 'RF00177']
/** Short names for the family selector (the page heading uses the Rfam description from family.json). */
export const FAMILY_NAMES: Record<string, string> = { RF00167: 'Purine riboswitch', RF00162: 'SAM-I riboswitch', RF00059: 'TPP riboswitch', RF00177: 'Bacterial 16S rRNA' }

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
  RF00167: 'Bacterial mRNA element that binds guanine or adenine and regulates genes of purine metabolism and transport.',
  RF00162: 'Bacterial mRNA element that binds S-adenosylmethionine (SAM) and regulates genes of methionine and sulfur metabolism.',
  RF00059: 'mRNA element that binds thiamine pyrophosphate (TPP); the only riboswitch class known in bacteria, archaea and eukaryotes.',
  RF00177: 'The RNA of the bacterial small ribosomal subunit; this page uses a curated set of chains rather than every deposited ribosome.',
}
