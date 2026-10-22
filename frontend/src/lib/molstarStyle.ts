/*
 * Thin round backbone trace for the Structures tab 3D viewer (user decision 2026-10-22: Structures tab only; the Family Overview and the
 * Interactions comparison keep the Mol* cartoon), as in the PDBe-KB superposition viewer, instead of the flat nucleic-acid ribbon. Mol* cartoon with only the trace and gap
 * visuals, aspect ratio 1 (round) and an elliptical nucleic profile; it renders like a uniform putty but keeps the cartoon type, so the
 * emphasis and transparency logic is unchanged. Selected nucleotides and loops stay atom-level (ball-and-stick) on top of it.
 * Chosen 2026-10-22 from cartoon / trace-only cartoon / uniform putty / backbone on 4FE5 B + 1Y27 X and 16S 4V9D Y + 4U27 A.
 */
export const TRACE_SIZE = 0.35
export const TRACE_SIZE_EMPHASIS = 0.5
export const traceTypeParams = (sizeFactor = TRACE_SIZE) => ({
  sizeFactor, aspectRatio: 1, nucleicProfile: 'elliptical' as const, visuals: ['polymer-trace', 'polymer-gap'] as ('polymer-trace' | 'polymer-gap')[],
})
