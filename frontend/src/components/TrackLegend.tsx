import type { ConsensusStats } from '../types/backend'

const Sw = ({ color, dashed, hatch, bar }: { color: string; dashed?: boolean; hatch?: boolean; bar?: boolean }) => (
  <svg width="16" height="12" className="mr-1 inline-block align-middle" aria-hidden="true">
    {hatch ? <><defs><pattern id="lg-hatch" width="4" height="4" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="4" stroke={color} strokeWidth="1.2" /></pattern></defs><rect x="1" y="1" width="14" height="10" fill="url(#lg-hatch)" /></>
      : bar ? <rect x="1" y="4" width="14" height="4" fill={color} />
      : dashed ? <rect x="1.5" y="1.5" width="13" height="9" fill="none" stroke={color} strokeDasharray="2 1.5" /> : <rect x="1" y="1" width="14" height="10" fill={color} />}
  </svg>
)

/** Explains bar heights, colours and patterns of the consensus tracks. */
export function TrackLegend({ stats, instanceLabel }: { stats: ConsensusStats; instanceLabel?: string }) {
  const t = stats.pp_display_thresholds
  return (
    <div className="grid gap-x-6 gap-y-2 rounded border border-stone-200 bg-white p-2 text-xs text-stone-700 md:grid-cols-2" data-testid="track-legend">
      <div><Sw color="#2a78d6" /><b>Sequence coverage</b>: bar height = fraction of the {stats.sequence_level.n_unique_sequences} unique mapped RNA sequences that have a nucleotide (match state) at the consensus position. Hover a bar for the numerator/denominator and the number of sequences in which the position is deleted or uncovered.</div>
      <div><Sw color="#1baf7a" /><b>Structural coverage</b>: bar height = fraction of the {stats.structure_level.n_instances} RNA structural instances with experimentally observed coordinates for the nucleotide at that position. Hover for numerator/denominator.</div>
      <div><Sw color="#1baf7a" /><Sw color="#eda100" /><Sw color="#e34948" /><b>Alignment confidence</b>: bar height = mean approximate Infernal posterior probability (bin midpoints of cmalign's PP annotation) over the unique sequences aligned at the position; hover shows that number of sequences. Colour bands are provisional display classifications: green ≥ {t.high}, amber ≥ {t.medium}, red below.</div>
      <div><Sw color="#1baf7a" /> matched, coordinates present · <Sw color="#1baf7a" dashed /> matched, no experimental coordinates · <Sw color="#c3c2b7" bar /> deleted in this chain (inside its aligned span) · <Sw color="#c3c2b7" hatch /> uncovered (outside the aligned span) · ▲ insertion / local-end residues after a position (blue = insert, red = local end). <b>{instanceLabel ?? 'Selected instance'}</b> track; selected consensus positions are shaded violet.</div>
    </div>
  )
}
