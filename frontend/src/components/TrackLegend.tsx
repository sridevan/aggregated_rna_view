import type { ConsensusStats } from '../types/backend'
import type { ConsensusReference } from './ConsensusTracks'

const Sw = ({ color, dashed, hatch, bar }: { color: string; dashed?: boolean; hatch?: boolean; bar?: boolean }) => (
  <svg width="16" height="12" className="mr-1 inline-block align-middle" aria-hidden="true">
    {hatch ? <><defs><pattern id="lg-hatch" width="4" height="4" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="4" stroke={color} strokeWidth="1.2" /></pattern></defs><rect x="1" y="1" width="14" height="10" fill="url(#lg-hatch)" /></>
      : bar ? <rect x="1" y="4" width="14" height="4" fill={color} />
      : dashed ? <rect x="1.5" y="1.5" width="13" height="9" fill="none" stroke={color} strokeDasharray="2 1.5" /> : <rect x="1" y="1" width="14" height="10" fill={color} />}
  </svg>
)

const InsSw = () => (
  <svg width="18" height="12" className="mr-1 inline-block align-middle" aria-hidden="true"><rect x="1" y="1" width="16" height="9" rx="2" fill="#fef3c7" stroke="#b45309" strokeWidth="0.8" /><text x="9" y="8.3" fontSize="7.5" textAnchor="middle" fontWeight="600" fill="#92400e">+N</text></svg>
)

/** Explains bar heights, colours and patterns of the consensus tracks. */
export function TrackLegend({ stats, instanceLabel, reference }: { stats: ConsensusStats; instanceLabel?: string; reference?: ConsensusReference }) {
  const t = stats.pp_display_thresholds
  return (
    <details className="rounded border border-stone-200 bg-white text-xs text-stone-700" data-testid="track-legend-details">
      <summary className="cursor-pointer select-none px-2 py-1 font-medium text-stone-800" data-testid="legend-toggle">Legend &amp; help</summary>
      <div className="grid gap-x-6 gap-y-1.5 border-t border-stone-200 p-2 md:grid-cols-2" data-testid="track-legend">
        <div><Sw color="#2a78d6" /><b>Sequence coverage</b>: bar height = fraction of the {stats.sequence_level.n_unique_sequences} unique mapped RNA sequences that have a nucleotide (match state) at the consensus position. Hover a bar for the numerator/denominator and the number of sequences in which the position is deleted or uncovered.</div>
        <div><Sw color="#1baf7a" /><b>Structural coverage</b>: bar height = fraction of the {stats.structure_level.n_instances} RNA structural instances with experimentally observed coordinates for the nucleotide at that position. Hover for numerator/denominator.</div>
        <div><Sw color="#1baf7a" /><Sw color="#eda100" /><Sw color="#e34948" /><b>Alignment confidence</b>: bar height = mean approximate Infernal posterior probability (bin midpoints of cmalign's PP annotation) over the unique sequences aligned at the position; hover shows that number of sequences. Colour bands are provisional display classifications: green ≥ {t.high}, amber ≥ {t.medium}, red below.</div>
        <div data-testid="legend-instance"><b>{instanceLabel ?? 'Selected instance'}</b> track: <Sw color="#1baf7a" /> nucleotide mapped to the position, with coordinates · <Sw color="#1baf7a" dashed /> mapped nucleotide without experimental coordinates (in the sequence, not in the model) · <Sw color="#c3c2b7" bar /> deleted: the position has no nucleotide in this chain (inside its aligned span) · <Sw color="#c3c2b7" hatch /> uncovered: outside the part of the model aligned to this chain · <InsSw /> insertion: N nucleotides between two consensus positions (Rfam insert states; no consensus number, no extra column). Hover a marker for its sequence, click it for residue numbers. Selected consensus positions are shaded violet.</div>
        <div data-testid="legend-loops"><b>RNA loop annotations</b>: hairpin (HL), internal (IL) and junction (JL: J3, J4) loops of the selected chain from the BGSU RNA 3D Motif Atlas. Each block is one loop segment on the consensus positions; a loop with several segments has several blocks. Overlapping loops are drawn in separate lanes. The selected loop is outlined in red; clicking it again clears the selection. A loop's membership can include nucleotides without a consensus position (insertions, flanks), which are selected with it but have no column.</div>
        <div data-testid="legend-reference"><b>Reference and alignment details</b>{reference ? <> (Rfam {reference.rfamAcc}{reference.rfamId ? `, ${reference.rfamId}` : ''}, {stats.clen} consensus match positions{reference.clanMembers.length ? `; clan members ${reference.clanMembers.join(', ')}` : ''})</> : null}: <i>Consensus sequence (RF)</i> is the reference annotation of the covariance model; upper case marks a conserved consensus nucleotide. <i>Secondary structure (SS_cons)</i> is the consensus structure of the Rfam seed in WUSS notation: matching brackets are base pairs; <code>: , _ -</code> are unpaired. <i>Alignment confidence</i> is the mean Infernal posterior probability (PP) of the unique sequences aligned at the position.</div>
        <div className="md:col-span-2 text-stone-500">Axis: Rfam {stats.rfam_acc} consensus positions 1–{stats.clen}. Click a column to select that nucleotide (click it again to clear); one selection at a time, shared with the 2D and 3D viewers, which highlight and centre on it.</div>
      </div>
    </details>
  )
}
