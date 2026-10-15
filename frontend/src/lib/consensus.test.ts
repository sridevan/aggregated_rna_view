import { describe, expect, it } from 'vitest'
import { readFileSync } from 'node:fs'
import { insertionsByAnchor, positionStates, residueToConsensus, translateConsensus } from './consensus'
import type { ChainMapping } from '../types/backend'

const load = (name: string) => JSON.parse(readFileSync(new URL(`../../../data/families/RF00167/mappings/${name}.json`, import.meta.url), 'utf8')) as ChainMapping

describe('translateConsensus', () => {
  it('maps consensus 50 to polymer 64 in 9ndd/A and to polymer 34 in 1y26/X', () => {
    const t = translateConsensus(load('9ndd_1_A_A'), [50])
    expect(t.highlightable).toHaveLength(1)
    expect(t.highlightable[0]).toMatchObject({ label_asym_id: 'A', label_seq_id: 64, auth_asym_id: 'A', author_residue_number: 64, nucleotide_parent: 'G' })
    const u = translateConsensus(load('1y26_1_X_A'), [50])
    expect(u.highlightable[0]).toMatchObject({ label_asym_id: 'A', label_seq_id: 34, auth_asym_id: 'X', author_residue_number: 46, nucleotide_parent: 'G' })
  })
  it('reports deletion, uncovered and unobserved positions without highlighting them', () => {
    const t = translateConsensus(load('9nbh_1_A_A'), [1, 16, 39, 50, 64, 102])
    expect(t.deleted).toEqual([39, 64])
    expect(t.highlightable.map((r) => r.rfam_consensus_pos)).toEqual([1, 102])
    expect(t.unobserved.map((r) => r.rfam_consensus_pos)).toEqual([16, 50])
    const u = translateConsensus(load('6p2h_1_A_A'), [1, 25, 26, 80, 81])
    expect(u.uncovered).toEqual([1, 25, 81])
    expect(u.highlightable.map((r) => r.label_seq_id)).toEqual([6, 60])
  })
  it('inverse mapping and instance-level residues', () => {
    const m = load('1y26_1_X_A')
    expect(residueToConsensus(m, 34)).toMatchObject({ rfam_consensus_pos: 50, alignment_state: 'match' })
    expect(residueToConsensus(m, 1)).toMatchObject({ rfam_consensus_pos: null, alignment_state: 'not_submitted_flank' })
    expect(residueToConsensus(m, 999)).toBeNull()
    const s = positionStates(load('9nbh_1_A_A'), 102)
    expect(s[0]).toBe('match_observed'); expect(s[15]).toBe('match_unobserved'); expect(s[38]).toBe('deletion')
    expect(positionStates(load('6p2h_1_A_A'), 102)[0]).toBe('uncovered')
    expect(insertionsByAnchor(load('6p2h_1_A_A')).size).toBe(0)
  })
})
