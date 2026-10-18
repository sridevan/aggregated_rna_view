import { describe, expect, it } from 'vitest'
import { reducer, type GroupSelection, type SelectionState } from './selection'

const group = (id: string, sn = 'a_1_A_A'): GroupSelection => ({ kind: 'group', id, label: id, type: 'IL', source: { pdb_id: 'a', seq_name: sn }, members: [], segments: [] })
const s0: SelectionState = { active: null, instance: 'a_1_A_A', viewReset: 0 }

describe('selection toggling and view resets', () => {
  it('a loop clicked again clears the selection and requests the default views', () => {
    const s1 = reducer(s0, { type: 'toggleGroup', group: group('IL_A_001') })
    expect(s1.active).toMatchObject({ kind: 'group', id: 'IL_A_001' }); expect(s1.viewReset).toBe(0)
    const s2 = reducer(s1, { type: 'toggleGroup', group: group('IL_A_001') })
    expect(s2.active).toBeNull(); expect(s2.viewReset).toBe(1)
  })
  it('another loop replaces the selection without a view reset', () => {
    const s1 = reducer(reducer(s0, { type: 'toggleGroup', group: group('IL_A_001') }), { type: 'toggleGroup', group: group('HL_A_001') })
    expect(s1.active).toMatchObject({ id: 'HL_A_001' }); expect(s1.viewReset).toBe(0)
  })
  it('the same nucleotide clicked again clears it and requests the default views; another one replaces it', () => {
    const s1 = reducer(s0, { type: 'togglePosition', pos: 76 })
    const s2 = reducer(s1, { type: 'togglePosition', pos: 77 }); expect(s2.active).toEqual({ kind: 'position', pos: 77 }); expect(s2.viewReset).toBe(0)
    const s3 = reducer(s2, { type: 'togglePosition', pos: 77 }); expect(s3.active).toBeNull(); expect(s3.viewReset).toBe(1)
    const r1 = reducer(s0, { type: 'toggleResidue', seqName: 'a_1_A_A', labelSeqId: 5 }), r2 = reducer(r1, { type: 'toggleResidue', seqName: 'a_1_A_A', labelSeqId: 5 })
    expect(r2.active).toBeNull(); expect(r2.viewReset).toBe(1)
  })
  it('a nucleotide click replaces a loop (no toggle across kinds); Clear always requests the default views', () => {
    const s1 = reducer(reducer(s0, { type: 'toggleGroup', group: group('IL_A_001') }), { type: 'togglePosition', pos: 30 })
    expect(s1.active).toEqual({ kind: 'position', pos: 30 })
    expect(reducer(s1, { type: 'clear' })).toMatchObject({ active: null, viewReset: 1 })
  })
  it('an instance change clears a loop of the previous instance, keeps a consensus position', () => {
    const s1 = reducer(reducer(s0, { type: 'toggleGroup', group: group('IL_A_001') }), { type: 'instance', seqName: 'b_1_B_B' })
    expect(s1.active).toBeNull(); expect(s1.instance).toBe('b_1_B_B')
    expect(reducer(reducer(s0, { type: 'togglePosition', pos: 76 }), { type: 'instance', seqName: 'b_1_B_B' }).active).toEqual({ kind: 'position', pos: 76 })
  })
})
