import { expect, test, type Page } from '@playwright/test'
import { colX, selectColumn, stripGeom } from './helpers'

/* Family Overview Phase 2: per-instance BGSU loop tracks (HL / IL / JL) with whole-loop selection in 2D and Mol*. Expected values come from
   loops.json (original BGSU membership and strands), not from the rendering. */
type Nt = { seq_name: string; label_seq_id: number; alignment_state: string; rfam_consensus_pos: number | null; has_coordinates: boolean }
type Loop = { loop_id: string; loop_type: string; strands: { nucleotide_indices: [number, number] }[]; nucleotides: Nt[] }
type H = { __rnaExplorer: { viewerIds: number[]; rna2dSelected?: number[]; lastMolstarSelection?: { labelSeqIds: number[] }; focusRepresentation?: number[]; lastFocus?: number[]; selection: { positions: number[]; active: { kind: string; id?: string } | null } } }
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }) }
const hook = (page: Page) => page.evaluate(() => { const h = (window as unknown as H).__rnaExplorer; return { v: h.viewerIds, d2: h.rna2dSelected, mol: h.lastMolstarSelection?.labelSeqIds, focus: h.focusRepresentation, cam: h.lastFocus, pos: h.selection.positions, active: h.selection.active } })
const sorted = (a: number[]) => [...a].sort((x, y) => x - y)
let LOOPS: Loop[] = []
const loopOf = (id: string) => LOOPS.find((l) => l.loop_id === id)!
const members = (id: string, sn: string, coords = true) => sorted(loopOf(id).nucleotides.filter((n) => n.seq_name === sn && (!coords || n.has_coordinates)).map((n) => n.label_seq_id))
const open = async (page: Page) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  LOOPS = (await (await page.request.get('/families/RF00167/loops.json')).json()).loops
}

test('default collapsed; tracks of the selected instance with blocks per original segment, distinct IDs and lanes for overlaps', async ({ page }) => {
  await open(page)
  await expect(page.getByTestId('loops-section')).toHaveAttribute('data-open', 'false'); await expect(page.getByTestId('loops-body')).toHaveCount(0)
  await expect(page.getByTestId('reference-section')).toHaveAttribute('data-open', 'false')
  const y = async (id: string) => (await page.getByTestId(id).boundingBox())!.y
  expect(await y('consensus-tracks')).toBeLessThan(await y('loops-section')); expect(await y('loops-section')).toBeLessThan(await y('reference-section')); expect(await y('reference-section')).toBeLessThan(await y('consensus-axis'))
  await page.getByTestId('loops-toggle').click()
  await expect(page.getByTestId('loop-tracks')).toBeVisible()
  for (const t of ['HL', 'IL', 'JL']) await expect(page.getByTestId(`loop-track-${t}`)).toBeVisible()
  // every loop of 4FE5 B, with one block per original strand spanning that strand's match positions
  const sn = '4fe5_1_B_A', mine = LOOPS.filter((l) => l.nucleotides.some((n) => n.seq_name === sn))
  const blocks = await page.locator('[data-testid^="loop-block-"]').evaluateAll((es) => es.map((e) => ({ id: e.getAttribute('data-loop-id'), from: +e.getAttribute('data-from')!, to: +e.getAttribute('data-to')!, lane: +e.getAttribute('data-lane')!, fill: e.getAttribute('fill') })))
  for (const l of mine) {
    const exp = l.strands.map((s) => l.nucleotides.slice(s.nucleotide_indices[0], s.nucleotide_indices[1] + 1).filter((n) => n.seq_name === sn && n.alignment_state === 'match').map((n) => n.rfam_consensus_pos!)).filter((p) => p.length).map((p) => [Math.min(...p), Math.max(...p)])
    expect(blocks.filter((b) => b.id === l.loop_id).map((b) => [b.from, b.to])).toEqual(exp)
  }
  expect(new Set(blocks.map((b) => b.id))).toEqual(new Set(mine.map((l) => l.loop_id)))                            // no loop dropped
  expect(blocks.filter((b) => b.id === 'IL_4FE5_001').length).toBe(2)                                               // discontinuous internal loop
  expect(blocks.filter((b) => b.id === 'J3_4FE5_001').length).toBe(3); expect(blocks.filter((b) => b.id === 'J4_4FE5_001').length).toBe(4)
  // overlap: J3 and J4 share nucleotides -> separate lanes, both selectable; colours per category only
  expect(new Set(blocks.filter((b) => b.id?.startsWith('J')).map((b) => b.lane))).toEqual(new Set([0, 1]))
  await expect(page.getByTestId('loop-track-JL')).toHaveAttribute('data-lanes', '2')
  expect(new Set(blocks.filter((b) => b.id?.startsWith('HL')).map((b) => b.fill)).size).toBe(1)
  // all blocks on the shared grid: block x of consensus 29 equals the column centre of the primary tracks
  const g = await stripGeom(page, 'loop-tracks'), half = g.cell / 2
  const il = (await page.getByTestId('loop-block-IL_4FE5_001-0').boundingBox())!
  expect(Math.abs(il.x - ((await colX(page, 'consensus-tracks', 29)) - (half - 0.75) * g.scale))).toBeLessThan(1)                // left edge on the column of consensus 29
  expect(Math.abs(il.x + il.width - ((await colX(page, 'consensus-tracks', 31)) + (half - 0.75) * g.scale))).toBeLessThan(1)     // right edge on consensus 31
  // hover emphasises all segments of that loop and fades the others; the tooltip names the loop
  await page.getByTestId('loop-block-IL_4FE5_001-0').hover()
  await expect(page.getByTestId('loop-block-IL_4FE5_001-1')).toHaveAttribute('opacity', '1'); await expect(page.getByTestId('loop-block-HL_4FE5_001-0')).toHaveAttribute('opacity', '0.45')
  await expect(page.getByTestId('loop-hover-IL_4FE5_001-1')).toHaveCount(1)
  await expect(page.getByTestId('loop-tooltip-id')).toHaveText('IL_4FE5_001'); await expect(page.getByTestId('loop-tooltip-counts')).toHaveText('6 nucleotides · 2 segments')
  expect((await hook(page)).active).toBeNull()                                                                       // hover never selects
})

test('whole-loop selection: exact membership in 2D and 3D, combined focus, 2D zoom kept, replacement, clear, single nucleotides intact', async ({ page }) => {
  await open(page)
  await page.getByTestId('loops-toggle').click()
  // user zooms the 2D view: the selection must not change it
  const rna = page.getByTestId('rna2d-viewer'); await expect(rna).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await page.locator('[data-testid="rna2d-viewer"] [id^="rnaTopologyZoomIn-"]').first().click(); await page.waitForTimeout(500)   // the viewer's own zoom (wheel zoom is disabled by the viewer)
  const view2d = () => page.evaluate(() => { const s = document.querySelector('[data-testid="rna2d-viewer"] svg.rnaTopoSvg')!; return [s.getAttribute('viewBox'), ...[...s.querySelectorAll('g[transform]')].map((g) => g.getAttribute('transform'))].join('|') })
  const before = await view2d()
  expect(before).toMatch(/scale\(1\.2/)                                                                                // the wheel zoom actually changed the 2D view
  // internal loop selected from its second segment
  await page.getByTestId('loop-block-IL_4FE5_001-1').click()
  const sn = '4fe5_1_B_A', il = members('IL_4FE5_001', sn)
  await expect.poll(async () => sorted((await hook(page)).v)).toEqual(il)
  await expect.poll(async () => sorted((await hook(page)).d2 ?? [])).toEqual(il)
  await expect.poll(async () => sorted((await hook(page)).mol ?? [])).toEqual(il)
  await expect.poll(async () => sorted((await hook(page)).focus ?? [])).toEqual(il)
  await expect.poll(async () => sorted((await hook(page)).cam ?? [])).toEqual(il)                                       // one camera focus over both segments
  const h = await hook(page)
  expect(h.v.some((id) => id > 13 && id < 29)).toBe(false)                                                         // intervening stem excluded
  expect(sorted(h.pos)).toEqual([29, 30, 31, 47, 48, 49])
  await expect(page.getByTestId('selection-status')).toHaveText('Internal loop · IL_4FE5_001 · 2 segments · 6 nucleotides')
  await expect(page.getByTestId('selection-status')).toHaveAttribute('title', /4FE5\|1\|B\|U\|25.*4FE5\|1\|B\|A\|45/)
  await expect(page.getByTestId('loop-outline-IL_4FE5_001-0')).toHaveAttribute('stroke', '#d1202f')
  await expect(page.getByTestId('sel-highlightable-details')).not.toHaveAttribute('open', '')                       // member list folded for a loop
  await expect(page.getByTestId('sel-highlightable-details').locator('summary')).toHaveText('Highlighted (6)')
  await expect(rna).toContainText('Selected Residues 11 - 13, 29 - 31')
  expect(await view2d()).toBe(before)                                                                               // 2D zoom and pan unchanged
  // hairpin replaces it; junction (3 segments) replaces that
  await page.getByTestId('loop-block-HL_4FE5_002-0').click()
  await expect.poll(async () => sorted((await hook(page)).v)).toEqual(members('HL_4FE5_002', sn))
  await expect(page.getByTestId('selection-status')).toContainText('Hairpin loop · HL_4FE5_002 · 1 segment · 8 nucleotides')
  await page.getByTestId('loop-block-J3_4FE5_001-2').click()
  await expect.poll(async () => sorted((await hook(page)).mol ?? [])).toEqual(members('J3_4FE5_001', sn))
  await expect(page.getByTestId('selection-status')).toContainText('Three-way junction loop · J3_4FE5_001 · 3 segments')
  // the overlapping J4 in the other lane is its own loop
  await page.getByTestId('loop-block-J4_4FE5_001-0').click()
  await expect.poll(async () => sorted((await hook(page)).v)).toEqual(members('J4_4FE5_001', sn))
  // collapsing the section keeps the selection; clear removes it everywhere
  await page.getByTestId('loops-toggle').click(); await expect(page.getByTestId('loops-body')).toHaveCount(0)
  expect((await hook(page)).active?.id).toBe('J4_4FE5_001')
  await page.getByTestId('clear-selection').click()
  await expect.poll(async () => (await hook(page)).v).toEqual([]); await expect.poll(async () => (await hook(page)).focus ?? []).toEqual([])
  // ordinary single-nucleotide selection is unchanged
  await selectColumn(page, 76)
  await expect(page.getByTestId('sel-highlightable')).toContainText('76 → C56'); await expect.poll(async () => (await hook(page)).focus).toEqual([56])
})

test('instance switching clears the loop and shows the new chain\'s loops; inserts, off-grid and inter-chain loops keep their full membership', async ({ page }) => {
  await open(page)
  await page.getByTestId('loops-toggle').click(); await page.getByTestId('loop-block-HL_4FE5_001-0').click()
  await expect.poll(async () => (await hook(page)).active?.id).toBe('HL_4FE5_001')
  // 9IWF A: the internal loop includes two inserted nucleotides (no consensus position): selected, sent to the viewers, no extra column
  await page.getByTestId('instance-select').selectOption('9iwf_1_A_A'); await ready(page, '9iwf')
  expect((await hook(page)).active).toBeNull()                                                                       // incompatible loop cleared
  await expect(page.getByTestId('loops-section')).toHaveAttribute('data-open', 'true')                              // section state kept
  await expect(page.getByTestId('loop-block-HL_4FE5_001-0')).toHaveCount(0)
  const ilIns = LOOPS.find((l) => l.loop_id === 'IL_9IWF_001')!
  expect(ilIns.nucleotides.filter((n) => n.alignment_state === 'insert').map((n) => n.label_seq_id)).toEqual([59, 60])
  await page.locator('[data-testid^="loop-block-IL_9IWF_001-"]').first().click()
  await expect.poll(async () => sorted((await hook(page)).v)).toEqual(members('IL_9IWF_001', '9iwf_1_A_A'))
  expect((await hook(page)).v).toEqual(expect.arrayContaining([59, 60]))
  // 9NBH A: loops without any consensus position are listed and selectable; an inter-chain loop reports its other-chain members
  await page.getByTestId('instance-select').selectOption('9nbh_1_A_A'); await ready(page, '9nbh')
  await expect(page.getByTestId('loops-off-grid')).toContainText('IL_9NBH_001')
  await page.getByTestId('loops-off-grid').getByRole('button', { name: 'IL_9NBH_003' }).click()
  await expect(page.getByTestId('selection-status')).toContainText('IL_9NBH_003'); await expect(page.getByTestId('selection-status')).toContainText('4 in another chain (not shown)')
  await expect.poll(async () => sorted((await hook(page)).v)).toEqual(members('IL_9NBH_003', '9nbh_1_A_A'))
})
