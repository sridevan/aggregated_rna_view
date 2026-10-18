import { expect, test, type Page } from '@playwright/test'
import { colX } from './helpers'

/* Family Overview Phase 3: Rfam insertion markers (+N) between consensus columns of the selected chain. Expected values come from the
   per-chain mapping JSON (alignment_state 'insert', insertion_anchor, offsets), never from the rendering. */
type Res = { polymer_seq_position: number; author_residue_number: number; author_insertion_code: string; nucleotide_parent: string; alignment_state: string; insertion_anchor: number | null; has_coordinates: boolean }
type H = { __rnaExplorer: { viewerIds: number[]; rna2dSelected?: number[]; lastMolstarSelection?: { labelSeqIds: number[] }; focusRepresentation?: number[]; selection: { positions: number[]; active: { kind: string; id?: string; type?: string } | null; viewReset: number }
  insertions: { id: string; boundaries: number[]; length: number; runs: { id: string; left: number | null; right: number | null; sequence: string; labelSeqIds: number[] }[] }[] } }
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }) }
const hook = (page: Page) => page.evaluate(() => { const h = (window as unknown as H).__rnaExplorer; return { v: h.viewerIds, d2: h.rna2dSelected ?? [], mol: h.lastMolstarSelection?.labelSeqIds ?? [], focus: h.focusRepresentation ?? [], pos: h.selection.positions, active: h.selection.active, reset: h.selection.viewReset, ins: h.insertions } })
const inserts = async (page: Page, sn: string) => ((await (await page.request.get(`/families/RF00167/mappings/${sn}.json`)).json()).residues as Res[]).filter((r) => r.alignment_state === 'insert')
const choose = async (page: Page, sn: string) => { await page.getByTestId('instance-select').selectOption(sn); await ready(page, sn.split('_')[0]) }
/** screen x of the boundary after consensus position `pos` in a strip (shared geometry read from the strip: follows the zoom) */
const boundaryX = (page: Page, testId: string, pos: number) => colX(page, testId, pos, true)
const centreX = async (page: Page, testId: string) => { const b = (await page.getByTestId(testId).locator('rect').nth(1).boundingBox())!; return b.x + b.width / 2 }
/** geometry of the family-level coverage bars and the axis: must not depend on the selected instance */
const gridSignature = (page: Page) => page.evaluate(() => {
  // x in consensus columns ((x - LEFT) / cell): independent of the page width; heights and y in track units
  const s = document.querySelector('[data-testid="consensus-tracks"]')!, a = document.querySelector('[data-testid="consensus-axis"]')!, cell = +s.getAttribute('data-cell')!
  const col = (v: string | null) => ((+v! - 140) / cell).toFixed(3)
  const bars = [...s.querySelectorAll(':scope > g:not([data-testid]) > rect')].map((r) => [col(r.getAttribute('x')), r.getAttribute('y'), r.getAttribute('height')].join(','))
  return { columns: ((+s.getAttribute('viewBox')!.split(' ')[2] - 152) / cell).toFixed(3), axis: [...a.querySelectorAll('text')].map((t) => col(t.getAttribute('x'))).join(','), bars: bars.join(';'), nBars: bars.length, viewBox: s.getAttribute('viewBox') }
})

test('9IWF A: one +2 marker between 77 and 78 from the mapping; none for chains without insert states (4FE5 B, 1Y26 X flanks)', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  await expect(page.locator('[data-testid^="insertion-marker-"]')).toHaveCount(0)
  expect((await hook(page)).ins).toEqual([])
  const g0 = await gridSignature(page)
  const cell = +(await page.getByTestId('consensus-tracks').getAttribute('data-cell'))!
  expect(+g0.viewBox!.split(' ')[2]).toBeCloseTo(140 + 102 * cell + 12, 6); expect(g0.nBars).toBeGreaterThan(150)                     // 102 columns, unchanged
  await choose(page, '1y26_1_X_A')                                                                                                    // has terminal flanks only
  await expect(page.locator('[data-testid^="insertion-marker-"]')).toHaveCount(0)
  await choose(page, '9iwf_1_A_A')
  const exp = await inserts(page, '9iwf_1_A_A')
  expect(exp.map((r) => r.insertion_anchor)).toEqual([77, 77])
  const mk = page.getByTestId('insertion-marker-77')
  await expect(mk).toHaveAttribute('data-length', String(exp.length)); await expect(mk).toHaveText('+2')
  await expect(page.locator('[data-testid^="insertion-marker-"]')).toHaveCount(1)
  const h = await hook(page)
  expect(h.ins).toEqual([{ id: '77', boundaries: [77], length: 2, runs: [{ id: '77:59', kind: 'insert', left: 77, right: 78, sequence: exp.map((r) => r.nucleotide_parent).join(''), labelSeqIds: exp.map((r) => r.polymer_seq_position) }] }])
  // centred on the boundary between columns 77 and 78 of the shared grid; the grid, coverage bars and axis are unchanged
  expect(Math.abs((await centreX(page, 'insertion-marker-77')) - (await boundaryX(page, 'consensus-tracks', 77)))).toBeLessThan(1)
  { const g1 = await gridSignature(page); expect({ ...g1, viewBox: '' }).toEqual({ ...g0, viewBox: '' }); expect(g1.columns).toBe('102.000') }
  // accessible: a named button
  await expect(mk).toHaveAttribute('role', 'button'); await expect(mk).toHaveAttribute('aria-label', /Insertion, 2 nucleotides, between rfam positions 77 and 78/)
})

test('tooltip and click-open details with deposited identifiers; never a consensus selection; click again closes; Escape closes', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5'); await choose(page, '9iwf_1_A_A')
  const exp = await inserts(page, '9iwf_1_A_A')
  await page.getByTestId('insertion-marker-77').hover()
  await expect(page.getByTestId('insertion-tooltip-title')).toHaveText('Insertion · +2 nt')
  await expect(page.getByTestId('insertion-tooltip-boundary')).toHaveText('Between Rfam positions 77 and 78')
  await expect(page.getByTestId('insertion-tooltip-sequence')).toHaveText(`Sequence: ${exp.map((r) => r.nucleotide_parent).join('')}`)
  const tb = (await page.getByTestId('insertion-tooltip').boundingBox())!, vw = page.viewportSize()!.width
  expect(tb.x).toBeGreaterThanOrEqual(0); expect(tb.x + tb.width).toBeLessThanOrEqual(vw)
  // click: details (residue identifiers, coordinates, provenance); no consensus position, no viewer change
  await page.getByTestId('insertion-marker-77').click()
  const det = page.getByTestId('insertion-details'); await expect(det).toBeVisible(); await expect(page.getByTestId('insertion-tooltip')).toHaveCount(0)
  await expect(page.getByTestId('insertion-marker-77')).toHaveAttribute('data-open', 'true')
  const rows = await page.getByTestId('insertion-run-residues').locator('tbody tr').allInnerTexts()
  expect(rows.map((r) => r.split('\t').map((c) => c.trim()))).toEqual(exp.map((r) => [r.nucleotide_parent, String(r.polymer_seq_position), `${r.author_residue_number}${r.author_insertion_code}`, r.has_coordinates ? 'yes' : 'no']))
  await expect(page.getByTestId('insertion-run-boundary')).toHaveText('Between Rfam positions 77 and 78')
  await expect(page.getByTestId('insertion-details-provenance')).toContainText('RF00167'); await expect(page.getByTestId('insertion-details-provenance')).toContainText('cmalign')
  let h = await hook(page); expect(h.active).toBeNull(); expect(h.pos).toEqual([]); expect(h.v).toEqual([])
  // the same marker again closes; Escape and the close button too
  await page.getByTestId('insertion-marker-77').click(); await expect(det).toHaveCount(0)
  await page.getByTestId('insertion-marker-77').click(); await expect(det).toBeVisible(); await page.keyboard.press('Escape'); await expect(det).toHaveCount(0)
  await page.getByTestId('insertion-marker-77').focus(); await page.keyboard.press('Enter'); await expect(det).toBeVisible()                // keyboard
  await page.getByTestId('insertion-details-close').click(); await expect(det).toHaveCount(0)
  h = await hook(page); expect(h.active).toBeNull()
})

test('show the inserted nucleotides in 2D and 3D from the details; the marker clicked again closes and clears with the default views', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5'); await choose(page, '9iwf_1_A_A')
  const ids = (await inserts(page, '9iwf_1_A_A')).map((r) => r.polymer_seq_position)
  await page.getByTestId('insertion-marker-77').click(); await page.getByTestId('insertion-view-77:59').click()
  await expect.poll(async () => (await hook(page)).v).toEqual(ids)
  await expect.poll(async () => (await hook(page)).d2).toEqual(ids); await expect.poll(async () => (await hook(page)).focus).toEqual(ids)
  const h = await hook(page); expect(h.pos).toEqual([]); expect(h.active?.type).toBe('INS')                                            // no consensus position involved
  await expect(page.getByTestId('selection-status')).toContainText('Insertion +2 · between rfam positions 77 and 78 · C59 A60')
  await expect(page.getByTestId('insertion-marker-77')).toHaveAttribute('data-selected', 'true')
  await expect(page.getByTestId('insertion-view-77:59')).toHaveText('Clear from 2D and 3D')
  const r0 = h.reset
  await page.getByTestId('insertion-marker-77').click()
  await expect(page.getByTestId('insertion-details')).toHaveCount(0)
  await expect.poll(async () => (await hook(page)).active).toBeNull(); await expect.poll(async () => (await hook(page)).focus).toEqual([])
  expect((await hook(page)).reset).toBe(r0 + 1)                                                                                         // default 2D/3D views restored
  // a consensus column right next to the marker selects that position only, never the insertion
  const b = (await page.getByTestId('consensus-tracks').boundingBox())!, x78 = await colX(page, 'consensus-tracks', 78)
  await page.mouse.click(x78, b.y + 12)
  await expect.poll(async () => (await hook(page)).active).toEqual({ kind: 'position', pos: 78 })
  await page.mouse.click(x78, b.y + 12)
  await expect.poll(async () => (await hook(page)).active).toBeNull()
})

test('instance switching replaces the markers and clears details, tooltip and insertion selection; loop tracks unchanged and selectable', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5'); await choose(page, '9iwf_1_A_A')
  await page.getByTestId('loops-toggle').click()
  // block geometry in consensus columns (x and width divided by the shared cell width), so a page-width change cannot masquerade as a move
  const loopGeom = () => page.locator('[data-testid^="loop-block-"]').evaluateAll((es) => { const cell = +(document.querySelector('[data-testid="loop-tracks"]')!.getAttribute('data-cell'))!
    return es.map((e) => [e.getAttribute('data-testid'), ((+e.getAttribute('x')! - 140) / cell).toFixed(2), (+e.getAttribute('width')! / cell).toFixed(2)].join(':')) })
  const before = await loopGeom()
  await page.getByTestId('insertion-marker-77').click(); await page.getByTestId('insertion-view-77:59').click()
  await expect.poll(async () => (await hook(page)).active?.type).toBe('INS')
  // the internal loop containing the inserted nucleotides keeps its full membership (59 and 60 included)
  await page.locator('[data-testid^="loop-block-IL_9IWF_001-"]').first().click()
  await expect.poll(async () => (await hook(page)).active?.id).toBe('IL_9IWF_001')
  expect((await hook(page)).v).toEqual(expect.arrayContaining([59, 60]))
  expect(await loopGeom()).toEqual(before)
  // open details + hovered marker, then a chain without insertions: no stale marker, tooltip, details or selection
  await page.getByTestId('insertion-marker-77').click(); await expect(page.getByTestId('insertion-details')).toBeVisible()
  await page.getByTestId('insertion-view-77:59').click(); await expect.poll(async () => (await hook(page)).active?.type).toBe('INS')
  await page.getByTestId('insertion-marker-77').hover(); await page.waitForTimeout(150)
  await choose(page, '4fe5_1_B_A')
  await expect(page.locator('[data-testid^="insertion-marker-"]')).toHaveCount(0)
  await expect(page.getByTestId('insertion-details')).toHaveCount(0); await expect(page.getByTestId('insertion-tooltip')).toHaveCount(0)
  expect((await hook(page)).active).toBeNull(); await expect.poll(async () => (await hook(page)).mol).toEqual([])
  await expect(page.getByTestId('loops-section')).toHaveAttribute('data-open', 'true')
  // 9LJN A: same boundary, its own residues; details from 9IWF never carried over
  await page.mouse.move(5, 5)
  await choose(page, '9ljn_1_A_A')
  const exp = await inserts(page, '9ljn_1_A_A')
  await expect.poll(async () => (await hook(page)).ins[0]?.runs[0].labelSeqIds).toEqual(exp.map((r) => r.polymer_seq_position))
  await page.getByTestId('insertion-marker-77').click()
  expect((await page.getByTestId('insertion-run-residues').locator('tbody tr').allInnerTexts()).map((r) => r.split('\t')[1].trim())).toEqual(exp.map((r) => String(r.polymer_seq_position)))
})

test('responsive: the marker stays on its boundary at narrow widths; legend and track note explain the coordinate system', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5'); await choose(page, '9iwf_1_A_A')
  for (const w of [1400, 900, 640]) {
    await page.setViewportSize({ width: w, height: 1000 }); await page.waitForTimeout(300)
    // at overview scales the +N markers are hidden to avoid clutter: zoom in until they show, then check the anchoring
    while ((await page.getByTestId('instance-track').getAttribute('data-markers-shown')) !== 'true') await page.getByTestId('zoom-in').click()
    await page.getByTestId('insertion-marker-77').scrollIntoViewIfNeeded()
    expect(Math.abs((await centreX(page, 'insertion-marker-77')) - (await boundaryX(page, 'consensus-tracks', 77))), `width ${w}`).toBeLessThan(1)
    expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)
  }
  await page.setViewportSize({ width: 1400, height: 1000 })
  await page.getByTestId('instance-track-info').hover()
  await expect(page.getByTestId('instance-track-info-text')).toHaveText('Shows nucleotide availability at Rfam match positions. +N markers indicate insertions between positions.')
  await page.getByTestId('legend-toggle').click()
  await expect(page.getByTestId('legend-instance')).toContainText('insertion: N nucleotides between two consensus positions')
  await expect(page.getByTestId('legend-instance')).toContainText('deleted'); await expect(page.getByTestId('legend-instance')).toContainText('without experimental coordinates'); await expect(page.getByTestId('legend-instance')).toContainText('uncovered')
})
