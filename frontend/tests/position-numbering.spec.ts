import { expect, test, type Page } from '@playwright/test'
import { selectColumn } from './helpers'

/* Position numbering (display only): Rfam consensus by default; any chain's author numbering on request, independent of the active instance.
   Every expected value is read from the real per-chain mapping files (data/families/<acc>/mappings), never from the rendering. */
type Res = { polymer_seq_position: number; author_residue_number: number | string | null; author_insertion_code: string; residue_name: string; has_coordinates: boolean; alignment_state: string; rfam_consensus_pos: number | null }
type Mapping = { residues: Res[]; consensus_columns: { rfam_consensus_pos: number; state: string; polymer_seq_position: number | null }[] }
const J = async (page: Page, acc: string, sn: string): Promise<Mapping> => (await page.request.get(`/families/${acc}/mappings/${sn}.json`)).json()
/** the value the UI must show for a chain at a position, computed from the mapping file */
const expected = (m: Mapping, pos: number) => {
  const c = m.consensus_columns[pos - 1]
  if (c.state === 'deletion') return { value: 'Not present', status: 'Deletion' }
  if (c.state === 'uncovered') return { value: 'Unmapped', status: 'Outside aligned region' }
  const r = m.residues.find((x) => x.polymer_seq_position === c.polymer_seq_position)!
  return { value: r.author_residue_number == null ? 'Unavailable' : `${r.residue_name}${r.author_residue_number}${r.author_insertion_code}`, status: r.has_coordinates ? 'Coordinates available' : 'Coordinates unavailable' }
}
const authorNo = (m: Mapping, pos: number) => { const c = m.consensus_columns[pos - 1]; if (c.state !== 'match') return null; const r = m.residues.find((x) => x.polymer_seq_position === c.polymer_seq_position)!; return r.author_residue_number == null ? null : `${r.author_residue_number}${r.author_insertion_code}` }
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 }) }
const chainCell = (page: Page, pos: number) => page.getByTestId('consensus-tracks').locator(`[data-row="chain"][data-tip-pos="${pos}"]`)
const row = (page: Page, key: string) => page.getByTestId(`track-tooltip-value-${key}`)
const rowLabel = async (page: Page, key: string) => (await page.getByTestId(`track-tooltip-row-${key}`).locator('dt').textContent())
const ticks = async (page: Page) => page.getByTestId('consensus-axis').locator('[data-ref-tick]').evaluateAll((gs) => gs.map((g) => ({ pos: +g.getAttribute('data-ref-tick')!, text: g.textContent })))

test('RF00167: default Rfam numbering; a reference chain relabels the axis and tooltips only; active instance and reference are independent', async ({ page }) => {
  const urls: string[] = []; page.on('request', (r) => urls.push(r.url()))
  await page.goto('/RF00167'); await ready(page, '4fe5')
  const fe5 = await J(page, 'RF00167', '4fe5_1_B_A'), y26 = await J(page, 'RF00167', '1y26_1_X_A'), xnw = await J(page, 'RF00167', '2xnw_1_A_A')
  // 1. default: Rfam consensus numbering
  await expect(page.getByTestId('numbering-select')).toHaveValue('')
  await expect(page.getByTestId('axis-convention')).toHaveText('Rfam RF00167 position')
  await chainCell(page, 50).hover()
  await expect(page.getByTestId('track-tooltip-title')).toHaveText('Selected RNA chain · 4FE5 B')
  await expect(row(page, 'rfam')).toHaveText('50'); expect(await rowLabel(page, 'author')).toBe('Author residue')
  await expect(row(page, 'author')).toHaveText(expected(fe5, 50).value); await expect(row(page, 'status')).toHaveText('Coordinates available')
  await expect(row(page, 'reference')).toHaveCount(0)

  // zoom in and scroll a little: numbering changes must keep both
  await page.getByTestId('zoom-in').click(); await page.getByTestId('zoom-in').click()
  await page.getByTestId('tracks-scroll').evaluate((e) => { e.scrollLeft = 120 })
  const zoom = await page.getByTestId('zoom-controls').getAttribute('data-cell-px'), scroll = await page.getByTestId('tracks-scroll').evaluate((e) => e.scrollLeft)
  const xBefore = (await chainCell(page, 50).boundingBox())!.x

  // 2. reference 1Y26 X: axis labels are 1Y26 author numbers at the right Rfam columns; columns do not move; viewers untouched
  const nReq = urls.length
  await page.getByTestId('numbering-select').selectOption('1y26_1_X_A')
  await expect(page.getByTestId('numbering-control')).toHaveAttribute('data-state', 'ready')
  await expect(page.getByTestId('axis-convention')).toHaveText('1Y26 X · author numbering')
  const t = await ticks(page); expect(t.length).toBeGreaterThan(3)
  for (const k of t) expect(k.text).toBe(authorNo(y26, k.pos))
  expect(await page.getByTestId('zoom-controls').getAttribute('data-cell-px')).toBe(zoom)
  expect(await page.getByTestId('tracks-scroll').evaluate((e) => e.scrollLeft)).toBe(scroll)
  expect((await chainCell(page, 50).boundingBox())!.x).toBeCloseTo(xBefore, 1)
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', '4fe5'); await expect(page.getByTestId('instance-select')).toHaveValue('4fe5_1_B_A')
  const newReq = urls.slice(nReq)
  expect(newReq.filter((u) => u.includes('/mappings/')).every((u) => u.endsWith('/1y26_1_X_A.json'))).toBe(true)   // only the reference mapping
  expect(newReq.some((u) => /pdbe\/static\/entry|model-server|ModelServer|\.bcif|\.cif/.test(u))).toBe(false)          // no 2D / 3D reload

  // 3. tooltip: active chain's own residue, the reference residue separately, status
  await page.getByTestId('tracks-scroll').evaluate((e) => { e.scrollLeft = 0 })
  for (const pos of [50, 64, 5]) {
    await chainCell(page, pos).scrollIntoViewIfNeeded(); await chainCell(page, pos).hover()
    await expect(row(page, 'rfam')).toHaveText(String(pos))
    expect(await rowLabel(page, 'author')).toBe('4FE5 B author residue'); await expect(row(page, 'author')).toHaveText(expected(fe5, pos).value)
    expect(await rowLabel(page, 'reference')).toBe('1Y26 X reference residue'); await expect(row(page, 'reference')).toHaveText(expected(y26, pos).value)
  }
  expect(expected(y26, 64).value).toBe('Not present'); expect(expected(y26, 5).value).toBe('Unmapped')   // deletion and outside-alignment in the reference
  // coverage tooltips name the Rfam position and add the reference residue
  await page.getByTestId('consensus-tracks').locator('[data-row="seq"][data-tip-pos="50"]').hover()
  await expect(page.getByTestId('track-tooltip-title')).toHaveText('Sequence coverage · Rfam 50'); await expect(row(page, 'reference')).toHaveText(expected(y26, 50).value)

  // 4. switch the active instance: the numbering reference stays; mapped-without-coordinates chain state
  await page.getByTestId('instance-select').selectOption('2xnw_1_A_A'); await ready(page, '2xnw')
  await expect(page.getByTestId('numbering-select')).toHaveValue('1y26_1_X_A'); await expect(page.getByTestId('axis-convention')).toHaveText('1Y26 X · author numbering')
  await chainCell(page, 19).hover()
  await expect(page.getByTestId('track-tooltip-title')).toHaveText('Selected RNA chain · 2XNW A')
  await expect(row(page, 'author')).toHaveText(expected(xnw, 19).value); await expect(row(page, 'status')).toHaveText('Coordinates unavailable')
  await expect(row(page, 'reference')).toHaveText(expected(y26, 19).value)
  await chainCell(page, 64).hover(); await expect(row(page, 'author')).toHaveText('Not present'); await expect(row(page, 'status')).toHaveText('Deletion')
  await chainCell(page, 5).hover(); await expect(row(page, 'author')).toHaveText('Unmapped'); await expect(row(page, 'status')).toHaveText('Outside aligned region')

  // 5. selection summary: Rfam position, active chain author residue, reference residue; still drives 2D/3D by the Rfam position
  await page.mouse.move(5, 5)
  await selectColumn(page, 50)
  await expect(page.getByTestId('selection-status')).toHaveText(`Rfam position 50 → 2XNW A ${expected(xnw, 50).value} · 1Y26 X reference ${expected(y26, 50).value} · highlighted and centred in 2D/3D`)
  await expect.poll(() => page.evaluate(() => (window.__rnaExplorer?.selection as { positions: number[] }).positions)).toEqual([50])
  await expect.poll(() => page.evaluate(() => (window.__rnaExplorer?.rna2dSelected as number[] | undefined)?.length ?? 0)).toBe(1)

  // 6. reference = active chain: one author row only
  await page.getByTestId('numbering-select').selectOption('2xnw_1_A_A')
  await chainCell(page, 50).hover(); expect(await rowLabel(page, 'author')).toBe('Author residue'); await expect(row(page, 'reference')).toHaveCount(0)
  await expect(page.getByTestId('axis-convention')).toHaveText('2XNW A · author numbering')
  // 7. back to Rfam: original axis; selection unchanged
  await page.getByTestId('numbering-select').selectOption('')
  await expect(page.getByTestId('axis-convention')).toHaveText('Rfam RF00167 position'); await expect(page.getByTestId('consensus-axis').locator('[data-ref-tick]')).toHaveCount(0)
  expect(await page.evaluate(() => (window.__rnaExplorer?.selection as { positions: number[] }).positions)).toEqual([50])
})

test('RF00059: author insertion codes are shown as residue identifiers (9C6J A, 306A ...)', async ({ page }) => {
  await page.goto('/RF00059'); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 })
  const m = await J(page, 'RF00059', '9c6j_1_A_A')
  expect(authorNo(m, 4)).toBe('306A')
  await page.getByTestId('numbering-select').selectOption('9c6j_1_A_A')
  await expect(page.getByTestId('numbering-control')).toHaveAttribute('data-state', 'ready')
  while (await page.getByTestId('zoom-in').isEnabled()) await page.getByTestId('zoom-in').click()   // nucleotide detail
  await page.getByTestId('tracks-scroll').evaluate((e) => { e.scrollLeft = 0 })
  const t = await ticks(page); for (const k of t) expect(k.text).toBe(authorNo(m, k.pos))
  await chainCell(page, 4).hover()
  await expect(row(page, 'reference')).toHaveText(expected(m, 4).value); expect(expected(m, 4).value).toMatch(/306A$/)
})

test('RF00177: 16S reference numbering stays sparse at Fit all and follows zoom', async ({ page }) => {
  await page.goto('/RF00177'); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 120_000 })
  const m = await J(page, 'RF00177', '4v9d_1_BA_Y')
  await page.getByTestId('numbering-select').selectOption('4v9d_1_BA_Y')
  await expect(page.getByTestId('numbering-control')).toHaveAttribute('data-state', 'ready', { timeout: 60_000 })
  let t = await ticks(page); expect(t.length).toBeGreaterThan(5); expect(t.length).toBeLessThan(25)
  for (const k of t) expect(k.text).toBe(authorNo(m, k.pos))
  for (let i = 0; i < 6; i++) await page.getByTestId('zoom-in').click()
  t = await ticks(page); for (const k of t) expect(k.text).toBe(authorNo(m, k.pos))
  const xs = await page.getByTestId('consensus-axis').locator('[data-ref-tick] text').evaluateAll((es) => es.map((e) => { const r = e.getBoundingClientRect(); return [r.left, r.right] }))
  for (let i = 1; i < xs.length; i++) expect(xs[i][0]).toBeGreaterThan(xs[i - 1][1])   // no overlapping labels
})

test('shared with the Annotations tab: the same reference relabels the annotation axis and tooltips, survives tab switches, Rfam restores', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  const y26 = await J(page, 'RF00167', '1y26_1_X_A')
  await page.getByTestId('numbering-select').selectOption('1y26_1_X_A')
  await page.getByTestId('tab-annotations').click()
  await expect(page.getByTestId('numbering-select')).toHaveValue('1y26_1_X_A')                    // shared choice
  await expect(page.getByTestId('annotation-axis-convention')).toHaveText('1Y26 X · author numbering')
  const axis = page.getByTestId('annotation-axis')
  const t = await axis.locator('[data-ref-tick]').evaluateAll((ds) => ds.map((d) => ({ pos: +d.getAttribute('data-ref-tick')!, text: d.querySelector('span')!.textContent })))
  expect(t.length).toBeGreaterThan(3); for (const k of t) expect(k.text).toBe(authorNo(y26, k.pos))
  // tick centred on its Rfam column's marks (columns unchanged)
  const k0 = t[Math.floor(t.length / 2)]
  const tb = (await axis.locator(`[data-ref-tick="${k0.pos}"] span`).first().boundingBox())!, mk = page.getByTestId('track-sasa-median').locator(`rect[data-pos="${k0.pos}"]`)
  if (await mk.count()) { const mb = (await mk.boundingBox())!; expect(Math.abs(tb.x + tb.width / 2 - (mb.x + mb.width / 2))).toBeLessThan(3) }
  // tooltip and pinned details: Rfam residue + reference residue
  await page.getByTestId('track-JL').locator('rect[data-pos="50"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toContainText('Rfam residue 50'); await expect(page.getByTestId('annotation-tooltip')).toContainText(`1Y26 X reference residue: ${expected(y26, 50).value}`)
  await page.getByTestId('track-JL').locator('rect[data-pos="50"]').click()
  await expect(page.getByTestId('annotation-inspect')).toContainText(`1Y26 X reference residue: ${expected(y26, 50).value}`)
  // zoom keeps working with reference labels
  await page.getByTestId('annotations-zoom-in').click(); await page.getByTestId('annotations-zoom-in').click()
  for (const k of await axis.locator('[data-ref-tick]').evaluateAll((ds) => ds.map((d) => ({ pos: +d.getAttribute('data-ref-tick')!, text: d.querySelector('span')!.textContent })))) expect(k.text).toBe(authorNo(y26, k.pos))
  // switch back to Rfam here: the Family Overview follows
  await page.getByTestId('numbering-select').selectOption('')
  await expect(page.getByTestId('annotation-axis-convention')).toHaveText('Rfam RF00167 consensus position'); await expect(axis.locator('[data-ref-tick]')).toHaveCount(0)
  await expect(page.getByTestId('annotation-inspect')).not.toContainText('reference residue')
  await page.getByTestId('tab-overview').click()
  await expect(page.getByTestId('numbering-select')).toHaveValue(''); await expect(page.getByTestId('axis-convention')).toHaveText('Rfam RF00167 position')
})
