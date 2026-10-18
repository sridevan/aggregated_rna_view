import { expect, test, type Page } from '@playwright/test'
import { colX, selectColumn, stripGeom } from './helpers'

/* Family Overview "Sequence and annotations": collapsible section, one horizontal scroll container, zoom (−, +, Fit all) on one shared scale. */
type H = { __rnaExplorer: { molstar: unknown; viewerIds: number[]; rna2dSelected?: number[]; selection: { positions: number[]; active: { kind: string; id?: string } | null } } }
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 }) }
const cellPx = async (page: Page) => +(await page.getByTestId('zoom-controls').getAttribute('data-cell-px'))!
const sel = (page: Page) => page.evaluate(() => (window as unknown as H).__rnaExplorer.selection)
const STRIPS = ['consensus-tracks', 'loop-tracks', 'reference-tracks', 'consensus-axis']
/** the centre of column `pos` coincides in every strip */
const aligned = async (page: Page, pos: number) => { const xs = await Promise.all(STRIPS.map((id) => colX(page, id, pos))); return Math.max(...xs) - Math.min(...xs) }
const open = async (page: Page) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  await page.evaluate(() => { (window as unknown as { __mol0: unknown }).__mol0 = (window as unknown as H).__rnaExplorer.molstar })
  await page.getByTestId('loops-toggle').click(); await page.getByTestId('reference-toggle').click(); await expect(page.getByTestId('loop-tracks')).toBeVisible()
}
const sameViewer = (page: Page) => page.evaluate(() => (window as unknown as { __mol0: unknown }).__mol0 === (window as unknown as H).__rnaExplorer.molstar)

test('default: expanded, Fit all (102 columns in view, no scrolling), numbering 1..102; zoom steps; alignment at min, default and max zoom', async ({ page }) => {
  await open(page)
  await expect(page.getByTestId('tracks-toggle')).toHaveAttribute('aria-expanded', 'true')
  const sc = page.getByTestId('tracks-scroll')
  expect(await sc.evaluate((e) => e.scrollWidth <= e.clientWidth + 1)).toBe(true)
  await expect(page.getByTestId('zoom-out')).toBeDisabled(); await expect(page.getByTestId('zoom-fit')).toBeDisabled()
  expect((await page.getByTestId('consensus-axis').locator('text').allTextContents())).toEqual(['1', ...Array.from({ length: 10 }, (_, k) => String((k + 1) * 10))])
  const fit = await cellPx(page); expect(fit).toBeGreaterThan(7)                                                   // letters readable at the default scale
  for (const pos of [1, 50, 102]) expect(await aligned(page, pos)).toBeLessThan(0.5)
  // zoom in to the maximum: wider columns, horizontal scrolling, still aligned, numbering unchanged
  let prev = fit
  while (await page.getByTestId('zoom-in').isEnabled()) { await page.getByTestId('zoom-in').click(); const c = await cellPx(page); expect(c).toBeGreaterThan(prev); prev = c }
  expect(prev).toBe(32); expect(await sc.evaluate((e) => e.scrollWidth > 2 * e.clientWidth)).toBe(true)
  await sc.evaluate((e) => { e.scrollLeft = e.scrollWidth })                                                       // the far end
  expect(await aligned(page, 100)).toBeLessThan(0.5)
  const ticks = (await page.getByTestId('consensus-axis').locator('text').allTextContents()).map(Number); expect(ticks.at(-1)).toBeLessThanOrEqual(102)
  expect(ticks.every((t) => t === 1 || t % 2 === 0)).toBe(true)                                                     // step 2 at 32 px per position
  // zoom out one step, then Fit all
  await page.getByTestId('zoom-out').click(); expect(await cellPx(page)).toBe(24)
  await page.getByTestId('zoom-fit').click(); expect(Math.abs((await cellPx(page)) - fit)).toBeLessThan(0.01)
  expect(await sameViewer(page)).toBe(true); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready')   // viewers never rebuilt
})

test('collapse keeps the selection and the viewers; zoomed selection: clicks select the column under the pointer, highlight kept, 3D selection scrolled into view', async ({ page }) => {
  await open(page)
  await selectColumn(page, 76); await expect.poll(async () => (await sel(page)).positions).toEqual([76])
  await page.getByTestId('tracks-toggle').click()
  await expect(page.getByTestId('tracks-scroll')).toHaveCount(0); await expect(page.getByTestId('zoom-controls')).toHaveCount(0)
  await expect(page.getByTestId('tracks-heading')).toBeVisible()
  expect((await sel(page)).positions).toEqual([76]); await expect(page.getByTestId('selection-status')).toContainText('Rfam position 76')
  await page.getByTestId('tracks-toggle').press('Enter')                                                           // keyboard
  await expect(page.getByTestId('tracks-scroll')).toBeVisible(); await expect(page.getByTestId('consensus-tracks').locator('rect[opacity="0.14"]')).toHaveCount(1)
  // zoom in: the selected column stays in view and highlighted
  for (let i = 0; i < 4; i++) await page.getByTestId('zoom-in').click()
  const sc = page.getByTestId('tracks-scroll'), vis = async (pos: number) => { const x = await colX(page, 'consensus-tracks', pos), b = (await sc.boundingBox())!; return x > b.x + 140 * 1.2 && x < b.x + b.width }
  expect(await vis(76)).toBe(true)
  // click a visible column while scrolled: the consensus position under the pointer
  const g = await stripGeom(page, 'consensus-tracks')
  const target = 70; expect(await vis(target)).toBe(true)
  await page.mouse.click(await colX(page, 'consensus-tracks', target), g.b.y + 3)
  await expect.poll(async () => (await sel(page)).positions).toEqual([target])
  // a 3D click on a residue outside the visible range scrolls its column into view (no re-zoom)
  const z = await cellPx(page)
  await sc.evaluate((e) => { e.scrollLeft = e.scrollWidth })                                                       // look at the far end
  expect(await vis(20)).toBe(false)
  await page.evaluate(() => { const ev = new MouseEvent('PDB.molstar.click', { bubbles: true }) as MouseEvent & { eventData?: unknown }; ev.eventData = { entry_id: '4fe5', entity_id: '1', auth_asym_id: 'B', seq_id: 2, auth_seq_id: 16, ins_code: '' }; document.querySelector('[data-testid="molstar-viewer"]')!.dispatchEvent(ev) })
  await expect.poll(async () => (await sel(page)).positions).toEqual([20]); await expect.poll(() => vis(20)).toBe(true)
  expect(await cellPx(page)).toBe(z)
  expect(await sameViewer(page)).toBe(true)
})

test('zoomed loops: blocks per original segment, partial blocks, the whole discontinuous loop selected; insertion marker anchored at max zoom', async ({ page }) => {
  await open(page)
  for (let i = 0; i < 3; i++) await page.getByTestId('zoom-in').click()
  await page.getByTestId('loop-block-IL_4FE5_001-0').scrollIntoViewIfNeeded()
  await page.getByTestId('loop-block-IL_4FE5_001-0').click()
  await expect.poll(async () => (await page.evaluate(() => (window as unknown as H).__rnaExplorer.viewerIds)).length).toBe(6)   // both segments, 6 nucleotides
  for (const k of [0, 1]) await expect(page.getByTestId(`loop-outline-IL_4FE5_001-${k}`)).toHaveCount(1)
  const g = await stripGeom(page, 'loop-tracks')
  for (const k of [0, 1]) { const b = page.getByTestId(`loop-block-IL_4FE5_001-${k}`), from = +(await b.getAttribute('data-from'))!, to = +(await b.getAttribute('data-to'))!
    expect(Math.abs(+(await b.getAttribute('width'))! - ((to - from + 1) * g.cell - 1.5))).toBeLessThan(0.01) }        // never merged into one interval
  // insertion marker: 9IWF A at max zoom, centred on the 77|78 boundary
  await page.getByTestId('instance-select').selectOption('9iwf_1_A_A'); await ready(page, '9iwf')
  while (await page.getByTestId('zoom-in').isEnabled()) await page.getByTestId('zoom-in').click()
  const mk = page.getByTestId('insertion-marker-77'); await mk.scrollIntoViewIfNeeded()
  const pill = (await mk.locator('rect').nth(1).boundingBox())!
  expect(Math.abs(pill.x + pill.width / 2 - await colX(page, 'consensus-tracks', 77, true))).toBeLessThan(1)
})

test('narrow window: controls usable, no page overflow, tracks scroll inside their own container', async ({ page }) => {
  await page.setViewportSize({ width: 420, height: 900 })
  await page.goto('/RF00167'); await ready(page, '4fe5')
  await expect(page.getByTestId('zoom-in')).toBeVisible()
  await page.getByTestId('zoom-in').click(); await page.getByTestId('zoom-in').click()
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)
  expect(await page.getByTestId('tracks-scroll').evaluate((e) => e.scrollWidth > e.clientWidth)).toBe(true)
  const lb = (await page.getByTestId('consensus-tracks-labels').boundingBox())!, sb = (await page.getByTestId('tracks-scroll').boundingBox())!
  await page.getByTestId('tracks-scroll').evaluate((e) => { e.scrollLeft = 300 })
  const lb2 = (await page.getByTestId('consensus-tracks-labels').boundingBox())!; expect(Math.abs(lb2.x - lb.x)).toBeLessThan(1); expect(lb.x).toBeGreaterThanOrEqual(sb.x - 1)   // labels stay in view
})
