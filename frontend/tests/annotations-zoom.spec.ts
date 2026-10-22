import { expect, test, type Page } from '@playwright/test'

/** centre x of a mark (or null when the position is outside the visible window) */
const markX = async (page: Page, track: string, pos: number) => {
  const r = page.getByTestId(track).locator(`rect[data-pos="${pos}"]`)
  if (!(await r.count())) return null
  const b = (await r.boundingBox())!; return b.x + b.width / 2
}
const cellPx = async (page: Page) => +(await page.getByTestId('annotations-zoom-controls').getAttribute('data-cell-px'))!

test('RF00167 Annotations: Fit all default, zoom in/out keeps rows and axis aligned, labels fixed, Fit all restores', async ({ page }) => {
  await page.goto('/RF00167'); await page.getByTestId('tab-annotations').click()
  const ctl = page.getByTestId('annotations-zoom-controls')
  await expect(ctl).toHaveAttribute('data-fit', 'true')
  await expect(page.getByTestId('annotations-zoom-out')).toBeDisabled(); await expect(page.getByTestId('annotations-zoom-fit')).toBeDisabled()
  const fit = await cellPx(page)
  const axis = page.getByTestId('annotation-axis')
  await expect(axis).toHaveAttribute('data-first-pos', '1')
  const clen = +(await axis.getAttribute('data-last-pos'))!; expect(clen).toBeGreaterThan(60)
  const labelX0 = (await page.getByTestId('track-HL').locator('div').first().boundingBox())!.x

  // inspect position 50, then zoom in twice: position 50 stays in view and all rows + axis tick stay aligned
  await page.getByTestId('track-JL').locator('rect[data-pos="50"]').click()
  for (let i = 0; i < 2; i++) await page.getByTestId('annotations-zoom-in').click()
  await expect(ctl).toHaveAttribute('data-fit', 'false')
  expect(await cellPx(page)).toBeGreaterThan(fit)
  const xs = (await Promise.all(['track-sasa-median', 'track-HL', 'track-IL', 'track-JL'].map((t) => markX(page, t, 50)))).filter((v): v is number => v != null)
  expect(xs.length).toBeGreaterThanOrEqual(2)                                  // tracks with a mark (or no-data mark) at 50
  for (const v of xs) expect(Math.abs(v - xs[0])).toBeLessThan(1.5)
  const tick = page.getByTestId('annotation-axis-plot').locator('span', { hasText: /^50$/ })
  const tb = (await tick.boundingBox())!; expect(Math.abs(tb.x + tb.width / 2 - xs[0]!)).toBeLessThan(3)
  expect((await page.getByTestId('track-HL').locator('div').first().boundingBox())!.x).toBeCloseTo(labelX0, 0)   // labels do not move

  // zoom to the maximum: only part of the model is visible; the shared scrollbar pans every row together
  while (await page.getByTestId('annotations-zoom-in').isEnabled()) await page.getByTestId('annotations-zoom-in').click()
  expect(await cellPx(page)).toBe(32)
  const bar = page.getByTestId('annotations-hscroll')
  await expect(bar).toBeVisible()
  const first0 = +(await axis.getAttribute('data-first-pos'))!
  await bar.evaluate((e) => { e.scrollLeft = 0 })
  await expect(axis).toHaveAttribute('data-first-pos', '1')
  expect(first0).toBeGreaterThan(1)
  expect(await markX(page, 'track-HL', clen)).toBeNull()                      // far end not rendered while scrolled to the start
  await bar.evaluate((e) => { e.scrollLeft = e.scrollWidth })
  await expect(axis).toHaveAttribute('data-last-pos', String(clen))
  const a = await markX(page, 'track-sasa-median', clen - 3), b = await markX(page, 'track-JL', clen - 3)
  if (a != null && b != null) expect(Math.abs(a - b)).toBeLessThan(1.5)

  // horizontal wheel pans too
  const before = +(await axis.getAttribute('data-first-pos'))!
  await page.getByTestId('track-HL').hover(); await page.mouse.wheel(-400, 0)
  await expect.poll(async () => +(await axis.getAttribute('data-first-pos'))!).toBeLessThan(before)

  // tooltip and inspect still work while zoomed
  const vis = +(await axis.getAttribute('data-first-pos'))! + 2
  await page.getByTestId('track-sasa-median').locator(`rect[data-pos="${vis}"]`).hover()
  await expect(page.getByTestId('annotation-tooltip')).toContainText(`Rfam residue ${vis}`)

  await page.getByTestId('annotations-zoom-out').click()
  expect(await cellPx(page)).toBe(24)
  await page.getByTestId('annotations-zoom-fit').click()
  await expect(ctl).toHaveAttribute('data-fit', 'true'); await expect(axis).toHaveAttribute('data-first-pos', '1'); await expect(axis).toHaveAttribute('data-last-pos', String(clen))
  expect(await bar.evaluate((e) => e.clientHeight)).toBe(0)
  // category collapse and vertical scroll are unchanged
  await page.getByTestId('category-toggle-loops').click(); await expect(page.getByTestId('category-loops')).toHaveAttribute('data-open', 'false')
})

test('RF00177 Annotations: 1,533 positions zoom to readable marks with thinned axis labels', async ({ page }) => {
  await page.goto('/RF00177'); await page.getByTestId('tab-annotations').click()
  await expect(page.getByTestId('category-loops')).toBeVisible({ timeout: 60_000 })
  const axis = page.getByTestId('annotation-axis')
  await expect(axis).toHaveAttribute('data-last-pos', '1533')
  expect(await cellPx(page)).toBeLessThan(1)
  for (let i = 0; i < 7; i++) await page.getByTestId('annotations-zoom-in').click()
  const px = await cellPx(page); expect(px).toBeGreaterThanOrEqual(8)
  const first = +(await axis.getAttribute('data-first-pos'))!, last = +(await axis.getAttribute('data-last-pos'))!
  expect(first).toBeGreaterThan(500); expect(last).toBeLessThan(1100); expect(last - first).toBeLessThan(250)   // zoomed about the centre
  // rendered marks limited to the window; mark width follows the zoom
  const n = await page.getByTestId('track-HL').locator('rect[data-pos]').count(); expect(n).toBeLessThanOrEqual(last - first + 1)
  const w = (await page.getByTestId('track-HL').locator('rect[data-pos]').first().boundingBox())!.width; expect(w).toBeGreaterThan(px - 3)
  const labels = (await page.getByTestId('annotation-axis-plot').locator('span').allTextContents()).filter((t) => /^\d+$/.test(t)).map(Number)
  const step = labels[2] - labels[1]; expect(step * px).toBeGreaterThanOrEqual(60); expect(step).toBeLessThan(100)
})
