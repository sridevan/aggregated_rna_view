import { expect, test, type Page } from '@playwright/test'
import { selectColumn } from './helpers'

/* Long chains (> 200 nt): a selection brings its region (helix and surroundings) into focus in the 2D viewer; Clear restores the whole topology.
 * Short chains (RF00167) keep the whole fitted topology. */
type Focus = { ids: number[]; target: { x: number; y: number; w: number; h: number }; full: { x: number; y: number; w: number; h: number } }
const vb = async (page: Page) => (await page.getByTestId('rna2d-viewer').locator('svg.rnaTopoSvg').getAttribute('viewBox'))!.split(' ').map(Number)
const ready2d = async (page: Page) => expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 120_000 })

test('RF00177: selecting a position zooms the 2D onto its region; the selected nucleotide is inside the view; Clear restores the fitted view', async ({ page }) => {
  await page.goto('/RF00177'); await ready2d(page)
  const full = await vb(page)
  await selectColumn(page, 700)
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-focused', 'true')
  const f = (await page.evaluate(() => window.__rnaExplorer?.rna2dFocus)) as Focus
  expect(f.ids.length).toBeGreaterThan(0)
  await expect.poll(async () => (await vb(page))[2]).toBeCloseTo(f.target.w, 0)    // animation finished on the target window
  const v = await vb(page)
  expect(v[2]).toBeLessThan(full[2] * 0.5)                                            // zoomed in
  expect(v[2]).toBeGreaterThan(full[2] * 0.03)                                        // but with context, not one nucleotide
  // the selected nucleotide is drawn inside the visible window
  const svg = page.getByTestId('rna2d-viewer').locator('svg.rnaTopoSvg')
  const box = (await svg.boundingBox())!
  const nt = await page.evaluate((id) => { const el = document.querySelector(`svg.rnaTopoSvg [class*="_${id}"]`) as SVGGraphicsElement | null; if (!el) return null; const r = el.getBoundingClientRect(); return { x: r.x + r.width / 2, y: r.y + r.height / 2 } }, f.ids[0])
  expect(nt).not.toBeNull(); expect(nt!.x).toBeGreaterThan(box.x); expect(nt!.x).toBeLessThan(box.x + box.width); expect(nt!.y).toBeGreaterThan(box.y); expect(nt!.y).toBeLessThan(box.y + box.height)
  // a second selection elsewhere moves the focus
  await selectColumn(page, 1400)
  await expect.poll(async () => ((await page.evaluate(() => window.__rnaExplorer?.rna2dFocus)) as Focus).ids[0]).not.toBe(f.ids[0])
  // Clear selection: whole topology again
  await page.getByTestId('clear-selection').first().click()
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-focused', 'false')
  await expect.poll(async () => (await vb(page))[2]).toBeCloseTo(full[2], 0)
})

test('RF00167 (short chain): a selection keeps the whole fitted 2D topology', async ({ page }) => {
  await page.goto('/RF00167'); await ready2d(page)
  const full = await vb(page)
  await selectColumn(page, 50)
  await expect.poll(async () => ((await page.evaluate(() => window.__rnaExplorer?.rna2dSelected)) as number[] | undefined)?.length ?? 0).toBeGreaterThan(0)
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-focused', 'false')
  expect(await vb(page)).toEqual(full)
})
