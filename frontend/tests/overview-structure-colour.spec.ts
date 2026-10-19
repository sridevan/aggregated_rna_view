import { expect, test, type Page } from '@playwright/test'
import { selectColumn } from './helpers'

/* Family Overview Mol*: one default RNA colour (slate #64748B) for every instance, in both view modes and after switching, resetting and
   clearing; selection (nucleotide or whole loop) in violet #7C3AED. Checked on the representation parameters and on the rendered pixels. */
const SLATE = 0x64748b, VIOLET = 0x7c3aed
type Rep = { component: string; theme: string | null; value: number | null; hidden: boolean }
type H = { __rnaExplorer: { molstarColours?: Rep[]; focusRepresentation?: number[]; lastMolstarSelection?: { color?: number } } }
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }) }
const reps = (page: Page) => page.evaluate(() => (window as unknown as H).__rnaExplorer.molstarColours ?? [])
/** every visible polymer / region representation is uniform slate */
const visibleRnaColours = async (page: Page) => (await reps(page)).filter((r) => !r.hidden && /Rfam-mapped region|Polymer|Selected RNA chain/.test(r.component)).map((r) => [r.theme, r.value])

/** hue classes of the rendered molecule (central crop of the canvas: no overlay buttons, badge or axes) */
const pixelClasses = async (page: Page) => {
  await page.waitForTimeout(700)
  const png = (await page.getByTestId('molstar-viewer').screenshot()).toString('base64')
  return page.evaluate(async (b64) => {
    const img = new Image(); img.src = `data:image/png;base64,${b64}`; await img.decode()
    const c = document.createElement('canvas'); c.width = img.width; c.height = img.height
    const g = c.getContext('2d')!; g.drawImage(img, 0, 0)
    const x0 = Math.round(img.width * 0.1), y0 = Math.round(img.height * 0.15), w = Math.round(img.width * 0.8), h = Math.round(img.height * 0.7)
    const d = g.getImageData(x0, y0, w, h).data
    const n = { slate: 0, violet: 0, other: 0 }
    for (let i = 0; i < d.length; i += 4) {
      const r = d[i] / 255, gg = d[i + 1] / 255, b = d[i + 2] / 255, mx = Math.max(r, gg, b), mn = Math.min(r, gg, b), s = mx ? (mx - mn) / mx : 0
      if (mx > 0.93 && s < 0.08) continue                              // white background / specular
      if (s < 0.1) continue                                            // achromatic (outlines, text): no hue information
      let hue = mx === r ? 60 * (((gg - b) / (mx - mn)) % 6) : mx === gg ? 60 * ((b - r) / (mx - mn) + 2) : 60 * ((r - gg) / (mx - mn) + 4); if (hue < 0) hue += 360
      if (hue >= 195 && hue <= 235 && s < 0.45) n.slate++
      else if (hue >= 245 && hue <= 290 && s >= 0.4) n.violet++
      else n.other++
    }
    return n
  }, png)
}
const frac = (n: { slate: number; violet: number; other: number }, k: keyof typeof n) => n[k] / Math.max(1, n.slate + n.violet + n.other)

test('same slate RNA colour for 1Y27 X and 2EET A (and 4FE5 B) in both modes; switching, reset and clear keep it', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  for (const [sn, pdb] of [['4fe5_1_B_A', '4fe5'], ['1y27_1_X_A', '1y27'], ['2eet_1_A_A', '2eet']] as const) {
    if (pdb !== '4fe5') { await page.getByTestId('instance-select').selectOption(sn); await ready(page, pdb) }
    await expect.poll(() => visibleRnaColours(page)).toEqual([['uniform', SLATE]])                   // region mode: the region cartoon only
    const px = await pixelClasses(page)
    expect(frac(px, 'slate'), `${pdb} region pixels ${JSON.stringify(px)}`).toBeGreaterThan(0.95); expect(px.violet).toBe(0)
    await page.getByTestId('mode-complete').click(); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-mode', 'complete')
    await expect.poll(async () => (await visibleRnaColours(page)).length).toBeGreaterThan(0)
    expect(new Set((await visibleRnaColours(page)).map((x) => JSON.stringify(x)))).toEqual(new Set([JSON.stringify(['uniform', SLATE])]))
    expect(frac(await pixelClasses(page), 'slate'), `${pdb} complete`).toBeGreaterThan(0.85)            // ligands/ions keep their element colours
    await page.getByTestId('mode-region').click(); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-mode', 'region')
    await page.getByTestId('reset-orientation').click()
    await expect.poll(() => visibleRnaColours(page)).toEqual([['uniform', SLATE]])
  }
})

test('selected nucleotide and selected loop in violet over slate; clearing restores the slate-only display', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  await page.getByTestId('instance-select').selectOption('1y27_1_X_A'); await ready(page, '1y27')
  // single nucleotide
  await selectColumn(page, 50)
  await expect.poll(() => page.evaluate(() => ((window as unknown as H).__rnaExplorer.focusRepresentation ?? []).length)).toBe(1)
  expect(await page.evaluate(() => (window as unknown as H).__rnaExplorer.lastMolstarSelection?.color)).toBe(VIOLET)
  let px = await pixelClasses(page); expect(px.violet, JSON.stringify(px)).toBeGreaterThan(300); expect(px.slate).toBeGreaterThan(px.violet)   // surrounding RNA stays slate
  // whole loop: every resolved member (all segments) violet, the same colour as a nucleotide
  await page.getByTestId('loops-toggle').click()
  await page.locator('[data-testid^="loop-block-J3_1Y27_001-"]').first().click()
  await expect.poll(() => page.evaluate(() => ((window as unknown as H).__rnaExplorer.focusRepresentation ?? []).length)).toBeGreaterThan(10)
  expect(await page.evaluate(() => (window as unknown as H).__rnaExplorer.lastMolstarSelection?.color)).toBe(VIOLET)
  px = await pixelClasses(page); expect(px.violet, JSON.stringify(px)).toBeGreaterThan(2000)
  expect((await reps(page)).filter((r) => r.component === 'Rfam-mapped region').map((r) => [r.theme, r.value])).toEqual([['uniform', SLATE]])   // base colour unchanged underneath
  // clearing: no violet left
  await page.getByTestId('clear-selection').click()
  await expect.poll(() => page.evaluate(() => ((window as unknown as H).__rnaExplorer.focusRepresentation ?? []).length)).toBe(0)
  px = await pixelClasses(page); expect(px.violet, JSON.stringify(px)).toBe(0); expect(frac(px, 'slate')).toBeGreaterThan(0.95)
})
