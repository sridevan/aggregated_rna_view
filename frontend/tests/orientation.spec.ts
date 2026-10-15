import { expect, test, type Page } from '@playwright/test'
import { mkdirSync } from 'node:fs'

type Hook = { __rnaExplorer: { bgColor?: number; structuresLoaded?: number; orientationApplied?: { pdbId: string; rmsd: number; n_pairs: number; fit_quality: string } | null; lastMolstarSelection?: { labelSeqIds: number[] }; lastView?: { dir: number[]; up: number[] }; getCamera: () => { position: number[]; target: number[]; up: number[] } | null; selection: { positions: number[] } } }
const hook = (page: Page) => page.evaluate(() => { const h = (window as unknown as Hook).__rnaExplorer; return { bg: h.bgColor, n: h.structuresLoaded, orient: h.orientationApplied, view: h.lastView, cam: h.getCamera() } })
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready') }
const dirOf = (cam: { position: number[]; target: number[] }) => { const d = cam.target.map((t, i) => t - cam.position[i]); const l = Math.hypot(...d); return d.map((x) => x / l) }

test('single representative card; switching structures keeps a white canvas and a white overlay while loading', async ({ page }) => {
  await page.goto('/')
  await expect(page.getByTestId('display-representative')).toContainText('4FE5 chain B')
  await expect(page.getByTestId('representative')).toHaveCount(0)
  await ready(page, '4fe5')
  expect((await hook(page)).bg).toBe(0xffffff)
  await page.getByTestId('tab-structures').click()
  const seq = ['1y26_1_X_A', '1y27_1_X_A', '3rkf_1_A_A']
  for (const s of seq) {
    await page.getByTestId('instance-search').fill(s.slice(0, 4))
    await page.getByTestId(`row-${s}`).click()
    // while loading: overlay present and container white; after: renderer background white, exactly one structure loaded
    await expect(page.getByTestId('molstar-loading-overlay')).toBeVisible()
    expect(await page.getByTestId('molstar-loading-overlay').evaluate((e) => getComputedStyle(e).backgroundColor)).toBe('rgb(255, 255, 255)')
    await ready(page, s.slice(0, 4))
    await expect(page.getByTestId('molstar-loading-overlay')).toHaveCount(0)
    const h = await hook(page)
    expect(h.bg).toBe(0xffffff); expect(h.n).toBe(1)
  }
})

test('precomputed orientation is applied (not for insufficient fits), camera direction is fixed, selection still maps after the transform', async ({ page }) => {
  await page.goto('/')
  await page.getByTestId('tab-annotations').click()
  await page.getByTestId('instance-select').selectOption('1y26_1_X_A')
  await ready(page, '1y26')
  let h = await hook(page)
  expect(h.orient).toMatchObject({ pdbId: '1y26', fit_quality: 'good', n_pairs: 67 }); expect(h.orient!.rmsd).toBeCloseTo(1.8984, 3)
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-oriented', 'true')
  await expect(page.getByTestId('orientation-badge')).toContainText('RMSD 1.90 Å over 67 C4′ pairs (good)')
  const d1 = dirOf(h.cam!)
  // selection after the transform: consensus 50 -> label_seq_id 34 highlighted, camera unchanged
  await page.getByTestId('goto-input').fill('50'); await page.getByTestId('goto-button').click()
  await expect(page.getByTestId('sel-highlightable')).toContainText('50 → G34 (auth 46')
  await expect.poll(() => page.evaluate(() => (window as unknown as Hook).__rnaExplorer.lastMolstarSelection?.labelSeqIds)).toEqual([34])
  await page.waitForTimeout(500)
  const h2 = await hook(page); expect(dirOf(h2.cam!)).toEqual(d1)
  // another instance: same view direction in the common frame, selection preserved
  await page.getByTestId('instance-select').selectOption('1y27_1_X_A')
  await ready(page, '1y27')
  h = await hook(page)
  expect(h.orient).toMatchObject({ pdbId: '1y27', fit_quality: 'good' })
  const d2 = dirOf(h.cam!); for (let i = 0; i < 3; i++) expect(Math.abs(d2[i] - d1[i])).toBeLessThan(1e-6)
  expect(await page.evaluate(() => (window as unknown as Hook).__rnaExplorer.selection.positions)).toEqual([50])
  await expect.poll(() => page.evaluate(() => (window as unknown as Hook).__rnaExplorer.lastMolstarSelection?.labelSeqIds)).toEqual([33])
  // reference: identity; insufficient fit (9nbh): shown as deposited
  await page.getByTestId('instance-select').selectOption('4fe5_1_B_A')
  await ready(page, '4fe5')
  await expect(page.getByTestId('orientation-badge')).toContainText('reference structure')
  await page.getByTestId('instance-select').selectOption('9nbh_1_A_A')
  await ready(page, '9nbh')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-oriented', 'false')
  await expect(page.getByTestId('orientation-badge')).toContainText('as deposited (fit insufficient')
  // reset orientation restores the fixed direction after the user rotates the camera
  await page.getByTestId('instance-select').selectOption('1y26_1_X_A')
  await ready(page, '1y26')
  await page.getByTestId('molstar-viewer').scrollIntoViewIfNeeded()
  const vb = (await page.getByTestId('molstar-viewer').boundingBox())!   // drag inside the 3D canvas wherever the tab layout puts it
  await page.mouse.move(vb.x + vb.width * 0.5, vb.y + vb.height * 0.6); await page.mouse.down(); await page.mouse.move(vb.x + vb.width * 0.65, vb.y + vb.height * 0.75, { steps: 10 }); await page.mouse.up()
  await page.waitForTimeout(400)
  const rotated = dirOf((await hook(page)).cam!)
  expect(rotated.some((x, i) => Math.abs(x - d1[i]) > 1e-3)).toBe(true)
  await page.getByTestId('reset-orientation').click()
  await expect.poll(async () => dirOf((await hook(page)).cam!).map((x) => +x.toFixed(6))).toEqual(d1.map((x) => +x.toFixed(6)))
})

test('screenshots of several structures in the common reference orientation', async ({ page }) => {
  mkdirSync('../docs/screenshots', { recursive: true })
  await page.setViewportSize({ width: 1400, height: 1000 })
  await page.goto('/')
  await page.getByTestId('tab-annotations').click()
  for (const s of ['4fe5_1_B_A', '1y26_1_X_A', '1y27_1_X_A', '3rkf_1_A_A', '2xnw_1_A_A']) {
    await page.getByTestId('instance-select').selectOption(s)
    await ready(page, s.slice(0, 4))
    await page.waitForTimeout(1200)
    await page.getByTestId('molstar-viewer').screenshot({ path: `../docs/screenshots/m2_2_orientation_${s.slice(0, 4)}.png` })
  }
  await page.setViewportSize({ width: 1400, height: 1000 })
  await page.getByTestId('tab-overview').click()
  await page.getByTestId('tab-structures').click(); await page.getByTestId('instance-search').fill('4fe5'); await page.getByTestId('row-4fe5_1_B_A').click()
  await page.getByTestId('tab-overview').click()
  await ready(page, '4fe5')
  await page.waitForTimeout(800)
  await page.screenshot({ path: '../docs/screenshots/m2_2_family_overview.png', fullPage: true })
})
