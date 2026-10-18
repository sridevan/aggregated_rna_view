import { expect, test, type Page } from '@playwright/test'
import { selectColumn } from './helpers'

type View = { occupancy: { width: number; height: number }; aspect: number; atoms: number; dir: number[] }
type Hook = { __rnaExplorer: { lastView?: View; getCamera: () => { position: number[]; target: number[]; up: number[]; radius: number } } }
const view = (page: Page) => page.evaluate(() => (window as unknown as Hook).__rnaExplorer.lastView!)
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }) }
const dirOf = (c: { position: number[]; target: number[] }) => { const d = c.target.map((t, i) => t - c.position[i]); const l = Math.hypot(...d); return d.map((x) => +(x / l).toFixed(6)) }

test('desktop layout: rows in order, compact family summary, 40/60 top-aligned metadata row, chips and headings', async ({ page }) => {
  await page.setViewportSize({ width: 1400, height: 1000 })
  await page.goto('/RF00167')
  await ready(page, '4fe5')
  const box = async (id: string) => (await page.getByTestId(id).boundingBox())!
  const about = await box('family-description'), stats = await box('stat-pdb-entries'), viewer = await box('viewer-card'), meta = await box('metadata-panel'), species = await box('species-chart')
  expect(about.y).toBeLessThan(stats.y); expect(stats.y).toBeLessThan(viewer.y); expect(viewer.y).toBeLessThan(meta.y)
  await expect(page.getByTestId('consensus-model')).toHaveCount(0)                                  // consensus model now lives (collapsed) in the track area, not as its own card
  // compact family summary (one card): three primary counts in one row, then the other statistics in a lightweight row, all values unchanged
  const fam = await (await page.request.get('/families/RF00167/family.json')).json()
  const primary = ['stat-pdb-entries', 'stat-rna-instances', 'stat-model-length']
  const pBoxes = await Promise.all(primary.map(async (id) => (await page.getByTestId(id).boundingBox())!))
  expect(new Set(pBoxes.map((b) => Math.round(b.y))).size).toBe(1)                                // one row
  const pFonts = await Promise.all(primary.map((id) => page.getByTestId(id).evaluate((e) => getComputedStyle(e).fontSize)))
  expect(new Set(pFonts).size).toBe(1)
  await expect(page.getByTestId('stat-model-length-card')).toContainText('Rfam consensus positions')
  await expect(page.getByTestId('stat-rna-instances-card')).toContainText(`RNA chain instances · ${fam.counts.instances_validation_ok} validated`)
  for (const id of ['stat-unique-sequences', 'stat-sequence-occupancy', 'stat-structural-occupancy']) {
    await expect(page.getByTestId(id)).toBeVisible()
    expect(await page.getByTestId(id).evaluate((e) => parseFloat(getComputedStyle(e).fontSize))).toBeLessThan(parseFloat(pFonts[0]))   // secondary row is lighter
  }
  await page.getByTestId('stat-secondary').getByTestId('info-tip').hover()
  await expect(page.getByRole('tooltip')).toContainText('Mean structural occupancy: mean over 102 positions of the fraction of 84 instances with coordinates')
  await page.mouse.move(0, 0)
  // the molecular view starts much higher than with the previous About card + six stat cards (~500 px at this width)
  expect(viewer.y - about.y).toBeLessThan(260)
  // viewer: full width, 500-600 px tall
  expect(viewer.width).toBeGreaterThan(about.width * 0.98)
  const canvasBox = await box('molstar-viewer'); expect(canvasBox.height).toBeGreaterThanOrEqual(500); expect(canvasBox.height).toBeLessThanOrEqual(600)
  // 40/60 split, same y, equal overall height
  expect(Math.abs(meta.width / (meta.width + species.width) - 0.4)).toBeLessThan(0.03)
  expect(Math.abs(meta.y - species.y)).toBeLessThan(2)                             // top-aligned; heights follow content
  const rep = await box('display-representative'), methods = await box('methods-card')
  expect(Math.abs(rep.height + methods.height + 16 - meta.height)).toBeLessThan(2)   // gap-4 = 16 px, no stretching
  expect(methods.y).toBeGreaterThan(rep.y + rep.height - 1)
  await expect(page.getByTestId('methods-list').locator('li span').first()).toHaveText('X-ray diffraction')
  // chips and heading sizes
  const chips = await page.getByTestId('section-eyebrow').allTextContents()
  for (const c of ['Molecular view', 'Structure', 'Statistics', 'Diversity']) expect(chips).toContain(c)   // the Reference card moved into the track area
  const h = page.getByTestId('viewer-card').locator('h2').first()
  await expect(h).toContainText('Representative RNA — 4FE5 chain B')
  const size = await h.evaluate((e) => parseFloat(getComputedStyle(e).fontSize)); expect(size).toBeGreaterThanOrEqual(17); expect(size).toBeLessThanOrEqual(21)
  const chipSize = await page.getByTestId('section-eyebrow').first().evaluate((e) => parseFloat(getComputedStyle(e).fontSize)); expect(chipSize).toBeLessThan(size)
  // representative metadata
  await expect(page.getByTestId('rep-id')).toContainText('4FE5 chain B')
  await expect(page.getByTestId('display-representative')).toContainText('X-ray diffraction'); await expect(page.getByTestId('display-representative')).toContainText('1.32 Å')
  await expect(page.getByTestId('display-representative')).toContainText('67 / 67'); await expect(page.getByTestId('display-representative')).toContainText('available')
  await page.getByTestId('display-representative').getByTestId('info-tip').hover()   // selection rule in the tooltip, not in the card body
  await expect(page.getByRole('tooltip')).toContainText('Selection rule')
  // no horizontal overflow
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true)
})

test('mobile layout: cards stack, no horizontal overflow, labels readable', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 })
  await page.goto('/RF00167')
  await ready(page, '4fe5')
  const meta = (await page.getByTestId('metadata-panel').boundingBox())!, species = (await page.getByTestId('species-chart').boundingBox())!
  expect(species.y).toBeGreaterThanOrEqual(meta.y + meta.height - 1)   // stacked
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)
  const label = page.locator('[data-testid^="species-row-tax:"]').first().getByTestId('species-label')
  expect((await label.boundingBox())!.width).toBeLessThanOrEqual(390)
  await expect(label).toBeVisible()
})

test('camera framing: tight fit of the visible region, fixed direction across instances, selection centres the camera, reset restores the fit', async ({ page }) => {
  await page.setViewportSize({ width: 1400, height: 1000 })
  await page.goto('/RF00167')
  await ready(page, '4fe5')
  const v0 = await view(page)
  expect(v0.atoms).toBeGreaterThan(1000)
  const occ0 = Math.max(v0.occupancy.width, v0.occupancy.height)
  expect(occ0).toBeGreaterThanOrEqual(0.6); expect(occ0).toBeLessThanOrEqual(0.95)   // the larger projected extent fills 60-95 % of the canvas
  expect(v0.occupancy.width).toBeLessThanOrEqual(0.95); expect(v0.occupancy.height).toBeLessThanOrEqual(0.95)
  const d0 = dirOf(await page.evaluate(() => (window as unknown as Hook).__rnaExplorer.getCamera()))
  for (const s of ['1y26_1_X_A', '1y27_1_X_A']) {
    await page.getByTestId('instance-select').selectOption(s)
    await ready(page, s.slice(0, 4))
    const v = await view(page); const occ = Math.max(v.occupancy.width, v.occupancy.height)
    expect(occ).toBeGreaterThanOrEqual(0.6); expect(occ).toBeLessThanOrEqual(0.95)
    expect(dirOf(await page.evaluate(() => (window as unknown as Hook).__rnaExplorer.getCamera()))).toEqual(d0)
  }
  // selecting a nucleotide centres the camera on it (target = residue centre, direction kept); reset restores the fit after a user rotation
  const camBefore = JSON.stringify(await page.evaluate(() => (window as unknown as Hook).__rnaExplorer.getCamera()))
  await selectColumn(page, 50)
  await expect(page.getByTestId('sel-highlightable')).toContainText('50 →')
  await expect.poll(async () => JSON.stringify(await page.evaluate(() => (window as unknown as Hook).__rnaExplorer.getCamera()))).not.toBe(camBefore)
  await page.waitForTimeout(500)
  const camSel = await page.evaluate(() => (window as unknown as Hook).__rnaExplorer.getCamera())
  const center = await page.evaluate(() => (window as unknown as { __rnaExplorer: { lastFocusCenter: number[] } }).__rnaExplorer.lastFocusCenter)
  expect(Math.hypot(...camSel.target.map((t, i) => t - center[i]))).toBeLessThan(1.5)
  expect(dirOf(camSel)).toEqual(d0)
  await page.getByTestId('molstar-viewer').scrollIntoViewIfNeeded()
  const vb = (await page.getByTestId('molstar-viewer').boundingBox())!
  await page.mouse.move(vb.x + vb.width * 0.5, vb.y + vb.height * 0.6); await page.mouse.down(); await page.mouse.move(vb.x + vb.width * 0.65, vb.y + vb.height * 0.75, { steps: 10 }); await page.mouse.up()
  await page.waitForTimeout(400)
  await page.getByTestId('reset-orientation').click()
  // position, target and up restored exactly; Mol* re-derives the clipping radius from the visible scene bounds (the selection's focus
  // representation is still shown), so the radius is compared within 1 %
  const pose = (c: { position: number[]; target: number[]; up: number[] }) => JSON.stringify({ position: c.position, target: c.target, up: c.up })
  const before = JSON.parse(camBefore)
  await expect.poll(async () => pose(await page.evaluate(() => (window as unknown as Hook).__rnaExplorer.getCamera()))).toBe(pose(before))
  expect(Math.abs((await page.evaluate(() => (window as unknown as Hook).__rnaExplorer.getCamera())).radius - before.radius)).toBeLessThan(0.01 * before.radius)
})
