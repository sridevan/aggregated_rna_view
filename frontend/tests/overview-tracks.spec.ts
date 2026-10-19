import { expect, test, type Page } from '@playwright/test'
import { colX, selectColumn } from './helpers'

/* Family Overview Phase 1: three primary tracks, a collapsed "Reference and alignment details" section (alignment confidence, RF, SS_cons),
   one shared Rfam position axis; every strip on the same 102-column grid; selection shared with the 2D and 3D viewers. */
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }) }
/** screen x of the centre of consensus column `pos` in a strip (viewBox: label column 140 units, 9 units per position) */
const selRect = (page: Page, id: string) => page.getByTestId(id).locator('rect[opacity="0.14"]')

test('default: three primary tracks, reference collapsed (heading only), one axis; consensus model only in the track area', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  const tracks = page.getByTestId('consensus-tracks')
  const labels = await page.getByTestId('consensus-tracks-labels').locator('text').allTextContents()
  expect(labels.filter((t) => /Sequence coverage|Structural coverage|state per position/.test(t)).length).toBe(3)
  expect(labels.join(' ')).toMatch(/Sequence coverage.*Structural coverage.*4FE5 B/)
  expect(labels.join(' ')).not.toContain('Alignment confidence')
  await expect(page.getByTestId('reference-section')).toHaveAttribute('data-open', 'false')
  await expect(page.getByTestId('reference-toggle')).toHaveAttribute('aria-expanded', 'false')
  await expect(page.getByTestId('reference-toggle')).toHaveText('▶Reference and alignment details')
  await expect(page.getByTestId('reference-body')).toHaveCount(0)
  await expect(page.getByTestId('consensus-axis')).toHaveCount(1); await expect(page.getByTestId('consensus-axis-labels')).toContainText('Rfam RF00167 position')
  // order: tracks, reference heading, axis
  const y = async (id: string) => (await page.getByTestId(id).boundingBox())!.y
  expect(await y('consensus-tracks')).toBeLessThan(await y('reference-section')); expect(await y('reference-section')).toBeLessThan(await y('consensus-axis'))
  // structured per-position tooltips kept (wording covered in overview-tooltips.spec.ts)
  await tracks.locator('rect[data-row="seq"][data-tip-pos="50"]').hover()
  await expect(page.getByTestId('track-tooltip-title')).toHaveText('Sequence coverage · Rfam 50'); await expect(page.getByTestId('track-tooltip-main')).toHaveText('29/29 unique sequences (100%)')
  await tracks.locator('rect[data-row="chain"][data-tip-pos="50"]').hover()
  await expect(page.getByTestId('track-tooltip-value-status')).toHaveText('Coordinates available')
})

test('expanded: alignment confidence, RF and SS_cons on the same grid; collapse again; selection unaffected and shared with 2D/3D', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  const fam = await (await page.request.get('/families/RF00167/family.json')).json()
  await selectColumn(page, 76)
  await expect(page.getByTestId('sel-highlightable')).toContainText('76 → C56')
  await expect.poll(() => page.evaluate(() => (window as unknown as { __rnaExplorer: { focusRepresentation?: number[] } }).__rnaExplorer.focusRepresentation)).toEqual([56])
  await page.getByTestId('reference-toggle').press('Enter')                                                  // keyboard
  await expect(page.getByTestId('reference-toggle')).toHaveAttribute('aria-expanded', 'true')
  await expect(page.getByTestId('reference-toggle')).toHaveAttribute('aria-controls', (await page.getByTestId('reference-body').getAttribute('id'))!)
  const ref = page.getByTestId('reference-tracks')
  await expect(page.getByTestId('reference-tracks-labels')).toContainText('Alignment confidence')
  await ref.locator('rect[data-row="pp"][data-tip-pos="50"]').hover()
  await expect(page.getByTestId('track-tooltip-title')).toHaveText('Alignment confidence · Rfam 50'); await expect(page.getByTestId('track-tooltip-secondary')).toHaveText('29 aligned unique sequences')
  // RF and SS_cons: exactly the family strings, one character per consensus position
  expect((await page.getByTestId('reference-row-rf').locator('text').allTextContents()).map((t) => t[0]).join('')).toBe(fam.consensus_sequence)
  expect((await page.getByTestId('reference-row-ss').locator('text').allTextContents()).map((t) => t[0]).join('')).toBe(fam.ss_cons)
  await expect(page.getByTestId('consensus-model')).toHaveCount(0)                                                // model line removed; kept in the info popup
  await page.getByRole('button', { name: 'About the reference and alignment details' }).click()
  await expect(page.getByTestId('info-popover')).toContainText('Rfam consensus sequence, reference secondary structure and alignment confidence.'); await page.keyboard.press('Escape')
  await page.getByTestId('legend-toggle').click(); await expect(page.getByTestId('legend-reference')).toContainText('Rfam RF00167, Purine, 102 consensus match positions'); await page.getByTestId('legend-toggle').click()   // model metadata in Legend & help
  await expect(page.getByTestId('reference-tracks-labels')).toContainText('Consensus sequence (RF)'); await expect(page.getByTestId('reference-tracks-labels')).toContainText('Secondary structure (SS_cons)')
  // same grid: column centres of every strip coincide (selected column too) and the axis tick 50 lines up
  for (const pos of [1, 50, 76, 102]) {
    const xs = await Promise.all(['consensus-tracks', 'reference-tracks', 'consensus-axis'].map((id) => colX(page, id, pos)))
    expect(Math.max(...xs) - Math.min(...xs)).toBeLessThan(0.5)
  }
  for (const id of ['consensus-tracks', 'reference-tracks', 'consensus-axis']) {
    const r = (await selRect(page, id).first().boundingBox())!; expect(Math.abs(r.x + r.width / 2 - await colX(page, id, 76))).toBeLessThan(1)
  }
  const tick = (await page.getByTestId('consensus-axis').locator('text', { hasText: /^50$/ }).boundingBox())!
  expect(Math.abs(tick.x + tick.width / 2 - await colX(page, 'consensus-tracks', 50))).toBeLessThan(2)
  // clicking a column of the reference strip selects it everywhere
  const rb = (await ref.boundingBox())!
  await page.mouse.click(await colX(page, 'reference-tracks', 30), rb.y + 5)
  await expect(page.getByTestId('sel-highlightable')).toContainText('30 →')
  // collapse: selection unchanged, axis still aligned
  await page.getByTestId('reference-toggle').click()
  await expect(page.getByTestId('reference-body')).toHaveCount(0); await expect(page.getByTestId('sel-highlightable')).toContainText('30 →')
  const tick2 = (await page.getByTestId('consensus-axis').locator('text', { hasText: /^50$/ }).boundingBox())!
  expect(Math.abs(tick2.x + tick2.width / 2 - await colX(page, 'consensus-tracks', 50))).toBeLessThan(2)
  // instance change updates the selected-chain track
  await page.getByTestId('instance-select').selectOption('1y26_1_X_A'); await ready(page, '1y26')
  await expect(page.getByTestId('consensus-tracks-labels')).toContainText('1Y26 X')
  await expect(page.getByTestId('instance-track').locator('rect[stroke-dasharray], rect[fill="#c3c2b7"]').first()).toBeAttached()
})

test('narrow viewport: no horizontal page overflow, strips stay aligned', async ({ page }) => {
  await page.setViewportSize({ width: 600, height: 900 })
  await page.goto('/RF00167'); await ready(page, '4fe5')
  await page.getByTestId('reference-toggle').click()
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)
  const xs = await Promise.all(['consensus-tracks', 'reference-tracks', 'consensus-axis'].map((id) => colX(page, id, 60)))
  expect(Math.max(...xs) - Math.min(...xs)).toBeLessThan(0.5)
})
