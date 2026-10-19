import { expect, test, type Page } from '@playwright/test'

/* Show ligands in the shared Mol* superposition viewer (Structures tab, Interactions Compare two): off by default; ligand instances within
   4 Å of each displayed chain (ModelServer surroundings), ions and additives left out; sticks with carbons in the structure's colour; the
   control sits in the toolbar above the canvas. */
type Lig = { show: boolean; built: string[]; visible: string[]; shown: Record<string, string[]> }
const hook = (page: Page, key: string) => page.evaluate((k) => (window.__rnaExplorer as Record<string, unknown> | undefined)?.[k] as Lig | undefined, key)
const ready = (page: Page, id: string) => expect(page.getByTestId(id)).toHaveAttribute('data-status', 'ready', { timeout: 120_000 })

test('Structures tab: TPP bound vs apo representatives, toggle, toolbar above the canvas', async ({ page }) => {
  await page.goto('/RF00059'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click(); await ready(page, 'alignment-viewer')
  const v = page.getByTestId('alignment-viewer'), ctl = v.getByTestId('viewer-show-ligands-control'), box = v.getByTestId('viewer-show-ligands')
  await expect(box).not.toBeChecked(); expect((await hook(page, 'alignment-viewerLigands'))?.built ?? []).toEqual([])      // off by default, nothing fetched
  const bar = (await v.getByTestId('viewer-toolbar').boundingBox())!, c = (await v.locator('canvas').first().boundingBox())!
  expect(c.y).toBeGreaterThanOrEqual(bar.y + bar.height - 1)
  await box.check()
  await expect.poll(async () => (await hook(page, 'alignment-viewerLigands'))?.shown['9c6k_1_X_A']).toEqual(['TPP X203'])
  await expect(ctl).toHaveAttribute('data-shown', '9c6k_1_X_A:TPP X203'); await expect(ctl).toHaveText('Show ligands')   // names only in the tooltip
  await expect(ctl).toHaveAttribute('title', /Shown: TPP/)
  // compare with the apo representative: its ligands are fetched, none nearby; the bound one keeps TPP
  await page.getByTestId('compare-on').click()
  await expect.poll(async () => (await v.getAttribute('data-structures')) ?? '').toContain(','); await ready(page, 'alignment-viewer')
  await expect.poll(async () => (await hook(page, 'alignment-viewerLigands'))?.shown).toEqual({ '9c6k_1_X_A': ['TPP X203'], '8f4o_1_A_A': [] })
  await expect(ctl).toHaveAttribute('data-shown', '9c6k_1_X_A:TPP X203;8f4o_1_A_A:')
  await box.uncheck()
  await expect.poll(async () => (await hook(page, 'alignment-viewerLigands'))?.visible).toEqual([])
  await expect(ctl).not.toHaveAttribute('title', /Shown:/)
})

test('Interactions Compare two: ligands follow the local base fit; apo chain has none', async ({ page }) => {
  await page.goto('/RF00167'); await page.getByTestId('tab-interactions').click(); await page.getByTestId('imode-compare').click()
  await page.getByTestId('icompare-a').selectOption('5e54_1_B_B'); await page.getByTestId('icompare-b').selectOption('5swe_1_X_A'); await ready(page, 'interaction-pair-viewer')
  const v = page.getByTestId('interaction-pair-viewer')
  await v.getByTestId('viewer-show-ligands').check()
  await expect.poll(async () => (await hook(page, 'interaction-pair-viewerLigands'))?.shown).toEqual({ '5e54_1_B_B': [], '5swe_1_X_A': ['ADE X101'] })
  await page.getByTestId('ipage-size').selectOption('0'); await page.getByTestId('icmp-row-basepair:26-56').click()
  await expect(v).toHaveAttribute('data-fit', 'local')
  expect((await hook(page, 'interaction-pair-viewerLigands'))?.visible).toEqual(['5e54_1_B_B', '5swe_1_X_A'])
  await expect(page.getByTestId('pair-viewer-label-a')).toBeVisible()
})
