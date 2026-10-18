import { expect, test, type Page } from '@playwright/test'

type Inst = { seq_name: string; pdb_id: string }
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }) }

test('first visit: display representative is selected by full identity, its page is shown and its row highlighted, Mol* loads it', async ({ page }) => {
  await page.goto('/RF00167')
  const fam = await (await page.request.get('/families/RF00167/family.json')).json()
  const inst: Inst[] = await (await page.request.get('/families/RF00167/instances.json')).json()
  const rep = fam.display_representative
  expect(rep.seq_name).toBe(`${rep.pdb_id}_${rep.entity_id}_${rep.auth_asym_id}_${rep.label_asym_id}`)
  await page.getByTestId('tab-structures').click()
  await expect(page.getByTestId('page-info')).toHaveText(`Page 1 of ${Math.ceil(inst.length / 10)}`)       // default order: page 1, the active group's representative on top
  await page.getByTestId('instance-search').fill(rep.pdb_id)
  await expect(page.getByTestId(`row-${rep.seq_name}`)).toHaveAttribute('aria-selected', 'true')            // the shared selection (4FE5 B) is unchanged
  await expect(page.locator('[data-testid="instance-table"] tbody tr[data-selected="true"]')).toHaveCount(1)
  await expect(page.getByTestId('structures-current')).toContainText(`${rep.pdb_id.toUpperCase()} chain ${rep.auth_asym_id}`)
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })                       // Structures tab: 3D alignment of the selected structural group
  await page.getByTestId('tab-overview').click()                                      // the single-instance Mol* view lives in the Family Overview
  await ready(page, rep.pdb_id)
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-mode', 'region')
  await page.getByTestId('viewer-details-toggle').click(); await expect(page.getByTestId('alignment-reference-structure')).toContainText('(this structure)')
})

test('manual selection persists across pages and tabs; hidden-by-search indication; no reset to the representative', async ({ page }) => {
  await page.goto('/RF00167')
  await page.getByTestId('tab-structures').click()
  const inst: Inst[] = await (await page.request.get('/families/RF00167/instances.json')).json()
  await page.getByTestId('instance-search').fill('1y26')
  await page.getByTestId('row-1y26_1_X_A').click()
  await page.getByTestId('instance-search').fill('')
  await expect(page.getByTestId('row-1y26_1_X_A')).toHaveAttribute('aria-selected', 'true')            // clearing the search keeps the selected row's page
  await page.getByTestId('page-next').click(); await page.getByTestId('page-next').click()
  await expect(page.getByTestId('structures-current')).toContainText('1Y26 chain X')        // paging never changes the selection
  await page.getByTestId('tab-overview').click()
  await expect(page.getByTestId('instance-select')).toHaveValue('1y26_1_X_A')
  await expect(page.getByTestId('overview-current')).toContainText('1Y26 chain X')
  await page.getByTestId('tab-structures').click()
  await expect(page.getByTestId('page-info')).toHaveText(/^Page 1 of/)                             // the tab opens on page 1 (active group representative on top)
  await expect(page.getByTestId('structures-current')).toContainText('1Y26 chain X')                // still selected (not reset to 4FE5)
  await page.getByTestId('instance-search').fill('1y26'); await expect(page.getByTestId('row-1y26_1_X_A')).toHaveAttribute('aria-selected', 'true'); await page.getByTestId('instance-search').fill('')
  // search that excludes the selection: 3D unchanged, clear indication, clear-search restores the row
  await page.getByTestId('instance-search').fill('thermo')
  await expect(page.getByTestId('selected-hidden')).toContainText('1Y26 chain X is hidden by the current search')
  await expect(page.getByTestId('row-1y26_1_X_A')).toHaveCount(0)
  await expect(page.getByTestId('structures-current')).toContainText('1Y26 chain X')
  await page.getByTestId('selected-hidden').getByRole('button').click()
  await expect(page.getByTestId('selected-hidden')).toHaveCount(0)
  await expect(page.getByTestId('row-1y26_1_X_A')).toHaveAttribute('aria-selected', 'true')
})

test('viewer never stays blank: repeated tab switching always reaches ready, with the white background and orientation applied', async ({ page }) => {
  await page.goto('/RF00167')
  for (let i = 0; i < 3; i++) {
    await page.getByTestId('tab-structures').click(); await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
    await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-ligands')).toBeVisible()
    await page.getByTestId('tab-overview').click(); await ready(page, '4fe5')
    const h = await page.evaluate(() => { const x = (window as unknown as { __rnaExplorer: { bgColor?: number; orientationApplied?: { fit_quality: string } | null } }).__rnaExplorer; return { bg: x.bgColor, q: x.orientationApplied?.fit_quality } })
    expect(h.bg).toBe(0xffffff); expect(h.q).toBe('good')
    await expect(page.getByTestId('molstar-loading-overlay')).toHaveCount(0)
    await expect(page.getByTestId('molstar-error')).toHaveCount(0)
  }
})
