import { expect, test, type Page } from '@playwright/test'

/* RF00177 (bacterial 16S rRNA, 60 curated E. coli chains): fourth family, served at /RF00177. Expected values from data/families/RF00177. */
const J = async (page: Page, f: string) => (await page.request.get(`/families/RF00177/${f}`)).json()
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 120_000 }) }

test('Family Overview: 60 chains, reference 4V9D BA, 1,533-position overview, loops loaded per chain (not the family-wide file)', async ({ page }) => {
  const fam = await J(page, 'family.json'), inst = await J(page, 'instances.json')
  expect(fam.display_representative.seq_name).toBe('4v9d_1_BA_Y'); expect(fam.display_representative.automatic_choice).toBe('4v9o_35_HA_PG')
  const urls: string[] = []; page.on('request', (r) => urls.push(r.url()))
  await page.goto('/RF00177'); await ready(page, '4v9d')
  await expect(page.getByTestId('family-heading')).toHaveText(fam.description); await expect(page.getByText(`${inst.length} validated`)).toBeVisible()
  await expect(page.getByTestId('instance-select')).toHaveValue('4v9d_1_BA_Y')
  // Fit all: every position present, no horizontal scrolling, axis labels sparse
  expect(await page.getByTestId('consensus-tracks').locator('[data-row="seq"]').count()).toBeGreaterThan(1500)
  expect(await page.getByTestId('tracks-scroll').evaluate((e) => e.scrollWidth <= e.clientWidth + 1)).toBe(true)
  const ticks = (await page.getByTestId('consensus-axis').locator('text').allTextContents()).map(Number); expect(ticks.length).toBeLessThan(20); expect(ticks.at(-1)).toBe(1500)
  // loops: only the selected chain's file is fetched
  await page.getByTestId('loops-toggle').click(); await expect(page.getByTestId('loop-tracks')).toBeVisible({ timeout: 60_000 })
  expect(urls.some((u) => u.endsWith('/families/RF00177/loops_by_instance/4v9d_1_BA_Y.json'))).toBe(true)
  expect(urls.some((u) => u.endsWith('/families/RF00177/loops.json'))).toBe(false)
  const nBlocks = await page.locator('[data-testid^="loop-block-"]').count(); expect(nBlocks).toBeGreaterThan(100)
  // zoom to nucleotide detail still selects the right consensus position
  while (await page.getByTestId('zoom-in').isEnabled()) await page.getByTestId('zoom-in').click()
  await expect(page.getByTestId('reference-section')).toBeVisible()
})

test('Structures: head-rotation groups 50 + 10 (provisional), no floor line; table opens with the representative of group 1', async ({ page }) => {
  const sc = await J(page, 'structure_clustering.json'), c = sc.metrics.discrepancy.clustering
  expect(c.groups.map((g: { n_instances: number }) => g.n_instances)).toEqual([50, 10]); expect(c.validation).toBe('provisional'); expect(c.floor.value).toBeNull()
  await page.goto('/RF00177'); await page.getByTestId('tab-structures').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 120_000 })
  await expect(page.getByTestId('group-chip-2')).toContainText('10 RNA instances')
  await expect.poll(() => page.locator('[data-testid="instance-table"] tbody tr').first().getAttribute('data-testid')).toBe(`row-${c.groups[0].representative}`)
  await expect(page.getByTestId('dendro-floor')).toHaveCount(0); await expect(page.getByTestId('overview-floor')).toHaveCount(0)
})

test('Annotations: readable axis (labels every 100), antibiotic rows without putrescine, no pocket category; pocket downloads unavailable', async ({ page }) => {
  await page.goto('/RF00177'); await page.getByTestId('tab-annotations').click()
  await expect(page.getByText('Ligand binding sites').first()).toBeVisible({ timeout: 60_000 })
  const labels = (await page.getByTestId('annotation-axis').locator('span').allTextContents()).filter((t) => /^\d+$/.test(t)).map(Number)
  expect(labels.length).toBeLessThanOrEqual(17); expect(labels.slice(1).every((v) => v % 100 === 0)).toBe(true)
  await expect(page.getByText('Apramycin · 3 entries')).toBeVisible()
  await expect(page.getByText('1,4-DIAMINOBUTANE', { exact: false })).toHaveCount(0)                     // putrescine: additive, not a ligand row
  // no pocket category (user decision), but the missing data are named rather than dropped silently
  await expect(page.getByTestId('category-pockets')).toHaveCount(0)
  await expect(page.getByTestId('categories-unavailable')).toHaveText('Not available for RF00177: Predicted binding pockets.')
  await page.getByTestId('download-button').click()
  await expect(page.getByTestId('download-pockets')).toHaveAttribute('aria-disabled', 'true')
  await expect(page.getByTestId('download-loops')).toHaveAttribute('href', /\/families\/RF00177\/loops\.json$/)
})
