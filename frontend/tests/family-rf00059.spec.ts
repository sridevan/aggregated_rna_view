import { expect, test, type Page } from '@playwright/test'

/* RF00059 (TPP riboswitch): third family, served at /RF00059. Expected values come from data/families/RF00059/*.json. */
const J = async (page: Page, f: string) => (await page.request.get(`/families/RF00059/${f}`)).json()
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 }) }

test('Family Overview: heading, counts, representative 2GDI X, the 3-nt insertion after 24, family selector', async ({ page }) => {
  const fam = await J(page, 'family.json'), inst = await J(page, 'instances.json')
  await page.goto('/RF00059'); await ready(page, fam.display_representative.pdb_id)
  await expect(page.getByTestId('family-heading')).toHaveText(fam.description)
  await expect(page.getByTestId('family-select')).toHaveValue('RF00059')
  await expect(page.getByText(`${inst.length} validated`)).toBeVisible()
  expect(fam.display_representative.seq_name).toBe('2gdi_1_X_A')
  const ins = (await J(page, 'mappings/2gdi_1_X_A.json')).residues.filter((r: { alignment_state: string }) => r.alignment_state === 'insert')
  expect(ins.map((r: { insertion_anchor: number }) => r.insertion_anchor)).toEqual([24, 24, 24])
  await page.getByTestId('insertion-marker-24').hover()
  await expect(page.getByTestId('insertion-tooltip-title')).toHaveText('Insertion · +3 nt')
  await expect(page.getByTestId('insertion-tooltip-sequence')).toHaveText(`Sequence: ${ins.map((r: { nucleotide_parent: string }) => r.nucleotide_parent).join('')}`)
  await page.mouse.move(5, 5)
  await page.getByTestId('family-select').selectOption('RF00162'); await expect(page).toHaveURL(/\/RF00162$/)
})

test('Structures: apo vs bound groups (35 + 2), duplicate floor line drawn, 9C6J not clustered', async ({ page }) => {
  const sc = await J(page, 'structure_clustering.json'), disc = sc.metrics.discrepancy.clustering
  expect(disc.groups.map((g: { n_instances: number }) => g.n_instances)).toEqual([35, 2]); expect(disc.floor.value).toBeCloseTo(0.1268, 4)
  await page.goto('/RF00059'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 })
  await expect(page.getByTestId('group-chip-2')).toContainText('2 RNA instances')
  await expect(page.getByTestId('dendro-floor').or(page.getByTestId('overview-floor')).first()).toBeAttached()   // a vertical SVG line has zero width, so 'attached', not 'visible'
  await page.getByTestId('group-chip-2').click()
  await expect.poll(async () => (await page.getByTestId('alignment-viewer').getAttribute('data-structures'))?.split(',')).toHaveLength(1)   // the apo group's representative alone
  await page.getByTestId('instance-search').fill('8f4o'); await page.getByTestId('row-8f4o_1_A_A').locator('[data-testid="cell-title"]').click(); await page.getByTestId('instance-search').fill('')
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-structures', '8f4o_1_A_A')
  await expect(page.getByTestId('fit-warning-8f4o_1_A_A')).toContainText('poor fit')                    // apo P5/P3 displaced: placement flagged
  await page.getByTestId('instance-search').fill('9c6j'); await page.getByTestId('row-9c6j_1_A_A').locator('[data-testid="cell-title"]').click()
  await expect(page.getByTestId('selected-not-clustered')).toContainText('13 of the 64 required core positions observed')
})

test('Annotations: TPP and analogue rows; additives not shown as ligands', async ({ page }) => {
  await page.goto('/RF00059'); await page.getByTestId('tab-annotations').click()
  await expect(page.getByText('Thiamine diphosphate · 9 entries', { exact: false }).first()).toBeVisible({ timeout: 30_000 })
  for (const t of ['Median SASA', 'hairpin loops', 'Pocket occurrence', 'GTP']) await expect(page.getByText(t, { exact: false }).first()).toBeVisible()
  await expect(page.getByText('TETRAETHYLENE GLYCOL', { exact: false })).toHaveCount(0)
})
