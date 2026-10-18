import { expect, test, type Page } from '@playwright/test'

/* RF00162 (SAM-I riboswitch): the second family served at /RF00162. Expected values come from data/families/RF00162/*.json. */
const J = async (page: Page, f: string) => (await page.request.get(`/families/RF00162/${f}`)).json()
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 }) }

test('Family Overview: heading, links, counts, representative, local-end marker, loops; family selector', async ({ page }) => {
  const fam = await J(page, 'family.json'), inst = await J(page, 'instances.json')
  await page.goto('/RF00162'); await ready(page, fam.display_representative.pdb_id)
  await expect(page.getByTestId('family-heading')).toHaveText(fam.description)
  await expect(page.getByTestId('family-badge')).toHaveAttribute('href', 'https://rfam.org/family/RF00162')
  await expect(page.getByTestId('family-select')).toHaveValue('RF00162')
  await expect(page.getByText(`${inst.length} validated`)).toBeVisible()
  await expect(page.getByTestId('instance-select')).toHaveValue(fam.display_representative.seq_name)
  // 5FJC A: the 2-nt local-end run (GA) after consensus 60 is a dashed marker, never a consensus column
  const m = (await J(page, `mappings/${fam.display_representative.seq_name}.json`)).residues.filter((r: { alignment_state: string }) => r.alignment_state === 'local_end')
  expect(m.map((r: { insertion_anchor: number }) => r.insertion_anchor)).toEqual([60, 60])
  const mk = page.getByTestId('insertion-marker-60'); await expect(mk).toHaveText('+2')
  await expect(mk.locator('rect').nth(1)).toHaveAttribute('stroke-dasharray', /\d/)
  await mk.hover()
  await expect(page.getByTestId('insertion-tooltip-title')).toHaveText('Local-end residues · +2 nt')
  await expect(page.getByTestId('insertion-tooltip-boundary')).toHaveText('Between Rfam positions 60 and 61')
  await expect(page.getByTestId('insertion-tooltip-sequence')).toHaveText(`Sequence: ${m.map((r: { nucleotide_parent: string }) => r.nucleotide_parent).join('')}`)
  // loops of the selected chain and whole-loop selection
  await page.mouse.move(5, 5); await page.getByTestId('loops-toggle').click()
  const first = page.locator('[data-testid^="loop-block-HL_5FJC_"]').first(); await first.click()
  await expect(page.getByTestId('selection-status')).toContainText('Hairpin loop · HL_5FJC_')
  // the family selector loads the other family's URL
  await page.getByTestId('family-select').selectOption('RF00167')
  await expect(page).toHaveURL(/\/RF00167$/); await expect(page.getByTestId('family-badge')).toHaveText('RF00167')
})

test('Structures: groups from structure_clustering.json, representative on top, no floor line, TM groups', async ({ page }) => {
  const sc = await J(page, 'structure_clustering.json')
  const disc = sc.metrics.discrepancy.clustering, tm = sc.metrics.tm.clustering
  await page.goto('/RF00162'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 })
  for (const g of disc.groups) await expect(page.getByTestId(`group-chip-${g.group_id}`)).toContainText(`${g.n_instances} RNA instance`)
  expect(disc.groups.map((g: { n_instances: number }) => g.n_instances)).toEqual([26, 2, 1])
  await expect.poll(() => page.locator('[data-testid="instance-table"] tbody tr').first().getAttribute('data-testid')).toBe(`row-${disc.groups[0].representative}`)
  expect(disc.floor.value).toBeNull()
  await expect(page.getByTestId('dendro-floor')).toHaveCount(0); await expect(page.getByTestId('overview-floor')).toHaveCount(0)
  // the singleton group: its representative alone in 3D
  await page.getByTestId('group-chip-3').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-structures', disc.groups[2].representative)
  await page.getByTestId('metric-tm').click()
  for (const g of tm.groups) await expect(page.getByTestId(`group-chip-${g.group_id}`)).toContainText(`${g.n_instances} RNA instance`)
})

test('Annotations: SAM / SAH / SFG ligand rows, the 5\' GTP, loops and pockets present', async ({ page }) => {
  await page.goto('/RF00162'); await page.getByTestId('tab-annotations').click()
  for (const t of ['SAM', 'SAH', 'SFG', 'GTP', 'Median SASA', 'hairpin loops', 'Pocket occurrence']) await expect(page.getByText(t, { exact: false }).first()).toBeVisible({ timeout: 30_000 })
  await expect(page.getByText('S-adenosylmethionine · 25 entries')).toBeVisible()
})
