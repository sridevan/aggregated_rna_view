import { expect, test, type Page } from '@playwright/test'

/* RF00005 (tRNA, curated 60 chains): fifth family, served at /RF00005. Expected values come from data/families/RF00005/*.json. */
const J = async (page: Page, f: string) => (await page.request.get(`/families/RF00005/${f}`)).json()
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 }) }
type Res = { alignment_state: string; insertion_anchor: number | null; nucleotide_parent: string }

test('Family Overview: 60 chains, default 1EHZ A, curated selection, type II variable-arm insertion of 9QRP C', async ({ page }) => {
  const fam = await J(page, 'family.json'), inst = await J(page, 'instances.json')
  expect(fam.clen).toBe(71); expect(inst.length).toBe(60); expect(fam.provenance.selection.n_chains).toBe(60)
  expect(fam.display_representative.seq_name).toBe('1ehz_1_A_A')
  await page.goto('/RF00005'); await ready(page, '1ehz')
  await expect(page.getByTestId('family-heading')).toHaveText(fam.description)
  await expect(page.getByTestId('family-select')).toHaveValue('RF00005')
  await expect(page.getByText('60 validated')).toBeVisible()
  // the long variable arm of a type II tRNA is an insertion after consensus 44 (never extra columns)
  const ins = ((await J(page, 'mappings/9qrp_2_T_C.json')).residues as Res[]).filter((r) => r.alignment_state === 'insert' && r.insertion_anchor === 44)
  expect(ins.length).toBeGreaterThan(5)
  await page.getByTestId('instance-select').selectOption('9qrp_2_T_C'); await ready(page, '9qrp')
  await expect(page.getByTestId('consensus-tracks').locator('[data-row="seq"]')).toHaveCount(71)
  await page.getByTestId('insertion-marker-44').hover()
  await expect(page.getByTestId('insertion-tooltip-title')).toHaveText(`Insertion · +${ins.length} nt`)
  await expect(page.getByTestId('insertion-tooltip-sequence')).toHaveText(`Sequence: ${ins.map((r) => r.nucleotide_parent).join('')}`)
})

test('Structures: provisional groups 56 + 2 + 1 + 1 (discrepancy), no floor; the lambda-form 1J2B is its own group', async ({ page }) => {
  const sc = await J(page, 'structure_clustering.json'), disc = sc.metrics.discrepancy.clustering, tm = sc.metrics.tm.clustering
  expect(disc.groups.map((g: { n_instances: number }) => g.n_instances)).toEqual([56, 2, 1, 1]); expect(disc.validation).toBe('provisional'); expect(disc.floor.value).toBeNull()
  expect(tm.groups.map((g: { n_instances: number }) => g.n_instances)).toEqual([55, 2, 1, 1, 1]); expect(tm.validation).toBe('provisional')
  expect(disc.groups[3].members).toEqual(['1j2b_1_C_A'])
  await page.goto('/RF00005'); await page.getByTestId('tab-structures').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 })
  await expect(page.getByTestId('dendro-floor')).toHaveCount(0); await expect(page.getByTestId('overview-floor')).toHaveCount(0)
  await page.getByTestId('group-chip-4').click()
  await expect.poll(async () => (await page.getByTestId('alignment-viewer').getAttribute('data-structures'))?.split(',')).toEqual(['1j2b_1_C_A'])
})

test('Annotations: J5 junctions in the junction track, nonstandard nucleotides, few ligands, no pockets', async ({ page }) => {
  const tr = await J(page, 'tracks.json')
  expect(tr.tracks.junction_loops_union.loop_types).toEqual(['J3', 'J4', 'J5'])
  const p7 = tr.positions[6].junction_loops_union
  expect(p7.by_type_n_instances.J5).toBeGreaterThan(0)
  await page.goto('/RF00005'); await page.getByTestId('tab-annotations').click()
  await expect(page.getByText('Ligand binding sites').first()).toBeVisible({ timeout: 30_000 })
  for (const t of ['Median SASA', 'hairpin loops', 'junction loops', 'PSU', 'Bottromycin A2']) await expect(page.getByText(t, { exact: false }).first()).toBeVisible()
  await expect(page.getByText('Predicted binding pockets')).toHaveCount(0)
  // the junction mark at consensus 7 counts the union and names J5 in its details
  await page.getByTestId('track-JL').locator('rect[data-pos="7"]').click()
  await expect(page.getByTestId('annotation-inspect')).toContainText(`J5: ${p7.by_type_n_instances.J5} instances`)
  await expect(page.getByTestId('annotation-inspect')).toContainText(`Observed in ${p7.denominators.n_positive_observable_and_evaluated} of ${p7.denominators.n_observable_and_evaluated} RNA instances`)
  await page.getByTestId('download-button').click()
  await expect(page.getByTestId('download-pockets')).toContainText('not available for RF00005')
})

test('Structures: group 1 representative has no internal gaps (6UFG B); the gapped overall medoid 9TNL Y is named', async ({ page }) => {
  const g1 = (await J(page, 'structure_clustering.json')).metrics.discrepancy.clustering.groups[0]
  expect(g1.representative).toBe('6ufg_1_B_A'); expect(g1.medoid).toBe('9tnl_6_Y_F'); expect(g1.representative_internal_gaps).toEqual([])
  await page.goto('/RF00005'); await page.getByTestId('tab-structures').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 })
  await expect(page.getByTestId('group-representative')).toHaveText('6UFG B')
  await expect(page.getByTestId('group-medoid-note')).toContainText('9TNL Y')
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-structures', '6ufg_1_B_A')
})
