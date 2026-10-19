import { expect, test, type Page } from '@playwright/test'

/* RF00011 (bacterial RNase P RNA class B, curated 50 chains of G. stearothermophilus P RNA, cryo-EM): sixth family, served at /RF00011.
   Expected values come from data/families/RF00011/*.json. */
const J = async (page: Page, f: string) => (await page.request.get(`/families/RF00011/${f}`)).json()
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 120_000 }) }
type Group = { n_instances: number; members: string[]; representative: string }

test('Family Overview: 50 chains, default 9OWQ A by override (no X-ray chain), 2D topology, viewer toolbar on one line', async ({ page }) => {
  const fam = await J(page, 'family.json'), inst = await J(page, 'instances.json')
  expect(fam.clen).toBe(365); expect(inst.length).toBe(50); expect(fam.provenance.selection.n_chains).toBe(50)
  expect(fam.display_representative.seq_name).toBe('9owq_1_A_A'); expect(fam.display_representative.selection_rule).toMatch(/^override \(user, 2026-10-19\)/)
  expect(fam.display_representative.automatic_choice).toBeNull()                       // all cryo-EM: the automatic rule finds no chain
  expect(new Set(inst.map((r: { experimental_method: string }) => r.experimental_method))).toEqual(new Set(['Electron Microscopy']))
  await page.goto('/RF00011'); await ready(page, '9owq')
  await expect(page.getByTestId('family-heading')).toHaveText(fam.description)
  await expect(page.getByTestId('family-select')).toHaveValue('RF00011')
  await expect(page.getByText('Explore aggregated RNA structures, conformational diversity, RNA–RNA interactions and residue-level annotations.')).toBeVisible()
  await expect(page.getByText('50 validated')).toBeVisible()
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 120_000 })
  // no small-molecule ligand near 9OWQ A (only Mg2+/Ca2+): Show ligands disabled, the note in the tooltip, every control on one row
  const ctl = page.getByTestId('show-ligands-control')
  await expect(page.getByTestId('show-ligands')).toBeDisabled(); await expect(ctl).toHaveText('Show ligands')
  const bar = page.getByTestId('molstar-toolbar'), top = (id: string) => bar.getByTestId(id).boundingBox().then((r) => Math.round(r!.y))
  expect(Math.abs(await top('molstar-clear-selection') - await top('mode-region'))).toBeLessThanOrEqual(2)
})

test('Structures: discrepancy k = 2 stable (closed 29 / open 21), TM k = 2 provisional (42 / 8), no floor', async ({ page }) => {
  const sc = await J(page, 'structure_clustering.json'), disc = sc.metrics.discrepancy.clustering, tm = sc.metrics.tm.clustering
  expect(disc.groups.map((g: Group) => g.n_instances)).toEqual([29, 21]); expect(disc.validation).toBe('stable'); expect(disc.floor.value).toBeNull()
  expect(tm.groups.map((g: Group) => g.n_instances)).toEqual([42, 8]); expect(tm.validation).toBe('provisional')
  expect(disc.groups[0].members).toContain('9owq_1_A_A')                                 // the closed default is in the closed group
  expect(disc.groups[1].members).toEqual(expect.arrayContaining(['9oy2_1_A_A', '9oy3_1_A_A', '9oy4_1_A_A']))   // tetraloop mutants: open
  expect(sc.dataset.n_clustered_instances).toBe(50); expect(sc.dataset.n_excluded_instances).toBe(0)
  await page.goto('/RF00011'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 120_000 })
  await expect(page.getByTestId('dendro-floor')).toHaveCount(0); await expect(page.getByTestId('overview-floor')).toHaveCount(0)
  await page.getByTestId('group-chip-2').click()
  await expect.poll(async () => (await page.getByTestId('alignment-viewer').getAttribute('data-structures'))?.split(',')).toEqual([disc.groups[1].representative])
})

test('Annotations and Interactions: no ligand tracks (metal ions only), no pockets, no nonstandard nucleotides, J6 junctions; FR3D for all 50 chains', async ({ page }) => {
  const tr = await J(page, 'tracks.json'), lig = await J(page, 'ligands.json'), it = await J(page, 'interactions.json')
  expect(tr.tracks.junction_loops_union.loop_types).toContain('J6')
  expect(Object.values(lig.by_ccd).every((c) => (c as { ligand_class: string }).ligand_class === 'metal_ion')).toBe(true)
  expect(it.default_numbering_reference.seq_name).toBe('9owq_1_A_A'); expect(it.summary.n_instances_annotated).toBe(50)
  await page.goto('/RF00011'); await page.getByTestId('tab-annotations').click()
  await expect(page.getByTestId('categories-unavailable')).toContainText('Predicted binding pockets')
  await expect(page.getByTestId('modifications-none')).toBeVisible()
  await page.getByTestId('tab-interactions').click()
  await expect(page.getByTestId('interactions-tab')).toHaveAttribute('data-circular', 'false')   // 365 positions: table and 3D only
  await expect(page.getByTestId('icircular-off')).toBeVisible()
})
