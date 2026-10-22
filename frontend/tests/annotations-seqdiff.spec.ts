import { expect, test, type Page } from '@playwright/test'

/* Reported sequence differences (mmCIF _struct_ref_seq_dif of the family RNA chains; the PDBe API mutated_AA_or_NA covers proteins only):
   one presence track per deposited classification, verbatim; per-chain observations open their own chain in the Family Overview.
   Expected values come from sequence_difference_tracks.json, never from the rendering. */
type Obs = { pdb_id: string; auth_asym_id: string; seq_name: string; reference: string | null; deposited: string | null; author_residue_number: string; author_insertion_code: string; change_kind: string }
type SD = { classifications: Record<string, { n_on_consensus_axis: number; n_instances: number; n_entries: number }>; positions: { classifications: Record<string, { n_instances: number; n_entries: number; observation_ids: string[]; changes: Record<string, number> }> }[]
  observations: Record<string, Obs>; off_axis: Obs[]; counts: { n_off_axis: number } }
const open = async (page: Page, acc = 'RF00167') => { await page.goto(`/${acc}`); await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-loops')).toBeVisible() }
const J = async (page: Page, acc: string, f: string) => (await page.request.get(`/families/${acc}/${f}`)).json()
const slug = (c: string) => `seqdiff-${c.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`

test('category between loops and nonstandard nucleotides; one track per classification on the axis, labels verbatim; marks only where reported', async ({ page }) => {
  await open(page)
  const sd: SD = await J(page, 'RF00167', 'sequence_difference_tracks.json')
  expect(await page.locator('section[data-testid^="category-"]').evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')))).toEqual(['category-sasa', 'category-loops', 'category-seqdiff', 'category-modifications'])
  await expect(page.getByTestId('category-toggle-seqdiff')).toContainText('Reported sequence differences')
  const onAxis = Object.entries(sd.classifications).filter(([, c]) => c.n_on_consensus_axis > 0).sort((a, b) => b[1].n_instances - a[1].n_instances || a[0].localeCompare(b[0])).map(([c]) => c)
  expect(onAxis).toEqual(['Engineered mutation', 'variant', 'Conflict'])                                        // SEE REMARK 999 rows are all in flanks
  expect(await page.getByTestId('category-body-seqdiff').locator('[data-testid^="track-"]').evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')))).toEqual(onAxis.map((c) => `track-${slug(c)}`))
  for (const c of onAxis) {
    const want = sd.positions.flatMap((p, i) => (p.classifications[c] ? [i + 1] : []))
    expect(await page.getByTestId(`track-${slug(c)}`).locator('rect[data-pos]').evaluateAll((rs) => rs.map((r) => +r.getAttribute('data-pos')!))).toEqual(want)
  }
  for (const c of onAxis) await expect(page.getByTestId(`legend-${slug(c)}`)).toHaveText(c)                  // one legend colour per classification
  await page.getByTestId('category-header-seqdiff').getByTestId('info-tip').click()
  await expect(page.getByTestId('info-popover')).toContainText('an empty position does not mean the nucleotide is unchanged'); await page.keyboard.press('Escape')
  await expect(page.getByTestId('seqdiff-off-axis')).toContainText(`${sd.counts.n_off_axis} reported sequence differences have no Rfam position`)
  await expect(page.getByTestId('seqdiff-off-axis')).toContainText('4LX5 A G13 added')
  await expect(page.getByTestId('category-modifications')).toBeVisible()                                         // the nonstandard-nucleotide tracks are unchanged and separate
})

test('consensus 80 (C74U): tooltip, pinned observations and highlighting in the right RNA chain', async ({ page }) => {
  await open(page)
  const sd: SD = await J(page, 'RF00167', 'sequence_difference_tracks.json')
  const c = sd.positions[79].classifications['Engineered mutation']
  expect(c.changes).toEqual({ 'C>U': 2, 'U>C': 1 }); expect(c.observation_ids.map((id) => sd.observations[id].pdb_id).sort()).toEqual(['2b57', '3fo4', '4lx6'])
  const mark = page.getByTestId('track-seqdiff-engineered-mutation').locator('rect[data-pos="80"]')
  await mark.hover()
  await expect(page.getByTestId('annotation-tooltip')).toHaveText(`Rfam residue 80Engineered mutation: C→U ×2, U→CReported in ${c.n_instances} RNA instances (${c.n_entries} PDB entries)2B57 A C74U3FO4 A C74U4LX6 A U74C`)
  await mark.click()
  const inspect = page.getByTestId('annotation-inspect')
  await expect(inspect).toContainText('ClassificationEngineered mutation (as recorded in the PDB entry)')
  await expect(inspect).toContainText('not the Rfam consensus')
  await expect(page.getByTestId('annotation-observations').getByRole('button')).toHaveText(['2B57 A C74U', '3FO4 A C74U', '4LX6 A U74C'])
  await page.getByTestId('annotation-observation-2').click()                                                      // 4LX6 A
  await expect(page.getByTestId('tab-overview')).toHaveAttribute('aria-selected', 'true')
  await expect(page.getByTestId('instance-select')).toHaveValue('4lx6_1_A_A')
  await expect(page.getByTestId('selection-status')).toContainText('Rfam position 80 → 4LX6 A C74')                // the deposited residue in that chain
  await expect(page.getByTestId('selection-origin')).toContainText('Recorded for this RNA chain (4LX6 A).')
})

test('two copies in one entry are two RNA instances, each opened in its own chain', async ({ page }) => {
  await open(page)
  const sd: SD = await J(page, 'RF00167', 'sequence_difference_tracks.json')
  const v = sd.positions[19].classifications.variant                                                             // 6VWT and 6VWV, chains A and B
  expect(v.n_instances).toBe(4); expect(v.n_entries).toBe(2)
  await page.getByTestId('track-seqdiff-variant').locator('rect[data-pos="20"]').click()
  const btns = page.getByTestId('annotation-observations').getByRole('button')
  await expect(btns).toHaveText(['6VWT A U16A', '6VWT B U16A', '6VWV A U16A', '6VWV B U16A'])
  await page.getByTestId('annotation-observation-1').click()
  await expect(page.getByTestId('instance-select')).toHaveValue('6vwt_1_B_B')
  await expect(page.getByTestId('selection-status')).toContainText('Rfam position 20 → 6VWT B A16')
})

test('downloads; families without differences on the axis say so; numbering changes labels only', async ({ page }) => {
  await open(page)
  await page.getByTestId('download-button').click()
  for (const [id, file] of [['seqdiff', 'sequence_differences.json'], ['seqdiff-tracks', 'sequence_difference_tracks.json']]) {
    const href = (await page.getByTestId(`download-${id}`).getAttribute('href'))!; expect(href.endsWith(`/families/RF00167/${file}`)).toBe(true)
    expect((await page.request.get(href)).ok()).toBe(true)
  }
  await page.keyboard.press('Escape')
  const before = await page.getByTestId('category-body-seqdiff').locator('rect[data-pos]').evaluateAll((rs) => rs.map((r) => r.getAttribute('data-pos')))
  await page.getByTestId('numbering-select').selectOption('1y26_1_X_A')
  expect(await page.getByTestId('category-body-seqdiff').locator('rect[data-pos]').evaluateAll((rs) => rs.map((r) => r.getAttribute('data-pos')))).toEqual(before)
  await page.getByTestId('track-seqdiff-engineered-mutation').locator('rect[data-pos="80"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toContainText('1Y26 X reference residue:')
  await page.getByTestId('numbering-select').selectOption('')
  // SAM-I: two rows, neither on a consensus position -> no track, a note instead
  await open(page, 'RF00162')
  await expect(page.getByTestId('category-seqdiff')).toHaveCount(0)
  await expect(page.getByTestId('seqdiff-none')).toHaveText('Reported sequence differences: none at Rfam positions (2 differences have no Rfam position: flanks, insertions or absent residues).')
})
