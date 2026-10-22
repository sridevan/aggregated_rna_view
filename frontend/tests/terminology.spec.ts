import { expect, test, type Page } from '@playwright/test'

/* Part A terminology: no user-facing "evaluable"/"unevaluable"/"denominator" wording in the Annotations tab (tooltips, pinned details, popovers,
   legend, download menu), and displayed counts/colour values unchanged (taken from the same JSON fields as before). */
const BANNED = /evaluable|unevaluable|evaluation denominator|source-evaluated|entry-entities/i

async function open(page: Page) {
  await page.goto('/RF00167')
  await page.getByTestId('tab-annotations').click()
  await expect(page.getByTestId('track-HL')).toBeVisible()
}

test('annotation wording uses "X of Y RNA instances" / "X of Y PDB entries"; no evaluable terminology anywhere visible', async ({ page }) => {
  await open(page)
  const tr = await (await page.request.get('/families/RF00167/tracks.json')).json()
  const pk = (await (await page.request.get('/families/RF00167/pocket_tracks.json')).json()).tracks.all_score_ge_0
  const sasa = await (await page.request.get('/families/RF00167/sasa_tracks.json')).json()
  const tip = page.getByTestId('annotation-tooltip'), inspect = page.getByTestId('annotation-inspect')
  // visible tab text (description, legend, axis, labels)
  expect(await page.getByTestId('annotations-tab').innerText()).not.toMatch(BANNED)
  // one representative mark per category: tooltip + pinned details (expanded)
  const marks: [string, number][] = [['track-HL', 53], ['track-JL', 53], ['track-IL', 30], ['track-lig-HPA', 26], ['track-pocket-occurrence', 53], ['track-sasa-median', 53], ['track-mod-CCC', 87]]
  for (const [id, pos] of marks) {
    const m = page.getByTestId(id).locator(`rect[data-pos="${pos}"]`)
    await m.scrollIntoViewIfNeeded(); await m.hover()
    await expect(tip).toBeVisible(); expect(await tip.innerText()).not.toMatch(BANNED)
    await m.click(); await inspect.locator('summary').click()
    expect(await inspect.innerText()).not.toMatch(BANNED)
    await page.getByTestId('annotations-clear').click()
  }
  // pocket at 53: tooltip numbers from the unchanged backend fields, details explain missing data in words
  const pm = page.getByTestId('track-pocket-occurrence').locator('rect[data-pos="53"]')
  await pm.hover()
  await expect(tip).toHaveText(`Rfam residue 53Pocket occurrence: ${pk[52].n_supporting_instances} of ${pk[52].n_evaluable_instances} RNA instances (16.3%)`)
  await pm.click(); await inspect.locator('summary').click()
  await expect(inspect).toContainText(`Not counted: ${pk[52].n_unevaluable_instances} RNA instances (no coordinates: ${pk[52].unevaluable_reasons.no_coordinates})`)
  await expect(inspect).toContainText('Frequencies are calculated from instances where annotation presence or absence can be determined.')
  await expect(inspect).toContainText('Successful predictions without a pocket here count as absent')
  await page.getByTestId('annotations-clear').click()
  // loop details: source-specific wording kept (not annotated by the source / source unavailable / overlap)
  const hl = page.getByTestId('track-HL').locator('rect[data-pos="53"]'); await hl.click(); await inspect.locator('summary').click()
  await expect(inspect).toContainText(`Not annotated as HL by the source in ${tr.positions[52].hairpin_loops.denominators.n_not_annotated_by_source} of 80 RNA instances`)
  await expect(inspect).toContainText('BGSU loop classifications can overlap.')
  await page.getByTestId('annotations-clear').click()
  // SASA keeps the quantitative wording, no positive/negative terms
  await page.getByTestId('track-sasa-median').locator('rect[data-pos="53"]').hover()
  const s53 = sasa.positions[52]
  await expect(tip).toHaveText(`Rfam residue 53Median SASA: ${s53.median_sasa_A2.toFixed(1)} Å²Based on ${s53.n_contributing_instances} RNA instances`)
  expect(await tip.innerText()).not.toMatch(/positive|negative|absent/i)
  // colour values unchanged: data-value still equals the backend frequency
  expect(await page.getByTestId('track-HL').locator('rect[data-pos="53"]').getAttribute('data-value')).toBe(tr.positions[52].hairpin_loops.frequency.toFixed(3))
  expect(await page.getByTestId('track-pocket-occurrence').locator('rect[data-pos="56"]').getAttribute('data-value')).toBe(pk[55].frequency.toFixed(3))
  // category popovers and download menu
  for (const c of ['sasa', 'loops', 'modifications', 'ligands', 'pockets']) {
    await page.getByTestId(`category-header-${c}`).getByTestId('info-tip').click()
    expect(await page.getByTestId('info-popover').innerText()).not.toMatch(BANNED)
    await page.keyboard.press('Escape')
  }
  await page.getByTestId('download-button').click()
  expect(await page.getByRole('menu').innerText()).not.toMatch(BANNED)
})
