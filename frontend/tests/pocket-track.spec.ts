import { expect, test, type Page } from '@playwright/test'

/** Predicted binding pockets category (Milestone 3E): one frequency track from pocket_tracks.json (default filter all Score >= 0), concise tooltip, downloads; other tracks unchanged. */
const open = async (page: Page) => { await page.goto('/RF00167'); await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-pockets')).toBeVisible() }
const lum = (rgb: string) => { const m = rgb.match(/\d+/g)!.map(Number); return 0.2126 * m[0] + 0.7152 * m[1] + 0.0722 * m[2] }

test('one pocket-occurrence track: marks, colours and tooltips follow pocket_tracks.json; category order and popover', async ({ page }) => {
  await open(page)
  const pt = await (await page.request.get('/families/RF00167/pocket_tracks.json')).json(); const tr = pt.tracks.all_score_ge_0
  await expect(page.getByTestId('category-toggle-pockets')).toContainText('Predicted binding pockets'); await expect(page.getByTestId('category-toggle-pockets')).toContainText('(1 track)')
  await expect(page.getByTestId('category-body-pockets').locator('[data-testid^="track-"]')).toHaveCount(1)
  const track = page.getByTestId('track-pocket-occurrence'); await expect(track).toContainText('Pocket occurrence')
  // order: loops, ligands, pockets, sasa
  const ys = await Promise.all(['sasa', 'loops', 'modifications', 'ligands', 'pockets'].map(async (id) => (await page.getByTestId(`category-${id}`).boundingBox())!.y))
  expect(ys).toEqual([...ys].sort((a, b) => a - b))
  // marks only where supporting instances exist; same geometry as the other tracks
  expect(await track.locator('rect[data-pos][data-value]').count()).toBe(tr.filter((p: { n_supporting_instances: number }) => p.n_supporting_instances > 0).length)
  const geo = await page.evaluate(() => [...document.querySelectorAll('[data-testid="track-pocket-occurrence"] rect[data-pos], [data-testid="track-HL"] rect[data-pos]')].map((r) => { const b = r.getBoundingClientRect(); return `${Math.round(b.width)}x${Math.round(b.height)}` }))
  expect(new Set(geo).size).toBe(1)
  // colour follows frequency: consensus 56 (0.76) darker than consensus 27 (~0.29)
  const fill = (pos: number) => track.locator(`rect[data-pos="${pos}"]`).getAttribute('fill')
  expect(lum((await fill(56))!)).toBeLessThan(lum((await fill(27))!))
  expect(await track.locator('rect[data-pos="56"]').getAttribute('data-value')).toBe(tr[55].frequency.toFixed(3))
  await expect(page.getByTestId('legend-occurrence-frequency')).toContainText('Occurrence frequency')
  // tooltip: concise, with the exact numerator / denominator
  await track.locator('rect[data-pos="56"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toHaveText('Rfam residue 56Pocket occurrence: 61 of 80 RNA instances (76.3%)')
  expect(tr[55].n_supporting_instances).toBe(61); expect(tr[55].n_evaluable_instances).toBe(80)
  // pinned details: entries, distinct pockets, score range, filter
  await track.locator('rect[data-pos="56"]').click()
  await page.getByTestId('annotation-inspect').locator('summary').click()
  await expect(page.getByTestId('annotation-inspect')).toContainText(`${tr[55].n_supporting_entries} PDB entries; ${tr[55].n_distinct_pockets} distinct predicted pockets`)
  await expect(page.getByTestId('annotation-inspect')).toContainText('fpocket Score')
  await track.locator('rect[data-pos="56"]').click()
  // category popover with the caveats; no track-level icon; clicking does not collapse
  await page.getByTestId('category-header-pockets').getByTestId('info-tip').hover()
  await expect(page.getByRole('tooltip')).toContainText('fpocketR-predicted RNA pockets aggregated by Rfam position')
  await expect(page.getByRole('tooltip')).toContainText('No predicted pocket does not establish absence of binding')
  await expect(page.getByTestId('category-pockets')).toHaveAttribute('data-open', 'true')
  await expect(page.locator('[data-testid="track-pocket-occurrence"] [data-testid="info-tip"]')).toHaveCount(0)
  await page.mouse.move(5, 5)
  // determined absences are empty, not grey: every position has a determinable result in RF00167, so no grey mark
  expect(await track.locator('rect[data-state="unevaluable"]').count()).toBe(tr.filter((p: { n_evaluable_instances: number }) => p.n_evaluable_instances === 0).length)
})

test('downloads for both pocket files; existing loop, ligand and SASA tracks unchanged', async ({ page }) => {
  await open(page)
  await page.getByTestId('download-button').click()
  for (const [id, file, label] of [['pockets', 'pockets.json', 'Predicted pockets · Per structure'], ['pocket-tracks', 'pocket_tracks.json', 'Predicted pockets · Aggregated']]) {
    const it = page.getByTestId(`download-${id}`); await expect(it).toContainText(label)
    const href = (await it.getAttribute('href'))!; expect(href.endsWith(`/families/RF00167/${file}`)).toBe(true); await expect(it).toHaveAttribute('download', `RF00167_${file}`)
    const j = await (await page.request.get(href)).json(); expect(j.rfam_acc).toBe('RF00167'); expect(j.method.fpocket_parameters.i).toBe(42)
  }
  const pk = await (await page.request.get('/families/RF00167/pockets.json')).json()
  expect(pk.counts.status_counts).toEqual({ pockets: 65, no_pockets: 19 }); expect(pk.instances[0].pockets.length === 0 || pk.instances[0].pockets[0]).toBeTruthy()
  await page.keyboard.press('Escape')
  await page.getByTestId('track-HL').locator('rect[data-pos="53"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toHaveText('Rfam residue 53Hairpin loop (HL)Observed in 66 of 80 RNA instances (82.5%)')
  await expect(page.locator('[data-testid^="track-lig-"]')).toHaveCount(29)
  await expect(page.getByTestId('track-sasa-median')).toBeVisible()
})
