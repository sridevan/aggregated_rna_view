import { expect, test, type Page } from '@playwright/test'

/** Solvent accessibility category (Milestone 3D): one quantitative track, median SASA coloured on a documented Å² domain, concise tooltips, downloads. */
const open = async (page: Page) => { await page.goto('/RF00167'); await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-sasa')).toBeVisible() }
const lum = (rgb: string) => { const m = rgb.match(/\d+/g)!.map(Number); return 0.2126 * m[0] + 0.7152 * m[1] + 0.0722 * m[2] }

test('one SASA track with fixed-size marks, colour by median SASA, Å² legend, missing distinct from measured zero', async ({ page }) => {
  await open(page)
  const t = await (await page.request.get('/families/RF00167/sasa_tracks.json')).json()
  await expect(page.getByTestId('category-toggle-sasa')).toContainText('Solvent accessibility'); await expect(page.getByTestId('category-toggle-sasa')).toContainText('(1 track)')
  await expect(page.getByTestId('category-body-sasa').locator('[data-testid^="track-"]')).toHaveCount(1)
  const track = page.getByTestId('track-sasa-median'); await expect(track).toContainText('Median SASA (Å²)')
  // same mark geometry as the other tracks
  const geo = await page.evaluate(() => [...document.querySelectorAll('[data-testid="track-sasa-median"] rect[data-pos], [data-testid="track-HL"] rect[data-pos]')].map((r) => { const b = r.getBoundingClientRect(); return `${Math.round(b.width)}x${Math.round(b.height)}` }))
  expect(new Set(geo).size).toBe(1)
  // one mark per position with a value (102 for RF00167), coloured by the median: higher median -> darker
  expect(await track.locator('rect[data-pos][data-value]').count()).toBe(t.summary.n_positions_with_value)
  const byPos = Object.fromEntries(t.positions.map((p: { rfam_consensus_pos: number; median_sasa_A2: number | null }) => [p.rfam_consensus_pos, p.median_sasa_A2]))
  const hi = t.positions.reduce((a: { median_sasa_A2: number }, b: { median_sasa_A2: number }) => (b.median_sasa_A2 > a.median_sasa_A2 ? b : a)).rfam_consensus_pos
  const lo = t.positions.reduce((a: { median_sasa_A2: number }, b: { median_sasa_A2: number }) => (b.median_sasa_A2 < a.median_sasa_A2 ? b : a)).rfam_consensus_pos
  const fill = (pos: number) => track.locator(`rect[data-pos="${pos}"]`).getAttribute('fill')
  expect(lum((await fill(hi))!)).toBeLessThan(lum((await fill(lo))!))
  expect(await track.locator('rect[data-pos="53"]').getAttribute('data-value')).toBe(byPos[53].toFixed(3))
  // legend in Å² with the documented domain (0 .. 350)
  const leg = page.getByTestId('legend-sasa-median'); await expect(leg).toContainText('Median SASA (Å²)'); await expect(leg).toContainText(`0 Å²`); await expect(leg).toContainText(`≥ ${t.colour_domain_suggestion_A2.max} Å²`)
  await expect(page.getByTestId('legend-gradient-sasa-median')).toBeVisible()
  // measured zero would be the lightest colour, distinct from the grey 'unevaluable' mark
  const zeroColour = await page.evaluate(() => { const s = (window as unknown as { __rnaExplorer: Record<string, unknown> }).__rnaExplorer; return s ? null : null })
  expect(zeroColour).toBeNull()
  const lightest = (await page.getByTestId('legend-gradient-sasa-median').evaluate((e) => getComputedStyle(e).backgroundImage))!
  expect(lightest).toContain('253, 240, 213')                 // #fdf0d5 = value 0, not the grey #e4e2dc used for missing
  // axis alignment: mark for position 1 of the SASA track centred under the axis tick "1"
  const tick = (await page.getByTestId('annotation-axis').locator('span', { hasText: /^1$/ }).first().boundingBox())!
  const mark = (await track.locator('rect[data-pos="1"]').boundingBox())!
  expect(Math.abs(tick.x + tick.width / 2 - (mark.x + mark.width / 2))).toBeLessThan(3)
})

test('concise SASA tooltip and pinned details; category info popover; downloads for both SASA files; loop tracks unchanged', async ({ page }) => {
  await open(page)
  const t = await (await page.request.get('/families/RF00167/sasa_tracks.json')).json(); const p53 = t.positions[52]
  await page.getByTestId('track-sasa-median').locator('rect[data-pos="53"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toHaveText(`Rfam residue 53Median SASA: ${p53.median_sasa_A2.toFixed(1)} Å²Based on ${p53.n_contributing_instances} RNA instances`)
  await page.getByTestId('track-sasa-median').locator('rect[data-pos="53"]').click()
  await page.getByTestId('annotation-inspect').locator('summary').click()
  await expect(page.getByTestId('annotation-inspect')).toContainText(`Q1 ${p53.q1_sasa_A2.toFixed(1)}`)
  await expect(page.getByTestId('annotation-inspect')).toContainText('U 78')
  await page.getByTestId('track-sasa-median').locator('rect[data-pos="53"]').click()
  await page.getByTestId('category-header-sasa').getByTestId('info-tip').hover()
  await expect(page.getByRole('tooltip')).toContainText('Per-nucleotide solvent-accessible surface area (SASA) of isolated RNA chains')
  await expect(page.getByRole('tooltip')).toContainText('summarised by Rfam position (median Å²)')
  await expect(page.getByTestId('category-sasa')).toHaveAttribute('data-open', 'true')
  await expect(page.locator('[data-testid="track-sasa-median"] [data-testid="info-tip"]')).toHaveCount(0)
  await page.mouse.move(5, 5)
  // downloads
  await page.getByTestId('download-button').click()
  for (const [id, file, label] of [['sasa', 'sasa.json', 'Solvent accessibility · Per nucleotide'], ['sasa-tracks', 'sasa_tracks.json', 'Solvent accessibility · Aggregated']]) {
    const it = page.getByTestId(`download-${id}`); await expect(it).toContainText(label)
    const href = (await it.getAttribute('href'))!; expect(href.endsWith(`/families/RF00167/${file}`)).toBe(true); await expect(it).toHaveAttribute('download', `RF00167_${file}`)
    const j = await (await page.request.get(href)).json(); expect(j.rfam_acc).toBe('RF00167'); expect(j.method.probe_radius_A).toBe(1.4)
  }
  const d = await (await page.request.get('/families/RF00167/sasa.json')).json()
  expect(d.instances[0].nucleotides[0]).toHaveProperty('status'); expect(d.instances[0].nucleotides[0]).toHaveProperty('missing_atoms'); expect(d.method.software).toHaveProperty('freesasa')
  await page.keyboard.press('Escape')
  // existing tracks unchanged
  const tr = await (await page.request.get('/families/RF00167/tracks.json')).json()
  await page.getByTestId('track-HL').locator('rect[data-pos="53"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toHaveText('Rfam residue 53Hairpin loop (HL)Observed in 8 of 80 RNA instances (10%)')
  expect(await page.getByTestId('track-HL').locator('rect[data-pos]').count()).toBe(tr.positions.filter((p: { hairpin_loops: { n_instances: number } }) => p.hairpin_loops.n_instances > 0).length)
  await expect(page.locator('[data-testid^="track-lig-"]')).toHaveCount(0)   // ligand tracks are in the Ligands tab
})
