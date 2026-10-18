import { expect, test, type Page } from '@playwright/test'

/* Final Annotations UI: compact header with dynamic counts, one shared teal occurrence-frequency scale (loops, ligands, pockets), separate SASA
   quantity scale, solid magenta nonstandard-nucleotide presence, one compact legend; no change to values, tracks or order. */
const open = async (page: Page) => { await page.goto('/RF00167'); await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-ligands')).toBeVisible() }
const TEAL_LIGHT = '#d2ecef', TEAL_DARK = '#0a4f5c'
const hex = (c: string) => [1, 3, 5].map((i) => parseInt(c.slice(i, i + 2), 16))
/** expected fill of the shared frequency scale for a value (same rule as scaleColor) */
const freqFill = (v: number) => { const bin = Math.min(4, Math.max(0, Math.floor(v * 5))); const a = hex(TEAL_LIGHT), b = hex(TEAL_DARK); return `rgb(${a.map((x, i) => Math.round(x + (b[i] - x) * (bin / 4))).join(',')})` }

test('compact header: "Annotations", dynamic family counts, Download still at the top right', async ({ page }) => {
  await open(page)
  const fam = await (await page.request.get('/families/RF00167/family.json')).json()
  const tab = page.getByTestId('annotations-tab')
  await expect(tab.locator('h2').first()).toHaveText('Annotations')
  await expect(page.getByTestId('annotations-summary')).toHaveText(`Residue-level annotations across ${fam.counts.instances} RNA instances from ${fam.counts.pdb_entries} PDB entries.`)
  expect(fam.counts.instances).toBe(84); expect(fam.counts.pdb_entries).toBe(66)
  for (const gone of ['Aggregated annotations on the Rfam consensus axis', 'Hover for a summary', 'click a mark to pin', 'independent of the Family Overview']) await expect(tab).not.toContainText(gone)
  await expect(tab.getByTestId('section-eyebrow')).toHaveCount(0)
  const h = (await tab.locator('h2').first().boundingBox())!, dl = (await page.getByTestId('download-button').boundingBox())!, card = (await tab.locator('section').first().boundingBox())!
  expect(Math.abs(dl.y + dl.height / 2 - (h.y + h.height / 2))).toBeLessThan(12); expect(card.x + card.width - (dl.x + dl.width)).toBeLessThan(24)
  await page.getByTestId('download-button').click(); await expect(page.getByTestId('download-items').getByRole('menuitem')).toHaveCount(9); await page.keyboard.press('Escape')
  // less vertical space above the tracks than before (header + legend + axis previously ended at ~380 px at 1400x1000)
  const axis = (await page.getByTestId('annotation-axis').boundingBox())!
  expect(axis.y - card.y).toBeLessThan(150)
})

test('one shared teal frequency scale for loops, ligands and pockets; SASA orange; nonstandard magenta; values unchanged', async ({ page }) => {
  await open(page)
  const tr = await (await page.request.get('/families/RF00167/tracks.json')).json()
  const pk = (await (await page.request.get('/families/RF00167/pocket_tracks.json')).json()).tracks.all_score_ge_0
  // every frequency mark in the three categories uses exactly the shared bin colour of its value
  const marks = await page.evaluate(() => [...document.querySelectorAll('[data-testid="category-loops"] rect[data-value], [data-testid="category-ligands"] rect[data-value], [data-testid="category-pockets"] rect[data-value]')]
    .map((r) => ({ cat: r.closest('section')!.getAttribute('data-testid'), v: Number(r.getAttribute('data-value')), fill: r.getAttribute('fill') })))
  expect(new Set(marks.map((m) => m.cat))).toEqual(new Set(['category-loops', 'category-ligands', 'category-pockets']))
  for (const m of marks) expect(m.fill).toBe(freqFill(m.v))
  // equal bins give equal colours across categories (e.g. >= 80 % in a loop, a ligand and the pocket track)
  const top = (cat: string) => new Set(marks.filter((m) => m.cat === cat && m.v >= 0.8).map((m) => m.fill))
  expect(top('category-loops')).toEqual(top('category-ligands')); expect(top('category-ligands').size).toBe(1)
  const low = (cat: string) => new Set(marks.filter((m) => m.cat === cat && m.v > 0 && m.v < 0.2).map((m) => m.fill))
  expect(low('category-loops')).toEqual(low('category-pockets'))
  // values unchanged (data-value is still the backend frequency)
  expect(await page.getByTestId('track-HL').locator('rect[data-pos="53"]').getAttribute('data-value')).toBe(tr.positions[52].hairpin_loops.frequency.toFixed(3))
  expect(await page.getByTestId('track-pocket-occurrence').locator('rect[data-pos="56"]').getAttribute('data-value')).toBe(pk[55].frequency.toFixed(3))
  // SASA keeps the orange quantity colours (not teal); nonstandard marks solid magenta
  const sasaFills = await page.getByTestId('track-sasa-median').locator('rect[data-value]').evaluateAll((rs) => rs.map((r) => r.getAttribute('fill')!))
  expect(sasaFills.length).toBe(102); expect(sasaFills.every((f) => { const [r, , b] = f.match(/\d+/g)!.map(Number); return r > b })).toBe(true)
  expect(new Set(await page.locator('[data-testid^="track-mod-"] rect[data-value]').evaluateAll((rs) => rs.map((r) => r.getAttribute('fill'))))).toEqual(new Set(['#c2185b']))
  // tooltips unchanged
  const tip = page.getByTestId('annotation-tooltip')
  await page.getByTestId('track-mod-50N').locator('rect[data-pos="22"]').hover()
  await expect(tip).toContainText('Rfam residue 2250N — '); await expect(tip).toContainText('Observed in 1 of 83 RNA instances (1.2%)')
  await page.getByTestId('track-pocket-occurrence').locator('rect[data-pos="56"]').hover()
  await expect(tip).toHaveText('Rfam residue 56Pocket occurrence: 61 of 80 RNA instances (76.3%)')
})

test('compact legend: shared frequency, SASA, nonstandard presence, no data vs no annotation; tracks and order intact; narrow screen', async ({ page }) => {
  await open(page)
  const legend = page.getByTestId('scale-legend')
  await expect(legend.locator('[data-testid^="legend-"]').evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')))).resolves.toEqual(
    ['legend-occurrence-frequency', 'legend-sasa-median', 'legend-gradient-sasa-median', 'legend-nonstandard-presence', 'legend-no-data', 'legend-no-annotation'])
  await expect(page.getByTestId('legend-occurrence-frequency')).toContainText('Occurrence frequency (0–100%)')
  for (const b of ['<20%', '20–40%', '40–60%', '60–80%', '≥80%']) await expect(page.getByTestId('legend-occurrence-frequency')).toContainText(b)
  for (const old of ['RNA loop frequency', 'Ligand contact frequency', 'Predicted pocket occurrence']) await expect(legend).not.toContainText(old)
  const t = await (await page.request.get('/families/RF00167/sasa_tracks.json')).json()
  await expect(page.getByTestId('legend-sasa-median')).toContainText(`Median SASA (Å²):0 Å²≥ ${t.colour_domain_suggestion_A2.max} Å²`)
  await expect(page.getByTestId('legend-nonstandard-presence')).toHaveText('Nonstandard nucleotide observed')
  await expect(page.getByTestId('legend-no-data')).toHaveText('No data'); await expect(page.getByTestId('legend-no-annotation')).toHaveText('Not observed')
  await expect(page.getByTestId('legend-no-annotation')).toHaveAttribute('title', /not a confirmed absence/)
  await expect(page.getByTestId('legend-no-data')).toHaveAttribute('title', /SASA: no qualifying measurement/)
  // legend is secondary: small text, at most two rows on desktop
  expect(await legend.evaluate((e) => parseFloat(getComputedStyle(e).fontSize))).toBeLessThanOrEqual(11)
  expect((await legend.boundingBox())!.height).toBeLessThanOrEqual(40)
  // categories, order, track counts
  const ids = await page.locator('section[data-testid^="category-"]').evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')))
  expect(ids).toEqual(['category-sasa', 'category-loops', 'category-modifications', 'category-ligands', 'category-pockets'])
  for (const [id, n] of [['sasa', 1], ['loops', 3], ['modifications', 4], ['ligands', 29], ['pockets', 1]] as const) await expect(page.getByTestId(`category-${id}`).locator('[data-testid^="track-"]')).toHaveCount(n)
  // axis alignment: tick 50 centred on position 50 of each category
  const tick = (await page.getByTestId('annotation-axis').locator('span', { hasText: /^50$/ }).boundingBox())!
  for (const tr of ['track-sasa-median', 'track-JL']) { const m = (await page.getByTestId(tr).locator('rect[data-pos="50"]').boundingBox())!; expect(Math.abs(tick.x + tick.width / 2 - (m.x + m.width / 2))).toBeLessThan(3) }
  // narrow: legend wraps without horizontal overflow; tooltip readable and inside
  await page.setViewportSize({ width: 390, height: 800 })
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)
  const lb = (await legend.boundingBox())!; expect(lb.x + lb.width).toBeLessThanOrEqual(390)
  const m = page.getByTestId('track-HL').locator('rect[data-pos="53"]'); await m.scrollIntoViewIfNeeded(); await m.hover()
  const tb = (await page.getByTestId('annotation-tooltip').boundingBox())!; expect(tb.x).toBeGreaterThanOrEqual(0); expect(tb.x + tb.width).toBeLessThanOrEqual(390)
})

test('no data (grey) and no annotation (empty) stay distinct with the shared scale', async ({ page }) => {
  // RF00167 has a determinable result at every position, so one pocket position is made undeterminable in the served copy (test-only; data on disk unchanged)
  await page.route('**/families/RF00167/pocket_tracks.json', async (route) => {
    const j = await (await route.fetch()).json()
    const p = j.tracks.all_score_ge_0[0]; p.n_unevaluable_instances += p.n_evaluable_instances; p.n_evaluable_instances = 0; p.n_supporting_instances = 0; p.frequency = null
    await route.fulfill({ json: j })
  })
  await open(page)
  const track = page.getByTestId('track-pocket-occurrence')
  await expect(track.locator('rect[data-pos="1"]')).toHaveAttribute('data-state', 'unevaluable')
  expect(await track.locator('rect[data-pos="1"]').getAttribute('fill')).toBe('#e4e2dc')
  // a determined position without a pocket draws no mark (empty = no annotation)
  const pk = (await (await page.request.get('/families/RF00167/pocket_tracks.json')).json()).tracks.all_score_ge_0
  const empty = pk.find((p: { rfam_consensus_pos: number; n_evaluable_instances: number; n_supporting_instances: number }) => p.rfam_consensus_pos > 1 && p.n_evaluable_instances > 0 && p.n_supporting_instances === 0).rfam_consensus_pos
  await expect(track.locator(`rect[data-pos="${empty}"]`)).toHaveCount(0)
  // legend swatches match the two renderings
  expect(await page.getByTestId('legend-no-data').locator('span').first().evaluate((e) => getComputedStyle(e).backgroundColor)).toBe('rgb(228, 226, 220)')
  expect(await page.getByTestId('legend-no-annotation').locator('span').first().evaluate((e) => getComputedStyle(e).backgroundColor)).toBe('rgb(255, 255, 255)')
})
