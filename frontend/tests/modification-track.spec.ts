import { expect, test, type Page } from '@playwright/test'

/* Nonstandard nucleotides category: one Occurrence track from modification_tracks.json, fixed-colour presence marks, exact counts in tooltips,
   pinned details, category popover, two downloads; optional dataset (missing file hides only this category). */
const open = async (page: Page) => { await page.goto('/RF00167'); await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-ligands')).toBeVisible() }

test('one row per on-axis component (like ligand rows): fixed-colour marks at the data positions, aligned; tooltips and pinned details from the data', async ({ page }) => {
  await open(page)
  const t = await (await page.request.get('/families/RF00167/modification_tracks.json')).json()
  const cat = page.getByTestId('category-modifications')
  // rows: components with at least one consensus-axis observation, most instances first; off-axis-only components (GTP, GDP, 23G) have no row
  const onAxis = Object.entries(t.components as Record<string, { n_on_consensus_axis: number; n_instances: number }>).filter(([, c]) => c.n_on_consensus_axis > 0)
    .sort((a, b) => b[1].n_instances - a[1].n_instances || a[0].localeCompare(b[0])).map(([k]) => k)
  expect(onAxis).toEqual(['CCC', '50L', '50N', 'IU'])
  await expect(page.getByTestId('category-toggle-modifications')).toContainText('Nonstandard nucleotides'); await expect(page.getByTestId('category-toggle-modifications')).toContainText('(4 tracks)')
  expect(await cat.locator('[data-testid^="track-"]').evaluateAll((rs) => rs.map((r) => r.getAttribute('data-testid')))).toEqual(onAxis.map((c) => `track-mod-${c}`))
  for (const c of ['GTP', 'GDP', '23G']) await expect(page.getByTestId(`track-mod-${c}`)).toHaveCount(0)
  await expect(page.getByTestId('track-mod-CCC')).toContainText('CCC'); await expect(page.getByTestId('track-mod-CCC')).toContainText('2 instances')
  // category order: SASA, loops, nonstandard nucleotides, ligands, pockets
  const ys = await Promise.all(['sasa', 'loops', 'modifications', 'ligands', 'pockets'].map(async (id) => (await page.getByTestId(`category-${id}`).boundingBox())!.y))
  expect(ys).toEqual([...ys].sort((a, b) => a - b))
  // marks per row exactly where the data has that component; one fixed colour; zero positions empty
  const expected: Record<string, number[]> = { CCC: [87], '50L': [84], '50N': [22], IU: [32, 35] }
  for (const [ccd, posList] of Object.entries(expected)) {
    const fromData = t.positions.filter((p: { ccd: { ccd_id: string }[] }) => p.ccd.some((x) => x.ccd_id === ccd)).map((p: { rfam_consensus_pos: number }) => p.rfam_consensus_pos)
    expect(fromData).toEqual(posList)
    const marks = page.getByTestId(`track-mod-${ccd}`).locator('rect[data-pos][data-value]')
    expect(await marks.evaluateAll((rs) => rs.map((r) => Number(r.getAttribute('data-pos'))))).toEqual(posList)
    expect(new Set(await marks.evaluateAll((rs) => rs.map((r) => r.getAttribute('fill'))))).toEqual(new Set(['#c2185b']))
  }
  expect(await cat.locator('rect[data-state="unevaluable"]').count()).toBe(4 * t.positions.filter((p: { denominator: number }) => p.denominator === 0).length)
  const ccc = page.getByTestId('track-mod-CCC')
  expect(await ccc.locator('rect[data-pos="87"]').getAttribute('data-value')).toBe((2 / 77).toFixed(3))
  // same geometry as other tracks; aligned with the SASA mark at the same position
  const geo = await page.evaluate(() => [...document.querySelectorAll('[data-testid^="track-mod-"] rect[data-pos], [data-testid="track-HL"] rect[data-pos]')].map((r) => { const b = r.getBoundingClientRect(); return `${Math.round(b.width)}x${Math.round(b.height)}` }))
  expect(new Set(geo).size).toBe(1)
  const a = (await ccc.locator('rect[data-pos="87"]').boundingBox())!; const ref = (await page.getByTestId('track-sasa-median').locator('rect[data-pos="87"]').boundingBox())!
  expect(Math.abs(a.x - ref.x)).toBeLessThan(1)
  await expect(page.getByTestId('legend-nonstandard-presence')).toHaveText('Nonstandard nucleotide observed')
  // tooltips
  const tip = page.getByTestId('annotation-tooltip')
  await ccc.locator('rect[data-pos="87"]').hover()
  await expect(tip).toHaveText("Rfam residue 87CCC — Cytidine-5'-phosphate-2',3'-cyclic phosphateObserved in 2 of 77 RNA instances (2.6%)")
  for (const [ccd, pos, k, n, pc] of [['50N', 22, 1, 83, '1.2%'], ['IU', 32, 1, 84, '1.2%'], ['IU', 35, 1, 84, '1.2%'], ['50L', 84, 1, 83, '1.2%']] as const) {
    await page.getByTestId(`track-mod-${ccd}`).locator(`rect[data-pos="${pos}"]`).hover()
    await expect(tip).toContainText(`Rfam residue ${pos}${ccd} — `); await expect(tip).toContainText(`Observed in ${k} of ${n} RNA instances (${pc})`)
  }
  // tooltip stays inside the viewport for the bottom-most row (flips above the pointer)
  await page.setViewportSize({ width: 1100, height: 640 })
  const last = page.getByTestId('track-mod-IU').locator('rect[data-pos="35"]')
  await last.evaluate((el) => { const r = el.getBoundingClientRect(); window.scrollBy(0, r.bottom - window.innerHeight + 12) })
  const mb = (await last.boundingBox())!; expect(mb.y + mb.height).toBeGreaterThan(640 - 20)
  await page.mouse.move(mb.x + mb.width / 2, mb.y + mb.height / 2)
  await expect(tip).toBeVisible()
  const box = (await tip.boundingBox())!; const vp = page.viewportSize()!
  expect(box.y + box.height).toBeLessThan(mb.y + mb.height)
  expect(box.x + box.width).toBeLessThanOrEqual(vp.width); expect(box.x).toBeGreaterThanOrEqual(0); expect(box.y).toBeGreaterThanOrEqual(0); expect(box.y + box.height).toBeLessThanOrEqual(vp.height)
  await page.setViewportSize({ width: 1400, height: 1000 })
  // pinned details
  await ccc.locator('rect[data-pos="87"]').click()
  const inspect = page.getByTestId('annotation-inspect')
  await expect(inspect).toContainText('CCC · Rfam residue 87'); await inspect.locator('summary').click()
  await expect(inspect).toContainText("CCC — CYTIDINE-5'-PHOSPHATE-2',3'-CYCLIC PHOSPHATE; parent C (CCD)")
  await expect(inspect).toContainText('Observed in 2 of 77 RNA instances (2.6%) in 2 PDB entries')
  await expect(inspect).toContainText('14 of the 77 RNA instances have no coordinates for this nucleotide')
  await page.getByTestId('track-mod-50L').locator('rect[data-pos="84"]').click()
  await expect(inspect).toContainText('50L — (5S)-2-amino-8-(5-O-phosphono-beta-D-ribofuranosyl)imidazo[1,2-a][1,3,5]triazin-4(8H)-one; parent G (NAKB)')
  await page.getByTestId('annotations-clear').click()
})

test('category popover, collapse, downloads; other categories unchanged', async ({ page }) => {
  await open(page)
  const icon = page.getByTestId('category-header-modifications').getByTestId('info-tip')
  await icon.scrollIntoViewIfNeeded(); await icon.click()
  const pop = page.getByTestId('info-popover')
  await expect(pop).toHaveText('Nonstandard nucleotide components in deposited RNA sequences, mapped to Rfam positions. These are not necessarily biological RNA modifications.')
  await expect(page.getByTestId('category-modifications')).toHaveAttribute('data-open', 'true')
  const b = (await pop.boundingBox())!; const vp = page.viewportSize()!
  expect(b.y + b.height).toBeLessThanOrEqual(vp.height); expect(b.x + b.width).toBeLessThanOrEqual(vp.width)
  await page.keyboard.press('Escape'); await expect(pop).toHaveCount(0)
  await page.getByTestId('category-toggle-modifications').click()
  await expect(page.getByTestId('category-modifications')).toHaveAttribute('data-open', 'false'); await expect(icon).toBeVisible()
  await icon.focus(); await expect(pop).toBeVisible(); await page.keyboard.press('Escape')
  await page.getByTestId('category-toggle-modifications').click()
  // downloads
  await page.getByTestId('download-button').click()
  await expect(page.locator('[data-testid^="download-"][role="menuitem"]')).toHaveCount(9)
  for (const [id, file, label, key] of [['modifications', 'modifications.json', 'Nonstandard nucleotides · Detailed', 'observations'], ['modification-tracks', 'modification_tracks.json', 'Nonstandard nucleotides · Aggregated', 'positions']]) {
    const it = page.getByTestId(`download-${id}`); await expect(it).toContainText(label)
    const href = (await it.getAttribute('href'))!; expect(href.endsWith(`/families/RF00167/${file}`)).toBe(true); await expect(it).toHaveAttribute('download', `RF00167_${file}`)
    const j = await (await page.request.get(href)).json(); expect(j.rfam_acc).toBe('RF00167'); expect(j).toHaveProperty(key); expect(j.provenance).toHaveProperty('pdbe_modified_AA_or_NA')
  }
  const d = await (await page.request.get('/families/RF00167/modifications.json')).json()
  expect(d.observations).toHaveLength(13); expect(d.observations.filter((o: { ccd_id: string; rfam_consensus_pos: number | null }) => o.ccd_id === '23G' && o.rfam_consensus_pos === null)).toHaveLength(2)
  await page.keyboard.press('Escape')
  // existing categories unchanged
  await page.getByTestId('track-HL').locator('rect[data-pos="53"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toHaveText('Rfam residue 53Hairpin loop (HL)Observed in 8 of 80 RNA instances (10%)')
  await expect(page.locator('[data-testid^="track-lig-"]')).toHaveCount(29)
  await expect(page.getByTestId('track-pocket-occurrence')).toBeVisible(); await expect(page.getByTestId('track-sasa-median')).toBeVisible()
  for (const id of ['legend-occurrence-frequency', 'legend-sasa-median', 'legend-no-data', 'legend-no-annotation']) await expect(page.getByTestId(id)).toBeVisible()
})

test('family without nonstandard-nucleotide data: category hidden, download disabled, other categories intact', async ({ page }) => {
  await page.route('**/families/*/modification_tracks.json', (r) => r.fulfill({ status: 404, body: 'not found' }))
  await page.route('**/families/*/modifications.json', (r) => r.fulfill({ status: 404, body: 'not found' }))
  await open(page)
  await expect(page.getByTestId('category-sasa')).toBeVisible()
  await expect(page.getByTestId('category-modifications')).toHaveCount(0)
  await expect(page.locator('section[data-testid^="category-"]')).toHaveCount(4)
  await expect(page.getByTestId('legend-nonstandard-presence')).toHaveCount(0)
  await page.getByTestId('download-button').click()
  await expect(page.getByTestId('download-modifications')).toHaveAttribute('aria-disabled', 'true')
  await expect(page.getByTestId('download-modification-tracks')).toHaveAttribute('aria-disabled', 'true')
  await expect(page.getByTestId('download-tracks')).not.toHaveAttribute('aria-disabled', 'true')
})
