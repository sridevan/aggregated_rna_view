import { expect, test, type Page } from '@playwright/test'
import { selectColumn } from './helpers'

/** Annotations tab (Milestone 3C refinement): PDBe-KB-style fixed-size marks coloured by frequency, scrollable categories, concise tooltips, info popovers, downloads; independent of Family Overview. */
const open = async (page: Page) => { await page.goto('/RF00167'); await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-ligands')).toBeVisible() }
const lum = (rgb: string) => { const m = rgb.match(/\d+/g)!.map(Number); return 0.2126 * m[0] + 0.7152 * m[1] + 0.0722 * m[2] }

test('two categories; HL/IL/JL and one track per default-track CCD; marks have identical size and colour follows frequency', async ({ page }) => {
  await open(page)
  await expect(page.locator('section[data-testid^="category-"]')).toHaveCount(5)   // RNA loop annotations, Ligand binding sites, Predicted binding pockets, Solvent accessibility, Nonstandard nucleotides
  await expect(page.getByTestId('category-toggle-loops')).toContainText('RNA loop annotations')
  await expect(page.getByTestId('category-toggle-ligands')).toContainText('Ligand binding sites')
  await expect(page.getByTestId('annotations-tab')).not.toContainText('RNA 2D structure')
  // information icons only in the two category headers, none on tracks
  await expect(page.getByTestId('annotations-tab').getByTestId('info-tip')).toHaveCount(5)
  await expect(page.locator('[data-testid^="track-"] [data-testid="info-tip"]')).toHaveCount(0)
  await expect(page.getByTestId('annotations-tab').getByTestId('instance-select')).toHaveCount(0)
  await expect(page.getByTestId('annotations-tab').getByTestId('molstar-viewer')).toHaveCount(0)
  for (const id of ['HL', 'IL', 'JL']) await expect(page.getByTestId(`track-${id}`)).toBeVisible()
  const lig = await (await page.request.get('/families/RF00167/ligands.json')).json()
  const ccds = Object.entries(lig.by_ccd as Record<string, { n_entries_in_default_track: number }>).filter(([, v]) => v.n_entries_in_default_track > 0).map(([k]) => k)
  await expect(page.locator('[data-testid^="track-lig-"]')).toHaveCount(ccds.length)
  // identical mark geometry across all tracks (viewBox units: width CELL-2, height MARK_H) and identical rendered heights
  const geo = await page.evaluate(() => [...document.querySelectorAll('[data-testid^="track-"] rect[data-pos]')].map((r) => `${r.getAttribute('width')}x${r.getAttribute('height')}`))
  expect(geo.length).toBeGreaterThan(100); expect(new Set(geo).size).toBe(1)
  const px = await page.evaluate(() => [...document.querySelectorAll('[data-testid="track-HL"] rect[data-pos], [data-testid="track-lig-HPA"] rect[data-pos]')].map((r) => { const b = r.getBoundingClientRect(); return `${Math.round(b.width)}x${Math.round(b.height)}` }))
  expect(new Set(px).size).toBe(1)
  const h = +px[0].split('x')[1]; expect(h).toBeGreaterThanOrEqual(14); expect(h).toBeLessThanOrEqual(18)
  // row height 28-32 px
  const rowH = (await page.getByTestId('track-HL').boundingBox())!.height; expect(rowH).toBeGreaterThanOrEqual(28); expect(rowH).toBeLessThanOrEqual(32)
  // colour intensity follows the existing frequency: consensus 37 (freq 1.0) is darker than a low-frequency HL position
  const tr = await (await page.request.get('/families/RF00167/tracks.json')).json()
  const low = tr.positions.find((p: { hairpin_loops: { frequency: number | null; n_instances: number } }) => p.hairpin_loops.n_instances > 0 && (p.hairpin_loops.frequency ?? 0) < 0.3).rfam_consensus_pos
  const fill = (pos: number) => page.locator(`[data-testid="track-HL"] rect[data-pos="${pos}"]`).getAttribute('fill')
  expect(lum((await fill(37))!)).toBeLessThan(lum((await fill(low))!))
  expect(await page.locator('[data-testid="track-HL"] rect[data-pos="37"]').getAttribute('data-value')).toBe(tr.positions[36].hairpin_loops.frequency.toFixed(3))
  // unevaluable positions (no observable instance) are grey marks, evaluated-but-unannotated positions are empty
  const unevalPos = tr.positions.find((p: { hairpin_loops: { denominators: { n_observable_and_evaluated: number } } }) => p.hairpin_loops.denominators.n_observable_and_evaluated === 0)
  if (unevalPos) await expect(page.locator(`[data-testid="track-HL"] rect[data-pos="${unevalPos.rfam_consensus_pos}"]`)).toHaveAttribute('data-state', 'unevaluable')
  await expect(page.getByTestId('scale-legend')).toContainText('Occurrence frequency (0–100%)')   // one shared scale for loops, ligands, pockets
  // axis aligned with the marks: first tick (position 1) centred on the first cell of a track
  const axisTick = await page.getByTestId('annotation-axis').locator('span', { hasText: /^1$/ }).first().boundingBox()
  const firstMark = await page.locator('[data-testid="track-JL"] rect[data-pos="1"]').boundingBox()
  expect(Math.abs(axisTick!.x + axisTick!.width / 2 - (firstMark!.x + firstMark!.width / 2))).toBeLessThan(3)
})

test('concise tooltips with correct counts (overlapping classes at 53 kept); header popovers do not toggle the category; pinned details keep the evidence', async ({ page }) => {
  await open(page)
  const tr = await (await page.request.get('/families/RF00167/tracks.json')).json(); const lig = await (await page.request.get('/families/RF00167/ligands.json')).json()
  const tip = page.getByTestId('annotation-tooltip')
  const p53 = tr.positions[52]
  expect(p53.hairpin_loops.n_instances).toBe(66); expect(p53.junction_loops_union.n_instances).toBe(76); expect(p53.internal_loops.n_instances).toBe(13)   // validated audit values
  await page.getByTestId('track-HL').locator('rect[data-pos="53"]').hover()
  await expect(tip).toHaveText('Rfam residue 53Hairpin loop (HL)Observed in 66 of 80 RNA instances (82.5%)')
  await page.getByTestId('track-JL').locator('rect[data-pos="53"]').hover()
  await expect(tip).toHaveText('Rfam residue 53Junction loop (JL)Observed in 76 of 80 RNA instances (95%)')
  await page.getByTestId('track-IL').locator('rect[data-pos="53"]').hover()
  await expect(tip).toHaveText('Rfam residue 53Internal loop (IL)Observed in 13 of 80 RNA instances (16.3%)')
  await page.getByTestId('track-HL').locator('rect[data-pos="37"]').hover()
  await expect(tip).toHaveText(`Rfam residue 37Hairpin loop (HL)Observed in ${tr.positions[36].hairpin_loops.denominators.n_positive_observable_and_evaluated} of ${tr.positions[36].hairpin_loops.denominators.n_observable_and_evaluated} RNA instances (${Math.round(1000 * tr.positions[36].hairpin_loops.frequency) / 10}%)`)
  const hpa = new Set(lig.bound_molecules.filter((b: { ccd_id: string; in_default_track: boolean; consensus_positions: number[] }) => b.ccd_id === 'HPA' && b.in_default_track && b.consensus_positions.includes(26)).map((b: { pdb_id: string }) => b.pdb_id))
  await page.getByTestId('track-lig-HPA').locator('rect[data-pos="26"]').hover()
  const e26 = tr.positions[25].ligand_binding_sites.by_ccd_entries.HPA   // denominator: entries with bound HPA and an observed nucleotide here
  expect(e26.n_positive_observable_and_evaluated).toBe(hpa.size)
  await expect(tip).toHaveText(`Rfam residue 26HPA — HypoxanthineRecorded contact in ${hpa.size} of ${e26.n_observable_and_evaluated} observable PDB entries (${Math.round(1000 * hpa.size / e26.n_observable_and_evaluated) / 10}%)`)
  for (const banned of ['HL_', 'J3', 'J4', 'motif', 'PP', '/ ']) expect(await tip.textContent()).not.toContain(banned)
  // header popovers: visible when collapsed, do not toggle, not clipped
  await page.getByTestId('category-header-loops').getByTestId('info-tip').click()
  await expect(page.getByRole('tooltip')).toContainText('BGSU hairpin (HL), internal (IL) and junction (JL) loops mapped to Rfam positions')
  await expect(page.getByRole('tooltip')).toContainText('Annotations may overlap or vary between structures')
  await expect(page.getByTestId('category-loops')).toHaveAttribute('data-open', 'true')   // click on ⓘ did not collapse
  let pop = (await page.getByRole('tooltip').boundingBox())!; const vp = page.viewportSize()!
  expect(pop.x).toBeGreaterThanOrEqual(0); expect(pop.x + pop.width).toBeLessThanOrEqual(vp.width); expect(pop.y + pop.height).toBeLessThanOrEqual(vp.height + 1)
  await page.mouse.move(5, 5); await page.keyboard.press('Escape')
  await page.getByTestId('category-toggle-ligands').click()
  await expect(page.getByTestId('category-ligands')).toHaveAttribute('data-open', 'false')
  await expect(page.getByTestId('category-header-ligands').getByTestId('info-tip')).toBeVisible()   // still visible when collapsed
  await page.getByTestId('category-header-ligands').getByTestId('info-tip').focus()                // keyboard reachable
  await expect(page.getByRole('tooltip')).toContainText('PDBe ligand-contact annotations aggregated by Rfam position and PDB entry')
  await expect(page.getByRole('tooltip')).toContainText('Entity-level annotations may not apply to every RNA chain')
  pop = (await page.getByRole('tooltip').boundingBox())!; expect(pop.x + pop.width).toBeLessThanOrEqual(vp.width)
  await expect(page.getByTestId('category-ligands')).toHaveAttribute('data-open', 'false')
  await page.getByTestId('category-toggle-ligands').click()
  await expect(page.getByTestId('category-ligands')).toHaveAttribute('data-open', 'true')
  // click-to-pin keeps the evidence (motif groups) out of the tooltip but available
  await page.getByTestId('track-HL').locator('rect[data-pos="53"]').click()
  await expect(page.getByTestId('annotation-inspect')).toContainText('Rfam residue 53')
  await page.getByTestId('annotation-inspect').locator('summary').click()
  await expect(page.getByTestId('annotation-inspect')).toContainText('HL_56046.1')
  await page.getByTestId('track-HL').locator('rect[data-pos="53"]').click()
  await expect(page.getByTestId('annotation-inspect')).toHaveCount(0)
})

test('ligand category scrolls independently with a visible header; tooltips are not clipped; collapse works; Download menu serves the detailed JSON', async ({ page }) => {
  await open(page)
  const body = page.getByTestId('category-body-ligands')
  const dims = await body.evaluate((e) => ({ client: e.clientHeight, scroll: e.scrollHeight }))
  expect(dims.client).toBeLessThanOrEqual(230); expect(dims.scroll).toBeGreaterThan(dims.client)
  const small = await page.getByTestId('category-body-loops').evaluate((e) => ({ client: e.clientHeight, scroll: e.scrollHeight }))
  expect(small.scroll).toBe(small.client); expect(small.client).toBeLessThan(120)   // natural height for a 3-track category
  const headerBefore = (await page.getByTestId('category-ligands').getByRole('button').first().boundingBox())!
  const firstRowBefore = (await page.getByTestId('track-lig-HPA').boundingBox())!
  await body.evaluate((e) => { e.scrollTop = 150 })
  await expect.poll(() => body.evaluate((e) => e.scrollTop)).toBe(150)
  const headerAfter = (await page.getByTestId('category-ligands').getByRole('button').first().boundingBox())!
  expect(headerAfter.y).toBe(headerBefore.y)                                        // header outside the scrolling area
  expect((await page.getByTestId('track-lig-HPA').boundingBox())!.y).toBeLessThan(firstRowBefore.y)
  expect(await page.getByTestId('category-body-loops').evaluate((e) => e.scrollTop)).toBe(0)   // other category unaffected
  // a tooltip for a row near the bottom edge of the scroll box is rendered fixed, fully inside the viewport
  const lastVisible = page.locator('[data-testid^="track-lig-"] rect[data-pos="25"]').nth(6)
  await lastVisible.hover()
  const tb = (await page.getByTestId('annotation-tooltip').boundingBox())!; const vp = page.viewportSize()!
  expect(tb.x).toBeGreaterThanOrEqual(0); expect(tb.x + tb.width).toBeLessThanOrEqual(vp.width); expect(tb.height).toBeGreaterThan(20)
  expect(await page.getByTestId('annotation-tooltip').evaluate((e) => getComputedStyle(e).position)).toBe('fixed')
  // collapse / expand
  await page.getByTestId('category-ligands').getByRole('button').first().click()
  await expect(page.getByTestId('category-ligands')).toHaveAttribute('data-open', 'false')
  await expect(page.locator('[data-testid^="track-lig-"]')).toHaveCount(0)
  await page.getByTestId('category-ligands').getByRole('button').first().click()
  await expect(page.locator('[data-testid^="track-lig-"]').first()).toBeVisible()
  // Download dropdown: three datasets for the selected family, served from the existing static files
  await expect(page.getByTestId('download-items')).toHaveCount(0)
  await page.getByTestId('download-button').click()
  await expect(page.getByTestId('download-items')).toBeVisible()
  await expect(page.getByTestId('download-items').getByRole('menuitem')).toHaveCount(9)   // loops, ligands, tracks, pockets x2, SASA x2, nonstandard nucleotides x2
  for (const [id, file, label] of [['loops', 'loops.json', 'RNA loop annotations'], ['ligands', 'ligands.json', 'Ligand binding annotations'], ['tracks', 'tracks.json', 'Aggregated annotation tracks']]) {
    const it = page.getByTestId(`download-${id}`)
    await expect(it).toContainText(label)
    const href = (await it.getAttribute('href'))!
    expect(href.endsWith(`/families/RF00167/${file}`)).toBe(true)
    await expect(it).toHaveAttribute('download', `RF00167_${file}`)
    const r = await page.request.get(href); expect(r.ok()).toBe(true); const j = await r.json(); expect(j.rfam_acc).toBe('RF00167')
  }
  await page.keyboard.press('Escape'); await expect(page.getByTestId('download-items')).toHaveCount(0)
  await page.getByTestId('download-button').click(); await page.getByTestId('download-items').getByRole('menuitem').first().focus(); await page.keyboard.press('ArrowDown')
  await expect(page.getByTestId('download-ligands')).toBeFocused()
  const dl = page.waitForEvent('download'); await page.getByTestId('download-ligands').click(); const d = await dl
  expect(d.suggestedFilename()).toBe('RF00167_ligands.json'); await expect(page.getByTestId('download-items')).toHaveCount(0)   // menu closes after a selection
  const loops = await (await page.request.get('/families/RF00167/loops.json')).json()
  const L = loops.loops[0]
  for (const k of ['loop_id', 'loop_type', 'motif_group', 'nucleotides', 'strands', 'consensus_positions', 'instances', 'pdb_id']) expect(L).toHaveProperty(k)
  expect(loops.counts.by_type).toHaveProperty('J3'); expect(loops.counts.by_type).toHaveProperty('J4')
  expect(loops.loops.some((x: { loop_type: string }) => x.loop_type === 'HL') && loops.loops.some((x: { loop_type: string }) => x.loop_type === 'J3')).toBe(true)   // overlapping classes preserved
  const tracks = await (await page.request.get('/families/RF00167/tracks.json')).json()
  expect(tracks.positions[52].hairpin_loops).toHaveProperty('denominators'); expect(tracks.positions[52].hairpin_loops.n_instances).toBe(66); expect(tracks).toHaveProperty('source_availability')
})

test('narrow screen: no horizontal overflow, tracks still aligned', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 })
  await open(page)
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)
  const a = (await page.locator('[data-testid="track-HL"] rect[data-pos="37"]').boundingBox())!, b = (await page.locator('[data-testid="track-lig-HPA"] rect[data-pos="25"]').boundingBox())!
  const cell = a.width + 2 * (a.width / 8)   // mark = CELL-2 of CELL
  expect(Math.abs((b.x - a.x) - (25 - 37) * cell)).toBeLessThan(1.5)
})

test('clicking a track pins its details locally and never touches the Family Overview selection, instance or Mol* camera', async ({ page }) => {
  await page.goto('/RF00167')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await selectColumn(page, 50)
  await expect(page.getByTestId('sel-highlightable')).toContainText('50 → ')
  type W = { __rnaExplorer: { selection?: { positions: number[]; instance: string }; lastFocus?: number[] } }
  await page.waitForTimeout(600)
  const before = await page.evaluate(() => { const h = (window as unknown as W).__rnaExplorer; return { sel: h.selection, focus: h.lastFocus } })
  await page.getByTestId('tab-annotations').click()
  await expect(page.getByTestId('track-HL')).toBeVisible()
  await page.getByTestId('track-HL').locator('rect[data-pos="37"]').click()
  await expect(page.getByTestId('annotation-inspect')).toContainText('Rfam residue 37')
  await expect(page.getByTestId('tab-annotations')).toHaveAttribute('aria-selected', 'true')
  await page.getByTestId('tab-overview').click()
  await expect(page.getByTestId('sel-highlightable')).toContainText('50 → ')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  const after = await page.evaluate(() => { const h = (window as unknown as W).__rnaExplorer; return { sel: h.selection, focus: h.lastFocus } })
  expect(after.sel).toEqual(before.sel); expect(after.focus).toEqual(before.focus)
})
