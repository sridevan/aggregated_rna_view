import { expect, test } from '@playwright/test'
import { selectColumn } from './helpers'
import { mkdirSync } from 'node:fs'

/** Captures review screenshots of the reorganised tabs (Milestone 3C): Family Overview with the structure explorer, Structures, Annotations tracks. */
const OUT = '../docs/screenshots'
type W = { __rnaExplorer: { molstar?: { plugin: { canvas3d?: { props: { renderer: { backgroundColor: number } } } } } } }
const bgHex = (page: import('@playwright/test').Page) => page.evaluate(() => { const c = (window as unknown as W).__rnaExplorer.molstar?.plugin.canvas3d?.props.renderer.backgroundColor; return c == null ? null : '#' + c.toString(16).padStart(6, '0') })

test('screenshots of the three tabs', async ({ page }) => {
  mkdirSync(OUT, { recursive: true })
  await page.setViewportSize({ width: 1400, height: 1000 })
  await page.goto('/RF00167')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready')
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready')
  await page.waitForTimeout(1500)
  expect(await bgHex(page)).toBe('#ffffff')
  await page.screenshot({ path: `${OUT}/m3c_family_overview_4fe5.png`, fullPage: true })
  await page.getByTestId('instance-select').selectOption('1y26_1_X_A')
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-pdb', '1y26')
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', '1y26')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready')
  await selectColumn(page, 50)
  await expect(page.getByTestId('sel-highlightable')).toContainText('50 → G34')
  await page.waitForTimeout(2000)
  expect(await bgHex(page)).toBe('#ffffff')
  await page.screenshot({ path: `${OUT}/m3c_family_overview_1y26_selection.png`, fullPage: true })
  await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await page.waitForTimeout(1000)
  await page.screenshot({ path: `${OUT}/m3c_structures.png`, fullPage: true })
  await page.getByTestId('tab-annotations').click()
  await expect(page.getByTestId('category-loops')).toBeVisible()
  await page.getByTestId('track-HL').locator('rect[data-pos="53"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toBeVisible()
  await page.screenshot({ path: `${OUT}/m3c_annotations.png`, fullPage: true })
  await page.mouse.move(5, 5)
  await page.getByTestId('track-sasa-median').locator('rect[data-pos="53"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toContainText('Median SASA')
  await page.getByTestId('category-sasa').scrollIntoViewIfNeeded()
  await page.screenshot({ path: `${OUT}/m3d_annotations_sasa.png`, fullPage: true })
  await page.mouse.move(5, 5)
  // predicted pockets are in the Ligands tab
  await page.getByTestId('tab-ligands').click(); await expect(page.getByTestId('category-pockets')).toBeVisible()
  await page.getByTestId('track-pocket-occurrence').locator('rect[data-pos="56"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toContainText('Pocket occurrence: 61 of 80 RNA instances')
  await page.screenshot({ path: `${OUT}/m3e_annotations_pockets.png`, fullPage: true })
  await page.mouse.move(5, 5)
  await page.getByTestId('category-header-pockets').getByTestId('info-tip').click()
  await expect(page.getByTestId('info-popover')).toBeVisible()
  await page.screenshot({ path: `${OUT}/m3e_popover_pockets.png`, fullPage: false })
  await page.keyboard.press('Escape')
  await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-sasa')).toBeVisible()
  await page.getByTestId('category-header-sasa').scrollIntoViewIfNeeded()
  await page.getByTestId('category-header-sasa').getByTestId('info-tip').click()
  await expect(page.getByTestId('info-popover')).toBeVisible()
  await page.screenshot({ path: `${OUT}/m3e_popover_sasa.png`, fullPage: false })
  await page.keyboard.press('Escape')
  // Nonstandard nucleotides (final Annotations tab: all five categories)
  await page.mouse.move(5, 5)
  await page.getByTestId('track-mod-CCC').locator('rect[data-pos="87"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toContainText('CCC — ')
  await page.screenshot({ path: `${OUT}/m3f_annotations_all_categories.png`, fullPage: true })
  await page.screenshot({ path: `${OUT}/m3g_annotations_final.png`, fullPage: true })
  await page.getByTestId('track-mod-CCC').locator('rect[data-pos="87"]').click()
  await page.getByTestId('annotation-inspect').locator('summary').click()
  await page.mouse.move(5, 5)
  await page.screenshot({ path: `${OUT}/m3f_nonstandard_pinned_details.png`, fullPage: true })
  await page.getByTestId('annotations-clear').click()
  await page.getByTestId('category-header-modifications').scrollIntoViewIfNeeded()
  await page.getByTestId('category-header-modifications').getByTestId('info-tip').click()
  await expect(page.getByTestId('info-popover')).toBeVisible()
  await page.screenshot({ path: `${OUT}/m3f_popover_nonstandard.png`, fullPage: false })
  await page.keyboard.press('Escape')
  await page.mouse.move(5, 5)
  await page.getByTestId('download-button').click()
  await expect(page.getByTestId('download-items')).toBeVisible()
  await page.screenshot({ path: `${OUT}/m3c_annotations_download_menu.png`, fullPage: false })
})

test('screenshots of the Structures tab (RF00167 structural clustering)', async ({ page }) => {
  const ready = () => expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await page.goto('/RF00167'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click(); await ready(); await page.waitForTimeout(800)
  await page.screenshot({ path: `${OUT}/m4_structures_default.png`, fullPage: true })
  await page.getByTestId('structures-section-b').screenshot({ path: `${OUT}/m4_structures_panels_aligned.png` })
  await page.getByTestId('table-mode-toggle').click(); await page.waitForTimeout(300)
  await page.screenshot({ path: `${OUT}/m4_structures_expanded.png`, fullPage: true })
  await page.getByTestId('table-mode-toggle').click()
  await page.getByTestId('group-chip-2').click(); await page.getByTestId('drilldown-toggle').click(); await ready(); await page.waitForTimeout(800)
  await page.screenshot({ path: `${OUT}/m4_structures_group2_drilldown.png`, fullPage: true })
  await page.getByTestId('alignment-panel').screenshot({ path: `${OUT}/m4_structures_3d_alignment.png` })
  await page.getByTestId('structures-section-b').screenshot({ path: `${OUT}/m4_structures_panels_aligned_drilldown.png` })
  await page.getByTestId('metric-tm').click(); await page.getByTestId('group-chip-2').click(); await page.getByTestId('drilldown-toggle').click(); await ready(); await page.waitForTimeout(800)
  await page.screenshot({ path: `${OUT}/m4_structures_tm.png`, fullPage: true })
})

test('screenshots: group colours in 3D and representative comparison paths', async ({ page }) => {
  const ready = () => expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  const shot = async (name: string) => { await ready(); await page.waitForTimeout(700); await page.getByTestId('structures-section-b').screenshot({ path: `${OUT}/${name}` }) }
  await page.goto('/RF00167'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click()
  await shot('m4_compare_g1_rep_alone.png')
  await page.getByTestId('group-chip-2').click(); await ready()
  await shot('m4_compare_g2_rep_alone.png')
  await page.getByTestId('group-chip-1').click(); await page.getByTestId('compare-on').click()
  await shot('m4_compare_g1_vs_g2.png')
  await page.getByTestId('compare-off').click(); await page.getByTestId('metric-tm').click()
  await page.getByTestId('group-chip-2').click(); await page.getByTestId('compare-select').selectOption('3'); await page.getByTestId('compare-on').click()
  await shot('m4_compare_tm_g2_vs_g3.png')
})

test('screenshots: associated features (both metrics)', async ({ page }) => {
  await page.goto('/RF00167'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click(); await expect(page.getByTestId('associated-features')).toBeVisible()
  const shot = async (name: string) => { await page.waitForTimeout(300); await page.getByTestId('associated-features').screenshot({ path: `${OUT}/${name}` }) }
  await shot('m4_features_collapsed.png')
  await page.getByTestId('features-toggle').click()
  await shot('m4_features_D_group1.png')
  await page.getByTestId('group-chip-2').click(); await shot('m4_features_D_group2.png')
  await page.getByTestId('group-chip-1').click(); await page.getByTestId('compare-select').selectOption('2'); await page.getByTestId('compare-on').click(); await shot('m4_features_D_compare_g1_g2.png')
  await page.getByTestId('compare-off').click(); await page.getByTestId('metric-tm').click()
  await page.getByTestId('group-chip-2').click(); await page.getByTestId('compare-select').selectOption('4'); await page.getByTestId('compare-on').click(); await shot('m4_features_TM_compare_g2_g4.png')
})

test('screenshots: Family Overview residue focus', async ({ page }) => {
  await page.goto('/RF00167')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }); await page.waitForTimeout(800)
  await page.getByTestId('instance-toolbar').locator('..').screenshot({ path: `${OUT}/m5_overview_default.png` })
  await selectColumn(page, 76)
  await expect.poll(() => page.evaluate(() => (window as unknown as { __rnaExplorer: { focusRepresentation?: number[] } }).__rnaExplorer.focusRepresentation)).toEqual([56])
  await page.waitForTimeout(900)
  await page.getByTestId('instance-toolbar').locator('..').screenshot({ path: `${OUT}/m5_overview_residue_focus.png` })
})

test('screenshots: Family Overview tracks (Phase 1)', async ({ page }) => {
  const ready = () => expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await page.goto('/RF00167'); await ready(); await page.waitForTimeout(800)
  const explorer = page.getByTestId('instance-toolbar').locator('..')
  await page.getByTestId('consensus-track-area').screenshot({ path: `${OUT}/m6_tracks_default.png` })
  await page.getByTestId('reference-toggle').click(); await page.waitForTimeout(300)
  await page.getByTestId('consensus-track-area').screenshot({ path: `${OUT}/m6_tracks_reference_expanded.png` })
  await page.getByTestId('reference-toggle').click()
  await selectColumn(page, 76)
  await expect.poll(() => page.evaluate(() => (window as unknown as { __rnaExplorer: { focusRepresentation?: number[] } }).__rnaExplorer.focusRepresentation)).toEqual([56]); await page.waitForTimeout(900)
  await explorer.screenshot({ path: `${OUT}/m6_overview_selected_position.png` })
  await page.setViewportSize({ width: 600, height: 1000 }); await page.waitForTimeout(500)
  await page.getByTestId('consensus-track-area').screenshot({ path: `${OUT}/m6_tracks_narrow.png` })
})

test('screenshots: Family Overview loops (Phase 2)', async ({ page }) => {
  const ready = () => expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  const explorer = page.getByTestId('instance-toolbar').locator('..')
  const settle = async (ids: number) => { await expect.poll(() => page.evaluate(() => ((window as unknown as { __rnaExplorer: { focusRepresentation?: number[] } }).__rnaExplorer.focusRepresentation ?? []).length)).toBe(ids); await page.waitForTimeout(900) }
  await page.goto('/RF00167'); await ready(); await page.waitForTimeout(800)
  await page.getByTestId('consensus-track-area').screenshot({ path: `${OUT}/m7_default_sections_collapsed.png` })
  await page.getByTestId('loops-toggle').click(); await expect(page.getByTestId('loop-tracks')).toBeVisible()
  await page.getByTestId('consensus-track-area').screenshot({ path: `${OUT}/m7_loop_tracks_expanded.png` })
  await page.getByTestId('loop-block-HL_4FE5_001-0').click(); await settle(9)
  await explorer.screenshot({ path: `${OUT}/m7_selected_hairpin.png` })
  await page.getByTestId('loop-block-IL_4FE5_001-1').click(); await settle(6)
  await explorer.screenshot({ path: `${OUT}/m7_selected_internal_loop.png` })
  await page.getByTestId('molstar-viewer').screenshot({ path: `${OUT}/m7_3d_focus_complete_loop.png` })
  await page.getByTestId('loop-block-J3_4FE5_001-0').click(); await settle(19)
  await explorer.screenshot({ path: `${OUT}/m7_selected_junction_loop.png` })
  await page.getByTestId('loop-block-IL_4FE5_001-1').hover(); await page.waitForTimeout(200)
  await page.getByTestId('loop-tracks').screenshot({ path: `${OUT}/m7_overlap_J3_J4_lanes.png` })
  const rna = page.getByTestId('rna2d-viewer')
  await rna.locator('[id^="rnaTopologyZoomIn-"]').first().click(); await rna.locator('[id^="rnaTopologyZoomIn-"]').first().click(); await page.waitForTimeout(500)   // the viewer's own zoom (wheel zoom is disabled)
  await page.getByTestId('loop-block-IL_4FE5_001-0').click(); await settle(6)
  await rna.screenshot({ path: `${OUT}/m7_2d_loop_zoom_kept.png` })
})

test('screenshots: Family Overview loop and viewer refinements', async ({ page }) => {
  const ready = () => expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  const explorer = page.getByTestId('instance-toolbar').locator('..')
  const settle = async (ids: number) => { await expect.poll(() => page.evaluate(() => ((window as unknown as { __rnaExplorer: { focusRepresentation?: number[] } }).__rnaExplorer.focusRepresentation ?? []).length)).toBe(ids); await page.waitForTimeout(900) }
  /** the loop tracks plus the room around them where the hover card appears */
  const tracksClip = async () => { const b = (await page.getByTestId('consensus-track-area').boundingBox())!; return { x: b.x, y: Math.max(0, b.y - 10), width: b.width, height: b.height + 120 } }
  await page.goto('/RF00167'); await ready(); await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }); await page.waitForTimeout(800)
  await page.getByTestId('loops-toggle').click(); await expect(page.getByTestId('loop-tracks')).toBeVisible()
  // 1. selected discontinuous internal loop: simplified Mol* (cartoon context, loop nucleotides only in ball-and-stick, no water, no contact lines)
  await page.getByTestId('loop-block-IL_4FE5_001-0').click(); await settle(6)
  await page.getByTestId('molstar-viewer').screenshot({ path: `${OUT}/m8_3d_simplified_internal_loop.png` })
  await explorer.screenshot({ path: `${OUT}/m8_selected_internal_loop.png` })
  // 2. structured tooltip (motif group with BGSU annotation) and 3. red outline on both segments of the selected loop
  await page.getByTestId('loop-block-HL_4FE5_001-0').hover(); await expect(page.getByTestId('loop-tooltip')).toBeVisible(); await page.waitForTimeout(150)
  await page.screenshot({ path: `${OUT}/m8_loop_tooltip.png`, clip: await tracksClip() })
  await page.mouse.move(5, 5); await page.waitForTimeout(150)
  await page.getByTestId('loop-tracks').screenshot({ path: `${OUT}/m8_red_outline_internal_loop.png` })
  // hairpin and junction selections
  await page.getByTestId('loop-block-HL_4FE5_001-0').click(); await settle(9)
  await explorer.screenshot({ path: `${OUT}/m8_selected_hairpin.png` })
  await page.getByTestId('loop-block-J3_4FE5_001-0').click(); await settle(19)
  await explorer.screenshot({ path: `${OUT}/m8_selected_junction_loop.png` })
  await page.getByTestId('molstar-viewer').screenshot({ path: `${OUT}/m8_3d_simplified_junction_loop.png` })
  // J3 selected (red, three segments), the internal loop sharing two of its nucleotides hovered (dashed) with its tooltip
  await page.getByTestId('loop-block-IL_4FE5_001-1').hover(); await expect(page.getByTestId('loop-tooltip')).toBeVisible(); await page.waitForTimeout(150)
  await page.screenshot({ path: `${OUT}/m8_overlap_selected_vs_hover.png`, clip: await tracksClip() })
  // the same loop clicked again: selection cleared, default 2D and 3D views restored
  await page.mouse.move(5, 5)
  await page.getByTestId('loop-block-J3_4FE5_001-1').click(); await page.mouse.move(5, 5); await settle(0)
  await explorer.screenshot({ path: `${OUT}/m8_deselected_default_views.png` })
})

test('screenshots: Family Overview consistent RNA colour', async ({ page }) => {
  const ready = async (pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }); await page.waitForTimeout(900) }
  const explorer = page.getByTestId('instance-toolbar').locator('..')
  await page.goto('/RF00167'); await ready('4fe5')
  await page.getByTestId('instance-select').selectOption('1y27_1_X_A'); await ready('1y27')
  await page.getByTestId('molstar-viewer').screenshot({ path: `${OUT}/m9_colour_1y27_X_default.png` })
  await page.getByTestId('instance-select').selectOption('2eet_1_A_A'); await ready('2eet')
  await page.getByTestId('molstar-viewer').screenshot({ path: `${OUT}/m9_colour_2eet_A_default.png` })
  await page.getByTestId('mode-complete').click(); await page.waitForTimeout(900)
  await page.getByTestId('molstar-viewer').screenshot({ path: `${OUT}/m9_colour_2eet_A_complete.png` })
  await page.getByTestId('mode-region').click(); await page.waitForTimeout(600)
  await page.getByTestId('loops-toggle').click()
  await page.locator('[data-testid^="loop-block-J3_2EET_001-"]').first().click(); await page.mouse.move(5, 5)
  await expect.poll(() => page.evaluate(() => ((window as unknown as { __rnaExplorer: { focusRepresentation?: number[] } }).__rnaExplorer.focusRepresentation ?? []).length)).toBeGreaterThan(10); await page.waitForTimeout(900)
  await page.getByTestId('molstar-viewer').screenshot({ path: `${OUT}/m9_colour_2eet_A_selected_junction.png` })
  await explorer.screenshot({ path: `${OUT}/m9_colour_2eet_A_selected_junction_overview.png` })
})

test('screenshots: Family Overview Phase 3 insertion markers', async ({ page, browser }) => {
  const ready = async (pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }); await page.waitForTimeout(500) }
  const area = page.getByTestId('consensus-track-area')
  await page.goto('/RF00167'); await ready('4fe5')
  await area.screenshot({ path: `${OUT}/m10_no_insertions_4fe5_B.png` })
  await page.getByTestId('instance-select').selectOption('9iwf_1_A_A'); await ready('9iwf')
  await area.screenshot({ path: `${OUT}/m10_insertion_marker_9iwf_A.png` })
  const ab = (await area.boundingBox())!
  await page.getByTestId('insertion-marker-77').hover(); await expect(page.getByTestId('insertion-tooltip')).toBeVisible(); await page.waitForTimeout(150)
  await page.screenshot({ path: `${OUT}/m10_insertion_tooltip.png`, clip: { x: ab.x + ab.width * 0.45, y: ab.y, width: ab.width * 0.55, height: ab.height + 70 } })
  await page.getByTestId('insertion-marker-77').click(); await expect(page.getByTestId('insertion-details')).toBeVisible(); await page.mouse.move(5, 5); await page.waitForTimeout(150)
  const db = (await page.getByTestId('insertion-details').boundingBox())!
  await page.screenshot({ path: `${OUT}/m10_insertion_details.png`, clip: { x: ab.x + ab.width * 0.4, y: ab.y, width: ab.width * 0.6, height: db.y + db.height - ab.y + 10 } })
  await page.getByTestId('insertion-details-close').click()
  // loops open: blocks, coverage and marker share the grid
  await page.getByTestId('loops-toggle').click(); await expect(page.getByTestId('loop-tracks')).toBeVisible(); await page.waitForTimeout(200)
  await area.screenshot({ path: `${OUT}/m10_alignment_with_loops_9iwf_A.png` })
  // close-up of the marker on the 77|78 boundary, rendered at 4x pixel density (same CSS layout)
  const ctx = await browser.newContext({ viewport: { width: 1400, height: 1000 }, deviceScaleFactor: 4, baseURL: 'http://127.0.0.1:5173' })
  const p2 = await ctx.newPage(); await p2.goto('/RF00167')
  await p2.getByTestId('instance-select').selectOption('9iwf_1_A_A'); await expect(p2.getByTestId('insertion-marker-77')).toBeVisible({ timeout: 60_000 })
  const mk = (await p2.getByTestId('insertion-marker-77').boundingBox())!, a2 = (await p2.getByTestId('consensus-track-area').boundingBox())!
  await p2.screenshot({ path: `${OUT}/m10_marker_closeup.png`, clip: { x: mk.x - 60, y: a2.y, width: 120, height: 66 } })
  await ctx.close()
})

test('screenshots: Family Overview tooltip refinement', async ({ page }) => {
  await page.goto('/RF00167'); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  const area = page.getByTestId('consensus-track-area')
  const clip = async (extra = 90) => { const b = (await area.boundingBox())!; return { x: b.x, y: b.y, width: b.width, height: b.height + extra } }
  const shot = async (name: string, extra?: number) => { await page.waitForTimeout(150); await page.screenshot({ path: `${OUT}/${name}`, clip: await clip(extra) }); await page.mouse.move(5, 5); await page.waitForTimeout(100) }
  const cell = (row: string, pos: number, strip = 'consensus-tracks') => page.getByTestId(strip).locator(`[data-row="${row}"][data-tip-pos="${pos}"]`)
  await cell('seq', 41).hover(); await shot('m11_tooltip_sequence_coverage.png')
  await cell('str', 40).hover(); await shot('m11_tooltip_structural_coverage.png')
  await page.getByTestId('instance-select').selectOption('2xnw_1_A_A'); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', '2xnw')
  await cell('chain', 19).hover(); await shot('m11_tooltip_selected_chain_unavailable.png')
  await page.getByTestId('instance-track-info').hover(); await shot('m11_tooltip_selected_chain_info.png')
  await page.getByRole('button', { name: 'About the RNA loop annotations' }).hover(); await shot('m11_tooltip_loop_annotations_info.png')
  await page.getByRole('button', { name: 'About the reference and alignment details' }).hover(); await shot('m11_tooltip_reference_info.png')
})

test('screenshots: Family Overview simplified controls and Viewer details', async ({ page }) => {
  await page.goto('/RF00167'); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }); await page.waitForTimeout(800)
  const explorer = page.getByTestId('instance-toolbar').locator('..')
  await explorer.screenshot({ path: `${OUT}/m12_default_no_selection.png` })
  await page.getByTestId('loops-toggle').click(); await page.getByTestId('loop-block-IL_4FE5_001-1').click(); await page.mouse.move(5, 5)
  await expect.poll(() => page.evaluate(() => ((window as unknown as { __rnaExplorer: { focusRepresentation?: number[] } }).__rnaExplorer.focusRepresentation ?? []).length)).toBe(6); await page.waitForTimeout(900)
  await explorer.screenshot({ path: `${OUT}/m12_active_loop_selection.png` })
  await page.getByTestId('viewer-details-toggle').click(); await page.getByTestId('viewer-details').scrollIntoViewIfNeeded(); await page.waitForTimeout(200)
  const v = (await page.getByTestId('molstar-viewer').boundingBox())!, d = (await page.getByTestId('viewer-details').boundingBox())!
  await page.screenshot({ path: `${OUT}/m12_viewer_details_expanded.png`, clip: { x: d.x, y: v.y + v.height - 120, width: d.width, height: d.y + d.height - (v.y + v.height - 120) + 6 } })
})

test('screenshots: Family Overview viewers without footer; Viewer details with 3D alignment', async ({ page }) => {
  const ready = async (pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }); await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }); await page.waitForTimeout(800) }
  await page.goto('/RF00167'); await ready('4fe5')
  const v = (await page.getByTestId('rna2d-viewer').boundingBox())!
  await page.getByTestId('viewer-details').scrollIntoViewIfNeeded()
  const vv = (await page.getByTestId('rna2d-viewer').boundingBox())!
  const d0 = (await page.getByTestId('viewer-details').boundingBox())!
  await page.screenshot({ path: `${OUT}/m13_viewers_without_footer.png`, clip: { x: vv.x, y: vv.y, width: d0.x + d0.width - vv.x, height: d0.y + d0.height - vv.y + 4 } })
  await page.getByTestId('viewer-details-toggle').click(); await page.waitForTimeout(200)
  const det = page.getByTestId('viewer-details')
  await det.screenshot({ path: `${OUT}/m13_viewer_details_4fe5_B.png` })
  await page.getByTestId('instance-select').selectOption('1y27_1_X_A'); await ready('1y27')
  await det.screenshot({ path: `${OUT}/m13_viewer_details_1y27_X.png` })
})

test('screenshots: Structures tab follows the selected instance into 3D', async ({ page }) => {
  const d = await (await page.request.get('/families/RF00167/structure_clustering.json')).json()
  const g1 = d.metrics.discrepancy.clustering.groups[0], a = g1.members.find((m: string) => m !== g1.representative)
  await page.goto('/RF00167'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 })
  await page.getByTestId('instance-search').fill(a.split('_')[0]); await page.getByTestId(`row-${a}`).locator('[data-testid="cell-title"]').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-structures', a); await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 })
  await page.waitForTimeout(800)
  await page.getByTestId('alignment-panel').screenshot({ path: `${OUT}/m14_structures_selected_instance_in_3d.png` })
})

test('screenshots: Structures table opens with the active group representative on top', async ({ page }) => {
  await page.goto('/RF00167'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 }); await page.waitForTimeout(500)
  await page.getByTestId('instance-table').screenshot({ path: `${OUT}/m15_table_default_group_order.png` })
  await page.getByTestId('group-chip-2').click(); await page.waitForTimeout(800)
  await page.getByTestId('instance-table').screenshot({ path: `${OUT}/m15_table_group2_order.png` })
})

test('screenshots: compact Download JSON menu with an info tooltip', async ({ page }) => {
  await page.goto('/RF00167'); await page.getByTestId('tab-annotations').click()
  await page.getByTestId('download-button').click(); await expect(page.getByTestId('download-items')).toBeVisible()
  const card = (await page.getByTestId('download-menu').boundingBox())!
  await page.screenshot({ path: `${OUT}/m16_download_menu.png`, clip: { x: card.x + card.width - 520, y: card.y - 10, width: 540, height: 360 } })
  await page.getByTestId('download-info-sasa').getByTestId('info-tip').hover(); await expect(page.getByTestId('info-popover')).toBeVisible(); await page.waitForTimeout(150)
  const vw = page.viewportSize()!.width
  await page.screenshot({ path: `${OUT}/m16_download_menu_tooltip.png`, clip: { x: card.x + card.width - 520, y: card.y - 10, width: Math.min(700, vw - (card.x + card.width - 520)), height: 360 } })
})

test('screenshots: Sequence and annotations section (default, zoomed, collapsed)', async ({ page }) => {
  await page.goto('/RF00167'); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await page.getByTestId('loops-toggle').click(); await page.getByTestId('reference-toggle').click(); await expect(page.getByTestId('loop-tracks')).toBeVisible(); await page.waitForTimeout(300)
  const area = page.getByTestId('consensus-track-area')
  await area.screenshot({ path: `${OUT}/m17_rf00167_default.png` })
  for (let i = 0; i < 3; i++) await page.getByTestId('zoom-in').click()
  await page.getByTestId('loop-block-IL_4FE5_001-0').click(); await page.mouse.move(5, 5); await page.waitForTimeout(400)
  await area.screenshot({ path: `${OUT}/m17_rf00167_zoomed.png` })
  await page.getByTestId('tracks-toggle').click(); await page.waitForTimeout(200)
  const explorer = page.getByTestId('instance-toolbar').locator('..'), b = (await explorer.boundingBox())!
  await page.screenshot({ path: `${OUT}/m17_rf00167_collapsed.png`, clip: { x: b.x, y: b.y, width: b.width, height: 360 } })
})
