import { expect, test, type Page } from '@playwright/test'
import { selectColumn } from './helpers'

/* Family Overview controls: no manual consensus-position input; selection only through tracks, 2D, 3D, loops and insertion markers;
   no idle instruction text; Clear selection disabled without a selection; collapsed "Viewer details" below the viewers. */
type H = { __rnaExplorer: { viewerIds: number[]; rna2dSelected?: number[]; lastMolstarSelection?: { labelSeqIds: number[] }; selection: { active: { kind: string; pos?: number; id?: string } | null } } }
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }); await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }) }
const hook = (page: Page) => page.evaluate(() => { const h = (window as unknown as H).__rnaExplorer; return { active: h.selection.active, v: h.viewerIds, d2: h.rna2dSelected ?? [], mol: h.lastMolstarSelection?.labelSeqIds ?? [] } })
const status = (page: Page) => page.getByTestId('selection-status')

test('toolbar: instance selector and Clear only; no idle text; selection from tracks, 2D and 3D; Clear works and is disabled when idle', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  const bar = page.getByTestId('instance-toolbar')
  await expect(page.getByTestId('goto-input')).toHaveCount(0); await expect(page.getByTestId('goto-button')).toHaveCount(0)
  await expect(bar).not.toContainText('Consensus position'); await expect(bar.getByRole('button', { name: 'Select', exact: true })).toHaveCount(0)
  await expect(status(page)).toHaveCount(0); await expect(page.getByText('No nucleotide selected')).toHaveCount(0)
  await expect(page.getByTestId('clear-selection')).toBeDisabled()
  // tracks: sequence coverage row, structural coverage row and selected-chain row all select the column
  await selectColumn(page, 76)
  await expect.poll(async () => (await hook(page)).active).toEqual({ kind: 'position', pos: 76 })
  await expect(status(page)).toHaveText('Rfam position 76 → 4FE5 B C70 · highlighted and centred in 2D/3D')   // Rfam position + author residue of the selected chain
  await expect(page.getByTestId('clear-selection')).toBeEnabled()
  await page.getByTestId('consensus-tracks').locator('[data-row="str"][data-tip-pos="50"]').click()
  await expect.poll(async () => (await hook(page)).active).toEqual({ kind: 'position', pos: 50 })
  await page.getByTestId('consensus-tracks').locator('[data-row="chain"][data-tip-pos="40"]').click()
  await expect.poll(async () => (await hook(page)).active).toEqual({ kind: 'position', pos: 40 })
  // 2D click (consensus 76 = label_seq_id 56 of 4FE5 B)
  await page.locator('[data-testid="rna2d-viewer"] text.rnaview_4fe5_56').dispatchEvent('click')
  await expect.poll(async () => (await hook(page)).active).toEqual({ kind: 'position', pos: 76 }); await expect.poll(async () => (await hook(page)).mol).toEqual([56])
  // 3D click (pdbe-molstar residue event): label_seq_id 37 -> its consensus position
  await page.evaluate(() => { const ev = new MouseEvent('PDB.molstar.click', { bubbles: true }) as MouseEvent & { eventData?: unknown }; ev.eventData = { entry_id: '4fe5', entity_id: '1', auth_asym_id: 'B', seq_id: 37, auth_seq_id: 57, ins_code: '' }; document.querySelector('[data-testid="molstar-viewer"]')!.dispatchEvent(ev) })
  await expect.poll(async () => (await hook(page)).d2).toEqual([37]); expect((await hook(page)).active?.kind).toBe('position')
  // Clear: nothing selected, no status line, button disabled again
  await page.getByTestId('clear-selection').click()
  await expect.poll(async () => (await hook(page)).active).toBeNull(); await expect.poll(async () => (await hook(page)).d2).toEqual([])
  await expect(status(page)).toHaveCount(0); await expect(page.getByTestId('clear-selection')).toBeDisabled()
  // a column clicked twice toggles off without leaving a status line
  await selectColumn(page, 30); await expect(status(page)).toBeVisible(); await selectColumn(page, 30); await expect(status(page)).toHaveCount(0)
})

test('loops and insertion markers still select and summarise; the loop summary keeps ID, type and membership', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  await page.getByTestId('loops-toggle').click(); await page.getByTestId('loop-block-IL_4FE5_001-1').click()
  await expect(status(page)).toHaveText('Internal loop · IL_4FE5_001 · 2 segments · 6 nucleotides')
  await expect.poll(async () => (await hook(page)).v.length).toBe(6)
  await page.getByTestId('loop-block-IL_4FE5_001-0').click(); await expect(status(page)).toHaveCount(0)
  await page.getByTestId('instance-select').selectOption('9iwf_1_A_A'); await ready(page, '9iwf')
  await page.getByTestId('insertion-marker-77').click(); await page.getByTestId('insertion-view-77:59').click()
  await expect(status(page)).toContainText('Insertion +2'); await expect.poll(async () => (await hook(page)).v).toEqual([59, 60])
  await page.getByTestId('clear-selection').click(); await expect(status(page)).toHaveCount(0)
})

test('Viewer details: collapsed below the viewers; molecular-viewer text without internal identifiers; 3D alignment per instance; no footer under Mol*', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  const det = page.getByTestId('viewer-details'), mol = page.getByTestId('molstar-viewer')
  // no alignment footer inside or below the 3D viewer; Mol*'s orientation axes stay on
  await expect(page.getByTestId('orientation-badge')).toHaveCount(0)
  await expect(mol).not.toContainText('RMSD'); await expect(mol).not.toContainText('Reference frame')
  expect(await page.evaluate(() => (window as unknown as { __rnaExplorer: { molstar: { plugin: { canvas3d?: { props: { camera: { helper: { axes: { name: string } } } } } } } } }).__rnaExplorer.molstar.plugin.canvas3d!.props.camera.helper.axes.name)).toBe('on')
  // collapsed, directly below the viewers
  await expect(det).not.toHaveAttribute('open', ''); await expect(page.getByTestId('viewer-details-body')).toBeHidden()
  const vb = (await mol.boundingBox())!, db = (await det.boundingBox())!
  expect(db.y).toBeGreaterThan(vb.y + vb.height - 1); expect(db.y - (vb.y + vb.height)).toBeLessThan(20); expect(db.height).toBeLessThan(32)
  await page.getByTestId('viewer-details-toggle').click(); await expect(page.getByTestId('viewer-details-body')).toBeVisible()
  const viewers = page.getByTestId('viewer-details-viewers')
  await expect(viewers).toHaveText('Molecular viewers2D viewer: Displays the RNA secondary structure and highlights selected nucleotides or loops.3D viewer: Displays the Rfam-mapped RNA region by default. Select Complete structure to view the entire deposited structure.Selecting a nucleotide or loop highlights it in both viewers and focuses it in 3D. Click the selection again, or use Clear selection, to restore the default views.')
  for (const t of ['label_asym_id', 'label_seq_id', 'author chain', 'water', 'contact line']) await expect(det).not.toContainText(t)
  // 3D alignment: dynamic per instance (orientation.json)
  const o = await (await page.request.get('/families/RF00167/orientation.json')).json()
  await expect(page.getByTestId('alignment-reference-structure')).toHaveText('4FE5 B (this structure)')
  await expect(page.getByTestId('alignment-rmsd')).toHaveText(`${o.transforms['4fe5_1_B_A'].rmsd.toFixed(2)} Å`)
  await expect(page.getByTestId('alignment-alignment-quality')).toHaveText('good')
  await page.getByTestId('instance-select').selectOption('1y27_1_X_A'); await ready(page, '1y27')
  const t = o.transforms['1y27_1_X_A']
  await expect(page.getByTestId('alignment-reference-structure')).toHaveText('4FE5 B'); await expect(page.getByTestId('alignment-superposition')).toHaveText('Superposed on the reference')
  await expect(page.getByTestId('alignment-rmsd')).toHaveText(`${t.rmsd.toFixed(2)} Å`); await expect(page.getByTestId('alignment-matched-atom-pairs')).toHaveText(new RegExp(`^${t.n_pairs} \\(`))
  await expect(page.getByTestId('alignment-alignment-quality')).toHaveText(t.fit_quality)
  await expect(page.getByTestId('viewer-details-body')).toBeVisible()                                         // stays open across instance changes
  await page.getByTestId('viewer-details-toggle').click(); await expect(page.getByTestId('viewer-details-body')).toBeHidden()
  // selection, focus and reset unchanged
  await selectColumn(page, 50); await expect.poll(async () => (await hook(page)).mol.length).toBe(1)
  await page.getByTestId('reset-orientation').click(); await page.getByTestId('clear-selection').click(); await expect.poll(async () => (await hook(page)).mol).toEqual([])
})

test('layout: instance selector and Clear in one row on desktop, wrapping without overflow on narrow screens', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  const sel = (await page.getByTestId('instance-select').boundingBox())!, clr = (await page.getByTestId('clear-selection').boundingBox())!
  expect(Math.abs((sel.y + sel.height / 2) - (clr.y + clr.height / 2))).toBeLessThan(4); expect(clr.x - (sel.x + sel.width)).toBeLessThan(40)   // adjacent, no gap left by the removed input
  await selectColumn(page, 76); await expect(status(page)).toBeVisible()
  for (const w of [700, 420]) {
    await page.setViewportSize({ width: w, height: 900 }); await page.waitForTimeout(200)
    expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)
    const tb = (await page.getByTestId('instance-toolbar').boundingBox())!, c = (await page.getByTestId('clear-selection').boundingBox())!
    expect(c.x + c.width).toBeLessThanOrEqual(tb.x + tb.width + 1)
  }
})

test('3D controls: Clear selection replaces Focus selection; clears the nucleotide or loop and restores the default 2D and 3D views', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  type C = { __rnaExplorer: { getCamera: () => { position: number[] } | null; selection: { viewReset: number } } }
  const cam = () => page.evaluate(() => (window as unknown as C).__rnaExplorer.getCamera()!.position)
  const reset = () => page.evaluate(() => (window as unknown as C).__rnaExplorer.selection.viewReset)
  const dist = (a: number[], b: number[]) => Math.hypot(...a.map((x, i) => x - b[i]))
  const view2d = () => page.evaluate(() => { const s = document.querySelector('[data-testid="rna2d-viewer"] svg.rnaTopoSvg')!; return [...s.querySelectorAll('g[transform]')].map((g) => g.getAttribute('transform')).join('|') })
  await page.waitForTimeout(500); const home = await cam()
  const btn = page.getByTestId('molstar-clear-selection')
  await expect(page.getByTestId('focus-selection')).toHaveCount(0)
  await expect(btn).toHaveText('Clear selection'); await expect(btn).toBeDisabled()
  // a whole loop: focused close-up, user zoom in 2D, then Clear in the 3D controls
  await page.getByTestId('loops-toggle').click(); await page.getByTestId('loop-block-J3_4FE5_001-0').click()
  await expect(btn).toBeEnabled()
  await expect.poll(async () => dist(await cam(), home)).toBeGreaterThan(5)
  await page.locator('[data-testid="rna2d-viewer"] [id^="rnaTopologyZoomIn-"]').first().click(); await expect.poll(view2d).toMatch(/scale\(1\.2/)
  const r0 = await reset()
  await btn.click()
  await expect.poll(async () => (await hook(page)).active).toBeNull()
  await expect.poll(async () => (await hook(page)).mol).toEqual([]); await expect.poll(async () => (await hook(page)).d2).toEqual([])
  expect(await reset()).toBe(r0 + 1)
  await expect.poll(async () => dist(await cam(), home)).toBeLessThan(0.5)                                  // default 3D view
  await expect.poll(view2d).not.toMatch(/scale\(1\.2/)                                                      // default 2D view
  await expect(btn).toBeDisabled(); await expect(page.getByTestId('clear-selection')).toBeDisabled(); await expect(page.getByTestId('selection-status')).toHaveCount(0)
  // a single nucleotide too
  await selectColumn(page, 76); await expect(btn).toBeEnabled(); await btn.click(); await expect.poll(async () => (await hook(page)).active).toBeNull()
})

test('3D viewer controls sit in a toolbar above the canvas, never over it (wide and narrow)', async ({ page }) => {
  for (const width of [1400, 700]) {
    await page.setViewportSize({ width, height: 1000 })
    await page.goto('/RF00167'); await ready(page, '4fe5')
    const viewer = page.getByTestId('molstar-viewer'), bar = page.getByTestId('molstar-toolbar')
    for (const id of ['mode-region', 'mode-complete', 'show-ligands-control', 'reset-orientation', 'molstar-clear-selection']) await expect(bar.getByTestId(id)).toBeVisible()
    const b = (await bar.boundingBox())!, c = (await viewer.locator('canvas').first().boundingBox())!, v = (await viewer.boundingBox())!
    expect(c.y).toBeGreaterThanOrEqual(b.y + b.height - 1)                   // canvas starts below the toolbar
    expect(c.y + c.height).toBeLessThanOrEqual(v.y + v.height + 1)           // and fills the rest of the viewer box
    expect(Math.abs((b.height + c.height) - v.height)).toBeLessThan(6)   // borders (viewer 2 px, toolbar 1 px) and rounding
    for (const id of ['mode-region', 'reset-orientation', 'molstar-clear-selection']) { const r = (await bar.getByTestId(id).boundingBox())!; expect(r.y + r.height).toBeLessThanOrEqual(c.y + 1) }
  }
  // desktop: ticking Show ligands (4FE5 B: hypoxanthine HPA nearby) keeps every control on one line
  await page.setViewportSize({ width: 1400, height: 1000 }); await page.goto('/RF00167'); await ready(page, '4fe5')
  const bar = page.getByTestId('molstar-toolbar'), h0 = (await bar.boundingBox())!.height
  await page.getByTestId('show-ligands').check(); await expect(page.getByTestId('show-ligands-control')).toHaveAttribute('data-shown', /HPA/)
  await expect(page.getByTestId('show-ligands-control')).toHaveText('Show ligands')                     // no ligand names next to the checkbox
  await expect(page.getByTestId('show-ligands-control')).toHaveAttribute('title', /Shown: HPA/)
  expect((await bar.boundingBox())!.height).toBeLessThanOrEqual(h0 + 1)
  const top = (id: string) => bar.getByTestId(id).boundingBox().then((r) => Math.round(r!.y))
  expect(Math.abs(await top('molstar-clear-selection') - await top('mode-region'))).toBeLessThanOrEqual(2)   // same row (the segmented control has its own border)
})
