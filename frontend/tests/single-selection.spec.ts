import { expect, test, type Page } from '@playwright/test'
import { colX, selectColumn, stripGeom } from './helpers'

/** Family Overview single-nucleotide selection: one active selection, replace on click, clear on re-click, no shift ranges; Mol* highlights and centres; 2D follows; instance switching remaps. */
type W = { __rnaExplorer: { selection: { positions: number[]; instance: string; active: unknown }; lastMolstarSelection?: { pdbId: string; labelSeqIds: number[] }; lastFocus?: number[]; lastFocusCenter?: number[]; getCamera: () => { position: number[]; target: number[] } } }
const hook = (page: Page) => page.evaluate(() => { const h = (window as unknown as W).__rnaExplorer; return { sel: h.selection, mol: h.lastMolstarSelection, focus: h.lastFocus, center: h.lastFocusCenter, cam: h.getCamera() } })
const ready = async (page: Page, pdb: string, need2d = true) => {
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', need2d ? 'ready' : /ready|unavailable/, { timeout: 60_000 })
}
const rna2dSelectedIds = (page: Page) => page.evaluate(() => [...document.querySelectorAll('[data-testid="rna2d-viewer"] text.rnaviewEle')]
  .filter((e) => { const f = (e as SVGElement).getAttribute('fill') || ''; return f !== '' && f !== 'none' && f !== 'black' && f !== '#000000' && f !== '#323232' }).map((e) => +(e.className as unknown as SVGAnimatedString).baseVal.split('_').pop()!))
/** click a consensus column of the tracks SVG (viewBox geometry: label column 140 units, 9 units per position) */
const clickColumn = async (page: Page, pos: number, modifiers: ('Shift')[] = []) => {
  const g = await stripGeom(page, 'consensus-tracks')
  await page.mouse.click(await colX(page, 'consensus-tracks', pos), g.b.y + 12 * g.scale, { modifiers })
}

test('4FE5 B: one nucleotide at a time, replace / clear semantics, no shift range, highlighted and centred in Mol*, mirrored in 2D', async ({ page }) => {
  await page.setViewportSize({ width: 1400, height: 1000 })
  await page.goto('/RF00167')
  await ready(page, '4fe5')
  await expect(page.getByTestId('selection-status')).toHaveCount(0)
  const cam0 = JSON.stringify((await hook(page)).cam)
  await clickColumn(page, 50)
  await expect.poll(async () => (await hook(page)).sel.positions).toEqual([50])
  await expect(page.getByTestId('sel-highlightable')).toContainText('50 → ')
  await expect(page.getByTestId('selection-status')).toHaveAttribute('data-ok', 'true')
  await expect.poll(async () => (await hook(page)).mol?.labelSeqIds.length).toBe(1)
  const h1 = await hook(page)
  expect(h1.focus).toEqual(h1.mol!.labelSeqIds)
  await expect.poll(async () => JSON.stringify((await hook(page)).cam)).not.toBe(cam0)
  await page.waitForTimeout(500)
  const h1b = await hook(page)
  expect(Math.hypot(...h1b.cam.target.map((t, i) => t - h1b.center![i]))).toBeLessThan(1.5)   // camera target = residue centre
  await expect.poll(() => rna2dSelectedIds(page)).toEqual(h1.mol!.labelSeqIds)
  // replace
  await clickColumn(page, 36)
  await expect.poll(async () => (await hook(page)).sel.positions).toEqual([36])
  await expect.poll(async () => (await hook(page)).mol?.labelSeqIds.length).toBe(1)
  const h2 = await hook(page); expect(h2.mol!.labelSeqIds).not.toEqual(h1.mol!.labelSeqIds)
  await expect.poll(() => rna2dSelectedIds(page)).toEqual(h2.mol!.labelSeqIds)
  // shift-click never creates a range
  await clickColumn(page, 60, ['Shift'])
  await expect.poll(async () => (await hook(page)).sel.positions).toEqual([60])
  // manual input follows the same rules: same position again clears
  await selectColumn(page, 60)
  await expect.poll(async () => (await hook(page)).sel.positions).toEqual([])
  await expect(page.getByTestId('selection-status')).toHaveCount(0)
  await expect.poll(async () => (await hook(page)).mol?.labelSeqIds).toEqual([])
  await expect.poll(() => rna2dSelectedIds(page)).toEqual([])
  // clearing did not move the camera (compare with the camera at the time of the last focus)
  const camAfterClear = JSON.stringify((await hook(page)).cam)
  await page.waitForTimeout(400)
  expect(JSON.stringify((await hook(page)).cam)).toBe(camAfterClear)
  // selected column shaded in the tracks; re-click clears
  await clickColumn(page, 25)
  await expect.poll(async () => (await hook(page)).sel.positions).toEqual([25])
  await clickColumn(page, 25)
  await expect.poll(async () => (await hook(page)).sel.positions).toEqual([])
})

test('instance switching keeps the consensus position and remaps it: 1Y26 X, 5E54 A and B; unavailable positions are reported, never mis-highlighted', async ({ page }) => {
  await page.setViewportSize({ width: 1400, height: 1000 })
  await page.goto('/RF00167')
  await ready(page, '4fe5')
  await selectColumn(page, 50)
  await expect(page.getByTestId('sel-highlightable')).toContainText('50 → ')
  const inst = await (await page.request.get('/families/RF00167/instances.json')).json() as { seq_name: string; pdb_id: string; auth_asym_id: string; pdbe_2d_available: boolean | null }[]
  for (const [pdb, ch] of [['1y26', 'X'], ['5e54', 'A'], ['5e54', 'B']] as const) {
    const i = inst.find((r) => r.pdb_id === pdb && r.auth_asym_id === ch)!
    await page.getByTestId('instance-select').selectOption(i.seq_name)
    await ready(page, pdb, i.pdbe_2d_available !== false)
    const h = await hook(page)
    expect(h.sel.positions).toEqual([50]); expect(h.sel.instance).toBe(i.seq_name)
    await expect(page.getByTestId('sel-highlightable')).toContainText('50 → ')
    await expect.poll(async () => (await hook(page)).mol).toMatchObject({ pdbId: pdb })
    await expect.poll(async () => (await hook(page)).mol?.labelSeqIds.length).toBe(1)
    await expect.poll(async () => (await hook(page)).focus).toEqual((await hook(page)).mol!.labelSeqIds)   // focused after the new structure loaded
    if (i.pdbe_2d_available !== false) await expect.poll(() => rna2dSelectedIds(page)).toEqual((await hook(page)).mol!.labelSeqIds)
  }
  // 6P2H covers consensus 26-80 only: consensus 1 is uncovered there; 9NBH A: consensus 39 deleted, 50 without coordinates
  await selectColumn(page, 1)
  await page.getByTestId('instance-select').selectOption('6p2h_1_A_A'); await ready(page, '6p2h', false)
  await expect(page.getByTestId('selection-status')).toHaveAttribute('data-ok', 'false')
  await expect(page.getByTestId('selection-status')).toContainText('not covered by the alignment')
  await expect.poll(async () => (await hook(page)).mol?.labelSeqIds).toEqual([])
  await selectColumn(page, 39)
  await page.getByTestId('instance-select').selectOption('9nbh_1_A_A'); await ready(page, '9nbh', false)
  await expect(page.getByTestId('selection-status')).toContainText('Rfam position 39 is deleted in 9NBH A')
  await expect.poll(async () => (await hook(page)).mol?.labelSeqIds).toEqual([])
  await selectColumn(page, 50)
  await expect(page.getByTestId('selection-status')).toContainText('no experimental coordinates')
  await expect.poll(async () => (await hook(page)).mol?.labelSeqIds).toEqual([])
})

test('switching the instance while the first structure is still loading never leaves an empty scene; selection then still highlights and centres', async ({ page }) => {
  await page.goto('/RF00167')
  await page.getByTestId('instance-select').selectOption('1y26_1_X_A')   // immediately, while 4FE5 is loading
  await ready(page, '1y26')
  type V = { __rnaExplorer: { molstar: { plugin: { canvas3d?: { boundingSphereVisible: { radius: number } } } }; molstarComponents?: { label: string; hidden: boolean | null; present?: boolean }[] } }
  await expect.poll(() => page.evaluate(() => (window as unknown as V).__rnaExplorer.molstar.plugin.canvas3d!.boundingSphereVisible.radius)).toBeGreaterThan(5)
  expect((await page.evaluate(() => (window as unknown as V).__rnaExplorer.molstarComponents))!.find((c) => c.label === 'Rfam-mapped region' && c.present !== undefined)?.present).toBe(true)
  await selectColumn(page, 50)
  await expect.poll(async () => (await hook(page)).mol?.labelSeqIds).toEqual([34])
  await expect.poll(async () => (await hook(page)).focus).toEqual([34])
  await expect.poll(async () => { const h = await hook(page); return Math.hypot(...h.cam.target.map((t, i) => t - h.center![i])) }).toBeLessThan(1.5)   // once the focus animation has finished
})
