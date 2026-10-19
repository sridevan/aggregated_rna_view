import { expect, test, type Page } from '@playwright/test'
import { colX, selectColumn } from './helpers'

/* Family Overview loop/viewer refinements: simplified Mol* display (no water, no contact lines, no surroundings), structured loop tooltip
   with the explicit BGSU motif group only, toggle-off of loops and nucleotides restoring the default 2D and 3D views, and a red inset outline
   on every block of the selected loop. Expected values come from loops.json. */
type Nt = { seq_name: string | null; label_seq_id: number; alignment_state: string; has_coordinates: boolean }
type Loop = { loop_id: string; loop_type: string; motif_group: string | null; annotation_text: string | null; strands: unknown[]; nucleotides: Nt[] }
type Cam = { position: number[]; target: number[] }
type H = { __rnaExplorer: { viewerIds: number[]; rna2dSelected?: number[]; rna2dReset?: number; molstarReset?: number; lastMolstarSelection?: { labelSeqIds: number[] }; focusRepresentation?: number[]; focusComponents?: string[] | null; focusIgnoreHydrogens?: boolean | null
  molstarComponents?: { label?: string; hidden: boolean | null }[]; selection: { active: { kind: string; id?: string } | null; viewReset: number }; getCamera: () => Cam | null
  molstar: { plugin: { state: { data: { cells: Map<string, { obj?: { label?: string } }> } } } } } }
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }) }
const hook = (page: Page) => page.evaluate(() => { const h = (window as unknown as H).__rnaExplorer; return { v: h.viewerIds, d2: h.rna2dSelected ?? [], d2Reset: h.rna2dReset ?? 0, molReset: h.molstarReset ?? 0, mol: h.lastMolstarSelection?.labelSeqIds ?? [], focus: h.focusRepresentation ?? [],
  focusComponents: h.focusComponents, noH: h.focusIgnoreHydrogens, comps: h.molstarComponents ?? [], active: h.selection.active, cam: h.getCamera()! } })
/** labels of every Mol* state cell (representations included): the focus contact lines are a cell labelled "[Focus] Non-covalent interactions" */
const cellLabels = (page: Page) => page.evaluate(() => [...(window as unknown as H).__rnaExplorer.molstar.plugin.state.data.cells.values()].map((c) => c.obj?.label ?? ''))
/** transform tags of every Mol* state cell: the focus behaviour tags its target, surroundings and contact (NCI) representations */
const cellTags = (page: Page) => page.evaluate(() => [...(window as unknown as H).__rnaExplorer.molstar.plugin.state.data.cells.values()].flatMap((c) => (c as unknown as { transform: { tags?: string[] } }).transform.tags ?? []))
const view2d = (page: Page) => page.evaluate(() => { const s = document.querySelector('[data-testid="rna2d-viewer"] svg.rnaTopoSvg')!; return [s.getAttribute('viewBox'), ...[...s.querySelectorAll('g[transform]')].map((g) => g.getAttribute('transform'))].join('|') })
const sorted = (a: number[]) => [...a].sort((x, y) => x - y)
const dist = (a: number[], b: number[]) => Math.hypot(...a.map((x, i) => x - b[i]))
let LOOPS: Loop[] = []
const members = (id: string, sn: string) => sorted(LOOPS.find((l) => l.loop_id === id)!.nucleotides.filter((n) => n.seq_name === sn && n.has_coordinates).map((n) => n.label_seq_id))
const SN = '4fe5_1_B_A'
/** load, wait for both viewers, record the default 3D camera and 2D view */
const open = async (page: Page) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  LOOPS = (await (await page.request.get('/families/RF00167/loops.json')).json()).loops
  await page.waitForTimeout(500)
  const h = await hook(page)
  return { cam: h.cam, view: await view2d(page) }
}
const atDefault = async (page: Page, d: { cam: Cam; view: string }) => {
  await expect.poll(async () => dist((await hook(page)).cam.position, d.cam.position)).toBeLessThan(0.5)
  await expect.poll(async () => dist((await hook(page)).cam.target, d.cam.target)).toBeLessThan(0.5)
  // 2D default: the whole topology fitted again (same view box within 1 %; it is re-measured) and no zoom/pan transform left
  const vb = (v: string) => v.split('|')[0].split(' ').map(Number)
  await expect.poll(async () => { const now = await view2d(page); const a = vb(now), b = vb(d.view); return a.every((x, i) => Math.abs(x - b[i]) <= 0.01 * Math.max(b[2], b[3])) && now.split('|').slice(1).every((t) => /^translate\(0,0\) scale\(1\)$/.test(t)) }).toBe(true)
}
/** user zoom in the 2D viewer (its own zoom-in button; wheel zoom is disabled by the viewer) */
const zoom2d = async (page: Page) => { await page.locator('[data-testid="rna2d-viewer"] [id^="rnaTopologyZoomIn-"]').first().click(); await page.waitForTimeout(500); expect(await view2d(page)).toMatch(/scale\(1\.2/) }
/** click the centre of a consensus column in the primary tracks */
const clickColumn = async (page: Page, pos: number) => {
  const b = (await page.getByTestId('consensus-tracks').boundingBox())!
  await page.mouse.click(await colX(page, 'consensus-tracks', pos), b.y + b.height / 2)
}

test('simplified Mol* display: water hidden, focus = selected loop nucleotides only, no contact lines or surroundings', async ({ page }) => {
  await open(page)
  const h0 = await hook(page)
  expect(h0.focusComponents).toEqual(['target'])                                                        // no 'surroundings', no 'interactions'
  expect(h0.noH).toBe(true)                                                                             // heavy atoms only
  // water exists in 4FE5 but is hidden by default, in both view modes (visibility only: the component is kept)
  await page.getByTestId('mode-complete').click(); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-mode', 'complete')
  await expect.poll(async () => (await hook(page)).comps.filter((c) => c.label === 'Water')).toEqual([{ label: 'Water', hidden: true }])
  // the rest of the complete structure stays visible; PDBe's default Polymer and Ligand components are replaced by the Complete-structure
  // layers and the "Show ligands" layers (off by default)
  await expect.poll(async () => (await hook(page)).comps.filter((c) => !['Water', 'Rfam-mapped region', 'Polymer', 'Ligand'].includes(c.label)).every((c) => c.hidden === false)).toBe(true)
  await page.getByTestId('mode-region').click(); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-mode', 'region')
  await expect.poll(async () => (await hook(page)).comps.find((c) => c.label === 'Water')?.hidden).toBe(true)
  // select the junction loop: the focus shows exactly its resolved members (no 5 Å surroundings) and no interaction lines
  await page.getByTestId('loops-toggle').click(); await page.getByTestId('loop-block-J3_4FE5_001-0').click()
  await expect.poll(async () => sorted((await hook(page)).focus)).toEqual(members('J3_4FE5_001', SN))
  const tags = await cellTags(page)
  expect(tags.filter((t) => /structure-focus-target-repr/.test(t)).length).toBe(1)                     // the selected loop nucleotides (ball-and-stick)
  expect(tags.filter((t) => /structure-focus-surr-repr|structure-focus-surr-nci-repr/.test(t))).toEqual([])   // no surroundings, no contact lines
  expect((await cellLabels(page)).filter((l) => /Non-covalent|Interactions/i.test(l))).toEqual([])
  // single nucleotides use the same simplified focus
  await selectColumn(page, 76)
  await expect.poll(async () => (await hook(page)).focus).toEqual([56])
  expect((await cellTags(page)).filter((t) => /structure-focus-surr-repr|structure-focus-surr-nci-repr/.test(t))).toEqual([])
})

test('structured loop tooltip: category, BGSU loop id, counts; motif group only when BGSU assigns one; kept inside the viewport', async ({ page }) => {
  await open(page)
  await page.getByTestId('loops-toggle').click()
  await page.getByTestId('loop-block-HL_4FE5_002-0').hover()
  const tip = page.getByTestId('loop-tooltip')
  await expect(page.getByTestId('loop-tooltip-title')).toHaveText('Hairpin loop (HL)')
  await expect(page.getByTestId('loop-tooltip-id')).toHaveText('HL_4FE5_002')
  await expect(page.getByTestId('loop-tooltip-counts')).toHaveText('8 nucleotides · 1 segment')
  const hl = LOOPS.find((l) => l.loop_id === 'HL_4FE5_002')!
  await expect(page.getByTestId('loop-tooltip-motif')).toContainText(`Motif: ${hl.motif_group} (${hl.annotation_text})`)
  // discontinuous internal loop: segment count from the original strands; motif group id only (no name in the source)
  await page.getByTestId('loop-block-IL_4FE5_001-1').hover()
  await expect(page.getByTestId('loop-tooltip-title')).toHaveText('Internal loop (IL)'); await expect(page.getByTestId('loop-tooltip-counts')).toHaveText('6 nucleotides · 2 segments')
  const il = LOOPS.find((l) => l.loop_id === 'IL_4FE5_001')!
  await expect(page.getByTestId('loop-tooltip-motif')).toContainText(`Motif: ${il.motif_group}${il.annotation_text ? ` (${il.annotation_text})` : ''}`)
  await expect(page.getByTestId('loop-tooltip-motif')).not.toContainText('IL_4FE5')                     // loop ids never shown as a motif
  await page.getByTestId('loop-block-J4_4FE5_001-0').hover(); await expect(page.getByTestId('loop-tooltip-title')).toHaveText('Four-way junction loop (J4)')
  // inside the viewport, also for the right-most block on a narrow window
  await page.setViewportSize({ width: 700, height: 900 })
  const last = page.locator('[data-testid^="loop-block-"]').last()
  await last.scrollIntoViewIfNeeded(); await last.hover()
  await expect(tip).toBeVisible()
  const tb = (await tip.boundingBox())!
  expect(tb.x).toBeGreaterThanOrEqual(0); expect(tb.x + tb.width).toBeLessThanOrEqual(700); expect(tb.y).toBeGreaterThanOrEqual(0); expect(tb.y + tb.height).toBeLessThanOrEqual(900)
  await page.mouse.move(5, 5); await expect(tip).toHaveCount(0)
  // a loop without an explicit BGSU motif-group assignment: no motif line
  const bare = LOOPS.find((l) => !l.motif_group && l.nucleotides.some((n) => n.seq_name && n.alignment_state === 'match'))!
  const sn = bare.nucleotides.find((n) => n.seq_name && n.alignment_state === 'match')!.seq_name!
  await page.setViewportSize({ width: 1400, height: 1000 })
  await page.getByTestId('instance-select').selectOption(sn); await ready(page, sn.split('_')[0])
  await page.locator(`[data-testid^="loop-block-${bare.loop_id}-"]`).first().hover()
  await expect(page.getByTestId('loop-tooltip-id')).toHaveText(bare.loop_id)
  await expect(page.getByTestId('loop-tooltip-motif')).toHaveCount(0)
})

test('loop toggle: complete membership; same loop again clears and restores default 2D and 3D views; another loop replaces without a reset', async ({ page }) => {
  const def = await open(page)
  await page.getByTestId('loops-toggle').click()
  await zoom2d(page)
  await page.getByTestId('loop-block-IL_4FE5_001-0').click()
  const il = members('IL_4FE5_001', SN)
  await expect.poll(async () => sorted((await hook(page)).v)).toEqual(il)
  await expect.poll(async () => sorted((await hook(page)).d2)).toEqual(il); await expect.poll(async () => sorted((await hook(page)).focus)).toEqual(il)
  await expect.poll(async () => dist((await hook(page)).cam.position, def.cam.position)).toBeGreaterThan(5)          // close-up
  // another loop replaces the selection: no intermediate reset of either view
  const r0 = await hook(page)
  await page.getByTestId('loop-block-HL_4FE5_002-0').click()
  await expect.poll(async () => sorted((await hook(page)).mol)).toEqual(members('HL_4FE5_002', SN))
  const r1 = await hook(page); expect(r1.molReset).toBe(r0.molReset); expect(r1.d2Reset).toBe(r0.d2Reset)
  expect(await view2d(page)).toMatch(/scale\(1\.2/)                                                       // the user's 2D zoom is kept
  // the same loop again (any of its segments): cleared everywhere, both default views back
  await page.getByTestId('loop-block-IL_4FE5_001-1').click(); await expect.poll(async () => (await hook(page)).active?.id).toBe('IL_4FE5_001')
  await page.getByTestId('loop-block-IL_4FE5_001-0').click()
  await expect.poll(async () => (await hook(page)).active).toBeNull()
  await expect.poll(async () => (await hook(page)).v).toEqual([]); await expect.poll(async () => (await hook(page)).d2).toEqual([])
  await expect.poll(async () => (await hook(page)).mol).toEqual([]); await expect.poll(async () => (await hook(page)).focus).toEqual([])
  await expect(page.locator('[data-testid^="loop-outline-"]')).toHaveCount(0)
  await expect(page.getByTestId('selection-status')).toHaveCount(0)
  await atDefault(page, def)
  // Clear selection does the same
  await page.getByTestId('loop-block-J3_4FE5_001-1').click(); await expect.poll(async () => sorted((await hook(page)).focus)).toEqual(members('J3_4FE5_001', SN))
  await zoom2d(page)
  await page.getByTestId('clear-selection').click()
  await expect.poll(async () => (await hook(page)).focus).toEqual([]); await atDefault(page, def)
})

test('nucleotide toggle: clicking the selected column again clears it and restores both default views', async ({ page }) => {
  const def = await open(page)
  await zoom2d(page)
  await clickColumn(page, 76)
  await expect.poll(async () => (await hook(page)).focus).toEqual([56]); expect((await hook(page)).active).toEqual({ kind: 'position', pos: 76 })
  await clickColumn(page, 77)                                                                           // replacement, not a toggle
  await expect.poll(async () => (await hook(page)).focus).toEqual([57])
  await clickColumn(page, 77)
  await expect.poll(async () => (await hook(page)).active).toBeNull()
  await expect.poll(async () => (await hook(page)).focus).toEqual([]); await expect.poll(async () => (await hook(page)).d2).toEqual([])
  await atDefault(page, def)
  // a single nucleotide never selects a loop
  await page.getByTestId('loops-toggle').click(); await clickColumn(page, 30)
  expect((await hook(page)).active).toEqual({ kind: 'position', pos: 30 }); await expect(page.locator('[data-testid^="loop-outline-"]')).toHaveCount(0)
})

test('red outline on every segment of the selected loop only; geometry unchanged; hover distinct; overlaps selectable; instance change clears', async ({ page }) => {
  await open(page)
  await page.getByTestId('loops-toggle').click()
  // in consensus columns (divided by the shared cell width): a page-width change is not a geometry change
  const geom = () => page.locator('[data-testid^="loop-block-"]').evaluateAll((es) => { const cell = +(document.querySelector('[data-testid="loop-tracks"]')!.getAttribute('data-cell'))!
    return es.map((e) => [e.getAttribute('data-testid'), ((+e.getAttribute('x')! - 140) / cell).toFixed(2), (+e.getAttribute('width')! / cell).toFixed(2), e.getAttribute('y')].join(':')) })
  await expect(page.getByTestId('loop-block-J3_4FE5_001-1')).toBeVisible()   // loops are fetched on opening: wait before the reference geometry
  const g0 = await geom()
  await page.getByTestId('loop-block-J3_4FE5_001-1').click()
  // all three segments of J3 outlined in solid red; the overlapping J4 and every other loop not outlined
  for (const k of [0, 1, 2]) { const o = page.getByTestId(`loop-outline-J3_4FE5_001-${k}`); await expect(o).toHaveAttribute('stroke', '#d1202f'); expect(await o.getAttribute('stroke-dasharray')).toBeNull() }
  await expect(page.locator('[data-testid^="loop-outline-"]')).toHaveCount(3)
  await expect(page.getByTestId('loop-block-J3_4FE5_001-0')).toHaveAttribute('fill', await page.getByTestId('loop-block-J4_4FE5_001-0').getAttribute('fill') as string)   // category fill kept
  expect(await geom()).toEqual(g0)                                                                       // widths and positions unchanged
  // outline drawn inside the block
  const blk = (await page.getByTestId('loop-block-J3_4FE5_001-0').boundingBox())!, out = (await page.getByTestId('loop-outline-J3_4FE5_001-0').boundingBox())!
  expect(out.x).toBeGreaterThanOrEqual(blk.x - 0.01); expect(out.x + out.width).toBeLessThanOrEqual(blk.x + blk.width + 0.01)
  // hovering the overlapping J4 shows a dashed dark outline (not red); J3 stays outlined in red
  await page.getByTestId('loop-block-J4_4FE5_001-0').hover()
  const hv = page.getByTestId('loop-hover-J4_4FE5_001-0'); await expect(hv).toHaveAttribute('stroke-dasharray', /\d/); expect(await hv.getAttribute('stroke')).not.toBe('#d1202f')
  await expect(page.locator('[data-testid^="loop-outline-"]')).toHaveCount(3)
  // the overlapping loop is selectable and takes over the outline
  await page.getByTestId('loop-block-J4_4FE5_001-0').click()
  await expect(page.locator('[data-testid^="loop-outline-J4_4FE5_001-"]')).toHaveCount(4); await expect(page.locator('[data-testid^="loop-outline-J3_"]')).toHaveCount(0)
  await expect.poll(async () => sorted((await hook(page)).v)).toEqual(members('J4_4FE5_001', SN))
  // instance change: incompatible loop cleared, no stale highlight in 2D or 3D, section stays open
  await page.getByTestId('instance-select').selectOption('1y26_1_X_A'); await ready(page, '1y26')
  await expect.poll(async () => (await hook(page)).active).toBeNull()
  await expect(page.locator('[data-testid^="loop-outline-"]')).toHaveCount(0)
  await expect.poll(async () => (await hook(page)).mol).toEqual([]); await expect.poll(async () => (await hook(page)).focus).toEqual([])
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-pdb', '1y26'); await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await expect.poll(async () => (await hook(page)).d2).toEqual([])
  await expect(page.getByTestId('loops-section')).toHaveAttribute('data-open', 'true')
})
