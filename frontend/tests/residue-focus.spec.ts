import { expect, test, type Page } from '@playwright/test'
import { selectColumn } from './helpers'

/* Family Overview residue selection: Mol*'s structure focus (the selected nucleotide in ball-and-stick; no surroundings or contact lines) and a
   close-up camera; same residue in tracks, 2D and 3D; unavailable positions never focus another residue. */
type Cam = { position: number[]; target: number[]; radius: number }
type H = { __rnaExplorer: { focusRepresentation?: number[]; lastMolstarSelection?: { pdbId: string; labelSeqIds: number[] }; lastFocus?: number[]; lastFocusCenter?: number[]; selection?: { positions: number[] }; getCamera: () => Cam | null } }
const hook = (page: Page) => page.evaluate(() => { const h = (window as unknown as H).__rnaExplorer; return { focus: h.focusRepresentation ?? [], sel: h.lastMolstarSelection, cam: h.getCamera(), center: h.lastFocusCenter, positions: h.selection?.positions } })
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }) }
const select = async (page: Page, pos: number) => { await selectColumn(page, pos) }
const dist = (a: number[], b: number[]) => Math.hypot(...a.map((x, i) => x - b[i]))

test('resolved nucleotide: focus representation, close-up camera, 2D/track/3D agree; neighbours, far jumps and reset', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  const whole = (await hook(page)).cam!
  expect((await hook(page)).focus).toEqual([])                                                       // default: whole region, no focus
  // consensus 76 -> C56 of 4FE5 B: same residue everywhere
  await select(page, 76)
  await expect(page.getByTestId('sel-highlightable')).toContainText('76 → C56')
  await expect.poll(async () => (await hook(page)).focus).toEqual([56])
  const h1 = await hook(page); expect(h1.sel!.labelSeqIds).toEqual([56]); expect(h1.positions).toEqual([76])
  await expect(page.getByTestId('rna2d-viewer')).toContainText('56')
  // close-up: the camera moved much closer than the whole-structure view and centres the residue
  await expect.poll(async () => dist((await hook(page)).cam!.target, h1.center!)).toBeLessThan(0.5)
  await expect.poll(async () => { const c = (await hook(page)).cam!; return dist(c.position, c.target) }).toBeLessThan(dist(whole.position, whole.target) * 0.6)
  // neighbour and far selections replace the focus (one target only, no accumulation)
  await select(page, 77); await expect.poll(async () => (await hook(page)).focus).toEqual([57])
  await select(page, 25); await expect.poll(async () => (await hook(page)).focus).toEqual([7])
  const far = await hook(page)
  await expect.poll(async () => dist((await hook(page)).cam!.target, far.center!)).toBeLessThan(0.5)
  for (let k = 0; k < 6; k++) await select(page, 30 + k)
  await expect.poll(async () => (await hook(page)).focus).toEqual([(await hook(page)).sel!.labelSeqIds[0]])
  expect((await hook(page)).focus.length).toBe(1)
  // reset returns to the whole structure; the selection (and its focus) stays
  await page.getByTestId('reset-orientation').click()
  await expect.poll(async () => { const c = (await hook(page)).cam!; return dist(c.position, c.target) }).toBeGreaterThan(dist(whole.position, whole.target) * 0.9)   // back to the whole-structure distance
  expect((await hook(page)).focus.length).toBe(1)
  // clearing the selection removes the focus representation and restores the default whole-structure camera
  await select(page, 76); await expect.poll(async () => (await hook(page)).focus).toEqual([56])
  await expect.poll(async () => { const c = (await hook(page)).cam!; return dist(c.position, c.target) }).toBeLessThan(dist(whole.position, whole.target) * 0.6)
  await page.getByTestId('clear-selection').click()
  await expect.poll(async () => (await hook(page)).focus).toEqual([])
  await expect.poll(async () => dist((await hook(page)).cam!.position, whole.position)).toBeLessThan(0.5)
})

test('unavailable positions never focus another residue; switching instance re-focuses its own residue', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  await select(page, 76); await expect.poll(async () => (await hook(page)).focus).toEqual([56])
  // switching instance: the same consensus position, the new chain's residue (1Y26 X: consensus 76 -> its own label_seq_id)
  await page.getByTestId('instance-select').selectOption('1y26_1_X_A'); await ready(page, '1y26')
  await expect(page.getByTestId('sel-highlightable')).toContainText('76 →')
  const id1y26 = Number((await page.getByTestId('sel-highlightable').textContent())!.match(/76 → [ACGU](\d+)/)![1])
  await expect.poll(async () => (await hook(page)).focus).toEqual([id1y26])
  await expect.poll(async () => (await hook(page)).sel?.pdbId).toBe('1y26')
  // deleted in 1Y26 (consensus 64): no 3D focus, no camera move, status shown
  await page.waitForTimeout(800)                                                                    // let the focus animation settle
  const cam = (await hook(page)).cam!
  await select(page, 64)
  await expect(page.getByTestId('sel-deleted')).toContainText('64')
  await expect.poll(async () => (await hook(page)).focus).toEqual([])
  expect((await hook(page)).sel!.labelSeqIds).toEqual([]); expect(dist((await hook(page)).cam!.target, cam.target)).toBeLessThan(0.01)
  // no coordinates (2XNW A, consensus 19) and uncovered (2B57 A, consensus 19)
  await page.getByTestId('instance-select').selectOption('2xnw_1_A_A'); await ready(page, '2xnw')
  await select(page, 19)
  await expect(page.getByTestId('sel-unobserved')).toContainText('19'); await expect.poll(async () => (await hook(page)).focus).toEqual([])
  await page.getByTestId('instance-select').selectOption('2b57_1_A_A'); await ready(page, '2b57')
  await expect(page.getByTestId('sel-uncovered')).toContainText('19'); await expect.poll(async () => (await hook(page)).focus).toEqual([])
  // back to a resolved position after an unavailable one
  await select(page, 76); await expect.poll(async () => (await hook(page)).focus.length).toBe(1)
})
