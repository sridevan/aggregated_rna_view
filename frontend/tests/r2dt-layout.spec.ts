import { expect, test, type Page } from '@playwright/test'
import { selectColumn } from './helpers'

/** Exploration (branch explore/r2dt-2d-3d-viewer): the 2D layout toggle switches the PDBe RNA viewer between PDBe topology files and
 *  R2DT-generated layouts (Rfam template), with the consensus selection still mirrored by label_seq_id. */
type W = { __rnaExplorer: { rna2dSelected?: number[] } }
const shot = (page: Page, name: string) => page.getByTestId('rna2d-viewer').screenshot({ path: `../experiments/r2dt_viewer/screenshots/${name}.png` })

test('R2DT layout for 4FE5 B: toggle, selection mirrored, PDBe restored', async ({ page }) => {
  await page.goto('/RF00167')
  const v = page.getByTestId('rna2d-viewer')
  await expect(v).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await expect(v).toHaveAttribute('data-layout', 'pdbe')
  await shot(page, '4fe5_pdbe')
  await page.getByTestId('layout-2d-r2dt').click()
  await expect(v).toHaveAttribute('data-layout', 'r2dt')
  await expect(v).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await expect(page.getByTestId('layout-2d-note')).toContainText('template RF00167')
  await shot(page, '4fe5_r2dt')
  // consensus 50 → the same label_seq_id in both sources (R2DT positions are validated against the mapping)
  await selectColumn(page, 50)
  await expect.poll(() => page.evaluate(() => (window as unknown as W).__rnaExplorer.rna2dSelected?.length)).toBe(1)
  const id = await page.evaluate(() => (window as unknown as W).__rnaExplorer.rna2dSelected![0])
  await expect(v.locator(`.rnaview_4fe5_${id}`).first()).toBeVisible()
  await shot(page, '4fe5_r2dt_sel50')
  await page.getByTestId('layout-2d-pdbe').click()
  await expect(v).toHaveAttribute('data-layout', 'pdbe')
  await expect(v).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
})

test('R2DT layouts share the Rfam template across chains; the toggle persists across instances', async ({ page }) => {
  await page.goto('/RF00167')
  await page.getByTestId('layout-2d-r2dt').click()
  for (const [seq, pdb] of [['1y26_1_X_A', '1y26'], ['5e54_1_A_A', '5e54']]) {
    await page.getByTestId('instance-select').selectOption(seq)
    const v = page.getByTestId('rna2d-viewer')
    await expect(v).toHaveAttribute('data-pdb', pdb)
    await expect(v).toHaveAttribute('data-layout', 'r2dt')
    await expect(v).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
    await shot(page, `${pdb}_r2dt`)
  }
})

/** Other families: an R2DT layout where PDBe draws a different one (TPP 2HOJ A), one where PDBe has no topology (RNase P 9OV3 A), and
 *  a tRNA on a GtRNAdb template (1EHZ A). Screenshots of both sources feed experiments/r2dt_viewer/REPORT.md. */
for (const [acc, seq, pdb, pdbeStatus] of [['RF00059', '2hoj_1_A_A', '2hoj', 'ready'], ['RF00011', '9ov3_1_A_A', '9ov3', 'unavailable'], ['RF00005', '1ehz_1_A_A', '1ehz', 'ready']] as const) {
  test(`${acc} ${seq}: PDBe ${pdbeStatus}, R2DT layout ready`, async ({ page }) => {
    await page.goto(`/${acc}`)
    await page.getByTestId('instance-select').selectOption(seq)
    const v = page.getByTestId('rna2d-viewer')
    await expect(v).toHaveAttribute('data-pdb', pdb)
    await expect(v).toHaveAttribute('data-status', pdbeStatus, { timeout: 60_000 })
    await shot(page, `${pdb}_pdbe`)
    await page.getByTestId('layout-2d-r2dt').click()
    await expect(v).toHaveAttribute('data-layout', 'r2dt')
    await expect(v).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
    await shot(page, `${pdb}_r2dt`)
  })
}
