import { expect, test, type Page } from '@playwright/test'
import { selectColumn } from './helpers'

/** Milestone 3C validation: the shared instance selector in Family Overview drives 2D, 3D and the instance track; the consensus selection survives switching. */
type Inst = { seq_name: string; pdb_id: string; auth_asym_id: string; label_asym_id: string; pdbe_2d_available: boolean | null }
const viewers = async (page: Page, i: Inst) => {
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', i.pdb_id)
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-pdb', i.pdb_id)
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', i.pdbe_2d_available === false ? 'unavailable' : 'ready', { timeout: 60_000 })
  await expect(page.getByTestId('overview-current')).toContainText(`${i.pdb_id.toUpperCase()} chain ${i.auth_asym_id}`)
  await expect(page.getByTestId('consensus-tracks-labels')).toContainText(`${i.pdb_id.toUpperCase()} ${i.auth_asym_id}`)   // chain name in the sticky label column
}

test('default 4FE5 B, then 1Y26 X, 5E54 A and B, and the four 3RKF chains as separate instances, selection re-translated each time', async ({ page }) => {
  await page.goto('/RF00167')
  const inst: Inst[] = await (await page.request.get('/families/RF00167/instances.json')).json()
  const by = (pdb: string, ch: string) => inst.find((i) => i.pdb_id === pdb && i.auth_asym_id === ch)!
  await expect(page.getByTestId('instance-select')).toHaveValue('4fe5_1_B_A')
  await viewers(page, by('4fe5', 'B'))
  await selectColumn(page, 50)
  await expect(page.getByTestId('sel-highlightable')).toContainText('50 → ')
  const targets = [by('1y26', 'X'), by('5e54', 'A'), by('5e54', 'B'), ...['A', 'B', 'C', 'D'].map((c) => by('3rkf', c))]
  expect(targets.every(Boolean)).toBe(true)
  expect(new Set(targets.map((t) => t.seq_name)).size).toBe(7)                                        // chains of one entry are distinct instances
  expect(inst.filter((i) => i.pdb_id === '3rkf').length).toBe(4); expect(inst.filter((i) => i.pdb_id === '5e54').length).toBe(2)
  for (const t of targets) {
    await page.getByTestId('instance-select').selectOption(t.seq_name)
    await viewers(page, t)
    await expect(page.getByTestId('sel-highlightable')).toContainText('50 → ')   // consensus 50 is a match residue with coordinates in all of them
    type W = { __rnaExplorer: { selection: { positions: number[]; instance: string }; lastMolstarSelection?: { labelSeqIds: number[] } } }
    const s = await page.evaluate(() => (window as unknown as W).__rnaExplorer.selection)
    expect(s.positions).toEqual([50]); expect(s.instance).toBe(t.seq_name)
    await expect.poll(() => page.evaluate(() => (window as unknown as W).__rnaExplorer.lastMolstarSelection?.labelSeqIds?.length)).toBe(1)
  }
  // options are one per chain instance, labelled with PDB id + author chain
  const options = await page.getByTestId('instance-select').locator('option').allTextContents()
  expect(options.length).toBe(inst.length)
  expect(options.filter((o) => o.startsWith('3RKF ')).length).toBe(4)
})
