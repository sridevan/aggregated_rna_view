import { expect, test, type Page } from '@playwright/test'

/* Structures tab 3D: hovering the backbone tube shows "PDBID.chain RESN<author number> (author)" for the residue under the pointer. */
const viewer = (page: Page) => page.getByTestId('alignment-viewer')
type Probe = { x: number; y: number; label: string } | null

async function hoverResidue(page: Page, acc: string, target?: string) {
  await page.goto(`/${acc}`); await page.getByTestId('tab-structures').click()
  await expect(viewer(page)).toHaveAttribute('data-status', 'ready', { timeout: 120_000 })
  if (target) {   // bring a specific chain into the view (table row click), then probe that chain
    await page.getByTestId('instance-search').fill(target.split('_')[0])
    await page.getByTestId(`row-${target}`).locator('[data-testid="cell-title"]').click()
    await expect.poll(async () => (await viewer(page).getAttribute('data-structures')) ?? '', { timeout: 120_000 }).toContain(target)
    await expect(viewer(page)).toHaveAttribute('data-status', 'ready', { timeout: 120_000 })
  }
  const sn = target ?? ((await viewer(page).getAttribute('data-structures')) ?? '').split(',')[0]
  const box = (await viewer(page).boundingBox())!
  const tip = page.getByTestId('viewer-residue-hover')
  // try a few residues along the chain (one may be hidden behind another part of the tube)
  for (const i of [10, 30, 50, 70, 5]) {
    const p = (await page.evaluate(([s, k]) => (window.__rnaExplorer?.alignmentProbe as (a: string, b: number) => Probe)(s as string, k as number), [sn, i])) as Probe
    if (!p) continue
    await page.mouse.move(p.x - 3, p.y - 3); await page.mouse.move(p.x, p.y, { steps: 3 })   // modelToScreen gives page coordinates
    try { await expect(tip).toBeVisible({ timeout: 2_000 }) } catch { continue }
    return { sn, tip, expected: p.label, box }
  }
  throw new Error('no hover label appeared')
}

test('RF00167: hover label names PDB id, author chain, residue and author number', async ({ page }) => {
  const { sn, tip, box } = await hoverResidue(page, 'RF00167')
  const text = (await tip.textContent())!
  expect(text.startsWith(sn.split('_')[0].toUpperCase() + '.')).toBe(true)
  expect(text).toMatch(new RegExp(`^${sn.split('_')[0].toUpperCase()}\\.[A-Za-z0-9]+ [A-Z0-9]{1,3}-?\\d+[A-Z]? \\(author\\)$`))
  // label stays inside the viewer and disappears when the pointer leaves
  const t = (await tip.boundingBox())!; expect(t.x).toBeGreaterThanOrEqual(box.x); expect(t.x + t.width).toBeLessThanOrEqual(box.x + box.width + 1)
  await page.mouse.move(box.x - 20, box.y - 20)
  await expect(tip).toHaveCount(0)
})

test('RF00177: multi-letter author chain ids are shown in full', async ({ page }) => {
  const { tip, expected } = await hoverResidue(page, 'RF00177', '4v9d_1_BA_Y')
  expect(expected).toMatch(/^4V9D\.BA /)
  await expect(tip).toHaveText(/^[0-9A-Z]{4}\.[A-Za-z0-9]+ [A-Z0-9]{1,3}\d+[A-Z]? \(author\)$/)
})
