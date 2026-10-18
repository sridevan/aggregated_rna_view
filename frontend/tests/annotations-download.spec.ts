import { expect, test, type Page } from '@playwright/test'

/* Annotations tab Download menu: one compact row per JSON file (short label + info icon); descriptions only in the info tooltip;
   the icon never downloads or closes the menu; filenames and contents unchanged. */
const ITEMS: [string, string, string, RegExp][] = [
  ['loops', 'loops.json', 'RNA loop annotations', /BGSU RNA 3D Motif Atlas loops/],
  ['ligands', 'ligands.json', 'Ligand binding annotations', /Bound molecules/],
  ['tracks', 'tracks.json', 'Aggregated annotation tracks', /Per-position annotation frequencies/],
  ['pockets', 'pockets.json', 'Predicted pockets · Per structure', /fpocketR pockets per RNA chain/],
  ['pocket-tracks', 'pocket_tracks.json', 'Predicted pockets · Aggregated', /three pocket filters/],
  ['sasa', 'sasa.json', 'Solvent accessibility · Per nucleotide', /nucleotide solvent accessibility \(Å²\)/],
  ['sasa-tracks', 'sasa_tracks.json', 'Solvent accessibility · Aggregated', /median, quartiles/],
  ['modifications', 'modifications.json', 'Nonstandard nucleotides · Detailed', /Nonstandard nucleotide identities/],
  ['modification-tracks', 'modification_tracks.json', 'Nonstandard nucleotides · Aggregated', /Per-position nonstandard-nucleotide frequencies/],
]
const open = async (page: Page) => {
  await page.goto('/RF00167'); await page.getByTestId('tab-annotations').click()
  await expect(page.getByTestId('download-button')).toHaveText('Download JSON')
  await page.getByTestId('download-button').click(); await expect(page.getByTestId('download-items')).toBeVisible()
}

test('compact rows in the same order: short labels, unchanged files, descriptions only behind the info icons', async ({ page }) => {
  await open(page)
  const menu = page.getByTestId('download-items')
  expect(await menu.getByRole('menuitem').evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')))).toEqual(ITEMS.map(([id]) => `download-${id}`))
  for (const [id, file, label, desc] of ITEMS) {
    const it = page.getByTestId(`download-${id}`)
    await expect(it).toHaveText(label)                                                              // no "(JSON)", no filename, no paragraph
    await expect(it).toHaveAttribute('download', `RF00167_${file}`); expect((await it.getAttribute('href'))!.endsWith(`/families/RF00167/${file}`)).toBe(true)
    const b = (await page.getByTestId(`download-row-${id}`).boundingBox())!; expect(b.height).toBeLessThan(34)   // one line per row
    await page.getByTestId(`download-info-${id}`).getByTestId('info-tip').hover()
    await expect(page.getByTestId('info-popover')).toHaveText(desc)
    const p = (await page.getByTestId('info-popover').boundingBox())!, v = page.viewportSize()!
    expect(p.width).toBeLessThanOrEqual(321); expect(p.x).toBeGreaterThanOrEqual(0); expect(p.x + p.width).toBeLessThanOrEqual(v.width); expect(p.y + p.height).toBeLessThanOrEqual(v.height)
  }
  await expect(menu).not.toContainText('.json')
  const w = (await menu.boundingBox())!.width; expect(w).toBeLessThanOrEqual(300)
})

test('info icon: no download, menu stays open; keyboard focus shows the description; the label still downloads the exact file', async ({ page }) => {
  await open(page)
  let downloads = 0; page.on('download', () => { downloads++ })
  await page.getByTestId('download-info-pockets').getByTestId('info-tip').click()
  await expect(page.getByTestId('info-popover')).toBeVisible(); await page.waitForTimeout(600)
  expect(downloads).toBe(0); await expect(page.getByTestId('download-items')).toBeVisible(); await expect(page).toHaveURL(/\/RF00167$/)
  await page.getByTestId('download-info-pockets').getByTestId('info-tip').click(); await expect(page.getByTestId('info-popover')).toHaveCount(0)   // click again: closes the tooltip only
  await expect(page.getByTestId('download-items')).toBeVisible()
  await page.getByTestId('download-info-sasa').getByTestId('info-tip').focus()
  await expect(page.getByTestId('info-popover')).toContainText('nucleotide solvent accessibility')
  const exp = await (await page.request.get('/families/RF00167/sasa_tracks.json')).text()
  const dl = page.waitForEvent('download'); await page.getByTestId('download-sasa-tracks').click(); const d = await dl
  expect(d.suggestedFilename()).toBe('RF00167_sasa_tracks.json')
  const { readFileSync } = await import('node:fs'); expect(readFileSync((await d.path())!, 'utf8')).toBe(exp)          // contents unchanged
  await expect(page.getByTestId('download-items')).toHaveCount(0)
})

test('narrow screen: menu inside the viewport and scrollable; tooltip not clipped', async ({ page }) => {
  await page.setViewportSize({ width: 380, height: 560 })
  await open(page)
  const menu = page.getByTestId('download-items'), m = (await menu.boundingBox())!
  expect(m.x).toBeGreaterThanOrEqual(0); expect(m.x + m.width).toBeLessThanOrEqual(380)
  expect(await menu.evaluate((e) => getComputedStyle(e).overflowY)).toBe('auto')
  await page.getByTestId('download-modification-tracks').scrollIntoViewIfNeeded()
  await page.getByTestId('download-info-modification-tracks').getByTestId('info-tip').click()
  const p = (await page.getByTestId('info-popover').boundingBox())!
  expect(p.x).toBeGreaterThanOrEqual(0); expect(p.x + p.width).toBeLessThanOrEqual(380); expect(p.y).toBeGreaterThanOrEqual(0); expect(p.y + p.height).toBeLessThanOrEqual(560)
})
