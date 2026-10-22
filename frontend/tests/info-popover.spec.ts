import { expect, test, type Page } from '@playwright/test'

/** Category information popovers: concise texts, always inside the viewport (flip, clamp, height limit), click pins / click or Escape / outside closes, never toggles the category. */
const PAD = 12
const inside = async (page: Page) => {
  const b = (await page.getByTestId('info-popover').boundingBox())!; const vp = page.viewportSize()!
  expect(b.x).toBeGreaterThanOrEqual(PAD - 0.5); expect(b.x + b.width).toBeLessThanOrEqual(vp.width - PAD + 0.5); expect(b.y).toBeGreaterThanOrEqual(PAD - 0.5); expect(b.y + b.height).toBeLessThanOrEqual(vp.height - PAD + 0.5)
  expect(b.width).toBeLessThanOrEqual(340)
  // content not cut: the popover's own scroll height fits or it scrolls internally
  const fits = await page.getByTestId('info-popover').evaluate((e) => e.scrollHeight <= e.clientHeight + 1 || getComputedStyle(e).overflowY === 'auto')
  expect(fits).toBe(true)
  return b
}
/** the popover must sit directly below or above its icon, as its data-side says */
const adjacent = async (page: Page, id: string) => {
  const icon = (await page.getByTestId(`category-header-${id}`).getByTestId('info-tip').boundingBox())!; const b = (await page.getByTestId('info-popover').boundingBox())!
  const side = await page.getByTestId('info-popover').getAttribute('data-side')
  if (side === 'below') expect(Math.abs(b.y - (icon.y + icon.height + 6))).toBeLessThan(2); else expect(Math.abs(b.y + b.height - (icon.y - 6))).toBeLessThan(2)
  return b
}
const openCat = async (page: Page, id: string) => { await page.getByTestId(`category-header-${id}`).getByTestId('info-tip').click(); await expect(page.getByTestId('info-popover')).toBeVisible(); await adjacent(page, id); return inside(page) }

test('popovers stay in the viewport: below by default, flipped above near the bottom, clamped at the edges, limited in height on a short viewport, repositioned on scroll; interactions', async ({ page }) => {
  await page.setViewportSize({ width: 1400, height: 1000 })
  await page.goto('/RF00167'); await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-sasa')).toBeVisible()
  // texts
  await page.getByTestId('category-header-sasa').getByTestId('info-tip').hover()
  await expect(page.getByTestId('info-popover')).toHaveText('Per-nucleotide solvent-accessible surface area (SASA) of isolated RNA chains, summarised by Rfam position (median Å²).')
  await page.mouse.move(5, 5); await expect(page.getByTestId('info-popover')).toHaveCount(0)   // hover is transient
  await page.getByTestId('category-header-modifications').getByTestId('info-tip').hover()
  await expect(page.getByTestId('info-popover')).toHaveText('Nonstandard nucleotide components in deposited RNA sequences, mapped to Rfam positions. These are not necessarily biological RNA modifications.')
  await page.mouse.move(5, 5)
  for (const id of ['loops']) { await page.getByTestId(`category-header-${id}`).getByTestId('info-tip').hover(); await expect(page.getByTestId('info-popover')).toBeVisible(); await page.mouse.move(5, 5) }
  // near the top: below the icon
  await page.evaluate(() => window.scrollTo(0, 0))
  let b = await openCat(page, 'loops'); const icon = (await page.getByTestId('category-header-loops').getByTestId('info-tip').boundingBox())!
  expect(b.y).toBeGreaterThan(icon.y + icon.height); await expect(page.getByTestId('info-popover')).toHaveAttribute('data-side', 'below')
  await page.keyboard.press('Escape'); await expect(page.getByTestId('info-popover')).toHaveCount(0)
  // pinned by click: stays after the mouse leaves; second click closes; outside click closes
  await page.getByTestId('category-header-loops').getByTestId('info-tip').click(); await page.mouse.move(5, 5); await expect(page.getByTestId('info-popover')).toBeVisible()
  await page.getByTestId('category-header-loops').getByTestId('info-tip').click(); await expect(page.getByTestId('info-popover')).toHaveCount(0)
  await page.getByTestId('category-header-loops').getByTestId('info-tip').click(); await page.mouse.click(700, 60); await expect(page.getByTestId('info-popover')).toHaveCount(0)
  await expect(page.getByTestId('category-loops')).toHaveAttribute('data-open', 'true')   // never toggled the category
  // near the bottom of the viewport: shrink the viewport so the SASA icon sits ~40 px above the bottom edge -> flips above
  await page.evaluate(() => window.scrollTo(0, 0))
  const iconAtTop = (await page.getByTestId('category-header-sasa').getByTestId('info-tip').boundingBox())!
  await page.setViewportSize({ width: 1400, height: Math.round(iconAtTop.y + iconAtTop.height + 40) })
  const sasaIcon = (await page.getByTestId('category-header-sasa').getByTestId('info-tip').boundingBox())!
  expect(page.viewportSize()!.height - sasaIcon.y).toBeLessThan(80)
  b = await openCat(page, 'sasa'); expect(b.y + b.height).toBeLessThanOrEqual(sasaIcon.y); await expect(page.getByTestId('info-popover')).toHaveAttribute('data-side', 'above')
  await page.keyboard.press('Escape'); await page.setViewportSize({ width: 1400, height: 1000 })
  // resizing repositions: a popover opened near the right edge moves left when the viewport narrows
  b = await openCat(page, 'modifications'); const xWide = b.x
  await page.setViewportSize({ width: 900, height: 1000 }); await page.waitForTimeout(100)
  const bNarrow = await inside(page); expect(bNarrow.x).toBeLessThan(xWide); await adjacent(page, 'modifications'); await page.keyboard.press('Escape')
  await page.setViewportSize({ width: 1400, height: 1000 })
  // scrolling repositions: in a short viewport the page scrolls; the popover follows its icon
  await page.setViewportSize({ width: 1400, height: 520 }); await page.evaluate(() => window.scrollTo(0, 0))
  b = await openCat(page, 'loops')
  await page.evaluate(() => window.scrollBy(0, 150)); await page.waitForTimeout(150)
  expect(await page.evaluate(() => window.scrollY)).toBeGreaterThan(100)
  await inside(page); await adjacent(page, 'loops')                 // still attached to the (moved) icon after the page scrolled
  await page.keyboard.press('Escape'); await page.setViewportSize({ width: 1400, height: 1000 }); await page.evaluate(() => window.scrollTo(0, 0))
  // collapsed category: icon still works
  await page.getByTestId('category-toggle-modifications').click(); await expect(page.getByTestId('category-modifications')).toHaveAttribute('data-open', 'false')
  await openCat(page, 'modifications'); await page.keyboard.press('Escape'); await page.getByTestId('category-toggle-modifications').click()
  // the Ligands tab uses the same category headers: pocket text, collapsed pocket category
  await page.getByTestId('tab-ligands').click(); await expect(page.getByTestId('category-pockets')).toBeVisible()
  await page.getByTestId('category-header-pockets').getByTestId('info-tip').hover()
  await expect(page.getByTestId('info-popover')).toHaveText('Computationally predicted RNA binding pockets identified using fpocketR. These annotations are independent of the PDBe ligand-contact data.')
  await page.mouse.move(5, 5)
  await page.getByTestId('category-toggle-pockets').click(); await expect(page.getByTestId('category-pockets')).toHaveAttribute('data-open', 'false')
  await openCat(page, 'pockets'); await page.keyboard.press('Escape'); await page.getByTestId('category-toggle-pockets').click()
  // right edge: the icons are at the right of the headers; popover is clamped inside
  await page.getByTestId('category-header-ligands').scrollIntoViewIfNeeded()
  b = await openCat(page, 'ligands'); const vp = page.viewportSize()!; expect(b.x + b.width).toBeLessThanOrEqual(vp.width - PAD + 0.5); await page.keyboard.press('Escape')
})

test('narrow and short viewport: popover fits with internal scrolling and viewport padding', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 500 })
  await page.goto('/RF00167'); await page.getByTestId('tab-ligands').click(); await expect(page.getByTestId('category-pockets')).toBeVisible()
  for (const [tab, ids] of [['ligands', ['ligands', 'pockets']], ['annotations', ['sasa', 'loops', 'modifications']]] as const) {
    await page.getByTestId(`tab-${tab}`).click(); await expect(page.getByTestId(`category-${ids[0]}`)).toBeVisible()
    for (const id of ids) {
      await page.getByTestId(`category-header-${id}`).scrollIntoViewIfNeeded()
      const b = await openCat(page, id); expect(b.width).toBeLessThanOrEqual(390 - 2 * PAD + 1)
      await page.keyboard.press('Escape'); await expect(page.getByTestId('info-popover')).toHaveCount(0)
    }
  }
  // keyboard: focus opens, Escape closes
  await page.getByTestId('category-header-loops').getByTestId('info-tip').focus(); await expect(page.getByTestId('info-popover')).toBeVisible(); await inside(page)
  await page.keyboard.press('Escape'); await expect(page.getByTestId('info-popover')).toHaveCount(0)
})
