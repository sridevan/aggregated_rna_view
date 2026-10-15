import { expect, test } from '@playwright/test'
import { mkdirSync } from 'node:fs'

/** Captures review screenshots of the three tabs (also asserts the white Mol* background through the renderer). */
const OUT = '../docs/screenshots'
type W = { __rnaExplorer: { molstar?: { plugin: { canvas3d?: { props: { renderer: { backgroundColor: number } } } } } } }
const bgHex = (page: import('@playwright/test').Page) => page.evaluate(() => { const c = (window as unknown as W).__rnaExplorer.molstar?.plugin.canvas3d?.props.renderer.backgroundColor; return c == null ? null : '#' + c.toString(16).padStart(6, '0') })

test('screenshots of the three tabs', async ({ page }) => {
  mkdirSync(OUT, { recursive: true })
  await page.setViewportSize({ width: 1400, height: 1000 })
  await page.goto('/')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready')
  await page.waitForTimeout(1500)
  expect(await bgHex(page)).toBe('#ffffff')
  await page.screenshot({ path: `${OUT}/m2_1_family_overview.png`, fullPage: true })
  await page.getByTestId('tab-structures').click()
  await page.getByTestId('instance-search').fill('1y26')
  await page.getByTestId('row-1y26_1_X_A').click()
  await page.getByTestId('instance-search').fill('')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', '1y26')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready')
  await page.waitForTimeout(1500)
  await page.locator('[data-testid="instance-table"] thead').getByTestId('info-tip').nth(4).hover()
  await page.screenshot({ path: `${OUT}/m2_1_structures.png`, fullPage: true })
  await page.getByTestId('tab-annotations').click()
  await page.getByTestId('instance-select').selectOption('1y26_1_X_A')
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready')
  await page.getByTestId('goto-input').fill('50'); await page.getByTestId('goto-button').click()
  await page.getByTestId('goto-input').fill('36'); await page.getByTestId('goto-button').click()
  await expect(page.getByTestId('sel-highlightable')).toContainText('50 → G34')
  await page.waitForTimeout(2000)
  expect(await bgHex(page)).toBe('#ffffff')
  await page.screenshot({ path: `${OUT}/m2_1_annotations.png`, fullPage: true })
})
