import { expect, test } from '@playwright/test'

/* The Rfam accession is in the URL: /{rfam_acc}. "/" is the home page; any other non-accession path is replaced by "/". An accession
   without generated data shows a clear message instead of an error. */
test('/ is the home page; other paths go home; /rf00167 is normalised; the family comes from the path', async ({ page }) => {
  await page.goto('/')
  await expect(page).toHaveURL(/127\.0\.0\.1:5173\/$/); await expect(page.getByTestId('home-page')).toBeVisible()
  await page.goto('/whatever/else'); await expect(page).toHaveURL(/:5173\/$/); await expect(page.getByTestId('home-page')).toBeVisible()
  await page.goto('/rf00167/?x=1#h'); await expect(page).toHaveURL(/\/RF00167\?x=1#h$/)
  await expect(page.getByTestId('family-badge')).toHaveText('RF00167')
  await page.goto('/RF00167'); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  expect(await page.evaluate(() => history.length)).toBeLessThanOrEqual(5)                       // replaceState: no extra history entries
})

test('an accession without data shows "not available" with a link to the available family', async ({ page }) => {
  await page.goto('/RF99999')
  await expect(page).toHaveURL(/\/RF99999$/)
  await expect(page.getByTestId('family-unavailable')).toContainText('Rfam family RF99999 is not available in this view.')
  await page.getByTestId('family-unavailable').getByRole('link', { name: 'RF00167' }).click()
  await expect(page).toHaveURL(/\/RF00167$/); await expect(page.getByTestId('family-badge')).toHaveText('RF00167')
})
