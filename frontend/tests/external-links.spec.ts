import { expect, test } from '@playwright/test'

/* External links: Rfam accession in the header; PDB ID (entry) and Chain (macromolecule entity) in the structures table. */
test('Rfam badge and structures-table links point to Rfam / PDBe and never select the row', async ({ page, context }) => {
  // never leave the machine: block the link targets (not the PDBe data the viewers load)
  await context.route(/rfam\.org\/family\/|pdbe\.org\/[0-9a-z]{4}$|\/pdbe\/entry\/pdb\//, (r) => r.abort())
  await page.goto('/RF00167')
  const badge = page.getByTestId('family-badge')
  await expect(badge).toHaveAttribute('href', 'https://rfam.org/family/RF00167'); await expect(badge).toHaveAttribute('target', '_blank')
  await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click()
  await expect(page.getByTestId('groups-card')).toBeVisible(); await expect(page.locator('[data-testid^="group-badge-"]').first()).toBeVisible()   // rows are in their final (group) order
  const instances = await (await page.request.get('/families/RF00167/instances.json')).json()
  const rows = page.locator('[data-testid="instance-table"] tbody tr[data-testid^="row-"]')
  await expect(rows.first()).toBeVisible()
  const n = await rows.count()
  for (let i = 0; i < n; i++) {
    const sn = (await rows.nth(i).getAttribute('data-testid'))!.slice(4)
    const inst = (instances.instances ?? instances).find((x: { seq_name: string }) => x.seq_name === sn)
    await expect(page.getByTestId(`link-pdb-${sn}`)).toHaveAttribute('href', `https://pdbe.org/${inst.pdb_id}`)
    await expect(page.getByTestId(`link-pdb-${sn}`)).toHaveText(inst.pdb_id.toUpperCase())
    await expect(page.getByTestId(`link-chain-${sn}`)).toHaveAttribute('href', `https://www.ebi.ac.uk/pdbe/entry/pdb/${inst.pdb_id}?activeTab=macromolecules&entity=${inst.entity_id}`)
    await expect(page.getByTestId(`link-chain-${sn}`)).toHaveText(inst.auth_asym_id)
  }
  // clicking a link opens a new tab and leaves the selected row unchanged
  // (the selected instance, 4FE5 B by default, need not be on the page shown: check it through the "Selected:" line)
  const current = page.getByTestId('structures-current'), before = await current.textContent()
  const sn = (await rows.nth(1).getAttribute('data-testid'))!.slice(4)
  expect(before).not.toContain(`${sn.split('_')[0].toUpperCase()} chain ${sn.split('_')[2]}`)
  const [popup] = await Promise.all([context.waitForEvent('page'), page.getByTestId(`link-chain-${sn}`).click()])
  await popup.close()
  await expect(current).toHaveText(before!); await expect(page.getByTestId(`row-${sn}`)).toHaveAttribute('data-selected', 'false')
  await page.getByTestId(`link-pdb-${sn}`).focus(); const [p2] = await Promise.all([context.waitForEvent('page'), page.keyboard.press('Enter')]); await p2.close()
  await expect(current).toHaveText(before!); await expect(page.getByTestId(`row-${sn}`)).toHaveAttribute('data-selected', 'false')
})
