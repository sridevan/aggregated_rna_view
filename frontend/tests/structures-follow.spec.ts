import { expect, test, type Page } from '@playwright/test'

/* Structures tab: selecting an instance (table row, heatmap/dendrogram leaf) brings it into the 3D alignment view and emphasises it.
   Groups come from structure_clustering.json (discrepancy: group 1 = 72 instances, rep shown alone by default; group 2 = 4 instances, all shown). */
type G = { group_id: number; representative: string; members: string[] }
const viewer = (page: Page) => page.getByTestId('alignment-viewer')
const ready = async (page: Page) => expect(viewer(page)).toHaveAttribute('data-status', 'ready', { timeout: 90_000 })
const shown = async (page: Page) => ((await viewer(page).getAttribute('data-structures')) ?? '').split(',').filter(Boolean)
const clickRow = async (page: Page, sn: string) => {
  await page.getByTestId('instance-search').fill(sn.split('_')[0])
  await page.getByTestId(`row-${sn}`).locator('[data-testid="cell-title"]').click()                 // not the PDB / chain links
  await page.getByTestId('instance-search').fill('')
}
let groups: G[] = [], excluded: string[] = []
const open = async (page: Page) => {
  const d = await (await page.request.get('/families/RF00167/structure_clustering.json')).json()
  groups = d.metrics.discrepancy.clustering.groups; excluded = d.dataset.excluded.map((e: { seq_name: string }) => e.seq_name)
  await page.goto('/RF00167'); await page.getByTestId('tab-structures').click(); await ready(page)
}

test('a member of the large group is shown with its representative and emphasised; a further member is added; a displayed one is only re-emphasised', async ({ page }) => {
  await open(page)
  const g1 = groups[0], [a, b] = g1.members.filter((m) => m !== g1.representative)
  expect(await shown(page)).toEqual([g1.representative])                                              // default: representative only (72 members)
  await clickRow(page, a)
  await expect(viewer(page)).toHaveAttribute('data-structures', `${g1.representative},${a}`); await ready(page)
  await expect(viewer(page)).toHaveAttribute('data-emphasis', a)
  await expect(page.getByTestId(`legend-${a}`)).toHaveAttribute('data-selected', 'true')
  await expect(page.getByTestId('structures-current')).toContainText(a.split('_')[0].toUpperCase())
  await clickRow(page, b)
  await expect(viewer(page)).toHaveAttribute('data-structures', `${g1.representative},${a},${b}`); await expect(viewer(page)).toHaveAttribute('data-emphasis', b)
  await clickRow(page, a)                                                                              // already displayed: same structures, emphasis moves
  await expect(viewer(page)).toHaveAttribute('data-emphasis', a); expect(await shown(page)).toEqual([g1.representative, a, b])
  // the emphasis is a restyle, not a reload: the viewer stays ready
  await expect(viewer(page)).toHaveAttribute('data-status', 'ready')
  // at most 5 per group: two more fill the overlay; the add control says so; a sixth is shown with the representative only
  const [c3, c4, c6] = g1.members.filter((x) => ![g1.representative, a, b].includes(x))
  await clickRow(page, c3); await clickRow(page, c4)
  await expect(viewer(page)).toHaveAttribute('data-structures', [g1.representative, a, b, c3, c4].join(','))
  await expect(page.getByTestId('add-member')).toBeDisabled(); await expect(page.getByTestId('add-member').locator('option').first()).toHaveText('Max 5 per group')
  await clickRow(page, c6)
  await expect(viewer(page)).toHaveAttribute('data-structures', `${g1.representative},${c6}`); await expect(viewer(page)).toHaveAttribute('data-emphasis', c6)
})

test('a member of another group switches to that group (all four shown), leaves a comparison; a non-clustered instance is shown on its own', async ({ page }) => {
  await open(page)
  const g2 = groups[1], m = g2.members.find((x) => x !== g2.representative)!
  await page.getByTestId('compare-on').click(); await ready(page)
  await clickRow(page, m)
  await expect(page.getByTestId('compare-off')).toHaveCount(0)                                         // comparison left
  await expect.poll(async () => (await shown(page)).sort()).toEqual([...g2.members].sort()); await ready(page)
  await expect(viewer(page)).toHaveAttribute('data-emphasis', m)
  await expect(page.getByTestId(`group-badge-${m}`)).toBeVisible()
  // not clustered: alone, explained, no emphasis needed
  const x = excluded[0]
  await clickRow(page, x)
  await expect(viewer(page)).toHaveAttribute('data-structures', x); await ready(page)
  await expect(page.getByTestId('alignment-caption')).toHaveText(`${x.split('_')[0].toUpperCase()} ${x.split('_')[2]} is not clustered; shown on its own`)
  await expect(page.getByTestId('selected-not-clustered')).toBeVisible()
  await expect(viewer(page)).toHaveAttribute('data-emphasis', '')
  // picking a group again restores that group's default view (Structural groups panel)
  await page.getByTestId('instance-search').fill(groups[0].representative.split('_')[0]); await page.getByTestId(`group-badge-${groups[0].representative}`).click()
  await expect(viewer(page)).toHaveAttribute('data-structures', groups[0].representative)
})
