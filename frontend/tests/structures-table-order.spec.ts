import { expect, test, type Page } from '@playwright/test'

/* Structures table default order: the active structural group first with its representative (medoid) on the top row of page 1, then
   the other groups in dendrogram order; an ordinary, re-sortable column default. Groups from structure_clustering.json. */
type G = { group_id: number; representative: string; members: string[] }
type M = { clustering: { groups: G[]; leaf_order?: number[] } }
const rowIds = (page: Page) => page.locator('[data-testid="instance-table"] tbody tr[data-testid^="row-"]').evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')!.slice(4)))
const groupHeader = (page: Page) => page.locator('[data-testid="instance-table"] thead th').filter({ hasText: 'Structural group' }).locator('button').first()
let D: { metrics: Record<'discrepancy' | 'tm', M> }
const open = async (page: Page) => {
  D = await (await page.request.get('/families/RF00167/structure_clustering.json')).json()
  await page.goto('/RF00167'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 })
}

test('on entering the tab: page 1, the active group representative on top, then its members; default sort shown; no auto-selection', async ({ page }) => {
  await open(page)
  const g1 = D.metrics.discrepancy.clustering.groups[0]
  await expect.poll(async () => (await rowIds(page))[0]).toBe(g1.representative)
  await expect(page.getByTestId('page-info')).toHaveText(/^Page 1 of/)
  expect((await rowIds(page)).every((sn) => g1.members.includes(sn))).toBe(true)                    // 72 members fill page 1
  await expect(page.getByTestId(`row-${g1.representative}`)).toContainText('medoid')
  await expect(groupHeader(page)).toHaveAttribute('aria-sort', 'ascending')
  await expect(page.getByTestId('structures-current')).toContainText('4FE5 chain B')                 // the shared selection is not changed
})

test('another group: its representative moves to the top (page 1); TM metric uses its own groups; user sorting and returning to the default', async ({ page }) => {
  await open(page)
  const [g1, g2] = D.metrics.discrepancy.clustering.groups
  await page.getByTestId('group-chip-2').click()
  await expect.poll(async () => (await rowIds(page))[0]).toBe(g2.representative)
  const ids = await rowIds(page)
  expect(ids.slice(0, g2.members.length).sort()).toEqual([...g2.members].sort())                  // the whole small group first
  expect(g1.members.includes(ids[g2.members.length])).toBe(true)                                  // then the other group
  await expect(page.getByTestId('page-info')).toHaveText(/^Page 1 of/)
  // TM-score: its own group 2 representative on top after selecting that group
  await page.getByTestId('metric-tm').click()
  const t2 = D.metrics.tm.clustering.groups[1]
  await page.getByTestId('group-chip-2').click()
  await expect.poll(async () => (await rowIds(page))[0]).toBe(t2.representative)
  // ordinary sorting still works; clicking Structural group restores the default order
  await page.locator('[data-testid="instance-table"] thead th').filter({ hasText: 'PDB ID' }).locator('button').first().click()
  await expect.poll(async () => { const r = await rowIds(page); return r.join() === [...r].sort().join() }).toBe(true)
  await groupHeader(page).click()
  await expect(groupHeader(page)).toHaveAttribute('aria-sort', 'ascending')
  await expect.poll(async () => (await rowIds(page))[0]).toBe(t2.representative)
})

test('a row selected deep in the active group stays visible; without clustering data the table keeps PDB order', async ({ page }) => {
  await open(page)
  const g1 = D.metrics.discrepancy.clustering.groups[0]
  await page.getByTestId('page-next').click(); await page.getByTestId('page-next').click()
  const target = (await rowIds(page))[3]
  await page.getByTestId(`row-${target}`).locator('[data-testid="cell-title"]').click()
  await expect(page.getByTestId(`row-${target}`)).toHaveAttribute('data-selected', 'true')            // still on screen
  await expect(page.getByTestId('page-info')).toHaveText(/^Page 3 of/)
  expect(g1.members.includes(target)).toBe(true)
  // missing clustering file: PDB ID order as before
  await page.route('**/families/*/structure_clustering.json', (r) => r.fulfill({ status: 404, body: 'not found' }))
  await page.goto('/RF00167'); await page.getByTestId('tab-structures').click(); await expect(page.getByTestId('no-clustering')).toBeVisible()
  const r = await rowIds(page); expect(r.join()).toBe([...r].sort().join())
})
