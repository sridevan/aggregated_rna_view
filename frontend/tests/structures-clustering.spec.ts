import { expect, test, type Page } from '@playwright/test'

/* RF00167 structural clustering in the Structures tab: totals and exclusions, fixed metric-specific groups and medoids, metric switching,
   the 50-instance rule with drill-down, compact/expanded table, coordinated group selection, 3D alignment and associated features. */
const open = async (page: Page) => { await page.goto('/RF00167'); await page.getByTestId('tab-structures').click(); await expect(page.getByTestId('groups-card')).toBeVisible() }
const viewerReady = (page: Page) => expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
type Data = { dataset: { instances: string[]; excluded: { seq_name: string; reason_text: string }[]; n_mapped_entries: number; n_clustered_entries: number }; metrics: Record<string, { clustering: { groups: { group_id: number; n_instances: number; n_entries: number; representative: string; members: string[]; within: { mean: number } }[]; assignments: Record<string, number>; leaf_order: number[] } }> }

test('totals, inclusion/exclusion, default discrepancy groups and medoids; all 84 rows kept', async ({ page }) => {
  await open(page)
  const d: Data = await (await page.request.get('/families/RF00167/structure_clustering.json')).json()
  const inst = await (await page.request.get('/families/RF00167/instances.json')).json()
  await expect(page.getByTestId('structures-summary')).toContainText(`84 RNA instances from 66 PDB entries · 76 clustered (61 PDB entries) · 8 not clustered`)
  expect(d.dataset.instances.length).toBe(76); expect(d.dataset.excluded.length).toBe(8); expect(inst.length).toBe(84)
  await expect(page.getByTestId('table-count')).toContainText('of 84 instances')
  // default: geometric discrepancy, largest group (72) selected, its medoid
  await expect(page.getByTestId('metric-discrepancy')).toHaveAttribute('aria-checked', 'true')
  const g = d.metrics.discrepancy.clustering.groups
  expect(g.map((x) => x.n_instances)).toEqual([72, 4])
  await expect(page.getByTestId('group-chip-1')).toHaveAttribute('aria-pressed', 'true')
  await expect(page.getByTestId('selected-group-counts')).toHaveText(`72 RNA instances · ${g[0].n_entries} PDB entries`)
  await expect(page.getByTestId('group-representative')).toHaveText('6UC8 B')
  await expect(page.getByTestId('within-mean')).toHaveText(g[0].within.mean.toFixed(2))
  // not-clustered rows stay visible with the reason, no group
  await page.getByTestId('instance-search').fill('9nbh')
  for (const sn of ['9nbh_1_A_A', '9nbh_1_B_B', '9nbh_1_C_C', '9nbh_1_D_D']) {
    await expect(page.getByTestId(`not-clustered-${sn}`)).toHaveText('Not clustered')
    await expect(page.getByTestId(`not-clustered-${sn}`)).toHaveAttribute('title', 'Not clustered: 26 of the 64 required core positions observed')
  }
  await page.getByTestId('instance-search').fill('6p2h'); await expect(page.getByTestId('not-clustered-6p2h_1_A_A')).toBeVisible()
  // badges carry the group of the assignment
  await page.getByTestId('instance-search').fill('5e54')
  await expect(page.getByTestId('group-badge-5e54_1_A_A')).toHaveAttribute('data-group', '2'); await expect(page.getByTestId('group-badge-5e54_1_B_B')).toHaveAttribute('data-group', '2')
})

test('50-instance rule: collapsed group view for 76 instances; drill-down only for groups of <= 50; group selection coordinates every panel', async ({ page }) => {
  await open(page)
  await expect(page.getByTestId('group-overview')).toBeVisible()
  await expect(page.getByTestId('pairwise-heatmap')).toHaveCount(0)                     // no 76 x 76 matrix by default
  await expect(page.getByTestId('overview-cut')).toHaveText('cut at 0.495')
  await expect(page.getByTestId('pairwise-note')).toHaveText('Groups with more than 50 instances are shown in collapsed form.')
  await expect(page.getByTestId('drilldown-toggle')).toHaveCount(0)                              // 72-member group: no instance-level matrix
  await viewerReady(page); await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-structures', '6uc8_1_B_A')   // large group: representative only
  // select group 2 from the collapsed dendrogram
  await page.getByTestId('overview-row-2').click()
  await expect(page.getByTestId('group-chip-2')).toHaveAttribute('aria-pressed', 'true')
  await expect(page.getByTestId('group-representative')).toHaveText('5E54 B')
  await page.getByTestId('instance-search').fill('5swd')
  await expect(page.getByTestId('row-5swd_1_A_A')).toHaveAttribute('data-in-group', 'true')
  await page.getByTestId('instance-search').fill('')
  await expect(page.getByTestId('associated-features')).toContainText('Structural group 2')
  await page.getByTestId('features-toggle').click()
  await expect(page.getByTestId('feature-col-g2').first()).toHaveAttribute('data-selected', 'true'); await expect(page.getByTestId('feature-col-g2').first()).toContainText('(n = 4 RNA instances)')
  // drill-down: 4 x 4 matrix with its own dendrogram, in dendrogram order
  await page.getByTestId('drilldown-toggle').click()
  const hm = page.getByTestId('drilldown-view')
  await expect(hm.locator('rect[data-testid^="hm-cell-"]')).toHaveCount(16)
  await expect(hm.getByTestId('instance-dendrogram')).toBeVisible()
  const labels = await hm.locator('text[data-testid^="hm-row-"]').evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')!.replace('hm-row-', '')))
  expect(labels.sort()).toEqual(['5e54_1_A_A', '5e54_1_B_B', '5swd_1_A_A', '5swd_1_B_B'])
  expect(await hm.getByTestId('hm-cell-0-0').getAttribute('data-value')).toBe('')                        // diagonal
  const d: Data = await (await page.request.get('/families/RF00167/structure_clustering.json')).json()
  const ids = d.dataset.instances
  const order = d.metrics.discrepancy.clustering.leaf_order.map((i) => ids[i]).filter((s) => d.metrics.discrepancy.clustering.assignments[s] === 2)
  expect(labels.length).toBe(4); expect(await hm.locator('text[data-testid^="hm-row-"]').evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')!.replace('hm-row-', '')))).toEqual(order)
  // 3D: all four members, representative first, distinct colours, poor-fit warnings from the existing alignment
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-structures', ['5e54_1_B_B', ...order.filter((s) => s !== '5e54_1_B_B')].join(','))
  await viewerReady(page)
  await expect(page.getByTestId('fit-warning-5e54_1_B_B')).toHaveText('poor fit'); await expect(page.getByTestId('fit-warning-5swd_1_B_B')).toHaveText('poor fit')
  await expect(page.getByTestId('fit-warning-5e54_1_A_A')).toHaveCount(0)
  // visibility toggle hides a structure without removing it from the legend
  await page.getByTestId('visible-5swd_1_A_A').uncheck()
  await expect(page.getByTestId('alignment-viewer')).not.toHaveAttribute('data-structures', /5swd_1_A_A/)
  await expect(page.getByTestId('legend-5swd_1_A_A')).toBeVisible()
  await viewerReady(page)
  const loaded = await page.evaluate(() => (window as unknown as { __rnaExplorer: { alignment: { loaded: string[]; transformed: boolean[] } } }).__rnaExplorer.alignment)
  expect(loaded.loaded.length).toBe(3); expect(loaded.transformed.every(Boolean)).toBe(true)
})

test('metric switching: TM-score uses its own groups, 1 − TM dendrogram and medoids; the selection falls back to its largest group', async ({ page }) => {
  await open(page)
  await page.getByTestId('group-chip-2').click()
  await page.getByTestId('metric-tm').click()
  await expect(page.getByTestId('pairwise-panel')).toHaveAttribute('data-metric', 'tm')
  await expect(page.getByTestId('group-chip-1')).toHaveAttribute('aria-pressed', 'true')                // largest TM group, not "group 2" carried over
  await expect(page.getByTestId('selected-group-counts')).toContainText('51 RNA instances')
  await expect(page.getByTestId('group-representative')).toHaveText('6UC8 B')
  await expect(page.getByTestId('pairwise-caption')).toHaveText('76 RNA instances · 4 structural groups · Fixed cut: 0.327provisional')
  await expect(page.getByTestId('group-overview')).toContainText('1 − TM-score at merge')
  await expect(page.getByTestId('within-heading')).toHaveText('Within-group 1 − TM-score')
  // group 2 under TM (17 instances): drill-down matrix shows TM-score values, legend says 1 − TM-score
  await page.getByTestId('group-chip-2').click()
  await expect(page.getByTestId('group-representative')).toHaveText('9V4V A')
  await page.getByTestId('drilldown-toggle').click()
  await expect(page.getByTestId('drilldown-view').locator('rect[data-testid^="hm-cell-"]')).toHaveCount(17 * 17)
  await expect(page.getByTestId('heatmap-legend')).toContainText('More similar ← 1 − TM-score → More different')
  const v = Number(await page.getByTestId('drilldown-view').getByTestId('hm-cell-0-1').getAttribute('data-value'))
  expect(v).toBeGreaterThan(0.5); expect(v).toBeLessThanOrEqual(1)                                          // a TM-score (similarity), not 1 − TM
  // switching back restores the discrepancy selection that was made there
  await page.getByTestId('metric-discrepancy').click()
  await expect(page.getByTestId('group-chip-2')).toHaveAttribute('aria-pressed', 'true')
  await expect(page.getByTestId('group-representative')).toHaveText('5E54 B')
})

test('compact/expanded table keeps metric, group, selection, search and sort; badges select groups; selected-group filter', async ({ page }) => {
  await open(page)
  await page.getByTestId('metric-tm').click(); await page.getByTestId('group-chip-3').click()
  await page.getByTestId('instance-search').fill('xpt')
  await page.locator('[data-testid="instance-table"] thead th > button', { hasText: 'Res. (Å)' }).click()
  await expect(page.getByTestId('structures-section-a')).toHaveAttribute('data-mode', 'compact')
  const before = await page.locator('[data-testid="instance-table"] tbody tr').evaluateAll((rs) => rs.map((r) => r.getAttribute('data-testid')))
  await page.getByTestId('table-mode-toggle').click()
  await expect(page.getByTestId('structures-section-a')).toHaveAttribute('data-mode', 'expanded')
  await expect(page.getByTestId('groups-card')).toHaveCount(0)
  await expect(page.locator('[data-testid="instance-table"] thead')).toContainText('Length (nt)'); await expect(page.locator('[data-testid="instance-table"] thead')).toContainText('Seq. cov. %')
  expect(await page.locator('[data-testid="instance-table"] tbody tr').evaluateAll((rs) => rs.map((r) => r.getAttribute('data-testid')))).toEqual(before)   // search + sort kept
  await expect(page.getByTestId('instance-search')).toHaveValue('xpt')
  await expect(page.getByTestId('pairwise-panel')).toBeVisible(); await expect(page.getByTestId('alignment-panel')).toBeVisible(); await expect(page.getByTestId('associated-features')).toBeVisible()
  // badge selects a group while the panel is hidden
  const badge = page.locator('[data-testid^="group-badge-"]').first(); const gid = await badge.getAttribute('data-group')
  await badge.click()
  await page.getByTestId('table-mode-toggle').click()
  await expect(page.getByTestId(`group-chip-${gid}`)).toHaveAttribute('aria-pressed', 'true')
  await expect(page.getByTestId('pairwise-panel')).toHaveAttribute('data-metric', 'tm')
  // selected group only (shown above 50 clustered instances)
  await page.getByTestId('instance-search').fill('')
  await page.getByTestId('metric-discrepancy').click(); await page.getByTestId('group-chip-2').click()
  await page.getByTestId('group-only').check()
  await expect(page.getByTestId('table-count')).toHaveText('Showing 1–4 of 4 instances (filtered from 84)')
  // selecting an individual instance highlights it and offers it to the 3D view
  await page.getByTestId('row-5swd_1_A_A').click()
  await expect(page.getByTestId('row-5swd_1_A_A')).toHaveAttribute('aria-selected', 'true')
  await expect(page.getByTestId('current-group')).toHaveText('Structural group 2')
  await expect(page.getByTestId('overview-instance-mark')).toHaveCount(1)
})

test('representative comparison, member adding with a cap, not-clustered selection, associated features with counts', async ({ page }) => {
  await open(page)
  await viewerReady(page)
  await page.getByTestId('compare-on').click()
  await expect(page.getByTestId('alignment-caption')).toHaveText('Representatives of Structural group 1 and Structural group 2')
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-structures', '6uc8_1_B_A,5e54_1_B_B')
  await viewerReady(page)
  await page.getByTestId('compare-off').click()
  for (let k = 0; k < 4; k++) { const opt = page.getByTestId('add-member').locator('option').nth(1); await page.getByTestId('add-member').selectOption(await opt.getAttribute('value')) }
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-structures', /^([^,]+,){4}[^,]+$/)   // 5 = cap per structural group
  await expect(page.getByTestId('add-member')).toBeDisabled()
  await page.getByTestId('reset-displayed').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-structures', '6uc8_1_B_A')
  // a not-clustered instance is explained, never added to a group
  await page.getByTestId('instance-search').fill('9nbh'); await page.getByTestId('row-9nbh_1_A_A').click()
  await expect(page.getByTestId('current-group')).toHaveText('not clustered')
  await expect(page.getByTestId('selected-not-clustered')).toContainText('9NBH A is not clustered (26 of the 64 required core positions observed)')
  // features of group 2: species and ligands with explicit k/n on both sides
  await page.getByTestId('group-chip-2').click()
  await page.getByTestId('features-toggle').click()
  await expect(page.getByTestId('cell-species-Vibrio vulnificus-g2')).toHaveText('4/4'); await expect(page.getByTestId('cell-species-Vibrio vulnificus-g1')).toHaveCount(0)   // single-group mode
  await page.getByTestId('compare-select').selectOption('1'); await page.getByTestId('compare-on').click()
  await expect(page.getByTestId('cell-species-Vibrio vulnificus-g1')).toHaveText('12/72')                                                   // comparison mode: both groups
  await page.getByTestId('compare-off').click()
  await expect(page.getByTestId('cell-method-X-ray diffraction-g2')).toHaveText('4/4')
})

test('missing clustering file: the table and the rest of the app still work', async ({ page }) => {
  await page.route('**/families/*/structure_clustering.json', (r) => r.fulfill({ status: 404, body: 'not found' }))
  await page.goto('/RF00167'); await page.getByTestId('tab-structures').click()
  await expect(page.getByTestId('no-clustering')).toBeVisible()
  await expect(page.locator('[data-testid="instance-table"] tbody tr')).toHaveCount(10)
  await expect(page.locator('[data-testid="instance-table"] thead')).toContainText('Rfam bits')
  await expect(page.getByTestId('groups-card')).toHaveCount(0)
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)
})

test('narrow screen: panels stack without horizontal page overflow', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 900 })
  await open(page)
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)
  const t = (await page.getByTestId('instance-table').boundingBox())!, g = (await page.getByTestId('groups-card').boundingBox())!
  expect(g.y).toBeGreaterThan(t.y)
})

test('panel layout: one-line caption, short note, details in the popover; equal-height panels in both metrics and modes; viewer resizes and goes fullscreen', async ({ page }) => {
  await open(page)
  const left = page.getByTestId('pairwise-panel'), right = page.getByTestId('alignment-panel'), viewer = page.getByTestId('alignment-viewer')
  const equal = async () => { const a = (await left.boundingBox())!, b = (await right.boundingBox())!; expect(Math.abs(a.height - b.height)).toBeLessThan(1); expect(Math.abs(a.y - b.y)).toBeLessThan(1) }
  // concise text: summary above, one sentence below; methodology, floor, threshold and caveat in the popover
  await expect(page.getByTestId('pairwise-caption')).toHaveText('76 RNA instances · 2 structural groups · Fixed cut: 0.495')
  await expect(left).not.toContainText('do not necessarily correspond'); await expect(left).not.toContainText('Dotted orange line')
  await left.getByTestId('info-tip').first().click()
  const pop = page.getByTestId('info-popover')
  for (const t of ['Geometric discrepancy compares base positions and orientations', 'dotted orange line marks the spread between duplicate structures']) await expect(pop).toContainText(t)
  await expect(pop).not.toContainText('UPGMA')                                                        // methods live in Clustering details
  await page.keyboard.press('Escape')
  await viewerReady(page); await equal()
  await expect.poll(async () => { const c = await viewer.locator('canvas').first().boundingBox(); return !!c && c.width >= 300 && c.height >= 300 }).toBe(true)   // 3D canvas at least 300 x 300
  // the note sits directly below the figure (no large gap)
  const fig = (await page.getByTestId('group-overview').boundingBox())!, note = (await page.getByTestId('pairwise-note').boundingBox())!
  expect(note.y - (fig.y + fig.height)).toBeLessThan(20)
  // drill-down (left grows): panels stay equal, the viewer grows with them and its canvas follows the box
  const h0 = (await viewer.boundingBox())!.height
  await page.getByTestId('group-chip-2').click(); await page.getByTestId('drilldown-toggle').click(); await viewerReady(page)
  await expect(page.getByTestId('pairwise-note')).toHaveText('Pairwise geometric discrepancy within Structural group 2, in dendrogram order.')
  await equal()
  const vb = (await viewer.boundingBox())!; expect(vb.height).toBeGreaterThan(h0)
  await expect.poll(async () => { const c = await viewer.locator('canvas').first().boundingBox(); return c ? Math.abs(c.height - vb.height) < 3 && Math.abs(c.width - vb.width) < 3 : false }).toBe(true)
  // TM metric, collapsed and drill-down
  await page.getByTestId('metric-tm').click(); await viewerReady(page); await equal()
  await page.getByTestId('group-chip-2').click(); await page.getByTestId('drilldown-toggle').click(); await viewerReady(page); await equal()
  await expect(page.getByTestId('pairwise-note')).toHaveText('Pairwise 1 − TM-score within Structural group 2, in dendrogram order.')
  // representative comparison and displayed-structure controls keep the layout
  await page.getByTestId('compare-on').click(); await viewerReady(page); await equal()
  await expect(page.getByTestId('viewer-legend')).toBeVisible()
  const legend = (await page.getByTestId('viewer-legend').boundingBox())!, panel = (await right.boundingBox())!
  expect(legend.y + legend.height).toBeLessThanOrEqual(panel.y + panel.height + 1)                             // controls never clipped
  // fullscreen and centre keep working
  await page.getByTestId('viewer-fullscreen').click()
  await expect(page.getByTestId('viewer-fullscreen')).toHaveText('Exit fullscreen')
  await page.getByTestId('viewer-fullscreen').click(); await expect(page.getByTestId('viewer-fullscreen')).toHaveText('Fullscreen')
  await page.getByTestId('viewer-centre').click(); await expect(viewer).toHaveAttribute('data-status', 'ready')
  // narrow: panels stack and size naturally, no horizontal overflow
  await page.setViewportSize({ width: 600, height: 900 })
  await expect.poll(async () => { const a = (await left.boundingBox())!, b = (await right.boundingBox())!; return b.y >= a.y + a.height - 0.5 }).toBe(true)   // stacked once the heatmap has reflowed
  expect((await viewer.boundingBox())!.height).toBeGreaterThanOrEqual(319)
  await expect.poll(() => page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)   // after the reflow settles
})

test('table: Length (nt) replaces Region (range in the tooltip, numeric sort); short method labels in both modes', async ({ page }) => {
  await open(page)
  const inst: { seq_name: string; pdb_id: string; map_start: number; map_end: number; mapped_length: number; experimental_method: string }[] = await (await page.request.get('/families/RF00167/instances.json')).json()
  const head = page.locator('[data-testid="instance-table"] thead')
  // compact: method shortened, full value on hover
  await page.getByTestId('instance-search').fill('1y26')
  await expect(page.getByTestId('row-1y26_1_X_A').getByTestId('cell-method')).toHaveText('X-ray')
  await expect(page.getByTestId('row-1y26_1_X_A').getByTestId('cell-method')).toHaveAttribute('title', 'X-ray diffraction')
  const em = inst.find((r) => r.experimental_method === 'Electron Microscopy')!
  await page.getByTestId('instance-search').fill(em.pdb_id)
  await expect(page.getByTestId(`row-${em.seq_name}`).getByTestId('cell-method')).toHaveText('EM')
  // expanded: Length (nt) with the region range in the tooltip; no Region column
  await page.getByTestId('table-mode-toggle').click()
  await expect(head).toContainText('Length (nt)'); await expect(head).not.toContainText('Region')
  await page.getByTestId('instance-search').fill('1y26')
  const r = inst.find((x) => x.seq_name === '1y26_1_X_A')!
  await expect(page.getByTestId('row-1y26_1_X_A').getByTestId('cell-length')).toHaveText(String(r.mapped_length))
  await expect(page.getByTestId('row-1y26_1_X_A').getByTestId('cell-length')).toHaveAttribute('title', `RNA region: ${r.map_start}–${r.map_end}`)
  await expect(page.getByTestId('row-1y26_1_X_A').getByTestId('cell-method')).toHaveText('X-ray')
  // numeric sort on the length
  await page.getByTestId('instance-search').fill('')
  // sorting first jumps to the selected row's page (a deferred effect): retry until page 1 holds
  const toFirstPage = async () => { await expect(async () => { while (await page.getByTestId('page-prev').isEnabled()) await page.getByTestId('page-prev').click(); await expect(page.getByTestId('page-info')).toHaveText(/^Page 1 of/, { timeout: 500 }) }).toPass({ timeout: 15_000 }) }
  const lenBtn = head.locator('th > button', { hasText: 'Length (nt)' })
  await lenBtn.click(); await expect(lenBtn).toHaveAttribute('aria-sort', 'ascending'); await toFirstPage()
  await expect.poll(async () => { await toFirstPage(); return Number((await page.getByTestId('cell-length').allTextContents())[0]) }, { timeout: 20_000 }).toBe(Math.min(...inst.map((x) => x.mapped_length)))
  const vals = (await page.getByTestId('cell-length').allTextContents()).map(Number)
  expect(vals).toEqual([...vals].sort((a, b) => a - b))
  await lenBtn.click(); await expect(lenBtn).toHaveAttribute('aria-sort', 'descending'); await toFirstPage()
  await expect.poll(async () => { await toFirstPage(); return Number((await page.getByTestId('cell-length').allTextContents())[0]) }, { timeout: 20_000 }).toBe(Math.max(...inst.map((x) => x.mapped_length)))
  // other columns kept
  for (const h of ['PDB ID', 'Chain', 'Structural group', 'Title', 'Organism', 'Method', 'Res. (Å)', 'Seq. cov. %', 'Struct. cov. %', 'Avg PP', 'Rfam bits']) await expect(head).toContainText(h)
})

test('structural groups panel: metric-named headings, Mean/Maximum to two decimals, ratio to one decimal; TM values are 1 − TM of the stored matrix', async ({ page }) => {
  await open(page)
  const d = await (await page.request.get('/families/RF00167/structure_clustering.json')).json()
  const ids: string[] = d.dataset.instances, n = ids.length
  const square = (upper: number[], diag: number) => { const M = ids.map(() => new Array(n).fill(diag)); let k = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { M[i][j] = M[j][i] = upper[k++] } return M }
  const check = async (metric: 'discrepancy' | 'tm', name: string) => {
    const m = d.metrics[metric], M = square(m.values_upper, metric === 'tm' ? 1 : 0), dist = (a: number, b: number) => (metric === 'tm' ? 1 - M[a][b] : M[a][b])
    for (const g of m.clustering.groups) {
      await page.getByTestId(`group-chip-${g.group_id}`).click()
      await expect(page.getByTestId('within-heading')).toHaveText(`Within-group ${name}`)
      await expect(page.getByTestId('across-heading')).toHaveText(`Across-group ${name}`)
      await expect(page.getByTestId('within-stats')).not.toContainText(' D'); await expect(page.getByTestId('within-stats')).toContainText('Mean'); await expect(page.getByTestId('within-stats')).toContainText('Maximum')
      // recompute from the stored matrix (1 − TM for TM): same values the panel shows, formatted to two decimals
      const idx = g.members.map((s: string) => ids.indexOf(s)), out = ids.map((_, i) => i).filter((i) => !idx.includes(i))
      const within: number[] = []; idx.forEach((a: number, x: number) => idx.slice(x + 1).forEach((b: number) => within.push(dist(a, b))))
      const across: number[] = []; idx.forEach((a: number) => out.forEach((b) => across.push(dist(a, b))))
      const wm = within.reduce((s, v) => s + v, 0) / within.length, am = across.reduce((s, v) => s + v, 0) / across.length
      expect(Math.abs(wm - g.within.mean)).toBeLessThan(1e-5); expect(Math.abs(am - g.across.mean)).toBeLessThan(1e-5)
      await expect(page.getByTestId('within-mean')).toHaveText(wm.toFixed(2)); await expect(page.getByTestId('within-max')).toHaveText(Math.max(...within).toFixed(2))
      await expect(page.getByTestId('across-mean')).toHaveText(am.toFixed(2)); await expect(page.getByTestId('across-ratio')).toHaveText(`${(am / wm).toFixed(1)}× the within-group mean`)
    }
  }
  await check('discrepancy', 'geometric discrepancy')
  await page.getByTestId('metric-tm').click()
  await check('tm', '1 − TM-score')
})

test('tooltips: short, metric-specific wording for both metrics; details moved to Clustering details; popovers inside the viewport', async ({ page }) => {
  await open(page)
  const tip = async (label: string) => { await page.getByRole('button', { name: label }).first().click(); const p = page.getByTestId('info-popover'); await expect(p).toBeVisible()
    const b = (await p.boundingBox())!, vp = page.viewportSize()!; expect(b.x).toBeGreaterThanOrEqual(0); expect(b.x + b.width).toBeLessThanOrEqual(vp.width); expect(b.y).toBeGreaterThanOrEqual(0); expect(b.y + b.height).toBeLessThanOrEqual(vp.height)
    const t = (await p.textContent())!.replace(/\s+/g, ' ').trim(); await page.keyboard.press('Escape'); return t }
  expect(await tip('About clustering coverage')).toBe('76 of 84 RNA instances were clustered. Inclusion requires at least 64 core positions with complete base, C4′ and C3′ coordinates. Comparisons use 57 positions shared by all included instances.')
  expect(await tip('About the representative')).toBe('The medoid: the RNA instance with the lowest average dissimilarity to other members of its structural group.')
  expect(await tip('About within-group dissimilarity')).toBe('Mean pairwise geometric discrepancy between RNA instances in this group. Lower values indicate greater structural similarity.')
  expect(await tip('About across-group dissimilarity')).toBe("Mean pairwise geometric discrepancy between this group's RNA instances and all other clustered instances.")
  expect(await tip('About structural groups')).toContain('do not necessarily correspond to biological conformational states')
  expect(await tip('About the 3D alignment')).toBe('RNA structures are superposed using C4′ atoms at shared Rfam positions. This visual alignment differs from that used to calculate geometric discrepancy. Poor-quality alignments are flagged.')
  await page.getByRole('button', { name: 'About the 3D alignment' }).click()
  const pb = (await page.getByTestId('info-popover').boundingBox())!, panel = (await page.getByTestId('alignment-panel').boundingBox())!
  expect(pb.width).toBeLessThanOrEqual(340); expect(pb.x).toBeGreaterThanOrEqual(panel.x - 1)                    // compact, kept over its own panel
  await page.keyboard.press('Escape')
  for (const label of ['About within-group dissimilarity', 'About across-group dissimilarity', 'About the comparison metric', 'About the representative']) {
    const t = await tip(label); expect(t).not.toMatch(/\bD\b/); expect(t).not.toMatch(/pairs|statistical test/)
  }
  await page.getByTestId('metric-tm').click()
  expect(await tip('About within-group dissimilarity')).toBe('Mean pairwise 1 − TM-score between RNA instances in this group. Lower values indicate greater structural similarity.')
  expect(await tip('About across-group dissimilarity')).toBe("Mean pairwise 1 − TM-score between this group's RNA instances and all other clustered instances.")
  expect(await tip('About the comparison metric')).toContain('Groups and the dendrogram use 1 − TM-score')
  // the detail is kept in Clustering details
  await page.getByTestId('clustering-details').locator('summary').click()
  const det = page.getByTestId('clustering-details')
  await expect(det).toContainText('Kabsch superposition of C4′ atoms'); await expect(det).toContainText('same for both metrics')
  for (const t of ['average (UPGMA) linkage on 1 − TM-score', 'provisional', 'Duplicate floor 0.131', 'pairs', 'at most 50 RNA instances', 'do not necessarily correspond']) await expect(det).toContainText(t)
  await expect(page.getByTestId('details-excluded')).toContainText('9NBH A (26 of the 64 required core positions observed)')
})

test('representative comparison highlights the stored dendrogram route; 3D colours follow the group colours in every mode', async ({ page }) => {
  await open(page)
  const d = await (await page.request.get('/families/RF00167/structure_clustering.json')).json()
  const hex = (rgb: string) => '#' + rgb.match(/\d+/g)!.slice(0, 3).map((v) => (+v).toString(16).padStart(2, '0')).join('')
  const groupColour = async (id: number) => hex(await page.getByTestId(`group-chip-${id}`).locator('span').first().evaluate((e) => getComputedStyle(e).backgroundColor))
  const colours = async () => ((await page.getByTestId('alignment-viewer').getAttribute('data-colors')) ?? '').toLowerCase().split(',')
  // expected route = path between the two group nodes through their lowest common ancestor in the stored linkage
  const route = (metric: string, a: number, b: number) => {
    const L: number[][] = d.metrics[metric].clustering.linkage, n = d.dataset.instances.length, parent = new Map<number, number>()
    L.forEach(([x, y], k) => { parent.set(x, n + k); parent.set(y, n + k) })
    const up = (v: number) => { const p = [v]; while (parent.has(p[p.length - 1])) p.push(parent.get(p[p.length - 1])!); return p }
    const node = (id: number) => d.metrics[metric].clustering.groups.find((g: { group_id: number }) => g.group_id === id).node
    const pa = up(node(a)), pb = up(node(b)); const lca = pa.find((v) => pb.includes(v))!
    return [...pa.slice(0, pa.indexOf(lca) + 1), ...pb.slice(0, pb.indexOf(lca))].sort((x, y) => x - y).join(',')
  }
  const g1 = await groupColour(1), g2 = await groupColour(2)
  // 1. Group 1 representative alone uses the Group 1 colour (also the table badge colour)
  await viewerReady(page)
  expect(await colours()).toEqual([g1]); expect((await page.getByTestId('swatch-6uc8_1_B_A').getAttribute('data-color'))!.toLowerCase()).toBe(g1)
  await page.getByTestId('instance-search').fill('6uc8')
  expect(hex(await page.getByTestId('group-badge-6uc8_1_B_A').evaluate((e) => getComputedStyle(e).borderColor))).toBe(g1)
  await page.getByTestId('instance-search').fill('')
  await expect(page.getByTestId('overview-route')).toHaveCount(0)
  // 2. Group 2 (4 members): representative exact colour, the others distinguishable shades of it
  await page.getByTestId('group-chip-2').click(); await viewerReady(page)
  const c2 = await colours(); expect(c2[0]).toBe(g2); expect(new Set(c2).size).toBe(4)
  expect((await page.getByTestId('alignment-viewer').getAttribute('data-structures'))!.split(',')[0]).toBe('5e54_1_B_B')
  // 3-4. compare group 2 with group 1: both representatives in their group colours, the stored route highlighted
  await page.getByTestId('compare-select').selectOption('1'); await page.getByTestId('compare-on').click(); await viewerReady(page)
  expect(await colours()).toEqual([g2, g1])
  await expect(page.getByTestId('overview-route')).toHaveAttribute('data-nodes', route('discrepancy', 2, 1))
  await expect(page.getByTestId('overview-compare-mark-1')).toBeVisible(); await expect(page.getByTestId('overview-compare-mark-2')).toBeVisible()
  await expect(page.getByTestId('group-chip-2')).toHaveAttribute('aria-pressed', 'true')                   // comparison does not change the selected group
  // 5. exit: route removed, group 2 and its default display restored
  await page.getByTestId('compare-off').click(); await viewerReady(page)
  await expect(page.getByTestId('overview-route')).toHaveCount(0); await expect(page.getByTestId('group-chip-2')).toHaveAttribute('aria-pressed', 'true')
  expect((await colours())[0]).toBe(g2); expect((await colours()).length).toBe(4)
  // 9. visibility, removal and camera still work with the new colours
  await page.getByTestId('visible-5swd_1_A_A').uncheck(); await expect(page.getByTestId('alignment-viewer')).not.toHaveAttribute('data-structures', /5swd_1_A_A/)
  await page.getByTestId('remove-5e54_1_A_A').click(); await viewerReady(page); expect((await colours())[0]).toBe(g2)
  await page.getByTestId('viewer-centre').click(); await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready')
  // 6-7. TM-score: four groups, arbitrary comparisons follow the TM tree and TM group colours
  await page.getByTestId('metric-tm').click(); await viewerReady(page)
  const t1 = await groupColour(1); expect(await colours()).toEqual([t1])
  for (const [a, b] of [[1, 4], [3, 2], [2, 4]]) {
    await page.getByTestId(`group-chip-${a}`).click(); await page.getByTestId('compare-select').selectOption(String(b)); await page.getByTestId('compare-on').click(); await viewerReady(page)
    expect(await colours()).toEqual([await groupColour(a), await groupColour(b)])
    await expect(page.getByTestId('overview-route')).toHaveAttribute('data-nodes', route('tm', a, b))
    await page.getByTestId('compare-off').click()
  }
  await expect(page.getByTestId('overview-route')).toHaveCount(0)
})

test('instance-level dendrogram marks the compared representatives (drill-down shows only its own group)', async ({ page }) => {
  await open(page)
  await page.getByTestId('group-chip-2').click(); await page.getByTestId('drilldown-toggle').click()
  await page.getByTestId('compare-on').click()
  const marks = page.getByTestId('drilldown-view').getByTestId('dendro-compare-mark')
  await expect(marks).toHaveCount(1)                                                // only group 2's representative is in this subtree
  await expect(page.getByTestId('drilldown-view').getByTestId('dendro-route')).toHaveCount(0)
})

test('associated features: counts-first tables per clustered RNA instance, recomputed from the source files for both metrics', async ({ page }) => {
  await open(page)
  const d = await (await page.request.get('/families/RF00167/structure_clustering.json')).json()
  const inst: { seq_name: string; pdb_id: string; entity_id: number; organism_species: string | null; organism: string | null; experimental_method: string; resolution: number | null; chain_length: number; n_match_with_coordinates: number }[] = await (await page.request.get('/families/RF00167/instances.json')).json()
  const lig = await (await page.request.get('/families/RF00167/ligands.json')).json()
  const mods = await (await page.request.get('/families/RF00167/modifications.json')).json()
  const by = new Map(inst.map((r) => [r.seq_name, r]))
  const ligOf = (sn: string) => new Set(lig.bound_molecules.filter((b: { in_default_track: boolean; entity_chains: string[] }) => b.in_default_track && b.entity_chains.includes(sn)).map((b: { ccd_id: string }) => b.ccd_id))
  const evaluated = (sn: string) => lig.entry_entity_sources[`${by.get(sn)!.pdb_id}:${by.get(sn)!.entity_id}`]?.status === 'evaluated'
  const modOf = (sn: string, outside: boolean) => new Set(mods.observations.filter((o: { seq_name: string; outside_mapped_region: boolean }) => o.seq_name === sn && o.outside_mapped_region === outside).map((o: { ccd_id: string }) => o.ccd_id))
  const count = (members: string[], known: (s: string) => boolean, has: (s: string) => boolean) => { const kn = members.filter(known); return `${kn.filter(has).length}/${kn.length}` }
  const median = (v: number[]) => { const s = [...v].sort((a, b) => a - b), m = Math.floor(s.length / 2); return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2 }
  // neutral presentation: no arrows, no enrichment colours, no 15-point filter, one descriptive note
  const panel = page.getByTestId('associated-features')
  await page.getByTestId('features-toggle').click()
  const showAll = async () => { const more = panel.locator('button[data-testid$="-toggle"][aria-expanded="false"]:not([data-testid="features-toggle"])'); while (await more.count()) await more.first().click() }
  await expect(page.getByTestId('features-summary')).toHaveText('Comparisons are descriptive; group sizes and repeated observations should be considered when interpreting differences.')
  expect(await panel.textContent()).not.toMatch(/[▲▼]|enrich|deplet|significan/i)
  expect(await panel.locator('[class*="emerald"], [class*="orange-"]').count()).toBe(0)
  await expect(page.getByTestId('features-show-all')).toHaveCount(0)
  type G = { group_id: number; members: string[] }
  const verify = async (cols: { id: string; members: string[]; sel: boolean }[]) => {
    await expect(panel.locator('th[data-testid^="feature-col-"]').first()).toBeVisible()
    expect(await panel.locator('th[data-testid^="feature-col-"]').evaluateAll((es) => [...new Set(es.map((e) => e.getAttribute('data-testid')))])).toEqual(cols.map((c) => `feature-col-${c.id}`))   // exactly these columns, in this order
    await showAll()                                                                                       // every row visible for the recomputation
    for (const col of cols) {
      await expect(page.getByTestId(`feature-col-${col.id}`).first()).toHaveAttribute('data-selected', String(col.sel))
      await expect(page.getByTestId(`feature-col-${col.id}`).first()).toContainText(`(n = ${col.members.length} RNA instances)`)
      // ligands: entity-level, deduplicated per instance, denominators = PDBe-evaluated instances
      for (const ccd of ['ADE', 'HPA', 'GUN']) { const cell = page.getByTestId(`cell-ligands-${ccd}-${col.id}`); if (await cell.count()) await expect(cell).toHaveText(count(col.members, evaluated, (s) => ligOf(s).has(ccd))) }
      // nonstandard nucleotides split by location: CCC is in the mapped region, GTP only outside it
      const ccc = page.getByTestId(`cell-nonstandard-CCC-${col.id}`); if (await ccc.count()) await expect(ccc).toHaveText(count(col.members, () => true, (s) => modOf(s, false).has('CCC')))
      const gtpOut = page.getByTestId(`cell-nonstandard-outside-GTP-${col.id}`); if (await gtpOut.count()) await expect(gtpOut).toHaveText(count(col.members, () => true, (s) => modOf(s, true).has('GTP')))
      await expect(page.getByTestId(`cell-nonstandard-GTP-${col.id}`)).toHaveCount(0)
      const sp = col.members.map((s) => by.get(s)!.organism_species ?? by.get(s)!.organism)[0]!
      await expect(page.getByTestId(`cell-species-${sp}-${col.id}`)).toHaveText(count(col.members, () => true, (s) => (by.get(s)!.organism_species ?? by.get(s)!.organism) === sp))
      const res = col.members.map((s) => by.get(s)!.resolution).filter((v): v is number => v != null)
      await expect(page.getByTestId(`median-resolution-${col.id}`)).toHaveText(median(res).toFixed(2)); await expect(page.getByTestId(`median-resolution-${col.id}`)).toHaveAttribute('data-n', String(res.length))
      await expect(page.getByTestId(`median-length-${col.id}`)).toHaveText(median(col.members.map((s) => by.get(s)!.chain_length)).toFixed(0))
    }
    // rows are the features observed in the displayed group(s): every row has at least one supporting instance
    const totals = await panel.locator('[data-testid^="feature-row-"]').evaluateAll((trs) => trs.map((tr) => [...tr.querySelectorAll('td')].slice(1).reduce((t, td) => t + Number(td.textContent!.split('/')[0]), 0)))
    expect(totals.every((t) => t > 0)).toBe(true)
  }
  const check = async (metric: 'discrepancy' | 'tm', pairs: [number, number][]) => {
    const groups = d.metrics[metric].clustering.groups as G[], gid = (id: number) => groups.find((g) => g.group_id === id)!
    for (const sel of groups.slice(0, 2)) {                                                               // single-group mode
      await page.getByTestId(`group-chip-${sel.group_id}`).click()
      await expect(page.getByTestId('features-compare-badge')).toHaveCount(0)
      await verify([{ id: `g${sel.group_id}`, members: sel.members, sel: true }])
    }
    for (const [a, b] of pairs) {                                                                         // comparison mode = the compared pair only
      await page.getByTestId(`group-chip-${a}`).click(); await page.getByTestId('compare-select').selectOption(String(b)); await page.getByTestId('compare-on').click()
      await expect(page.getByTestId('features-compare-badge')).toHaveText(gid(b).members.length ? `Structural group ${b}` : '')
      await verify([{ id: `g${a}`, members: gid(a).members, sel: true }, { id: `g${b}`, members: gid(b).members, sel: false }])
      await page.getByTestId('compare-off').click()
      await expect(page.getByTestId('features-compare-badge')).toHaveCount(0)
      await expect(panel.locator('th[data-testid^="feature-col-"]').first()).toHaveAttribute('data-testid', `feature-col-g${a}`)
    }
  }
  await check('discrepancy', [[1, 2]])
  await page.getByTestId('metric-tm').click()
  await expect(panel).toContainText('provisional groups')
  await check('tm', [[2, 4], [3, 1]])
  // tooltips: provenance wording, inside the viewport
  await page.getByRole('button', { name: 'About Ligand associations' }).click()
  await expect(page.getByTestId('info-popover')).toContainText('PDBe entity-level annotations and do not necessarily indicate direct contact')
  const b = (await page.getByTestId('info-popover').boundingBox())!; expect(b.x + b.width).toBeLessThanOrEqual(page.viewportSize()!.width)
})

test('associated features: collapsed by default, keyboard toggle, state kept across groups/metrics/tabs; five rows per category with independent Show all / Show fewer', async ({ page }) => {
  await open(page)
  const panel = page.getByTestId('associated-features'), toggle = page.getByTestId('features-toggle')
  // collapsed: heading, selected-group badge and chevron only
  await expect(panel).toHaveAttribute('data-open', 'false'); await expect(toggle).toHaveAttribute('aria-expanded', 'false')
  await expect(page.getByTestId('features-body')).toHaveCount(0)
  await expect(page.getByTestId('features-group-badge')).toHaveText('Structural group 1')
  expect((await panel.boundingBox())!.height).toBeLessThan(60)
  expect((await panel.innerText()).replace(/\s+/g, ' ').trim()).toBe('▶ Associated features Structural group 1')
  // group / metric changes update the badge without opening
  await page.getByTestId('group-chip-2').click(); await expect(page.getByTestId('features-group-badge')).toHaveText('Structural group 2'); await expect(panel).toHaveAttribute('data-open', 'false')
  await page.getByTestId('metric-tm').click(); await expect(page.getByTestId('features-group-badge')).toHaveText('Structural group 1'); await expect(panel).toHaveAttribute('data-open', 'false')
  // keyboard opens; state kept across group, metric and tab switches
  await toggle.focus(); await page.keyboard.press('Enter')
  await expect(toggle).toHaveAttribute('aria-expanded', 'true'); await expect(page.getByTestId('features-body')).toBeVisible()
  await expect(toggle).toHaveAttribute('aria-controls', (await page.getByTestId('features-body').getAttribute('id'))!)
  await page.getByTestId('group-chip-3').click(); await expect(panel).toHaveAttribute('data-open', 'true')
  await page.getByTestId('metric-discrepancy').click(); await expect(panel).toHaveAttribute('data-open', 'true')
  await page.getByTestId('tab-overview').click(); await page.getByTestId('tab-structures').click(); await expect(page.getByTestId('associated-features')).toHaveAttribute('data-open', 'true')
  // five rows by default; more only through the category's own control; others unaffected
  const rows = (id: string) => page.getByTestId(`features-${id}`).locator('tbody tr')
  await page.getByTestId('group-chip-1').click()                                                     // discrepancy group 1 (group 2 was selected earlier in this metric)
  await expect(rows('ligands')).toHaveCount(5); await expect(page.getByTestId('features-ligands-toggle')).toHaveText(/^Show all \d+$/)
  const nLig = Number((await page.getByTestId('features-ligands-toggle').textContent())!.match(/\d+/)![0]); expect(nLig).toBeGreaterThan(5)
  await expect(rows('species')).toHaveCount(5); await expect(page.getByTestId('features-species-toggle')).toHaveText('Show all 6')
  await expect(rows('nonstandard')).toHaveCount(4); await expect(page.getByTestId('features-nonstandard-toggle')).toHaveCount(0)       // <= 5: all rows, no control
  await expect(rows('method')).toHaveCount(2); await expect(page.getByTestId('features-numeric').locator('tbody tr')).toHaveCount(3)
  const order1 = await rows('ligands').locator('td:first-child').allTextContents()
  await page.getByTestId('features-ligands-toggle').press('Enter')
  await expect(rows('ligands')).toHaveCount(nLig); await expect(page.getByTestId('features-ligands-toggle')).toHaveText('Show fewer'); await expect(page.getByTestId('features-ligands-toggle')).toHaveAttribute('aria-expanded', 'true')
  await expect(rows('species')).toHaveCount(5)                                                                                          // independent
  // deterministic order: total supporting instances (desc), then identifier; unchanged by group selection
  const all1 = await rows('ligands').locator('td:first-child').allTextContents(); expect(all1.slice(0, 5)).toEqual(order1)
  const totals = await rows('ligands').evaluateAll((trs) => trs.map((tr) => [...tr.querySelectorAll('td')].slice(1).reduce((s, td) => s + Number(td.textContent!.split('/')[0]), 0)))
  expect(totals).toEqual([...totals].sort((a, b) => b - a))
  await page.getByTestId('features-ligands-toggle').click(); await expect(rows('ligands')).toHaveCount(5)          // Show fewer
  // another group: its own features only (single-group mode), category state and panel state unchanged
  await page.getByTestId('group-chip-2').click()
  await expect(page.getByTestId('features-group-badge')).toHaveText('Structural group 2'); await expect(page.getByTestId('associated-features')).toHaveAttribute('data-open', 'true')
  await expect(rows('ligands')).toHaveCount(1); await expect(page.getByTestId('cell-ligands-ADE-g2')).toHaveText('2/4'); await expect(page.getByTestId('features-ligands-toggle')).toHaveCount(0)
  await page.getByTestId('group-chip-1').click()
  // TM-score: same rules, order from the TM dataset (totals over all clustered instances)
  await page.getByTestId('metric-tm').click(); await page.getByTestId('group-chip-1').click(); await expect(rows('ligands')).toHaveCount(5)   // TM group 1 (group 3 was selected earlier in this metric)
  const tmTotals = await rows('ligands').evaluateAll((trs) => trs.map((tr) => [...tr.querySelectorAll('td')].slice(1).reduce((s, td) => s + Number(td.textContent!.split('/')[0]), 0)))
  expect(tmTotals).toEqual([...tmTotals].sort((a, b) => b - a))
  // collapse again
  await page.getByTestId('features-toggle').click(); await expect(page.getByTestId('features-body')).toHaveCount(0)
})
