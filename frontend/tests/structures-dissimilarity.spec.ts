import { expect, test, type Page } from '@playwright/test'

/* Structures tab, Structural Dissimilarity mode: table and heatmap in the stored TPPL order (the table can be re-sorted and restored),
   selection by seq_name shared by table, heatmap and 3D view (one RNA or exactly one pair), hover details, the 50-instance label rule,
   Reset, and mode switching without refetching or rebuilding. Expected values come from the data files. */
type Clu = { dataset: { instances: string[]; n_shared_positions: number }; metrics: { discrepancy: { values_upper: number[]; seriation: { order: string[]; seed: number; repetitions: number } } } }
type Row = { seq_name: string; title: string | null; experimental_method: string | null; resolution: number | null; organism: string | null; organism_species?: string | null }
const J = async (page: Page, acc: string, file: string) => (await page.request.get(`/families/${acc}/${file}`)).json()
const label = (sn: string) => { const [p, , a] = sn.split('_'); return `${p.toUpperCase()} ${a}` }
const viewer = (page: Page) => page.getByTestId('alignment-viewer')
const ready = (page: Page) => expect(viewer(page)).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
/** stored discrepancy of a pair of seq_names (upper triangle in dataset order) */
const valueOf = (d: Clu) => { const ids = d.dataset.instances, n = ids.length; return (a: string, b: string) => { let i = ids.indexOf(a), j = ids.indexOf(b); if (i > j) [i, j] = [j, i]; return d.metrics.discrepancy.values_upper[i * n - (i * (i + 1)) / 2 + (j - i - 1)] } }
async function openStructures(page: Page, acc: string) {
  await page.goto(`/${acc}`); await page.getByTestId('tab-structures').click(); await ready(page)
}

/** the first instance of the stored TPPL order: the default selection */
const firstOf = async (page: Page, acc: string) => ((await J(page, acc, 'structure_clustering.json')) as Clu).metrics.discrepancy.seriation.order[0]
const diag = (page: Page, sn: string, order: string[]) => page.getByTestId(`dm-cell-${order.indexOf(sn)}-${order.indexOf(sn)}`)
const allRows = async (page: Page) => { const out: string[] = []
  while (!(await page.getByTestId('page-prev').isDisabled())) await page.getByTestId('page-prev').click()
  for (;;) { out.push(...await rowIds(page)); if (await page.getByTestId('page-next').isDisabled()) break; await page.getByTestId('page-next').click() }
  return out }
const rowIds = (page: Page) => page.getByTestId('structures-section-a').locator('tbody tr').evaluateAll((rs) => rs.map((r) => r.getAttribute('data-testid')!.slice(4)))
const shown = (page: Page) => viewer(page).getAttribute('data-structures')
const checks = (page: Page) => page.locator('[data-testid^="dm-diag-check-"]').evaluateAll((gs) => gs.map((g) => Number(g.getAttribute('data-testid')!.slice(14))))

test('modes: Structural Dissimilarity first and default; each mode keeps its own one- or two-structure state; no stale structures, refetch or rebuild', async ({ page }) => {
  let fetches = 0
  page.on('request', (r) => { if (r.url().includes('structure_clustering.json')) fetches++ })
  const rep = await firstOf(page, 'RF00167')
  await openStructures(page, 'RF00167')
  await expect(page.getByTestId('view-dissimilarity')).toHaveAttribute('aria-checked', 'true')
  expect(await page.getByTestId('view-selector').getByRole('radio').allInnerTexts()).toEqual(['Structural Dissimilarity', 'Clustering'])
  await expect(page.getByTestId('dissimilarity-heatmap')).toBeVisible()
  for (const id of ['groups-card', 'pairwise-panel', 'metric-selector', 'table-mode-toggle', 'associated-features', 'compare-box', 'add-member', 'reset-displayed']) await expect(page.getByTestId(id)).toHaveCount(0)
  await expect(page.locator('[data-testid^="visible-"], [data-testid^="remove-"]')).toHaveCount(0)        // no visibility or remove controls
  await expect(page.getByTestId('structures-section-a').locator('thead')).not.toContainText('Structural group')
  const wideW = (await page.getByTestId('structures-section-a').locator('table').boundingBox())!.width
  const builds = async () => page.evaluate(() => (window.__rnaExplorer as Record<string, number>).dissimilarityBuilds)
  const b0 = await builds(), f0 = fetches

  // Clustering: its own state (here a representative comparison), unchanged by visits to Structural Dissimilarity
  await page.getByTestId('view-clustering').click()
  await expect(page.getByTestId('groups-card')).toBeVisible(); await expect(page.getByTestId('pairwise-panel')).toBeVisible(); await expect(page.getByTestId('compare-box')).toBeVisible()
  await expect(page.getByTestId('alignment-title')).toHaveText('3D structure')
  expect((await page.getByTestId('structures-section-a').locator('table').boundingBox())!.width).toBeLessThan(wideW - 200)   // compact table beside the groups panel
  await page.getByTestId('compare-on').click(); await ready(page); const pairC = (await shown(page))!
  expect(pairC.split(',').length).toBe(2)
  // back in Structural Dissimilarity: the representative alone (nothing carried over); back again: the comparison is intact
  await page.getByTestId('view-dissimilarity').click(); await expect(viewer(page)).toHaveAttribute('data-structures', rep)
  await expect(page.getByTestId('compare-box')).toHaveCount(0)
  await page.getByTestId('view-clustering').click(); await expect(viewer(page)).toHaveAttribute('data-structures', pairC)
  await page.getByTestId('compare-off').click(); await expect(viewer(page)).toHaveAttribute('data-structures', pairC.split(',')[0])
  await page.getByTestId('view-dissimilarity').click(); await expect(page.getByTestId('dissimilarity-heatmap')).toBeVisible()
  expect(fetches).toBe(f0); expect(await builds()).toBe(b0)
  // the mode survives a visit to another tab
  await page.getByTestId('tab-annotations').click(); await page.getByTestId('tab-structures').click()
  await expect(page.getByTestId('view-dissimilarity')).toHaveAttribute('aria-checked', 'true'); await expect(viewer(page)).toHaveAttribute('data-structures', rep)
})

test('default and sorting: table and heatmap share the TPPL order; the first instance is selected and shown; column sorts change only the table; Seriation restores it', async ({ page }) => {
  test.setTimeout(300_000)   // walks every table page several times
  const d: Clu = await J(page, 'RF00167', 'structure_clustering.json'), order = d.metrics.discrepancy.seriation.order
  const inst: (Row & { pdb_id: string; auth_asym_id: string })[] = await J(page, 'RF00167', 'instances.json')
  const builds = async () => page.evaluate(() => (window.__rnaExplorer as Record<string, number>).dissimilarityBuilds)
  await openStructures(page, 'RF00167')
  const b0 = await builds()
  // default: row 1 = diagonal 1 = the 3D content = order[0]
  await expect(page.getByTestId('order-by')).toHaveValue('seriation'); await expect(page.getByTestId('sort-note')).toHaveCount(0)
  await expect(page.getByTestId('page-info')).toHaveText(/^Page 1 of/)
  expect((await rowIds(page))[0]).toBe(order[0]); await expect(page.getByTestId(`row-${order[0]}`)).toHaveAttribute('data-selected', 'true')
  expect(await checks(page)).toEqual([0])
  await expect(viewer(page)).toHaveAttribute('data-structures', order[0]); await expect(page.getByTestId('alignment-title')).toHaveText('3D structure')
  await expect(page.getByTestId('structures-current')).toContainText(order[0].split('_')[0].toUpperCase())
  // the whole table in seriation order, then the instances without a seriation position (not clustered), by PDB ID and chain; none twice
  const all = await allRows(page)
  const rest = inst.map((r) => r).filter((r) => !order.includes(r.seq_name)).sort((a, b) => a.pdb_id.toUpperCase().localeCompare(b.pdb_id.toUpperCase()) || a.auth_asym_id.localeCompare(b.auth_asym_id)).map((r) => r.seq_name)
  expect(all).toEqual([...order, ...rest]); expect(new Set(all).size).toBe(inst.length)
  for (const sn of rest.slice(0, 2)) { await page.getByTestId('instance-search').fill(sn.split('_')[0]); await expect(page.getByTestId(`row-${sn}`).getByTestId('not-in-matrix')).toBeVisible() }
  await page.getByTestId('instance-search').fill('')
  // row 5 (seriation order) -> diagonal 5
  while (!(await page.getByTestId('page-prev').isDisabled())) await page.getByTestId('page-prev').click()
  const r5 = (await rowIds(page))[4]; expect(r5).toBe(order[4])
  await page.getByTestId(`row-${r5}`).locator('[data-testid="cell-title"]').click()
  expect(await checks(page)).toEqual([4]); await expect(viewer(page)).toHaveAttribute('data-structures', r5)
  // sort by PDB ID, then by resolution: only the table changes; the heatmap order, values and selection stay
  const heat = async () => page.locator('[data-testid^="dm-cell-"]').evaluateAll((rs) => rs.slice(0, 400).map((r) => r.getAttribute('data-value')).join())
  const h0 = await heat()
  await page.getByTestId('instance-table').getByRole('button', { name: /^PDB ID/ }).click()
  await expect(page.getByTestId('order-by')).toHaveValue('column'); await expect(page.getByTestId('order-by')).toContainText('PDB ID ▲')
  await expect(page.getByTestId('sort-note')).toHaveText('Table sorted by PDB ID · Heatmap ordered by structural similarity')
  await expect(page.getByTestId('instance-table').getByRole('button', { name: /^PDB ID ▲/ })).toBeVisible()
  const byPdb = await allRows(page)
  expect(byPdb.map((sn) => sn.split('_')[0])).toEqual([...byPdb.map((sn) => sn.split('_')[0])].sort())
  expect(await heat()).toBe(h0); expect(await checks(page)).toEqual([4]); await expect(viewer(page)).toHaveAttribute('data-structures', r5)
  await page.getByTestId('instance-table').getByRole('button', { name: /^Res\./ }).click()
  await expect(page.getByTestId('sort-note')).toHaveText('Table sorted by Res. (Å) · Heatmap ordered by structural similarity')
  expect(await heat()).toBe(h0)
  // selecting after a column sort: the right diagonal (by identity, not position)
  const pick = (await rowIds(page))[2]
  await page.getByTestId(`row-${pick}`).locator('[data-testid="cell-title"]').click()
  if (order.includes(pick)) expect(await checks(page)).toEqual([order.indexOf(pick)]); else expect(await checks(page)).toEqual([])
  await expect(viewer(page)).toHaveAttribute('data-structures', pick)
  // filtering keeps the sort mode; paging keeps the selection and the 3D content
  await page.getByTestId('instance-search').fill('9'); await expect(page.getByTestId('order-by')).toHaveValue('column'); await page.getByTestId('instance-search').fill('')
  await expect(page.getByTestId('order-by')).toHaveValue('column')
  await expect(page.getByTestId(`row-${pick}`)).toBeVisible()   // clearing the search returns to the selected row's page first
  await page.getByTestId(await page.getByTestId('page-next').isDisabled() ? 'page-prev' : 'page-next').click(); await expect(page.getByTestId(`row-${pick}`)).toHaveCount(0)
  await expect(viewer(page)).toHaveAttribute('data-structures', pick); expect(await checks(page)).toEqual(order.includes(pick) ? [order.indexOf(pick)] : [])
  // a pair, then Seriation (default): exact TPPL order again, the pair kept
  await page.getByTestId('dm-cell-40-10').click(); await expect(viewer(page)).toHaveAttribute('data-structures', `${order[40]},${order[10]}`)
  await page.getByTestId('order-by').selectOption('seriation')
  await expect(page.getByTestId('order-by')).toHaveValue('seriation'); await expect(page.getByTestId('sort-note')).toHaveCount(0)
  expect(await allRows(page)).toEqual([...order, ...rest])
  await expect(viewer(page)).toHaveAttribute('data-structures', `${order[40]},${order[10]}`); await expect(page.getByTestId('dm-selection')).toHaveAttribute('data-row', '40')
  expect(await heat()).toBe(h0)
  // no rebuild of the matrix (and no TPPL in the browser at all) through sorting, paging, selection
  expect(await builds()).toBe(b0)
})

test('matrix: stored TPPL order on rows and columns, stored values in every checked cell, symmetric, grey diagonal; labels hidden above 50 instances', async ({ page }) => {
  const d: Clu = await J(page, 'RF00167', 'structure_clustering.json'), v = valueOf(d)
  const order = d.metrics.discrepancy.seriation.order
  expect(d.metrics.discrepancy.seriation).toMatchObject({ seed: 39873, repetitions: 10 })
  expect([...order].sort()).toEqual([...d.dataset.instances].sort())
  await openStructures(page, 'RF00167'); await page.getByTestId('view-dissimilarity').click()
  const hm = page.getByTestId('dissimilarity-heatmap')
  await expect(hm).toHaveAttribute('data-n', '76'); await expect(hm).toHaveAttribute('data-labels', 'false')
  await expect(page.locator('[data-testid^="dm-row-"]')).toHaveCount(0)
  await expect(page.locator('[data-testid^="dm-cell-"]')).toHaveCount(76 * 76)
  const cells = await page.locator('[data-testid^="dm-cell-"]').evaluateAll((rs) => rs.map((r) => [r.getAttribute('data-testid'), r.getAttribute('data-value')]))
  const val = new Map(cells as [string, string][])
  for (let r = 0; r < 76; r += 5) for (let c = 0; c < 76; c += 3) {
    const got = val.get(`dm-cell-${r}-${c}`)!
    if (r === c) expect(got).toBe('')
    else { expect(Number(got)).toBeCloseTo(v(order[r], order[c]), 6); expect(got).toBe(val.get(`dm-cell-${c}-${r}`)) }
  }
  expect(await page.getByTestId('dm-cell-7-7').getAttribute('fill')).toBe('#d7dae0')
})

test('labels: shown for 50 or fewer instances (RF00162, in seriation order), hidden for 60 (RF00005)', async ({ page }) => {
  const d: Clu = await J(page, 'RF00162', 'structure_clustering.json')
  await openStructures(page, 'RF00162'); await page.getByTestId('view-dissimilarity').click()
  await expect(page.getByTestId('dissimilarity-heatmap')).toHaveAttribute('data-labels', 'true')
  const rows = await page.locator('[data-testid^="dm-row-"]').evaluateAll((ts) => ts.map((t) => t.getAttribute('data-testid')!.slice(7)))
  expect(rows).toEqual(d.metrics.discrepancy.seriation.order)
  await expect(page.getByTestId(`dm-row-${rows[0]}`)).toHaveText(label(rows[0]))
  await openStructures(page, 'RF00005'); await page.getByTestId('view-dissimilarity').click()
  await expect(page.getByTestId('dissimilarity-heatmap')).toHaveAttribute('data-n', '60')
  await expect(page.getByTestId('dissimilarity-heatmap')).toHaveAttribute('data-labels', 'false'); await expect(page.locator('[data-testid^="dm-row-"]')).toHaveCount(0)
})

test('selection sync: hover details; diagonal and table rows show one RNA; below the diagonal shows exactly that pair (replacing any previous one); upper triangle does nothing; Reset keeps sort and page', async ({ page }) => {
  const d: Clu = await J(page, 'RF00167', 'structure_clustering.json'), v = valueOf(d), order = d.metrics.discrepancy.seriation.order
  const inst: Row[] = await J(page, 'RF00167', 'instances.json')
  await openStructures(page, 'RF00167')
  const tip = page.getByTestId('dm-tooltip')

  // hover (below and above the diagonal; the diagonal gives the instance metadata)
  await page.getByTestId('dm-cell-40-10').hover()
  await expect(page.getByTestId('dm-tooltip-a')).toHaveText(label(order[40])); await expect(page.getByTestId('dm-tooltip-b')).toHaveText(label(order[10]))
  await expect(page.getByTestId('dm-tooltip-value')).toHaveText(`Geometric discrepancy ${v(order[40], order[10]).toFixed(3)}`)
  await expect(tip).toContainText(`${d.dataset.n_shared_positions} shared Rfam positions compared`); await expect(tip).toContainText('Click to superpose')
  await page.getByTestId('dm-cell-10-40').hover()
  await expect(page.getByTestId('dm-tooltip-value')).toHaveText(`Geometric discrepancy ${v(order[10], order[40]).toFixed(3)}`); await expect(tip).not.toContainText('Click to superpose')
  const row5 = inst.find((r) => r.seq_name === order[5])!
  await page.getByTestId('dm-cell-5-5').hover()
  await expect(page.getByTestId('dm-tooltip-id')).toHaveText(label(order[5]))
  if (row5.title) await expect(page.getByTestId('dm-tooltip-title')).toHaveText(row5.title)
  if (row5.experimental_method) await expect(page.getByTestId('dm-tooltip-method')).toContainText(row5.experimental_method)
  if (row5.resolution != null) await expect(page.getByTestId('dm-tooltip-resolution')).toHaveText(`Resolution: ${row5.resolution} Å`)
  await expect(page.getByTestId('dm-tooltip-organism')).toContainText((row5.organism_species ?? row5.organism)!)

  // diagonal click: that RNA alone, its table row selected
  await page.getByTestId('dm-cell-5-5').click()
  await expect(viewer(page)).toHaveAttribute('data-structures', order[5]); expect(await checks(page)).toEqual([5])
  await expect(page.getByTestId(`row-${order[5]}`)).toHaveAttribute('data-selected', 'true')
  // upper triangle: nothing
  await page.getByTestId('dm-cell-10-40').click()
  await expect(viewer(page)).toHaveAttribute('data-structures', order[5]); await expect(page.getByTestId('dm-selection')).toHaveCount(0)

  // below the diagonal: exactly that pair, superposed in the one viewer, with the stored value
  const pairState = async (r: number, c: number) => {
    const [a, b] = [order[r], order[c]]
    await expect(viewer(page)).toHaveAttribute('data-structures', `${a},${b}`); await ready(page)
    await expect(viewer(page)).toHaveCount(1)
    await expect(page.getByTestId('alignment-title')).toHaveText('3D superposition')
    await expect(page.getByTestId('viewer-labels')).toHaveAttribute('data-state', 'pair'); await expect(page.getByTestId('viewer-labels')).toContainText('Comparing:')
    await expect(page.getByTestId('viewer-labels')).toContainText('↔'); await expect(page.locator('[data-testid^="viewer-label-"]')).toHaveCount(2)
    await expect(page.getByTestId('alignment-caption')).toHaveText(`Comparing ${label(a)} and ${label(b)}`)
    await expect(page.getByTestId('structures-current')).toHaveText(`Selected: ${label(a)} ↔ ${label(b)}`)
    expect([await page.getByTestId(`swatch-${a}`).getAttribute('data-color'), await page.getByTestId(`swatch-${b}`).getAttribute('data-color')].join(',')).toBe(await viewer(page).getAttribute('data-colors'))
    await expect(page.getByTestId('pair-summary')).toHaveText(`Geometric discrepancy ${v(a, b).toFixed(3)} · ${d.dataset.n_shared_positions} shared Rfam positions`)
    await expect(page.getByTestId('dm-selection')).toHaveAttribute('data-row', String(r)); await expect(page.getByTestId('dm-selection')).toHaveAttribute('data-col', String(c))
    expect((await checks(page)).sort((x, y) => x - y)).toEqual([c, r])
    await expect(page.getByTestId('table-pair-note')).toContainText(`${label(a)} ↔ ${label(b)}`)
    for (const sn of [a, b]) {   // both rows are marked selected, whichever page they are on
      await page.getByTestId('instance-search').fill(sn.split('_')[0]); await expect(page.getByTestId(`row-${sn}`)).toHaveAttribute('data-selected', 'true')
    }
    await page.getByTestId('instance-search').fill('')
  }
  await page.getByTestId('dm-cell-40-10').click(); await pairState(40, 10)
  await page.getByTestId('dm-cell-70-2').click(); await pairState(70, 2)                 // replaces, never adds
  // a table row: that structure alone, the pair cleared
  await page.getByTestId('instance-search').fill(order[33].split('_')[0]); await page.getByTestId(`row-${order[33]}`).locator('[data-testid="cell-title"]').click(); await page.getByTestId('instance-search').fill('')
  await expect(viewer(page)).toHaveAttribute('data-structures', order[33]); await expect(page.getByTestId('alignment-title')).toHaveText('3D structure')
  await expect(page.getByTestId('pair-summary')).toHaveCount(0); await expect(page.getByTestId('dm-selection')).toHaveCount(0); expect(await checks(page)).toEqual([33])
  // a not-clustered row: shown alone, no diagonal
  const x = inst.map((r) => r.seq_name).find((sn) => !order.includes(sn))!
  await page.getByTestId('instance-search').fill(x.split('_')[0]); await page.getByTestId(`row-${x}`).locator('[data-testid="cell-title"]').click(); await page.getByTestId('instance-search').fill('')
  await expect(viewer(page)).toHaveAttribute('data-structures', x); expect(await checks(page)).toEqual([])

  // Reset: the first seriation instance alone, checked and selected; table sort and page untouched
  await page.getByTestId('dm-cell-70-2').click(); await expect(viewer(page)).toHaveAttribute('data-structures', `${order[70]},${order[2]}`)
  await page.getByTestId('instance-table').getByRole('button', { name: /^Res\./ }).click(); await page.getByTestId('page-next').click(); await page.getByTestId('page-next').click()
  const pageBefore = await page.getByTestId('page-info').innerText(), rowsBefore = await rowIds(page)
  await page.getByTestId('viewer-reset').click()
  await expect(viewer(page)).toHaveAttribute('data-structures', order[0]); await expect(page.getByTestId('alignment-title')).toHaveText('3D structure')
  await expect(page.getByTestId('pair-summary')).toHaveCount(0); await expect(page.getByTestId('dm-selection')).toHaveCount(0); expect(await checks(page)).toEqual([0])
  await expect(page.getByTestId('page-info')).toHaveText(pageBefore); expect(await rowIds(page)).toEqual(rowsBefore)
  await expect(page.getByTestId('order-by')).toHaveValue('column')
  await page.getByTestId('order-by').selectOption('seriation'); await expect(page.getByTestId(`row-${order[0]}`)).toHaveAttribute('data-selected', 'true')
})

test('narrow window: the matrix fits its panel, no horizontal page scroll', async ({ page }) => {
  await page.setViewportSize({ width: 760, height: 900 })
  await openStructures(page, 'RF00167'); await page.getByTestId('view-dissimilarity').click()
  const panel = (await page.getByTestId('dissimilarity-panel').boundingBox())!, svg = (await page.getByTestId('dm-svg').boundingBox())!
  expect(svg.x + svg.width).toBeLessThanOrEqual(panel.x + panel.width + 1)
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true)
})

test('info tip gives the short meaning; the ordering method (TPPL citation, repetitions, seed) is kept in "Ordering method"', async ({ page }) => {
  await openStructures(page, 'RF00167')
  const d = (await J(page, 'RF00167', 'structure_clustering.json')) as Clu
  await page.getByTestId('dissimilarity-panel').getByTestId('info-tip').first().click()
  await expect(page.getByTestId('info-popover')).toHaveText(`Pairwise geometric discrepancy for ${d.dataset.instances.length} RNA instances at ${d.dataset.n_shared_positions} shared Rfam positions; lower means more similar. Rows and columns are ordered by tree-penalized path length seriation (TPPL), independently of structural groups.`)
  await page.keyboard.press('Escape')
  const m = page.getByTestId('dissimilarity-methods')
  await expect(m).not.toHaveAttribute('open', '')
  await m.locator('summary').click()
  await expect(m).toContainText(`Tree-penalized path length seriation (TPPL; Aliyev & Zirbel 2023), ${d.metrics.discrepancy.seriation.repetitions} repetitions, seed ${d.metrics.discrepancy.seriation.seed}.`)
  await expect(m).toContainText('the order uses the dissimilarity matrix only, not the structural groups')
})
