import { expect, test, type Page } from '@playwright/test'

/** Milestone 2.1 refinements: representative display structure, table, tooltips, viewer framing, legends. Needs PDBe network access. */
const waitMolstar = async (page: Page, pdb?: string) => {
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready')
  if (pdb) await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb)
}

test('display structure follows the documented rule and loads first (4fe5/B, Rfam region, white background, no control panel)', async ({ page }) => {
  await page.goto('/')
  const fam = await (await page.request.get('/families/RF00167/family.json')).json()
  const d = fam.display_representative
  expect(d.experimental_method).toMatch(/x-ray/i); expect(d.pdbe_2d_available).toBe(true); expect(d.organism_species).not.toBe('synthetic construct')
  expect(d.n_match_with_coordinates / d.n_match).toBeGreaterThanOrEqual(0.95); expect(d.fraction_clen_matched).toBeGreaterThanOrEqual(d.min_clen_fraction_used)
  await expect(page.getByTestId('display-representative')).toContainText(d.pdb_id.toUpperCase())
  await expect(page.getByTestId('display-representative').locator('.italic')).toHaveText(d.organism_species)
  await expect(page.getByTestId('representative')).toHaveCount(0)   // the backend coverage representative is no longer shown as a second card
  await expect(page.getByTestId('family-description')).toContainText(fam.comment.slice(0, 60))
  // banner: subordinate app title, family description as heading, accession badge linking to Rfam, no repeated short id
  await expect(page.getByTestId('app-title')).toHaveText('Aggregated RNA 3D View')
  await expect(page.getByTestId('family-heading')).toHaveText(fam.description)
  await expect(page.getByTestId('family-badge')).toHaveText(fam.rfam_acc)
  await expect(page.getByTestId('family-badge')).toHaveAttribute('href', fam.rfam_url)
  await expect(page.getByTestId('app-header')).not.toContainText(`${fam.rfam_id} ·`)
  expect(await page.getByTestId('family-heading').evaluate((e) => parseFloat(getComputedStyle(e).fontSize))).toBeGreaterThan(await page.getByTestId('app-title').evaluate((e) => parseFloat(getComputedStyle(e).fontSize)))
  await expect(page.getByTestId('rfam-link')).toHaveAttribute('href', fam.rfam_url)
  await expect(page.getByTestId('stat-sequence-occupancy')).toBeVisible(); await expect(page.getByTestId('stat-structural-occupancy')).toBeVisible()
  // species chart: values match instances.json (entry deduplication), italics, synthetic apart, toggle switches the unit
  const inst = await (await page.request.get('/families/RF00167/instances.json')).json() as { pdb_id: string; organism_species: string | null; species_tax_id: string | null }[]
  const nSpecies = new Set(inst.filter((r) => r.species_tax_id && r.organism_species !== 'synthetic construct').map((r) => r.species_tax_id)).size
  const entriesWithSpecies = new Set(inst.filter((r) => r.species_tax_id && r.organism_species !== 'synthetic construct').map((r) => r.pdb_id)).size
  await expect(page.getByTestId('species-summary')).toContainText(`${nSpecies} distinct identified species`)
  await expect(page.getByTestId('species-summary')).toContainText(`${entriesWithSpecies} of ${new Set(inst.map((r) => r.pdb_id)).size} PDB entries have species-level assignments`)
  const first = page.locator('[data-testid^="species-row-tax:"]').first()
  await expect(first.getByTestId('species-label')).toHaveAttribute('data-italic', 'true')
  await expect(first.getByTestId('species-label')).toHaveText('Bacillus subtilis')
  const bsEntries = new Set(inst.filter((r) => r.organism_species === 'Bacillus subtilis').map((r) => r.pdb_id)).size
  await expect(first.getByTestId('species-count')).toHaveText(String(bsEntries))
  await expect(first).toHaveAttribute('title', new RegExp(`${bsEntries} PDB entries = \\d+\\.\\d % of 66`))
  expect(await first.getByTestId('species-bar').evaluate((e) => (e as HTMLElement).style.width)).toBe('100%')   // largest category spans the track
  await expect(page.getByTestId('species-row-synthetic').getByTestId('species-label')).toHaveAttribute('data-italic', 'false')
  const synthEntries = new Set(inst.filter((r) => r.organism_species === 'synthetic construct').map((r) => r.pdb_id)).size
  await expect(page.getByTestId('species-row-synthetic').getByTestId('species-count')).toHaveText(String(synthEntries))
  await page.getByTestId('unit-instances').click()
  const bsInst = inst.filter((r) => r.organism_species === 'Bacillus subtilis').length
  await expect(first.getByTestId('species-count')).toHaveText(String(bsInst))
  await expect(first).toHaveAttribute('title', new RegExp(`${bsInst} RNA chain instances = \\d+\\.\\d % of ${inst.length}`))
  await page.getByTestId('unit-entries').click()
  await waitMolstar(page, d.pdb_id)
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-mode', 'region')
  const bg = await page.evaluate(() => getComputedStyle(document.querySelector('[data-testid="molstar-viewer"] > div')!).backgroundColor)
  expect(bg).toBe('rgb(255, 255, 255)')
  expect(await page.locator('[data-testid="molstar-viewer"] .msp-layout-region-right, [data-testid="molstar-viewer"] .msp-layout-region-left').count()).toBe(0)
  // the viewer is fitted: the canvas fills the panel
  const box = await page.locator('[data-testid="molstar-viewer"] canvas').boundingBox(); const panel = await page.getByTestId('molstar-viewer').boundingBox()
  expect(box!.width).toBeGreaterThan(panel!.width * 0.95)
  await page.getByTestId('mode-complete').click()
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-mode', 'complete')
})

test('structures table: columns, titles, author chains, italic species, tooltips, search and pagination', async ({ page }) => {
  await page.goto('/')
  await page.getByTestId('tab-structures').click()
  const inst = await (await page.request.get('/families/RF00167/instances.json')).json()
  const headers = await page.locator('[data-testid="instance-table"] thead th > button').allTextContents()
  expect(headers.map((h) => h.replace(/[▲▼]/g, '').trim())).toEqual(['PDB ID', 'Chain', 'Title', 'Organism', 'Method', 'Res. (Å)', 'Region', 'Seq. cov. %', 'Struct. cov. %', 'Avg PP', 'Rfam bits'])
  await expect(page.locator('[data-testid="instance-table"] tbody tr')).toHaveCount(10)
  // first visit opens the page holding the display representative (4FE5/B)
  const fam = await (await page.request.get('/families/RF00167/family.json')).json()
  const sorted = inst.slice().sort((a: { pdb_id: string }, b: { pdb_id: string }) => a.pdb_id.toUpperCase().localeCompare(b.pdb_id.toUpperCase()))
  const repPage = Math.floor(sorted.findIndex((r: { seq_name: string }) => r.seq_name === fam.display_representative.seq_name) / 10) + 1
  await expect(page.getByTestId('table-count')).toHaveText(`Showing ${(repPage - 1) * 10 + 1}–${Math.min(repPage * 10, inst.length)} of ${inst.length} instances`)
  await expect(page.getByTestId('page-info')).toHaveText(`Page ${repPage} of ${Math.ceil(inst.length / 10)}`)
  // footer layout: count bottom-left, compact Previous / Page / Next group bottom-right, nothing beside the search field
  const footer = await page.getByTestId('table-footer').boundingBox(); const count = await page.getByTestId('table-count').boundingBox(); const nav = await page.getByTestId('pagination').boundingBox()
  expect(count!.x).toBeLessThan(footer!.x + 24); expect(nav!.x + nav!.width).toBeGreaterThan(footer!.x + footer!.width - 24)
  const [prev, info, next] = await Promise.all([page.getByTestId('page-prev').boundingBox(), page.getByTestId('page-info').boundingBox(), page.getByTestId('page-next').boundingBox()])
  const g1 = info!.x - (prev!.x + prev!.width), g2 = next!.x - (info!.x + info!.width)
  expect(Math.abs(g1 - g2)).toBeLessThan(2); expect(g1).toBeLessThan(24)
  expect(await page.locator('[data-testid="instance-table"]').locator('..').locator('..').getByTestId('table-count').count()).toBe(1)
  // author chain only, title truncated with full title on hover
  const row = page.getByTestId('row-1y26_1_X_A')
  await page.getByTestId('instance-search').fill('1y26')
  await expect(row.locator('td').nth(1)).toHaveText('X')
  const title = row.getByTestId('cell-title')
  await expect(title).toHaveAttribute('title', inst.find((r: { seq_name: string }) => r.seq_name === '1y26_1_X_A').title)
  expect(await title.evaluate((e) => getComputedStyle(e).textOverflow)).toBe('ellipsis')
  await expect(row.getByTestId('cell-organism')).toHaveClass(/italic/)
  await expect(row.getByTestId('cell-organism')).toHaveText('Vibrio vulnificus')
  await page.getByTestId('instance-search').fill('')
  // strain names collapsed to species, source name preserved on hover
  await page.getByTestId('instance-search').fill('str. 168')
  const strainRows = page.locator('[data-testid="instance-table"] tbody tr')
  await expect(strainRows.first().getByTestId('cell-organism')).toHaveText('Bacillus subtilis')
  await expect(strainRows.first().getByTestId('cell-organism')).toHaveAttribute('title', /str\. 168/)
  // synthetic construct not italic
  await page.getByTestId('instance-search').fill('synthetic')
  await expect(page.locator('[data-testid="instance-table"] tbody tr').first().getByTestId('cell-organism')).not.toHaveClass(/italic/)
  // tooltips by mouse and keyboard
  const tips = page.getByTestId('info-tip')
  expect(await tips.count()).toBeGreaterThanOrEqual(6)
  // whole table fits the 1400px viewport (no horizontal overflow hiding the last columns)
  const tbl = await page.getByTestId('instance-table').boundingBox(); expect(tbl!.x + tbl!.width).toBeLessThanOrEqual(1400)
  await tips.nth(3).hover()
  await expect(page.getByRole('tooltip')).toContainText(/model length \(102\)/)
  await page.mouse.move(0, 0)
  await tips.nth(5).focus()
  await expect(page.getByRole('tooltip')).toContainText(/posterior probability/i)
  await page.keyboard.press('Escape')
  await expect(page.getByRole('tooltip')).toHaveCount(0)
  // pagination with preserved selection and reset on search
  await page.getByTestId('instance-search').fill('1y26')
  await page.getByTestId('row-1y26_1_X_A').click()
  await page.getByTestId('instance-search').fill('')
  await expect(page.getByTestId('row-1y26_1_X_A')).toHaveAttribute('aria-selected', 'true')
  await expect(page.getByTestId('page-info')).toHaveText(`Page 1 of ${Math.ceil(inst.length / 10)}`)   // 1Y26 sorts onto page 1
  await page.getByTestId('page-next').click()
  await expect(page.getByTestId('page-info')).toHaveText(`Page 2 of ${Math.ceil(inst.length / 10)}`)
  await expect(page.getByTestId('table-count')).toHaveText(`Showing 11–20 of ${inst.length} instances`)
  await expect(page.getByTestId('structures-current')).toContainText('1Y26 chain X')
  await page.getByTestId('page-prev').click()
  await expect(page.getByTestId('row-1y26_1_X_A')).toHaveAttribute('aria-selected', 'true')
  await page.getByTestId('page-next').click()
  await page.getByTestId('instance-search').fill('vulnificus')
  await expect(page.getByTestId('page-info')).toHaveText('Page 1 of 2')
  await expect(page.getByTestId('table-count')).toContainText('(filtered from 84)')
  await waitMolstar(page, '1y26')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-mode', 'region')
})

test('annotations: legend, tooltips with numerators, 2D fully visible, 3D not re-focused on selection, selection preserved across instances', async ({ page }) => {
  await page.goto('/')
  await page.getByTestId('tab-annotations').click()
  await page.getByTestId('instance-select').selectOption('1y26_1_X_A')
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready')
  await waitMolstar(page, '1y26')
  await expect(page.getByTestId('track-legend')).toContainText('fraction of the 29 unique mapped RNA sequences')
  await expect(page.getByTestId('track-legend')).toContainText('fraction of the 84 RNA structural instances')
  await expect(page.getByTestId('track-legend')).toContainText('provisional display classifications')
  const titles = await page.locator('[data-testid="consensus-tracks"] rect title').allTextContents()
  expect(titles.find((t) => t.startsWith('Consensus 50: sequence coverage'))).toMatch(/29 \/ 29 unique sequences = 100\.0 %/)
  expect(titles.find((t) => t.startsWith('Consensus 50: structural coverage'))).toMatch(/\/ 84 instances with observed coordinates/)
  expect(titles.find((t) => t.startsWith('Consensus 50: mean Infernal posterior'))).toMatch(/over 29 aligned unique sequences/)
  // 2D topology inside its panel, not clipped
  const panel = await page.getByTestId('rna2d-viewer').boundingBox()
  const bb = await page.evaluate(() => { const rs = [...document.querySelectorAll('[data-testid="rna2d-viewer"] text.rnaviewEle')].map((e) => e.getBoundingClientRect()); return { l: Math.min(...rs.map((r) => r.left)), t: Math.min(...rs.map((r) => r.top)), r: Math.max(...rs.map((r) => r.right)), b: Math.max(...rs.map((r) => r.bottom)), n: rs.length } })
  expect(bb.n).toBe(71)
  expect(bb.l).toBeGreaterThanOrEqual(panel!.x - 1); expect(bb.r).toBeLessThanOrEqual(panel!.x + panel!.width + 1); expect(bb.t).toBeGreaterThanOrEqual(panel!.y - 1); expect(bb.b).toBeLessThanOrEqual(panel!.y + panel!.height + 1)
  // selection colours but does not move the camera; Focus selection does
  type W = { __rnaExplorer: { getCamera: () => unknown; lastMolstarSelection: { labelSeqIds: number[] }; selection: { positions: number[] } } }
  const camera = () => page.evaluate(() => JSON.stringify((window as unknown as W).__rnaExplorer.getCamera()))
  const camBefore = await camera()
  await page.getByTestId('goto-input').fill('50'); await page.getByTestId('goto-button').click()
  await expect(page.getByTestId('sel-highlightable')).toContainText('50 → G34 (auth 46')
  await expect.poll(() => page.evaluate(() => (window as unknown as W).__rnaExplorer.lastMolstarSelection.labelSeqIds)).toEqual([34])
  await page.waitForTimeout(800)
  expect(await camera()).toBe(camBefore)
  await expect(page.getByTestId('focus-selection')).toBeEnabled()
  await page.getByTestId('focus-selection').click()   // explicit focus action sends the selection to Mol*'s camera focus
  await expect.poll(() => page.evaluate(() => (window as unknown as { __rnaExplorer: { lastFocus?: number[] } }).__rnaExplorer.lastFocus)).toEqual([34])
  // switch instance: selection preserved
  await page.getByTestId('instance-select').selectOption('4fe5_1_B_A')
  await waitMolstar(page, '4fe5')
  await expect(page.getByTestId('sel-highlightable')).toContainText('50 → ')
  expect(await page.evaluate(() => (window as unknown as W).__rnaExplorer.selection.positions)).toEqual([50])
})
