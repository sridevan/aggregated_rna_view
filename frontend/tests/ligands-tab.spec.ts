import { expect, test, type Page } from '@playwright/test'

/* Ligands tab: RNA-contacting ligand gallery (default ligand set of ligands.json), the per-CCD ligand-contact residue tracks and the
   fpocketR pocket track moved from the Annotations tab (same track browser, unchanged values), numbering as labels only, provenance texts.
   Expected values come from the family's JSON files, never from the rendering. */
type ByCcd = Record<string, { name: string | null; n_entries_in_default_track: number }>
type Bm = { pdb_id: string; ccd_id: string; in_default_track: boolean; consensus_positions: number[] }
const open = async (page: Page, acc = 'RF00167') => { await page.goto(`/${acc}`); await page.getByTestId('tab-ligands').click(); await expect(page.getByTestId('ligands-tab')).toBeVisible() }
const defaultSet = (byCcd: ByCcd) => Object.entries(byCcd).filter(([, v]) => v.n_entries_in_default_track > 0)
  .sort((a, b) => b[1].n_entries_in_default_track - a[1].n_entries_in_default_track || a[0].localeCompare(b[0]))

test('tab order; gallery = default ligand set with names and PDB entry counts; PDBe depictions; tracks moved unchanged', async ({ page }) => {
  await open(page)
  await expect(page.getByRole('tab')).toHaveText(['Family Overview', 'Structures', 'Interactions', 'Ligands', 'Annotations'])
  const lig = await (await page.request.get('/families/RF00167/ligands.json')).json()
  const set = defaultSet(lig.by_ccd)
  await expect(page.getByTestId('ligand-gallery-count')).toHaveText(`(${set.length})`)
  // first 12 in track order, then "Show all"
  const cards = page.locator('[data-testid^="ligand-card-"]')
  await expect(cards).toHaveCount(12)
  expect(await cards.evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')!.slice('ligand-card-'.length)))).toEqual(set.slice(0, 12).map(([c]) => c))
  await page.getByTestId('ligand-gallery-toggle').click(); await expect(cards).toHaveCount(set.length)
  // HPA: name and distinct PDB entries with a recorded RNA contact (no ligand-copy counts anywhere)
  await expect(page.getByTestId('ligand-name-HPA')).toHaveText('Hypoxanthine')
  await expect(page.getByTestId('ligand-entries-HPA')).toHaveText(`Found in ${lig.by_ccd.HPA.n_entries_in_default_track} PDB entries`)
  const hpaEntries = [...new Set((lig.bound_molecules as Bm[]).filter((b) => b.ccd_id === 'HPA' && b.in_default_track).map((b) => b.pdb_id.toUpperCase()))].sort()
  await expect(page.getByTestId('ligand-entries-HPA')).toHaveAttribute('title', `RNA contact recorded by PDBe in: ${hpaEntries.join(', ')}`)
  await expect(page.getByTestId('ligand-gallery-section')).not.toContainText(/cop(y|ies)|instances? of the ligand/i)
  await expect(page.getByTestId('ligand-link-HPA')).toHaveAttribute('href', 'https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/HPA')
  // authoritative depiction: the PDBe chemistry SVG for the CCD id, actually loaded
  const img = page.getByTestId('depiction-ADE')
  await expect(img).toHaveAttribute('src', 'https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/ADE_200.svg'); await expect(img).toHaveAttribute('data-source', 'pdbe')
  await expect.poll(() => img.evaluate((e: HTMLImageElement) => e.complete && e.naturalWidth > 0), { timeout: 20_000 }).toBe(true)
  // one residue track per gallery ligand, in the same order; the pocket track follows in its own category
  expect(await page.locator('[data-testid^="track-lig-"]').evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')!.slice('track-lig-'.length)))).toEqual(set.map(([c]) => c))
  await expect(page.getByTestId('category-toggle-ligands')).toContainText('Ligand-contact residues')
  await expect(page.getByTestId('category-toggle-pockets')).toContainText('Predicted binding pockets')
  await expect(page.getByTestId('track-pocket-occurrence')).toBeVisible()
  await expect(page.getByTestId('track-lig-HPA')).toContainText('Hypoxanthine · 12 entries')
  // colour values = the unchanged backend frequencies
  const tr = await (await page.request.get('/families/RF00167/tracks.json')).json()
  await expect(page.locator('[data-testid="track-lig-HPA"] rect[data-pos="79"]')).toHaveAttribute('data-value', tr.positions[78].ligand_binding_sites.by_ccd_entries.HPA.frequency.toFixed(3))
  await expect(page.getByTestId('scale-legend')).toContainText('Occurrence frequency (0–100%)')
})

test('tooltips and pinned details keep the PDB-entry counting; provenance texts; methodology section', async ({ page }) => {
  await open(page)
  const tr = await (await page.request.get('/families/RF00167/tracks.json')).json(); const lig = await (await page.request.get('/families/RF00167/ligands.json')).json()
  const tip = page.getByTestId('annotation-tooltip')
  const hpa = new Set((lig.bound_molecules as Bm[]).filter((b) => b.ccd_id === 'HPA' && b.in_default_track && b.consensus_positions.includes(26)).map((b) => b.pdb_id))
  const e26 = tr.positions[25].ligand_binding_sites.by_ccd_entries.HPA
  expect(e26.n_positive_observable_and_evaluated).toBe(hpa.size)
  await page.getByTestId('track-lig-HPA').locator('rect[data-pos="26"]').hover()
  await expect(tip).toHaveText(`Rfam residue 26HPA — HypoxanthineRecorded contact in ${hpa.size} of ${e26.n_observable_and_evaluated} observable PDB entries (${Math.round(1000 * hpa.size / e26.n_observable_and_evaluated) / 10}%)`)
  const e79 = tr.positions[78].ligand_binding_sites.by_ccd_entries.HPA   // 3 of 12 (4FE5, 4FEL, 4FEN)
  await page.getByTestId('track-lig-HPA').locator('rect[data-pos="79"]').click()
  const inspect = page.getByTestId('annotation-inspect')
  await expect(inspect).toContainText(`Recorded contact${e79.n_positive_observable_and_evaluated} of ${e79.n_observable_and_evaluated} observable PDB entries with bound HPA`)
  await expect(inspect).toContainText('Counting unitPDB entry (chains of one entry counted once)')
  await inspect.locator('summary').click(); await expect(inspect).toContainText('Supporting entries: 4FE5, 4FEL, 4FEN')
  await expect(inspect).toContainText('An observed-contact frequency, not a binding probability')
  await page.getByTestId('annotations-clear').click()
  // pocket tooltip unchanged
  const pk = (await (await page.request.get('/families/RF00167/pocket_tracks.json')).json()).tracks.all_score_ge_0
  await page.getByTestId('track-pocket-occurrence').locator('rect[data-pos="53"]').hover()
  await expect(tip).toHaveText(`Rfam residue 53Pocket occurrence: ${pk[52].n_supporting_instances} of ${pk[52].n_evaluable_instances} RNA instances (16.3%)`)
  await page.mouse.move(5, 5)
  // information texts (§8): gallery, residues, pockets, frequency legend, numbering
  const popover = async (loc: ReturnType<Page['getByTestId']>, text: string) => { await loc.getByTestId('info-tip').click(); await expect(page.getByTestId('info-popover')).toContainText(text); await page.keyboard.press('Escape') }
  await popover(page.getByTestId('ligand-gallery-section').locator('h2'), 'Ligands shown here have contacts reported by PDBe for RNA entities belonging to this Rfam family. Contacts are derived from Arpeggio and may not represent every RNA chain or biological assembly in an entry.')
  await popover(page.getByTestId('category-header-ligands'), 'The annotations are aggregated at RNA entity level rather than individual chain level.')
  await popover(page.getByTestId('category-header-pockets'), 'Computationally predicted RNA binding pockets identified using fpocketR. These annotations are independent of the PDBe ligand-contact data.')
  await popover(page.getByTestId('legend-info'), 'Contact frequency is the proportion of eligible PDB entries containing a ligand in which PDBe reports a contact at this position. Entries with unobserved residues are excluded from the denominator.')
  await popover(page.getByTestId('numbering-info'), 'Ligand-contact annotations remain aggregated at RNA entity level and are not necessarily specific to the selected chain.')
  // methodology: collapsed by default, then lists source, level, coverage, unit and the limitation
  const methods = page.getByTestId('ligand-methods')
  await expect(methods).toHaveAttribute('data-open', 'false'); await expect(page.getByTestId('ligand-methods-list')).toHaveCount(0)
  await methods.getByRole('button', { name: /About these annotations/ }).click()
  for (const t of ['PDBe API v2 binding_sites', 'Arpeggio', 'RNA entity', 'preferred biological assembly', 'Distinct PDB entries', 'not always available', 'not proof of an experimentally validated binding event'])
    await expect(page.getByTestId('ligand-methods-list')).toContainText(t)
})

test('selecting a ligand highlights and reveals its track (even in a collapsed category) without changing any value; search filters the gallery', async ({ page }) => {
  await open(page)
  const values = () => page.locator('[data-testid^="track-lig-"] rect[data-pos]').evaluateAll((es) => es.map((e) => `${e.closest('[data-testid^="track-"]')!.getAttribute('data-testid')}:${e.getAttribute('data-pos')}:${e.getAttribute('data-value')}`))
  const before = await values()
  await page.getByTestId('category-toggle-ligands').click(); await expect(page.getByTestId('category-ligands')).toHaveAttribute('data-open', 'false')
  // a ligand far down the list: shown via "Show all", its card selected, the category re-opened, the row highlighted and scrolled into view
  await page.getByTestId('ligand-gallery-toggle').click()
  const lig = await (await page.request.get('/families/RF00167/ligands.json')).json(); const set = defaultSet(lig.by_ccd)
  const last = set[set.length - 1][0]
  await page.getByTestId(`ligand-select-${last}`).click()
  await expect(page.getByTestId(`ligand-card-${last}`)).toHaveAttribute('data-selected', 'true')
  await expect(page.getByTestId('category-ligands')).toHaveAttribute('data-open', 'true')
  const row = page.getByTestId(`track-lig-${last}`)
  await expect(row).toHaveAttribute('data-highlighted', 'true')
  await expect.poll(async () => { const r = (await row.boundingBox())!, b = (await page.getByTestId('category-body-ligands').boundingBox())!; return r.y >= b.y - 1 && r.y + r.height <= b.y + b.height + 1 }).toBe(true)
  await expect(page.locator('[data-highlighted="true"]')).toHaveCount(1)
  expect(await values()).toEqual(before)                                    // highlight only: every mark and value unchanged
  // selecting it again clears the highlight; selecting another moves it
  await page.getByTestId(`ligand-ccd-${last}`).click(); await expect(page.locator('[data-highlighted="true"]')).toHaveCount(0)
  await page.getByTestId('ligand-select-ADE').click(); await expect(page.getByTestId('track-lig-ADE')).toHaveAttribute('data-highlighted', 'true')
  // search: name or CCD id, case-insensitive; every match shown; no match message
  await page.getByTestId('ligand-search').fill('guanine')
  const shown = await page.locator('[data-testid^="ligand-card-"]').evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')!.slice('ligand-card-'.length)))
  expect(shown).toEqual(set.filter(([c, v]) => c.toLowerCase().includes('guanine') || (v.name ?? '').toLowerCase().includes('guanine')).map(([c]) => c))
  expect(shown).toContain('GUN'); await expect(page.getByTestId('ligand-gallery-toggle')).toHaveCount(0)
  await page.getByTestId('ligand-search').fill('hpa'); await expect(page.locator('[data-testid^="ligand-card-"]')).toHaveCount(1)
  await page.getByTestId('ligand-search').fill('zzzz-none'); await expect(page.getByTestId('ligand-gallery-no-match')).toBeVisible()
  await expect(page.locator('[data-testid^="track-lig-"]')).toHaveCount(set.length)   // the search filters the gallery only, never the tracks
})

test('depiction fallback: RCSB image when PDBe fails, a placeholder when both fail', async ({ page }) => {
  await page.route('https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/**', (r) => r.abort())
  await page.route('https://cdn.rcsb.org/images/ccd/unlabeled/G/**', (r) => r.abort())
  await open(page)
  await expect(page.getByTestId('depiction-HPA')).toHaveAttribute('data-source', 'rcsb')
  await expect(page.getByTestId('depiction-HPA')).toHaveAttribute('src', 'https://cdn.rcsb.org/images/ccd/unlabeled/H/HPA.svg')
  await expect(page.getByTestId('depiction-missing-GUN')).toHaveText('No depiction available')
})

test('numbering changes labels only: axis and tooltips follow the chosen chain, contacts and frequencies stay aggregated', async ({ page }) => {
  await open(page)
  const values = () => page.locator('[data-testid^="track-"] rect[data-pos]').evaluateAll((es) => es.map((e) => `${e.closest('[data-testid^="track-"]')!.getAttribute('data-testid')}:${e.getAttribute('data-pos')}:${e.getAttribute('data-value')}`))
  const before = await values()
  await page.getByTestId('numbering-select').selectOption('1y26_1_X_A')
  await expect(page.getByTestId('annotation-axis-convention')).toHaveText('1Y26 X · author numbering')
  expect(await values()).toEqual(before)                                   // same tracks, marks and values: no filtering to that chain
  await page.getByTestId('track-lig-HPA').locator('rect[data-pos="26"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toContainText('1Y26 X reference residue:')
  await page.getByTestId('track-pocket-occurrence').locator('rect[data-pos="53"]').hover()
  await expect(page.getByTestId('annotation-tooltip')).toContainText('1Y26 X reference residue:')
  // shared with the other tabs: Annotations shows the same numbering
  await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('annotation-axis-convention')).toHaveText('1Y26 X · author numbering')
  await page.getByTestId('numbering-select').selectOption('')
})

test('pinned details: View in Family Overview says the position came from the Ligands tab', async ({ page }) => {
  await open(page)
  await page.getByTestId('track-lig-HPA').locator('rect[data-pos="25"]').click()
  await page.getByTestId('view-in-overview').click()
  await expect(page.getByTestId('tab-overview')).toHaveAttribute('aria-selected', 'true')
  await expect(page.getByTestId('selection-origin')).toContainText('Opened from Ligands (HPA at Rfam position 25)')
})

test('ligand category scrolls independently with a visible header; Download menu serves the ligand and pocket files', async ({ page }) => {
  await open(page)
  const body = page.getByTestId('category-body-ligands')
  const dims = await body.evaluate((e) => ({ client: e.clientHeight, scroll: e.scrollHeight }))
  expect(dims.client).toBeLessThanOrEqual(230); expect(dims.scroll).toBeGreaterThan(dims.client)
  const headerBefore = (await page.getByTestId('category-ligands').getByRole('button').first().boundingBox())!
  await body.evaluate((e) => { e.scrollTop = 150 }); await expect.poll(() => body.evaluate((e) => e.scrollTop)).toBe(150)
  expect((await page.getByTestId('category-ligands').getByRole('button').first().boundingBox())!.y).toBe(headerBefore.y)
  expect(await page.getByTestId('category-body-pockets').evaluate((e) => e.scrollTop)).toBe(0)
  await page.getByTestId('download-button').click()
  await expect(page.getByTestId('download-items').getByRole('menuitem')).toHaveCount(4)
  for (const [id, file] of [['ligands', 'ligands.json'], ['tracks', 'tracks.json'], ['pockets', 'pockets.json'], ['pocket-tracks', 'pocket_tracks.json']]) {
    const href = (await page.getByTestId(`download-${id}`).getAttribute('href'))!
    expect(href.endsWith(`/families/RF00167/${file}`)).toBe(true)
    const r = await page.request.get(href); expect(r.ok()).toBe(true)
  }
})

test('other families: no hardcoded data; empty states when PDBe reports no ligand or fpocketR was not run', async ({ page }) => {
  for (const acc of ['RF00177', 'RF00011']) {
    await open(page, acc)
    const lig = await (await page.request.get(`/families/${acc}/ligands.json`)).json(); const n = defaultSet(lig.by_ccd).length
    await expect(page.getByTestId('ligand-gallery-count')).toHaveText(`(${n})`)
    await expect(page.locator('[data-testid^="track-lig-"]')).toHaveCount(n)
    await expect(page.getByTestId('category-pockets')).toHaveCount(0)
    await expect(page.getByTestId('categories-unavailable')).toHaveText(`Not available for ${acc}: Predicted binding pockets.`)
    if (!n) { await expect(page.getByTestId('ligand-gallery-empty')).toBeVisible(); await expect(page.getByTestId('category-empty-ligands')).toBeVisible() }
  }
})

test('narrow screen: no horizontal overflow; gallery stacks; tracks aligned with the pocket track', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 })
  await open(page)
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)
  const c = await page.locator('[data-testid^="ligand-card-"]').evaluateAll((es) => es.slice(0, 2).map((e) => e.getBoundingClientRect().left))
  expect(Math.abs(c[0] - c[1])).toBeLessThan(1)                             // one column
  const a = (await page.locator('[data-testid="track-lig-HPA"] rect[data-pos="25"]').boundingBox())!, b = (await page.locator('[data-testid="track-pocket-occurrence"] rect[data-pos="25"]').boundingBox())!
  expect(Math.abs(a.x - b.x)).toBeLessThan(1)
})
