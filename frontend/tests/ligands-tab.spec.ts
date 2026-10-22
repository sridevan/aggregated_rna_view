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
  await expect(page.getByRole('tab')).toHaveText(['Family Overview', 'Structures', 'Ligands', 'Interactions', 'Annotations'])
  // the tab bar shows no scrollbar (only its 1 px bottom border around the content box)
  expect(await page.getByTestId('tab-bar').evaluate((e) => [getComputedStyle(e).overflowY, e.offsetHeight - e.clientHeight, e.offsetWidth - e.clientWidth])).toEqual(['hidden', 1, 0])
  const lig = await (await page.request.get('/families/RF00167/ligands.json')).json()
  const set = defaultSet(lig.by_ccd)
  await expect(page.getByTestId('ligand-gallery-count')).toHaveText(`(${set.length})`)
  // pages of 9 in track order (3 x 3 on desktop); every ligand reachable through the pages
  const cards = page.locator('[data-testid^="ligand-card-"]')
  const ids = () => cards.evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')!.slice('ligand-card-'.length)))
  const pages = Math.ceil(set.length / 9)
  await expect(page.getByTestId('ligand-page-indicator')).toHaveText(`Page 1 of ${pages}`)
  await expect(page.getByTestId('ligand-page-summary')).toHaveText(`Showing 1–9 of ${set.length} ligands`)
  expect(await ids()).toEqual(set.slice(0, 9).map(([c]) => c))
  const x = await cards.evaluateAll((es) => es.slice(0, 4).map((e) => Math.round(e.getBoundingClientRect().left)))
  expect(new Set(x.slice(0, 3)).size).toBe(3); expect(x[3]).toBe(x[0])                                       // three per row
  const seen: string[] = [...await ids()]
  for (let k = 2; k <= pages; k++) { await page.getByTestId('ligand-page-next').click(); await expect(page.getByTestId('ligand-page-indicator')).toHaveText(`Page ${k} of ${pages}`); seen.push(...await ids()) }
  expect(seen).toEqual(set.map(([c]) => c)); await expect(page.getByTestId('ligand-page-next')).toBeDisabled()
  await expect(page.getByTestId('ligand-page-summary')).toHaveText(`Showing ${(pages - 1) * 9 + 1}–${set.length} of ${set.length} ligands`)
  await page.getByTestId('ligand-page-1').click(); await expect(page.getByTestId('ligand-page-1')).toHaveAttribute('aria-current', 'page'); await expect(page.getByTestId('ligand-page-prev')).toBeDisabled()
  // HPA: name and distinct PDB entries with a recorded RNA contact (no ligand-copy counts anywhere)
  await expect(page.getByTestId('ligand-name-HPA')).toHaveText('Hypoxanthine')
  await expect(page.getByTestId('ligand-entries-HPA')).toHaveText(`Found in ${lig.by_ccd.HPA.n_entries_in_default_track} PDB entries`)
  const hpaEntries = [...new Set((lig.bound_molecules as Bm[]).filter((b) => b.ccd_id === 'HPA' && b.in_default_track).map((b) => b.pdb_id.toUpperCase()))].sort()
  await expect(page.getByTestId('ligand-entries-HPA')).toHaveAttribute('title', `RNA contact recorded by PDBe in: ${hpaEntries.join(', ')}`)
  await expect(page.getByTestId('ligand-gallery-section')).not.toContainText(/cop(y|ies)|instances? of the ligand/i)
  // the CCD id is the external PDBe link (new tab); the old separate link is gone
  const link = page.getByTestId('ligand-ccd-HPA')
  await expect(link).toHaveAttribute('href', 'https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/HPA'); await expect(link).toHaveAttribute('target', '_blank')
  await expect(link).toContainText('HPA'); await expect(link).toContainText('↗')
  await expect(page.locator('[data-testid^="ligand-link-"]')).toHaveCount(0)
  const [popup] = await Promise.all([page.waitForEvent('popup'), link.click()])
  expect(popup.url()).toContain('/chemicalCompound/show/HPA'); await popup.close()
  await expect(page.getByTestId('ligand-card-HPA')).toHaveAttribute('data-selected', 'false')            // the link never selects the card
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
  // a ligand on the last gallery page: its card selected, the category re-opened, the row highlighted and scrolled into view
  const lig = await (await page.request.get('/families/RF00167/ligands.json')).json(); const set = defaultSet(lig.by_ccd)
  const last = set[set.length - 1][0]
  await page.getByTestId(`ligand-page-${Math.ceil(set.length / 9)}`).click()
  await page.getByTestId(`ligand-select-${last}`).click()
  await expect(page.getByTestId(`ligand-card-${last}`)).toHaveAttribute('data-selected', 'true')
  await expect(page.getByTestId('category-ligands')).toHaveAttribute('data-open', 'true')
  const row = page.getByTestId(`track-lig-${last}`)
  await expect(row).toHaveAttribute('data-highlighted', 'true')
  await expect.poll(async () => { const r = (await row.boundingBox())!, b = (await page.getByTestId('category-body-ligands').boundingBox())!; return r.y >= b.y - 1 && r.y + r.height <= b.y + b.height + 1 }).toBe(true)
  await expect(page.locator('[data-highlighted="true"]')).toHaveCount(1)
  expect(await values()).toEqual(before)                                    // highlight only: every mark and value unchanged
  // the highlight survives a page change (the chip names it) and is cleared from the chip or by selecting the card again
  await page.getByTestId('ligand-page-1').click()
  await expect(page.getByTestId('ligand-selected-chip')).toContainText(`Highlighted track: ${last}`); await expect(row).toHaveAttribute('data-highlighted', 'true')
  await page.getByTestId('ligand-selected-clear').click(); await expect(page.locator('[data-highlighted="true"]')).toHaveCount(0); await expect(page.getByTestId('ligand-selected-chip')).toHaveCount(0)
  await page.getByTestId('ligand-locate-ADE').click(); await expect(page.getByTestId('track-lig-ADE')).toHaveAttribute('data-highlighted', 'true')
  await expect(page.getByTestId('ligand-locate-ADE')).toHaveText('Clear track highlight')
  await page.getByTestId('ligand-select-ADE').click(); await expect(page.locator('[data-highlighted="true"]')).toHaveCount(0)
  await page.getByTestId('ligand-select-ADE').click()
  // search: name or CCD id, case-insensitive; pagination over the matches, back to page 1 on every new query
  await page.getByTestId('ligand-page-2').click()
  await page.getByTestId('ligand-search').fill('guanine')
  const match = set.filter(([c, v]) => c.toLowerCase().includes('guanine') || (v.name ?? '').toLowerCase().includes('guanine')).map(([c]) => c)
  const shown = await page.locator('[data-testid^="ligand-card-"]').evaluateAll((es) => es.map((e) => e.getAttribute('data-testid')!.slice('ligand-card-'.length)))
  expect(shown).toEqual(match.slice(0, 9)); expect(shown).toContain('GUN')
  await expect(page.getByTestId('ligand-page-summary')).toHaveText(`Showing 1–${Math.min(9, match.length)} of ${match.length} ${match.length === 1 ? 'ligand' : 'ligands'} matching “guanine”`)
  if (match.length <= 9) await expect(page.getByTestId('ligand-page-indicator')).toHaveCount(0)
  await page.getByTestId('ligand-search').fill('hpa'); await expect(page.locator('[data-testid^="ligand-card-"]')).toHaveCount(1)
  await page.getByTestId('ligand-search').fill('zzzz-none'); await expect(page.getByTestId('ligand-gallery-no-match')).toBeVisible(); await expect(page.getByTestId('ligand-pagination')).toHaveCount(0)
  await expect(page.getByTestId('track-lig-ADE')).toHaveAttribute('data-highlighted', 'true')          // the gallery search never changes the tracks or the highlight
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
  // tabs still scroll sideways on a phone, without a visible scrollbar
  const bar = page.getByTestId('tab-bar')
  expect(await bar.evaluate((e) => e.scrollWidth > e.clientWidth && e.offsetHeight - e.clientHeight === 1)).toBe(true)
  await bar.evaluate((e) => { e.scrollLeft = 150 }); expect(await bar.evaluate((e) => e.scrollLeft)).toBeGreaterThan(0)
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)
  const c = await page.locator('[data-testid^="ligand-card-"]').evaluateAll((es) => es.slice(0, 2).map((e) => e.getBoundingClientRect().left))
  expect(Math.abs(c[0] - c[1])).toBeLessThan(1)                             // one column
  const a = (await page.locator('[data-testid="track-lig-HPA"] rect[data-pos="25"]').boundingBox())!, b = (await page.locator('[data-testid="track-pocket-occurrence"] rect[data-pos="25"]').boundingBox())!
  expect(Math.abs(a.x - b.x)).toBeLessThan(1)
})

/** a valid V2000 molfile with the given coordinates (test fixture served in place of the remote SDF) */
const sdf = (atoms: [string, number, number, number][]) => ['FIX', '  test', '', `${String(atoms.length).padStart(3)}  0  0  0  0  0  0  0  0  0999 V2000`,
  ...atoms.map(([el, x, y, z]) => `${x.toFixed(4).padStart(10)}${y.toFixed(4).padStart(10)}${z.toFixed(4).padStart(10)} ${el.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`), 'M  END', '$$$$'].join('\n')
type LV = { __rnaExplorer: { ligandViewer?: { ccd: string; status: string; nAtoms?: number; source?: string; reason?: string } | null } }
const viewer = (page: Page) => page.evaluate(() => (window as unknown as LV).__rnaExplorer.ligandViewer)

test('View ideal 3D: modal with CCD id and name, PDBe ideal SDF in Mol* (ball-and-stick), provenance note, reset, Escape and focus return; no RNA viewer affected', async ({ page }) => {
  await open(page)
  const before = await page.locator('canvas').count()
  await page.getByTestId('ligand-view3d-ADE').click()
  const dlg = page.getByRole('dialog'); await expect(dlg).toBeVisible()
  await expect(page.getByTestId('ligand-card-ADE')).toHaveAttribute('data-selected', 'false')            // View 3D never selects the card
  await expect(page.getByTestId('ligand-viewer-title')).toHaveText('ADE — Adenine')
  await expect(page.getByTestId('ligand-viewer-subtitle')).toContainText('Ideal 3D coordinates from the PDB Chemical Component Dictionary')
  await expect(page.getByTestId('ligand-viewer-close')).toBeFocused()
  await expect.poll(() => viewer(page).then((v) => v?.status), { timeout: 30_000 }).toBe('ready')
  expect(await viewer(page)).toMatchObject({ ccd: 'ADE', nAtoms: 15, source: 'https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/ADE_ideal.sdf' })
  await expect(page.getByTestId('ligand-viewer')).toContainText('15 atoms · file from PDBe')
  await expect(page.getByTestId('ligand-viewer')).toHaveAttribute('data-status', 'ready')
  await expect(page.getByTestId('ligand-viewer-link')).toHaveAttribute('href', 'https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/ADE')
  await page.getByTestId('ligand-viewer-subtitle').getByTestId('info-tip').click()
  await expect(page.getByTestId('info-popover')).toHaveText('This viewer displays idealized coordinates from the PDB Chemical Component Dictionary. These coordinates represent a reference molecular geometry and may differ from the ligand conformation observed in an experimental RNA complex.')
  await page.getByTestId('ligand-viewer-subtitle').getByTestId('info-tip').click(); await expect(page.getByTestId('info-popover')).toHaveCount(0)
  await page.getByTestId('ligand-viewer-reset').click()
  await expect(page.locator('canvas')).toHaveCount(before + 1)                                             // one viewer for the modal
  await page.keyboard.press('Escape'); await expect(dlg).toHaveCount(0)
  await expect(page.getByTestId('ligand-view3d-ADE')).toBeFocused()
  await expect.poll(() => page.locator('canvas').count()).toBe(before)                                     // disposed on close
  // reopening and quickly switching: the last ligand is the one shown, still a single extra canvas
  for (const c of ['HPA', 'GUN']) { await page.getByTestId(`ligand-view3d-${c}`).click(); await page.getByTestId('ligand-viewer-close').click() }
  await page.getByTestId('ligand-view3d-6AP').click()
  await expect.poll(() => viewer(page).then((v) => v && `${v.ccd}:${v.status}`), { timeout: 30_000 }).toBe('6AP:ready')
  await expect(page.getByTestId('ligand-viewer')).toHaveAttribute('data-ccd', '6AP'); await expect.poll(() => page.locator('canvas').count()).toBe(before + 1)
  await page.getByTestId('ligand-viewer-backdrop').click({ position: { x: 5, y: 5 } }); await expect(dlg).toHaveCount(0)   // backdrop click closes
})

test('View ideal 3D: missing ideal coordinates, unreadable file and network failure are reported; RCSB is used only when PDBe fails; no other coordinates substituted', async ({ page }) => {
  await open(page)
  // all-zero coordinates = the CCD has no ideal coordinates (as for MIY or TAC)
  await page.route('https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/HPA_ideal.sdf', (r) => r.fulfill({ body: sdf([['C', 0, 0, 0], ['N', 0, 0, 0], ['O', 0, 0, 0]]) }))
  await page.getByTestId('ligand-view3d-HPA').click()
  await expect(page.getByTestId('ligand-viewer-error')).toHaveText('The Chemical Component Dictionary has no ideal coordinates for HPA, so no ideal 3D structure is shown.')
  expect(await viewer(page)).toMatchObject({ ccd: 'HPA', status: 'error', reason: 'no_ideal_coordinates' })
  await expect(page.getByTestId('ligand-viewer-reset')).toBeDisabled(); await page.getByTestId('ligand-viewer-close').click()
  // unreadable file
  await page.route('https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/GUN_ideal.sdf', (r) => r.fulfill({ body: '<html>not a molfile</html>' }))
  await page.getByTestId('ligand-view3d-GUN').click()
  await expect(page.getByTestId('ligand-viewer-error')).toHaveText('The ideal coordinate file for GUN could not be read.'); await page.getByTestId('ligand-viewer-close').click()
  // PDBe down: the identical RCSB ideal SDF is used, and the footer says so
  await page.route('https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/ADE_ideal.sdf', (r) => r.abort())
  await page.getByTestId('ligand-view3d-ADE').click()
  await expect.poll(() => viewer(page).then((v) => v?.status), { timeout: 30_000 }).toBe('ready')
  expect((await viewer(page))!.source).toBe('https://files.rcsb.org/ligands/download/ADE_ideal.sdf')
  await expect(page.getByTestId('ligand-viewer')).toContainText('15 atoms · file from RCSB'); await page.getByTestId('ligand-viewer-close').click()
  // both down: a clear message, the modal still closes
  await page.route('https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/6AP_ideal.sdf', (r) => r.abort())
  await page.route('https://files.rcsb.org/ligands/download/6AP_ideal.sdf', (r) => r.fulfill({ status: 404, body: 'not found' }))
  await page.getByTestId('ligand-view3d-6AP').click()
  await expect(page.getByTestId('ligand-viewer-error')).toHaveText('Ideal coordinates for 6AP could not be retrieved from PDBe or RCSB. Check the connection and try again.')
  await page.getByTestId('ligand-viewer-close').click(); await expect(page.getByRole('dialog')).toHaveCount(0)
})

test('names longer than 20 characters are shortened with Show full / Show less; DrugBank cross-references: PDBe, RCSB fallback, none', async ({ page }) => {
  await open(page)
  // 6AP "9h-purine-2,6-diamine" (21 characters)
  const n = page.getByTestId('ligand-name-text-6AP')
  await expect(n).toHaveText('9h-purine-2,6-diamin…'); await expect(page.getByTestId('ligand-name-6AP')).toHaveAttribute('title', '9h-purine-2,6-diamine')
  await page.getByTestId('ligand-name-toggle-6AP').click(); await expect(n).toHaveText('9h-purine-2,6-diamine'); await expect(page.getByTestId('ligand-name-toggle-6AP')).toHaveText('Show less')
  await expect(page.getByTestId('ligand-card-6AP')).toHaveAttribute('data-selected', 'false')                // the toggle never selects
  await expect(page.getByTestId('ligand-name-toggle-HPA')).toHaveCount(0); await expect(page.getByTestId('ligand-name-text-HPA')).toHaveText('Hypoxanthine')
  // PDBe mapping (ADE), none in PDBe or RCSB (6GU)
  const ade = page.getByTestId('ligand-drugbank-ADE')
  await expect(ade).toHaveAttribute('data-state', 'mapped', { timeout: 20_000 }); await expect(ade).toHaveAttribute('data-source', 'PDBe')
  await expect(page.getByTestId('ligand-drugbank-link-ADE')).toHaveAttribute('href', 'https://go.drugbank.com/drugs/DB00173')
  await expect(page.getByTestId('ligand-drugbank-6GU')).toHaveAttribute('data-state', 'none', { timeout: 20_000 }); await expect(page.getByTestId('ligand-drugbank-6GU')).toHaveText('No DrugBank mapping')
})

test('DrugBank: RCSB is asked only when the PDBe file has no DrugBank record; unreachable sources are not reported as "no mapping"', async ({ page }) => {
  const rcsb: string[] = []
  page.on('request', (r) => { if (r.url().includes('data.rcsb.org')) rcsb.push(r.url()) })
  await page.route('https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/GUN.cif', (r) => r.fulfill({ body: 'data_GUN\n_chem_comp.id GUN\n' }))
  await page.route('https://data.rcsb.org/rest/v1/core/chemcomp/GUN', (r) => r.fulfill({ json: { rcsb_chem_comp_related: [{ resource_name: 'DrugBank', resource_accession_code: 'DB99999' }] } }))
  await page.route('https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/6GU.cif', (r) => r.abort())
  await page.route('https://data.rcsb.org/rest/v1/core/chemcomp/6GU', (r) => r.abort())
  await open(page)
  const gun = page.getByTestId('ligand-drugbank-GUN')
  await expect(gun).toHaveAttribute('data-source', 'RCSB', { timeout: 20_000 }); await expect(gun).toContainText('DB99999'); await expect(gun).toContainText('via RCSB')
  await expect(page.getByTestId('ligand-drugbank-6GU')).toHaveAttribute('data-state', 'unavailable', { timeout: 20_000 })
  await expect(page.getByTestId('ligand-drugbank-ADE')).toHaveAttribute('data-source', 'PDBe', { timeout: 20_000 })
  expect(rcsb.filter((u) => u.endsWith('/ADE'))).toEqual([])                                                  // PDBe had it: RCSB never asked
})
