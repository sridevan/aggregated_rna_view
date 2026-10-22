import { expect, test, type Page } from '@playwright/test'
import { categoryColor, CATEGORY_STYLE, defaultComparePair } from '../src/lib/interactions'
import type { InteractionsDataset, StructureClusteringDataset } from '../src/types/backend'

/* Interactions tab (RF00167): family-wide conservation from interactions.json, status in the displayed structure, numbering reference vs
   displayed structure, lookup, supporting-structure inspection, and the two-instance comparison with its local superposition. */
type Pair = { id: string; category: string; pos_1: number; pos_2: number; n_observed: number; n_eligible: number; n_entries_observed: number; frequency: number | null
  /** schema 1.4: base pairs count true or near pairs ("paired"); Rfam SS_cons comparison */
  n_paired: number; frequency_paired: number | null; n_entries_paired: number; n_near_only: number; near_classes: Record<string, number>; consensus: 'both' | '2d' | '3d' | null }
type Inst = { seq_name: string; eligibility: string; observed: Record<string, number[]>; near: [string, string, number, number][] }
type DS = { pairs: Pair[]; instances: Inst[]; default_numbering_reference: { seq_name: string } }
let ds: DS
const REF = '4fe5_1_B_A'
const KL = 'basepair:41-66'   // G37-C61 (4FE5 B numbering), the L2-L3 kissing-loop pair
const open = async (page: Page, allRows = true) => {
  ds ??= await (await page.request.get('/families/RF00167/interactions.json')).json()
  await page.goto('/RF00167'); await page.getByTestId('tab-interactions').click()
  await expect(page.getByTestId('idiagram')).toBeVisible(); await expect(page.getByTestId('inumbering-control')).toHaveAttribute('data-state', 'ready')
  if (allRows) await page.getByTestId('ipage-size').selectOption('0')   // most checks look at every row; paging has its own test
}
const labelSeqs = async (page: Page, sn: string, positions: number[]) => {
  const m = await (await page.request.get(`/families/RF00167/mappings/${sn}.json`)).json() as { consensus_columns: { rfam_consensus_pos: number; polymer_seq_position: number | null }[] }
  return positions.map((p) => m.consensus_columns.find((c) => c.rfam_consensus_pos === p)!.polymer_seq_position)
}
/** independent expectation from the data file (base pairs: both bases eligible; observed index of the chain; else a near base pair) */
const statusIn = (sn: string, p: Pair) => {
  const inst = ds.instances.find((i) => i.seq_name === sn)!, e = inst.eligibility
  if (e[p.pos_1 - 1] !== 'E' || e[p.pos_2 - 1] !== 'E') return 'not_comparable'
  const k = ds.pairs.indexOf(p)
  if (Object.values(inst.observed).some((ks) => ks.includes(k))) return 'present'
  return p.category === 'basepair' && inst.near.some(([c, , i, j]) => c === 'basepair' && i === p.pos_1 && j === p.pos_2) ? 'near' : 'absent'
}
const famOrder = (a: Pair, b: Pair) => b.n_entries_paired - a.n_entries_paired || (b.frequency_paired ?? -1) - (a.frequency_paired ?? -1) || b.n_eligible - a.n_eligible || a.pos_1 - b.pos_1 || a.pos_2 - b.pos_2

test('family-wide table: every family interaction (not only the reference), stored numerators and denominators, status in the displayed structure', async ({ page }) => {
  await open(page)
  const bps = ds.pairs.filter((p) => p.category === 'basepair')
  await expect(page.locator('[data-testid^="irow-"]')).toHaveCount(bps.length)                    // base pairs by default
  await expect(page.getByTestId('idiagram')).toHaveAttribute('data-arc-count', String(bps.length))
  // interactions absent from or not comparable in the reference are listed too
  expect(bps.some((p) => statusIn(REF, p) !== 'present')).toBe(true)
  for (const p of bps) {
    await expect(page.getByTestId(`icons-${p.id}`)).toContainText(`${p.n_paired}/${p.n_eligible}`)
    await expect(page.getByTestId(`irow-${p.id}`)).toHaveAttribute('data-status', statusIn(REF, p))
    await expect(page.getByTestId(`iconsensus-${p.id}`)).toHaveAttribute('data-consensus', p.consensus ?? '')
  }
  // kissing loop: conserved in every eligible chain; labels in the default reference numbering
  await expect(page.getByTestId(`irow-${KL}`)).toContainText('G37–C61'); await expect(page.getByTestId(`icons-${KL}`)).toContainText('80/80')
  // missing nucleotide: 4FE5 B has no residue at Rfam 18 -> not comparable, never absent
  await expect(page.getByTestId('irow-basepair:18-88')).toHaveAttribute('data-status', 'not_comparable')
  await expect(page.getByTestId('istatus-basepair:18-88')).toHaveAttribute('title', /Rfam 18: no residue/)
  // conservation colour follows the numeric frequency (continuous scale)
  for (const p of bps.slice(0, 10)) await expect(page.getByTestId(`icons-${p.id}`)).toHaveAttribute('data-color', categoryColor('basepair', p.frequency_paired))
})

test('interaction types: Bp only by default; several categories together; multiple interactions grouped with their own frequencies', async ({ page }) => {
  await open(page)
  await expect(page.getByTestId('itypes-summary')).toHaveText('Bp')
  await page.getByTestId('itypes-toggle').click()
  await expect(page.getByTestId('ifilter-basepair')).toBeChecked()
  for (const c of ['stacking', 'baseribose', 'basephosphate']) await expect(page.getByTestId(`ifilter-${c}`)).not.toBeChecked()
  for (const c of ['stacking', 'baseribose', 'basephosphate']) await page.getByTestId(`ifilter-${c}`).check()
  await expect(page.getByTestId('itypes-summary')).toHaveText('All types')
  await page.keyboard.press('Escape')
  // one top-level row (and one arc) per pair of nucleotides; every annotation still present as its own member
  const key = (p: Pair) => (p.pos_1 < p.pos_2 ? `${p.pos_1}-${p.pos_2}` : `${p.pos_2}-${p.pos_1}`)
  const byKey = new Map<string, Pair[]>(); for (const p of ds.pairs) byKey.set(key(p), [...(byKey.get(key(p)) ?? []), p])
  const multi = [...byKey.entries()].filter(([, v]) => v.length > 1)
  expect(multi.length).toBeGreaterThan(0)
  await expect(page.getByTestId('idiagram')).toHaveAttribute('data-arc-count', String(byKey.size))
  await expect(page.locator('[data-testid^="igroup-multi:"]')).toHaveCount(multi.length)
  const [k, members] = multi[0]
  await expect(page.getByTestId(`idiagram-arc-multi:${k}`)).toHaveAttribute('data-multi', 'true')
  await page.getByTestId(`igroup-toggle-multi:${k}`).click()
  for (const m of members) await expect(page.getByTestId(`icons-${m.id}`)).toContainText(`${m.n_paired}/${m.n_eligible}`)
  // the group row states an explicit "any" occurrence, never a sum of member counts
  await expect(page.getByTestId(`icons-multi:${k}`)).toContainText('any')
  // grouped arc: neutral double line coloured by the explicit "any" occurrence, never by one category's colour
  const arc = page.getByTestId(`idiagram-arc-multi:${k}`)
  await expect(arc).toHaveAttribute('data-category', 'multiple'); await expect(arc).toHaveAttribute('data-dash', '')
  const any = await page.getByTestId(`icons-multi:${k}`).getAttribute('data-color'); expect(await arc.getAttribute('data-color')).toBe(any)
  expect(Object.values(CATEGORY_STYLE).flatMap((v) => v.ramp)).not.toContain(any)
  // clicking the grouped arc selects a member and lists every constituent category with its own frequency
  await page.locator(`[data-testid="idiagram"] [data-hit="multi:${k}"]`).dispatchEvent('click')
  await expect(page.getByTestId('idetails-multiple')).toBeVisible()
  await page.getByTestId(`irow-${members[1].id}`).click()
  await expect(page.getByTestId('idetails-multiple')).toBeVisible()
  for (const m of members) await expect(page.getByTestId(`idetails-member-${m.id}`)).toContainText(`${m.n_paired}/${m.n_eligible}`)
  await expect(page.getByTestId(`idetails-member-${members[1].id}`)).toHaveAttribute('aria-pressed', 'true')
})

test('category encoding: colour family and line pattern per type, intensity from frequency; dropdown Select all / Reset to Bp / none', async ({ page }) => {
  await open(page)
  await page.getByTestId('itypes-toggle').click(); await page.getByTestId('itypes-all').click()
  for (const c of ['basepair', 'stacking', 'baseribose', 'basephosphate']) await expect(page.getByTestId(`ifilter-${c}`)).toBeChecked()
  await expect(page.getByTestId('itypes-summary')).toHaveText('All types')
  // every single-type arc: its category, its dash pattern, and the category colour at its own frequency
  const singles = ds.pairs.filter((p) => ds.pairs.filter((q) => Math.min(q.pos_1, q.pos_2) === Math.min(p.pos_1, p.pos_2) && Math.max(q.pos_1, q.pos_2) === Math.max(p.pos_1, p.pos_2)).length === 1)
  for (const c of ['basepair', 'stacking', 'baseribose', 'basephosphate'] as const) {
    const p = singles.find((x) => x.category === c)!
    const arc = page.getByTestId(`idiagram-arc-${p.id}`)
    await expect(arc).toHaveAttribute('data-category', c); await expect(arc).toHaveAttribute('data-dash', CATEGORY_STYLE[c].dash ?? '')
    await expect(arc).toHaveAttribute('data-color', categoryColor(c, p.frequency_paired))
    await expect(page.getByTestId(`ilegend-${c}`)).toHaveAttribute('data-dash', CATEGORY_STYLE[c].dash ?? '')
    await expect(page.getByTestId(`ilegend-ramp-${c}`)).toHaveAttribute('data-ramp', CATEGORY_STYLE[c].ramp.join(','))   // one intensity scale per type
  }
  // the selected arc keeps its halo whatever its type
  const st = singles.find((x) => x.category === 'baseribose')!
  await page.locator(`[data-testid="idiagram"] [data-hit="${st.id}"]`).dispatchEvent('click')
  await expect(page.getByTestId('idiagram-selected')).toHaveAttribute('data-id', st.id)
  // Reset to Bp; then none selected: an explicit empty state (not "no interactions")
  await page.getByTestId('itypes-reset').click(); await expect(page.getByTestId('itypes-summary')).toHaveText('Bp')
  await page.getByTestId('ifilter-basepair').uncheck(); await expect(page.getByTestId('itypes-summary')).toHaveText('none selected')
  await expect(page.getByTestId('ino-types')).toBeVisible(); await expect(page.getByTestId('ifamily-table')).toContainText('No interaction types selected')
  await page.getByTestId('itypes-reset').click(); await expect(page.locator('[data-testid^="irow-"]')).toHaveCount(ds.pairs.filter((p) => p.category === 'basepair').length)
})

test('default order: PDB entries, then frequency, eligible chains, positions; header sorting and reset', async ({ page }) => {
  await open(page)
  const bps = ds.pairs.filter((p) => p.category === 'basepair')
  const expected = [...bps].sort(famOrder).map((p) => `irow-${p.id}`)
  const ids = () => page.getByTestId('ifamily-tbody').locator('tr').evaluateAll((rs) => rs.map((r) => r.getAttribute('data-testid')))
  expect(await ids()).toEqual(expected)
  await page.getByTestId('isort-pos').click()
  expect(await ids()).toEqual([...bps].sort((a, b) => a.pos_1 * 10000 + a.pos_2 - (b.pos_1 * 10000 + b.pos_2)).map((p) => `irow-${p.id}`))
  await page.getByTestId(`irow-${KL}`).click()
  await page.getByTestId('isort-reset').click()
  expect(await ids()).toEqual(expected)
  await expect(page.getByTestId(`irow-${KL}`)).toHaveAttribute('aria-selected', 'true')            // sorting keeps the selection
})

test('arcs connect the right positions; clicking an arc selects its row, details and 3D residues; numbering changes labels only', async ({ page }) => {
  await open(page)
  const arc = page.getByTestId(`idiagram-arc-${KL}`)
  await expect(arc).toHaveAttribute('data-i', '41'); await expect(arc).toHaveAttribute('data-j', '66')
  await page.locator(`[data-testid="idiagram"] [data-hit="${KL}"]`).dispatchEvent('click')
  await expect(page.getByTestId(`irow-${KL}`)).toHaveAttribute('aria-selected', 'true')
  await expect(page.getByTestId('idiagram')).toHaveAttribute('data-selected', KL)
  await expect(page.getByTestId('idetails-title')).toHaveText('G37–C61 · Bp cWW')
  await expect(page.getByTestId('idetails-conservation')).toContainText('80/80 eligible chains (100%)')
  await expect(page.getByTestId('idetails-status')).toHaveAttribute('data-state', 'present')
  await expect(page.getByTestId('iviewer-caption')).toHaveAttribute('data-highlighted', (await labelSeqs(page, REF, [41, 66])).join(','))
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 })
  // Rfam numbering: same selection, frequencies and geometry; labels change
  const d = await arc.getAttribute('d')
  await page.getByTestId('inumbering-select').selectOption('')
  await expect(page.getByTestId('idetails-title')).toHaveText('Rfam 41–Rfam 66 · Bp cWW')
  await expect(page.getByTestId(`icons-${KL}`)).toContainText('80/80'); await expect(page.getByTestId('idiagram')).toHaveAttribute('data-selected', KL)
  expect(await arc.getAttribute('d')).toBe(d)
  // another chain as numbering reference (5E54 B): labels follow it, status still refers to the displayed 4FE5 B
  await page.getByTestId('inumbering-select').selectOption('5e54_1_B_B')
  await expect(page.getByTestId('idetails-title')).toContainText('–'); await expect(page.getByTestId('idisplayed-select')).toHaveValue(REF)
})

test('absent interaction: residues highlighted without a contact; inspect a supporting structure keeps numbering and frequencies', async ({ page }) => {
  await open(page)
  // an interaction eligible but not annotated in 4FE5 B
  const p = ds.pairs.filter((q) => q.category === 'basepair' && statusIn(REF, q) === 'absent').sort((a, b) => b.n_observed - a.n_observed)[0]
  await page.getByTestId(`irow-${p.id}`).click()
  await expect(page.getByTestId('idetails-status')).toHaveAttribute('data-state', 'absent')
  await expect(page.getByTestId('iviewer-caption')).toContainText('not annotated'); await expect(page.getByTestId('iviewer-caption')).toContainText('no interaction is implied')
  await expect(page.getByTestId('iviewer-caption')).toHaveAttribute('data-highlighted', (await labelSeqs(page, REF, [p.pos_1, p.pos_2])).join(','))
  // supporting chains are the actual observations
  const support = ds.instances.filter((i) => statusIn(i.seq_name, p) === 'present').map((i) => i.seq_name)
  const opts = await page.getByTestId('isupport-select').locator('option').evaluateAll((os) => os.map((o) => (o as HTMLOptionElement).value).filter(Boolean))
  expect(opts.sort()).toEqual([...support].sort())
  const numbering = await page.getByTestId('inumbering-select').inputValue()
  await page.getByTestId('isupport-select').selectOption(opts[0])
  await expect(page.getByTestId('idisplayed-select')).toHaveValue(opts[0])
  await expect(page.getByTestId(`irow-${p.id}`)).toHaveAttribute('data-status', 'present')
  await expect(page.getByTestId('idetails-status')).toHaveAttribute('data-state', 'present')
  await expect(page.getByTestId('inumbering-select')).toHaveValue(numbering)
  await expect(page.getByTestId(`icons-${p.id}`)).toContainText(`${p.n_paired}/${p.n_eligible}`)
  await expect(page.getByTestId('interactions-tab')).toHaveAttribute('data-mode', 'family')   // not switched to comparison
  // every status column entry follows the new displayed structure
  for (const q of ds.pairs.filter((x) => x.category === 'basepair')) await expect(page.getByTestId(`irow-${q.id}`)).toHaveAttribute('data-status', statusIn(opts[0], q))
})

test('hovering an arc shows its residues, Rfam positions, FR3D classes, conservation and the displayed-structure status', async ({ page }) => {
  await open(page)
  await page.locator(`[data-testid="idiagram"] [data-hit="${KL}"]`).dispatchEvent('mousemove', { clientX: 300, clientY: 300 })
  const tip = page.getByTestId('iarc-tooltip')
  await expect(tip).toContainText('G37–C61'); await expect(tip).toContainText('Rfam 41–66'); await expect(tip).toContainText('Bp cWW · 80/80 (100%)')
  await expect(tip).toContainText('Present in 4FE5 B')
})

test('unified table search: residues, pairs, Rfam positions and class codes filter the table; only an explicit choice selects', async ({ page }) => {
  await open(page)
  // the toolbar lookup is gone
  await expect(page.getByTestId('ilookup-input')).toHaveCount(0); await expect(page.getByTestId('ilookup-go')).toHaveCount(0)
  const search = page.getByTestId('isearch'); const rows = page.locator('[data-testid^="irow-"]')
  await expect(search).toHaveAttribute('placeholder', 'Search interactions or residues')
  await search.fill('G37-C61')
  await expect(rows).toHaveCount(1); await expect(page.getByTestId(`irow-${KL}`)).toBeVisible()
  await expect(page.getByTestId('isearch-match')).toContainText('Between G37 and C61 (Rfam 41–66)')
  // typing never selects: no row, arc, details or 3D highlight yet
  await expect(page.getByTestId(`irow-${KL}`)).toHaveAttribute('aria-selected', 'false'); await expect(page.getByTestId('idiagram')).toHaveAttribute('data-selected', '')
  await expect(page.getByTestId('idetails-empty')).toBeVisible(); await expect(page.getByTestId('iviewer-caption')).toHaveAttribute('data-highlighted', '')
  // Enter (explicit) selects the first result: table, diagram, details and 3D follow
  await search.press('Enter')
  await expect(page.getByTestId(`irow-${KL}`)).toHaveAttribute('aria-selected', 'true'); await expect(page.getByTestId('idiagram')).toHaveAttribute('data-selected', KL)
  await expect(page.getByTestId('iviewer-caption')).toHaveAttribute('data-highlighted', (await labelSeqs(page, REF, [41, 66])).join(','))
  // a different search keeps that selection
  const bp41 = ds.pairs.filter((p) => p.category === 'basepair' && (p.pos_1 === 41 || p.pos_2 === 41)).length
  for (const q of ['G37', '37', 'Rfam 41']) { await search.fill(q); await expect(rows).toHaveCount(bp41) }
  await search.fill('Rfam 41-66'); await expect(rows).toHaveCount(1)
  await expect(page.getByTestId('idiagram')).toHaveAttribute('data-selected', KL)
  // class codes and regions are text
  // (text match: near "ncWW" classes match too)
  await search.fill('cWW'); await expect(rows).toHaveCount(ds.pairs.filter((p) => p.category === 'basepair' && [...Object.keys((p as unknown as { classes: Record<string, number> }).classes), ...Object.keys(p.near_classes)].some((c) => c.includes('cWW'))).length)
  // invalid / non-existent residues: explicit message, empty table
  await search.fill('A37'); await expect(page.getByTestId('isearch-error')).toContainText('is G37'); await expect(rows).toHaveCount(0)
  await search.fill('9999'); await expect(page.getByTestId('isearch-error')).toContainText('No residue 9999')
  // pairs with annotations in hidden types are reported, not silently dropped
  await search.fill('G32-A33'); await expect(page.getByTestId('isearch-hidden')).toContainText('more in types not shown')
  // Rfam numbering: the same pair is found by its Rfam positions
  await page.getByTestId('inumbering-select').selectOption(''); await search.fill('41-66'); await expect(page.getByTestId(`irow-${KL}`)).toBeVisible()
})

test('junction variation is visible through the region filter; kissing loop through the tertiary filter', async ({ page }) => {
  await open(page)
  await page.getByTestId('ifilter-more-toggle').click()
  await page.getByTestId('ifilter-region').selectOption('J2/3')
  const j23 = ds.pairs.filter((p) => p.category === 'basepair' && ((p as unknown as { region_1: string }).region_1 === 'J2/3' || (p as unknown as { region_2: string }).region_2 === 'J2/3'))
  await expect(page.locator('[data-testid^="irow-"]')).toHaveCount(j23.length)
  expect(j23.filter((p) => p.n_observed / p.n_eligible < 0.9).length).toBeGreaterThan(j23.length / 2)
  await page.getByTestId('ifilter-region').selectOption('')
  await page.getByTestId('ifilter-structure').selectOption('tertiary')
  await expect(page.getByTestId(`irow-${KL}`)).toBeVisible()
  await expect(page.getByTestId('irow-basepair:26-56')).toBeVisible()                              // junction U22-A52 (J1/2 x J2/3)
})

test('compare 5E54 B (apo) with 5SWE (adenine-bound): same layout, junction differences, local base-anchored superposition', async ({ page }) => {
  test.setTimeout(240_000)
  await open(page)
  await page.getByTestId('imode-compare').click()
  await expect(page.getByTestId('iemph-diff')).toHaveAttribute('aria-pressed', 'true')                    // differences emphasised by default
  await page.getByTestId('icompare-a').selectOption('5e54_1_B_B'); await page.getByTestId('icompare-b').selectOption('5swe_1_X_A')
  await expect(page.getByTestId('icompare-auto')).toHaveCount(0)
  // identical layout: the shared kissing-loop arc has the same path in both diagrams
  expect(await page.getByTestId(`idiagram-a-arc-${KL}`).getAttribute('d')).toBe(await page.getByTestId(`idiagram-b-arc-${KL}`).getAttribute('d'))
  for (const p of [1, 41, 102]) {
    const a = page.getByTestId(`idiagram-a-pos-${p}`), b = page.getByTestId(`idiagram-b-pos-${p}`)
    expect([await a.getAttribute('cx'), await a.getAttribute('cy')]).toEqual([await b.getAttribute('cx'), await b.getAttribute('cy')])
  }
  // Stagno et al. 2017: bound junction pairs only in 5SWE; apo pairs only in 5E54 B (Rfam ids; 4FE5 B labels U22-A52, A23-G46, A23-U74 ...)
  for (const id of ['basepair:26-56', 'basepair:27-50']) await expect(page.getByTestId(`icmp-row-${id}`)).toHaveAttribute('data-status', 'only_b')
  await expect(page.getByTestId('icmp-row-basepair:27-80')).toHaveAttribute('data-status', 'only_a')
  await expect(page.getByTestId(`icmp-row-${KL}`)).toHaveAttribute('data-status', 'shared')
  // colour = comparison status (not category); pattern = category
  const statusColor: Record<string, string> = { shared: '#a8a29e', only_a: '#D55E00', only_b: '#D55E00', reclassified: '#CC79A7', not_comparable: '#c9b98f' }
  await page.getByTestId('itypes-toggle').click(); await page.getByTestId('itypes-all').click(); await page.keyboard.press('Escape')
  for (const el of await page.locator('[data-testid="idiagram-a-arcs"] path[data-status]').all()) {
    const [stt, cat, color, dash] = await Promise.all(['data-status', 'data-category', 'data-color', 'data-dash'].map((x) => el.getAttribute(x)))
    expect(color).toBe(statusColor[stt!]); expect(dash).toBe(cat === 'multiple' ? '' : CATEGORY_STYLE[cat as 'basepair'].dash ?? '')
  }
  await page.getByTestId('itypes-toggle').click(); await page.getByTestId('itypes-reset').click(); await page.keyboard.press('Escape')
  // status filter
  await page.getByTestId('icmp-filter-shared').click()
  await expect(page.locator('[data-testid^="icmp-row-"][data-status="shared"]')).toHaveCount(0)
  await page.getByTestId('icmp-filter-shared').click()
  // 3D: exactly two chains; whole-RNA superposition until an interaction is selected
  const v = page.getByTestId('interaction-pair-viewer')
  await expect(v).toHaveAttribute('data-status', 'ready', { timeout: 90_000 }); await expect(v).toHaveAttribute('data-structures', '5e54_1_B_B,5swe_1_X_A')
  await expect(page.getByTestId('icompare-fit')).toHaveAttribute('data-state', 'none')
  const fetches = await page.evaluate(() => (window.__rnaExplorer?.interactionViewer as { fetches: number }).fetches)
  await page.getByTestId('icmp-row-basepair:26-56').click()
  await expect(page.getByTestId('icompare-fit')).toHaveAttribute('data-state', 'local')
  await expect(page.getByTestId('icompare-fit')).toContainText('Locally aligned on')
  const hook = await page.evaluate(() => window.__rnaExplorer?.interactionViewer as { fetches: number; fit: { rmsd: number; nAtoms: number }; geometry: Record<string, { raw: number; shown: number }>; sticks: Record<string, number[]> })
  expect(hook.fit.nAtoms).toBeGreaterThanOrEqual(6); expect(hook.fit.rmsd).toBeLessThan(0.5)
  // the pair's own geometry is preserved (rigid move of the whole chain fitted on ONE base)
  for (const g of Object.values(hook.geometry)) { expect(g.raw).not.toBeNull(); expect(Math.abs(g.shown - g.raw)).toBeLessThan(1e-6) }
  // the absent pair is still shown in A (both nucleotides resolved)
  expect(hook.sticks['5e54_1_B_B']).toEqual(await labelSeqs(page, '5e54_1_B_B', [26, 56]))
  expect(hook.sticks['5swe_1_X_A']).toEqual(await labelSeqs(page, '5swe_1_X_A', [26, 56]))
  // another interaction: no chain is fetched again
  await page.getByTestId(`icmp-row-${KL}`).click(); await expect(page.getByTestId('idiagram-a')).toHaveAttribute('data-selected', KL)
  await expect.poll(() => page.evaluate(() => (window.__rnaExplorer?.interactionViewer as { sticks: Record<string, number[]> }).sticks['5swe_1_X_A'].length)).toBe(2)
  expect(await page.evaluate(() => (window.__rnaExplorer?.interactionViewer as { fetches: number }).fetches)).toBe(fetches)
  // back to family conservation keeps the selection
  await page.getByTestId('icompare-back').click()
  await expect(page.getByTestId('interactions-tab')).toHaveAttribute('data-mode', 'family'); await expect(page.getByTestId(`irow-${KL}`)).toHaveAttribute('aria-selected', 'true')
})

test('families without interactions.json show a short notice', async ({ page }) => {
  await page.route('**/families/RF00167/interactions.json', (r) => r.fulfill({ status: 404, body: 'not found' }))
  await page.goto('/RF00167'); await page.getByTestId('tab-interactions').click()
  await expect(page.getByTestId('interactions-unavailable')).toBeVisible()
})

test('compare two: default pair with the highest valid discrepancy; a manual pair survives mode and tab switches', async ({ page }) => {
  await open(page)
  const sc = await (await page.request.get('/families/RF00167/structure_clustering.json')).json() as StructureClusteringDataset
  const rows = await (await page.request.get('/families/RF00167/instances.json')).json() as { seq_name: string; resolution: number | null }[]
  const want = defaultComparePair(ds as unknown as InteractionsDataset, sc, (sn) => rows.find((r) => r.seq_name === sn)?.resolution ?? null)!
  await page.getByTestId('imode-compare').click()
  await expect(page.getByTestId('icompare-a')).toHaveValue(want.a); await expect(page.getByTestId('icompare-b')).toHaveValue(want.b)
  await expect(page.getByTestId('icompare-auto')).toContainText('highest geometric discrepancy'); await expect(page.getByTestId('icompare-auto')).toContainText('≤ 3.0 Å')
  await expect(page.getByTestId('icompare-coverage-a')).toBeVisible(); await expect(page.getByTestId('icompare-coverage-b')).toBeVisible()
  // the comparison table lists differences first
  const st = await page.getByTestId('icompare-tbody').locator('tr').evaluateAll((rs) => rs.map((r) => r.getAttribute('data-status')))
  const rank = ['only_a', 'only_b', 'reclassified', 'not_comparable', 'shared']
  expect(st.map((x) => rank.indexOf(x!))).toEqual([...st.map((x) => rank.indexOf(x!))].sort((a, b) => a - b))
  // manual pair: kept across Family / Compare and across tabs
  await page.getByTestId('icompare-b').selectOption('5swe_1_X_A')
  await page.getByTestId('icompare-back').click(); await page.getByTestId('imode-compare').click()
  await expect(page.getByTestId('icompare-b')).toHaveValue('5swe_1_X_A'); await expect(page.getByTestId('icompare-auto')).toHaveCount(0)
  await page.getByTestId('tab-structures').click(); await page.getByTestId('tab-interactions').click(); await page.getByTestId('imode-compare').click()
  await expect(page.getByTestId('icompare-a')).toHaveValue(want.a); await expect(page.getByTestId('icompare-b')).toHaveValue('5swe_1_X_A')
})

for (const clen of [299, 300, 301]) {
  test(`circular diagrams only up to 300 consensus positions (clen ${clen})`, async ({ page }) => {
    const real = await (await page.request.get('/families/RF00167/interactions.json')).json() as InteractionsDataset & { clen: number }
    ds ??= real as unknown as DS
    await page.route('**/families/RF00167/interactions.json', (r) => r.fulfill({ json: { ...real, clen, position_regions: [...real.position_regions, ...Array(clen - real.position_regions.length).fill(null)],
      instances: real.instances.map((i) => ({ ...i, eligibility: i.eligibility.padEnd(clen, '-') })) } }))
    await page.goto('/RF00167'); await page.getByTestId('tab-interactions').click()
    await expect(page.getByTestId('ifamily-table')).toBeVisible(); await page.getByTestId('ipage-size').selectOption('0')
    const on = clen <= 300
    await expect(page.getByTestId('interactions-tab')).toHaveAttribute('data-circular', String(on))
    await expect(page.getByTestId('idiagram')).toHaveCount(on ? 1 : 0); await expect(page.getByTestId('icircular-off')).toHaveCount(on ? 0 : 1)
    // every interaction is still listed; the 3D viewer and comparison remain
    await expect(page.locator('[data-testid^="irow-"]')).toHaveCount(real.pairs.filter((p) => p.category === 'basepair').length)
    await expect(page.getByTestId('molstar-viewer')).toBeVisible()
    await page.getByTestId('imode-compare').click()
    await expect(page.getByTestId('icompare-table')).toBeVisible(); await expect(page.getByTestId('interaction-pair-viewer')).toBeVisible()
    await expect(page.getByTestId('icompare-diagrams')).toHaveCount(on ? 1 : 0)
  })
}

test('table and 3D viewer share width, height and top edge in both modes', async ({ page }) => {
  await open(page)
  const same = async (viewer: string) => {
    const t = (await page.getByTestId('itable-box').boundingBox())!, v = (await page.getByTestId(viewer).boundingBox())!
    for (const k of ['y', 'width', 'height'] as const) expect(Math.abs(t[k] - v[k])).toBeLessThan(1.5)
  }
  await same('molstar-viewer')
  await page.getByTestId('imode-compare').click(); await expect(page.getByTestId('interaction-pair-viewer')).toBeVisible()
  await same('interaction-pair-viewer')
  await expect(page.getByTestId('iemph-diff')).toHaveText('Show differences')
  await expect(page.getByTestId('icmp-filter-only_a')).toHaveText(/^Only in A \(\d+\)$/)
  await expect(page.getByTestId('icmp-filter-reclassified')).toHaveText(/^Different class \(\d+\)$/)
  await expect(page.getByTestId('icmp-filter-reclassified')).toHaveAttribute('title', /different FR3D class/)
})

test('pagination: 15 rows per page by default; first/previous/next/last and a typed page; paging keeps the selection; a selection elsewhere shows its page', async ({ page }) => {
  await open(page, false)
  const bps = ds.pairs.filter((p) => p.category === 'basepair'), n = Math.ceil(bps.length / 15)
  const order = [...bps].sort(famOrder)
  const rows = page.locator('[data-testid^="irow-"]'), cur = page.getByTestId('ipage-input')
  await expect(rows).toHaveCount(15); await expect(cur).toHaveValue('1'); await expect(page.getByTestId('ipage-count')).toHaveText(String(n))
  await expect(page.getByTestId('ipage-first')).toBeDisabled(); await expect(page.getByTestId('ipage-prev')).toBeDisabled()
  await page.getByTestId('ipage-next').click(); await expect(cur).toHaveValue('2')
  await expect(rows.first()).toHaveAttribute('data-testid', `irow-${order[15].id}`)
  // extremes
  await page.getByTestId('ipage-last').click(); await expect(cur).toHaveValue(String(n))
  await expect(page.getByTestId('ipage-last')).toBeDisabled(); await expect(page.getByTestId('ipage-next')).toBeDisabled()
  await expect(rows.last()).toHaveAttribute('data-testid', `irow-${order[order.length - 1].id}`)
  await page.getByTestId('ipage-first').click(); await expect(cur).toHaveValue('1')
  // typed page (clamped to the range)
  await cur.fill('2'); await cur.press('Enter'); await expect(cur).toHaveValue('2')
  await cur.fill('999'); await cur.press('Enter'); await expect(cur).toHaveValue(String(n))
  // select on page 2, page away: the selection stays
  await page.getByTestId('ipage-first').click(); await page.getByTestId('ipage-next').click()
  await rows.first().click(); const sel = order[15].id
  await page.getByTestId('ipage-prev').click(); await expect(page.getByTestId('idiagram')).toHaveAttribute('data-selected', sel)
  // a selection from the diagram shows its page
  const last = order[order.length - 1]
  await page.locator(`[data-testid="idiagram"] [data-hit="${last.id}"]`).dispatchEvent('click')
  await expect(page.getByTestId(`irow-${last.id}`)).toBeVisible(); await expect(cur).toHaveValue(String(n))
  // sorting after a selection starts at page 1 of the new order (not the selected row's page); the selection is kept
  await page.getByTestId('isort-entries').click(); await page.getByTestId('isort-entries').click()   // descending, then ascending
  await expect(cur).toHaveValue('1'); await expect(page.getByTestId('isort-entries')).toHaveText('Entries\u00a0▲')
  const minEntries = Math.min(...bps.map((p) => p.n_entries_observed))
  await expect(rows.first().locator('td').nth(4)).toHaveText(String(minEntries))
  await expect(page.getByTestId('idiagram')).toHaveAttribute('data-selected', last.id)
  await page.getByTestId('isort-reset').click(); await expect(cur).toHaveValue('1')
  // rows per page
  await page.getByTestId('ipage-size').selectOption('50'); await expect(rows).toHaveCount(Math.min(50, bps.length))
  await page.getByTestId('ipage-size').selectOption('15')
  // compare table is paged too
  await page.getByTestId('imode-compare').click()
  await expect(page.getByTestId('icompare-tbody').locator('tr')).toHaveCount(15)
})

test('16S (RF00177, 1,533 consensus positions): table and 3D only, paged; search and comparison work', async ({ page }) => {
  test.setTimeout(240_000)
  const d = await (await page.request.get('/families/RF00177/interactions.json')).json() as InteractionsDataset
  await page.goto('/RF00177'); await page.getByTestId('tab-interactions').click()
  await expect(page.getByTestId('itable-box')).toBeVisible({ timeout: 60_000 })
  await expect(page.getByTestId('interactions-tab')).toHaveAttribute('data-circular', 'false')
  await expect(page.getByTestId('idiagram')).toHaveCount(0); await expect(page.getByTestId('icircular-off')).toBeVisible()
  const bp = d.pairs.filter((p) => p.category === 'basepair').length
  await expect(page.getByTestId('ifamily-count')).toHaveText(`${bp} of ${bp}`)
  await expect(page.locator('[data-testid^="irow-"]')).toHaveCount(15); await expect(page.getByTestId('ipage-count')).toHaveText(String(Math.ceil(bp / 15)))
  await page.getByTestId('ipage-last').click(); await expect(page.getByTestId('ipage-input')).toHaveValue(String(Math.ceil(bp / 15)))
  // search through the 4V9D Y numbering: a residue resolves, and Enter selects its first interaction
  await page.getByTestId('isearch').fill('G1053'); await expect(page.getByTestId('isearch-match')).toContainText('Involving G1053')
  await page.getByTestId('isearch').press('Enter'); await expect(page.getByTestId('idetails')).toBeVisible()
  await expect(page.getByTestId('molstar-viewer')).toBeVisible()
  // comparison without diagrams; default pair from the stored discrepancy
  await page.getByTestId('isearch').fill(''); await page.getByTestId('imode-compare').click()
  await expect(page.getByTestId('icompare-auto')).toBeVisible(); await expect(page.getByTestId('icompare-diagrams')).toHaveCount(0)
  await expect(page.getByTestId('icompare-tbody').locator('tr')).toHaveCount(15)
  await expect(page.getByTestId('interaction-pair-viewer')).toHaveAttribute('data-structures', /^[^,]+,[^,]+$/)
  await expect(page.getByTestId('interaction-pair-viewer')).toHaveAttribute('data-polymer-style', 'cartoon')   // the thin trace is for the Structures tab only
})

test('annotation coverage counts near base pairs (9V50 A, 3.25 Å: 13 true + 6 near of 20) and flags a sparsely paired chain (tRNA 1J2B C)', async ({ page }) => {
  await open(page)
  await page.getByTestId('imode-compare').click()
  await page.getByTestId('icompare-a').selectOption('5e54_1_B_B'); await page.getByTestId('icompare-b').selectOption('9v50_1_A_A')
  await expect(page.getByTestId('icompare-coverage-a')).toContainText('19 of 19 eligible (100%)')
  await expect(page.getByTestId('icompare-coverage-b')).toContainText('19 of 20 eligible (95%)')
  await expect(page.getByTestId('icompare-coverage')).not.toContainText('sparse'); await expect(page.getByTestId('icompare-coverage')).toContainText('family median 100%')
  // a near pair is "Near" in that chain and counts as paired in the comparison (never as missing)
  await expect(page.getByTestId('icmp-row-basepair:20-86')).toHaveAttribute('data-status', 'reclassified')   // cWW in A, ncWW in B
  await page.goto('/RF00005'); await page.getByTestId('tab-interactions').click(); await expect(page.getByTestId('itable-box')).toBeVisible()
  await page.getByTestId('imode-compare').click()
  await page.getByTestId('icompare-a').selectOption('1ehz_1_A_A'); await page.getByTestId('icompare-b').selectOption('1j2b_1_C_A')
  await expect(page.getByTestId('icompare-coverage-b')).toContainText('16 of 21 eligible (76%)'); await expect(page.getByTestId('icompare-coverage-b')).toContainText('sparse annotation')
})

test('an interaction with several FR3D classes expands into a per-class breakdown', async ({ page }) => {
  await open(page)
  const p = ds.pairs.find((x) => x.category === 'basepair' && Object.keys((x as unknown as { classes: Record<string, number> }).classes).length > 1)!
  const classes = (p as unknown as { classes: Record<string, number> }).classes
  await expect(page.getByTestId(`iclass-row-${p.id}-${Object.keys(classes)[0]}`)).toHaveCount(0)
  await page.getByTestId(`iclass-toggle-${p.id}`).click()
  await expect(page.getByTestId(`irow-${p.id}`)).toHaveAttribute('aria-selected', 'false')            // the toggle does not select
  for (const [c, n] of Object.entries(classes)) await expect(page.getByTestId(`iclass-cons-${p.id}-${c}`)).toHaveText(new RegExp(`^${n}/${p.n_eligible} `))
  await page.getByTestId(`iclass-toggle-${p.id}`).click()
  await expect(page.locator(`[data-testid^="iclass-row-${p.id}-"]`)).toHaveCount(0)
  // single-class interactions have no toggle
  const one = ds.pairs.find((x) => x.category === 'basepair' && Object.keys((x as unknown as { classes: Record<string, number> }).classes).length === 1 && !Object.keys(x.near_classes).length)!
  await expect(page.getByTestId(`iclass-toggle-${one.id}`)).toHaveCount(0)
})

test('compare context: resolution next to coverage; same RNA in different constructs is flagged (TPP 2GDI vs 2HOJ)', async ({ page }) => {
  await page.goto('/RF00059'); await page.getByTestId('tab-interactions').click(); await expect(page.getByTestId('itable-box')).toBeVisible()
  await page.getByTestId('imode-compare').click()
  await page.getByTestId('icompare-a').selectOption('2gdi_1_X_A'); await page.getByTestId('icompare-b').selectOption('2hoj_1_A_A')
  await expect(page.getByTestId('icompare-resolution-a')).toHaveText(' (X-ray, 2.05 Å)'); await expect(page.getByTestId('icompare-resolution-b')).toHaveText(' (X-ray, 2.5 Å)')
  await expect(page.getByTestId('icompare-sequence')).toHaveAttribute('data-kind', 'construct')
  await expect(page.getByTestId('icompare-sequence')).toContainText('same RNA on the consensus'); await expect(page.getByTestId('icompare-sequence')).toContainText('different constructs')
  // short meaning in the info tip; the exact rules (terminal-construct exception, 1% threshold) in "Sequence classification rules"
  await page.getByTestId('icompare-sequence').getByTestId('info-tip').click()
  await expect(page.getByTestId('info-popover')).toContainText('Exact classification rules are in the methods.'); await page.keyboard.press('Escape')
  const rules = page.getByTestId('icompare-sequence-methods'); await rules.locator('summary').click()
  await expect(rules).toContainText('differ only at up to 3 positions within 5 positions of a chain end'); await expect(rules).toContainText('at most 1% of shared positions differ')
  await page.getByTestId('icompare-b').selectOption('8f4o_1_A_A'); await page.getByTestId('icompare-a').selectOption('2hoj_1_A_A')
  await expect(page.getByTestId('icompare-sequence')).toHaveAttribute('data-kind', 'identical')
  await page.getByTestId('icompare-a').selectOption('2cky_1_A_A')
  await expect(page.getByTestId('icompare-sequence')).toHaveAttribute('data-kind', 'different')
})

test('near base pairs: counted as paired, listed as their own class, "Near" in the displayed structure', async ({ page }) => {
  await open(page)
  const p = ds.pairs.find((x) => x.id === 'basepair:20-86')!
  expect(p.n_near_only).toBeGreaterThan(0)
  await expect(page.getByTestId(`icons-${p.id}`)).toContainText(`${p.n_observed + p.n_near_only}/${p.n_eligible}`)
  await page.getByTestId(`iclass-toggle-${p.id}`).click()
  await expect(page.getByTestId(`iclass-cons-${p.id}-ncWW`)).toHaveText(new RegExp(`^${p.near_classes.ncWW}/${p.n_eligible} `))
  await page.getByTestId('idisplayed-select').selectOption('9v50_1_A_A')
  await expect(page.getByTestId(`irow-${p.id}`)).toHaveAttribute('data-status', 'near')
  await page.getByTestId(`irow-${p.id}`).click()
  await expect(page.getByTestId('idetails-status')).toHaveAttribute('data-state', 'near'); await expect(page.getByTestId('idetails-status')).toContainText('ncWW')
  await expect(page.getByTestId('idetails-conservation')).toContainText(`${p.n_observed} true, ${p.n_near_only} near`)
  await expect(page.getByTestId('iviewer-caption')).toContainText('near base pair')
  await expect(page.getByTestId('interactions-provenance-near')).toContainText('fr3d-python 3feb9b9')
})

test('Rfam consensus column and filter: Both / 2D / 3D; 2D pairs get a row and a grey dashed arc (TPP 64-83)', async ({ page }) => {
  const tpp = await (await page.request.get('/families/RF00059/interactions.json')).json() as DS
  await page.goto('/RF00059'); await page.getByTestId('tab-interactions').click(); await expect(page.getByTestId('idiagram')).toBeVisible()
  await page.getByTestId('ipage-size').selectOption('0')
  const bps = tpp.pairs.filter((p) => p.category === 'basepair'), by = (c: string) => bps.filter((p) => p.consensus === c)
  expect(by('2d').length).toBeGreaterThan(0); expect(by('both').length + by('2d').length + by('3d').length).toBe(bps.length)
  await expect(page.getByTestId('iconsensus-basepair:64-83')).toHaveAttribute('data-consensus', '2d')
  await expect(page.getByTestId('icons-basepair:64-83')).toContainText(`0/${tpp.pairs.find((p) => p.id === 'basepair:64-83')!.n_eligible}`)
  const arc = page.getByTestId('idiagram-arc-basepair:64-83')
  await expect(arc).toHaveAttribute('data-consensus', '2d'); await expect(arc).toHaveAttribute('data-dash', '2 3'); await expect(arc).toHaveAttribute('data-color', '#a8a29e')
  await expect(page.getByTestId('ilegend-2d')).toBeVisible()
  await page.getByTestId('ifilter-more-toggle').click()
  for (const c of ['2d', 'both', '3d'] as const) {
    await page.getByTestId('ifilter-consensus').selectOption(c)
    await expect(page.locator('[data-testid^="irow-"]')).toHaveCount(by(c).length)
    await expect(page.getByTestId('idiagram')).toHaveAttribute('data-arc-count', String(by(c).length))
  }
  // the filter applies to base pairs only: other categories are hidden while it is set
  await page.getByTestId('itypes-toggle').click(); await page.getByTestId('itypes-all').click(); await page.getByTestId('itypes-toggle').click()
  await expect(page.locator('[data-testid^="irow-"]')).toHaveCount(by('3d').length)
  await page.getByTestId('ifilter-reset').click()   // the More filters panel is still open
  await page.getByTestId('irow-basepair:64-83').click()
  await expect(page.getByTestId('idetails-consensus')).toHaveAttribute('data-consensus', '2d')
  // sort by the column: Both, then 2D, then 3D
  await page.getByTestId('itypes-toggle').click(); await page.getByTestId('itypes-reset').click(); await page.getByTestId('itypes-toggle').click()
  await page.getByTestId('isort-consensus').click()
  const order = await page.locator('[data-testid^="iconsensus-"]').evaluateAll((es) => es.map((e) => e.getAttribute('data-consensus')))
  expect(order).toEqual([...order].sort((a, b) => ['both', '2d', '3d'].indexOf(a!) - ['both', '2d', '3d'].indexOf(b!)))
})
