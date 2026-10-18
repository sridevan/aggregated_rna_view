import { expect, test, type Page } from '@playwright/test'

/* Family Overview tooltip wording: structured per-position tooltips (heading · main value · optional muted detail) and short
   information tooltips; detail in Legend & help. Expected numbers come from consensus_stats.json, never from the rendering. */
type Stats = { sequence_level: { n_unique_sequences: number; n_match: number[]; coverage: (number | null)[]; n_deletion: number[]; n_uncovered: number[]; pp_mean: (number | null)[]; pp_n: number[] }
  structure_level: { n_instances: number; n_match_with_coordinates: number[]; n_match: number[]; structural_coverage: (number | null)[] } }
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 }) }
const cell = (page: Page, row: string, pos: number, strip = 'consensus-tracks') => page.getByTestId(strip).locator(`[data-row="${row}"][data-tip-pos="${pos}"]`)
const pct = (f: number) => `${(f * 100).toFixed(1).replace(/\.0$/, '')}%`
const tipText = async (page: Page) => ({ title: await page.getByTestId('track-tooltip-title').textContent(), main: await page.getByTestId('track-tooltip-main').textContent(),
  secondary: (await page.getByTestId('track-tooltip-secondary').count()) ? await page.getByTestId('track-tooltip-secondary').textContent() : null })
const inViewport = async (page: Page, testId: string) => { const b = (await page.getByTestId(testId).boundingBox())!, v = page.viewportSize()!; expect(b.x).toBeGreaterThanOrEqual(0); expect(b.y).toBeGreaterThanOrEqual(0); expect(b.x + b.width).toBeLessThanOrEqual(v.width); expect(b.y + b.height).toBeLessThanOrEqual(v.height); expect(b.width).toBeLessThanOrEqual(361) }
let S: Stats
const open = async (page: Page) => { await page.goto('/RF00167'); await ready(page, '4fe5'); S = await (await page.request.get('/families/RF00167/consensus_stats.json')).json() }

test('sequence and structural coverage: position, values and denominators unchanged; secondary detail only when relevant', async ({ page }) => {
  await open(page)
  const sl = S.sequence_level, st = S.structure_level
  await cell(page, 'seq', 41).hover()
  expect(await tipText(page)).toEqual({ title: 'Sequence coverage · Position 41', main: `${sl.n_match[40]}/${sl.n_unique_sequences} unique sequences (${pct(sl.coverage[40]!)})`, secondary: null })
  expect(sl.n_deletion[40] + sl.n_uncovered[40]).toBe(0)
  // a position with deletions or uncovered sequences gets the short secondary line
  const p = sl.coverage.findIndex((c, i) => c != null && c >= 0.5 && sl.n_deletion[i] + sl.n_uncovered[i] > 0) + 1
  await cell(page, 'seq', p).hover()
  const t = await tipText(page)
  expect(t.main).toBe(`${sl.n_match[p - 1]}/${sl.n_unique_sequences} unique sequences (${pct(sl.coverage[p - 1]!)})`)
  if (sl.n_deletion[p - 1]) expect(t.secondary).toMatch(new RegExp(`deleted in ${sl.n_deletion[p - 1]}`, "i"))
  if (sl.n_uncovered[p - 1]) expect(t.secondary).toMatch(new RegExp(`not covered by the alignment in ${sl.n_uncovered[p - 1]}`, "i"))
  // structural coverage 40: 80/84 (95.2%), 84 contain the nucleotide in sequence
  await cell(page, 'str', 40).hover()
  expect(await tipText(page)).toEqual({ title: 'Structural coverage · Position 40', main: `${st.n_match_with_coordinates[39]}/${st.n_instances} RNA instances (${pct(st.structural_coverage[39]!)})`, secondary: `${st.n_match[39]} instances contain the nucleotide in sequence` })
  expect(await tipText(page)).toEqual({ title: 'Structural coverage · Position 40', main: '80/84 RNA instances (95.2%)', secondary: '84 instances contain the nucleotide in sequence' })
  const same = st.n_match.findIndex((n, i) => st.structural_coverage[i] != null && n === st.n_match_with_coordinates[i]) + 1
  await cell(page, 'str', same).hover(); expect((await tipText(page)).secondary).toBeNull()                               // no redundant detail
  await inViewport(page, 'track-tooltip')
  await page.mouse.move(5, 5); await expect(page.getByTestId('track-tooltip')).toHaveCount(0)
  expect(await page.evaluate(() => (window as unknown as { __rnaExplorer: { selection: { active: unknown } } }).__rnaExplorer.selection.active)).toBeNull()   // hover never selects
})

test('selected chain: available, unavailable, deleted and outside the aligned region are labelled distinctly', async ({ page }) => {
  await open(page)
  // key-value rows: Rfam position, author residue of the selected chain (or why there is none), status (wording in position-numbering.spec.ts)
  const rows = async () => Object.fromEntries(await Promise.all(['rfam', 'author', 'status'].map(async (k) => [k, await page.getByTestId(`track-tooltip-value-${k}`).textContent()])))
  await cell(page, 'chain', 42).hover()
  await expect(page.getByTestId('track-tooltip-title')).toHaveText('Selected RNA chain · 4FE5 B')
  expect((await rows()).status).toBe('Coordinates available'); expect((await rows()).rfam).toBe('42')
  await cell(page, 'chain', 5).hover()
  expect(await rows()).toEqual({ rfam: '5', author: 'Unmapped', status: 'Outside aligned region' })
  for (const [sn, pdb, pos, status, author] of [['2xnw_1_A_A', '2xnw', 19, 'Coordinates unavailable', null], ['1y26_1_X_A', '1y26', 64, 'Deletion', 'Not present']] as const) {
    await page.getByTestId('instance-select').selectOption(sn); await ready(page, pdb)
    const m = await (await page.request.get(`/families/RF00167/mappings/${sn}.json`)).json()
    const col = m.consensus_columns[pos - 1]; expect(col.state).toBe(status === 'Deletion' ? 'deletion' : 'match')
    await cell(page, 'chain', pos).hover()
    await expect(page.getByTestId('track-tooltip-title')).toHaveText(`Selected RNA chain · ${pdb.toUpperCase()} ${sn.split('_')[2]}`)
    const r = await rows(); expect(r.status).toBe(status); expect(r.rfam).toBe(String(pos))
    if (author) expect(r.author).toBe(author)
    else { const res = m.residues.find((x: { polymer_seq_position: number }) => x.polymer_seq_position === col.polymer_seq_position); expect(r.author).toBe(`${res.residue_name}${res.author_residue_number}${res.author_insertion_code}`) }
  }
  // the information tooltip: short text; symbols in Legend & help
  await page.getByTestId('instance-track-info').hover()
  await expect(page.getByTestId('instance-track-info-text')).toHaveText('Shows nucleotide availability at Rfam match positions. +N markers indicate insertions between positions.')
  await inViewport(page, 'instance-track-info-card')
  // clicking a column still selects (and clicking again clears)
  await cell(page, 'chain', 42).click(); await expect(page.getByTestId('sel-highlightable')).toContainText('42 →')
  await cell(page, 'chain', 42).click(); await expect.poll(() => page.evaluate(() => (window as unknown as { __rnaExplorer: { selection: { active: unknown } } }).__rnaExplorer.selection.active)).toBeNull()
})

test('information tooltips are short with a heading; details in Legend & help; loop and reference rows keep structured tooltips', async ({ page }) => {
  await open(page)
  // RNA loop annotations (collapsed and expanded)
  await page.getByRole('button', { name: 'About the RNA loop annotations' }).hover()
  await expect(page.getByTestId('info-popover-title')).toHaveText('RNA loop annotations')
  await expect(page.getByTestId('info-popover')).toHaveText('RNA loop annotationsBGSU hairpin, internal and junction loops in the selected RNA chain. Click a block to select the entire loop, including discontinuous segments, in 2D and 3D.')
  await inViewport(page, 'info-popover'); await page.mouse.move(5, 5)
  await page.getByTestId('loops-toggle').click()
  await page.getByTestId('loop-block-HL_4FE5_002-0').hover()
  await expect(page.getByTestId('loop-tooltip-title')).toHaveText('Hairpin loop (HL)'); await expect(page.getByTestId('loop-tooltip-counts')).toHaveText('8 nucleotides · 1 segment')
  await expect(page.getByTestId('loop-tooltip-motif')).toContainText('Motif: HL_56046.1'); await inViewport(page, 'loop-tooltip')
  await page.mouse.move(5, 5)
  // Reference and alignment details: short popover in both states; RF / SS_cons / PP rows structured
  await page.getByRole('button', { name: 'About the reference and alignment details' }).click()
  await expect(page.getByTestId('info-popover')).toHaveText('Reference and alignment detailsRfam consensus sequence, secondary-structure annotation and alignment confidence.')
  await page.keyboard.press('Escape')
  await page.getByTestId('reference-toggle').click()
  const fam = await (await page.request.get('/families/RF00167/family.json')).json()
  await cell(page, 'rf', 30, 'reference-tracks').hover()
  expect(await tipText(page)).toEqual({ title: 'Consensus sequence (RF) · Position 30', main: fam.consensus_sequence[29], secondary: null })
  await cell(page, 'ss', 30, 'reference-tracks').hover()
  expect(await tipText(page)).toEqual({ title: 'Secondary structure (SS_cons) · Position 30', main: fam.ss_cons[29], secondary: null })
  await cell(page, 'pp', 30, 'reference-tracks').hover()
  const pp = S.sequence_level.pp_mean[29]!
  expect(await tipText(page)).toEqual({ title: 'Alignment confidence · Position 30', main: `Mean PP ${pp.toFixed(3)} (high)`, secondary: `${S.sequence_level.pp_n[29]} aligned unique sequences` })
  // the right-most column at a narrow width stays inside the viewport
  await page.setViewportSize({ width: 700, height: 900 }); await cell(page, 'seq', 87).scrollIntoViewIfNeeded(); await cell(page, 'seq', 87).hover(); await inViewport(page, 'track-tooltip')
  await page.setViewportSize({ width: 1400, height: 1000 }); await page.mouse.move(5, 5)
  // Legend & help keeps the methodological detail that left the popovers
  await page.getByTestId('legend-toggle').click()
  await expect(page.getByTestId('legend-loops')).toContainText('separate lanes'); await expect(page.getByTestId('legend-loops')).toContainText('outlined in red')
  await expect(page.getByTestId('legend-loops')).toContainText('clicking it again clears'); await expect(page.getByTestId('legend-loops')).toContainText('without a consensus position')
  await expect(page.getByTestId('legend-reference')).toContainText('WUSS'); await expect(page.getByTestId('legend-reference')).toContainText('posterior probability')
  await expect(page.getByTestId('legend-reference')).toContainText('clan members'); await expect(page.getByTestId('legend-reference')).toContainText('upper case')
})
