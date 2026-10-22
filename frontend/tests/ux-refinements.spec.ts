import { expect, test, type Page } from '@playwright/test'
import { selectColumn } from './helpers'

/* Targeted UX refinements: compact family summary, selection hint, annotation evidence (units, numerators/denominators, missing vs absent),
   representative comparison summary, and the explicit "View in Family Overview" action. Every expected value is read from the family JSON. */
const J = async (page: Page, acc: string, file: string) => (await page.request.get(`/families/${acc}/${file}`)).json()
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 }) }
const pct = (n: number, d: number) => `${d ? Math.round((1000 * n) / d) / 10 : 0}%`
const pinned = (page: Page) => page.getByTestId('annotation-inspect')
const row = (page: Page, key: string) => pinned(page).getByTestId(`annotation-summary-${key}`)
async function pin(page: Page, track: string, pos: number) {
  const m = page.getByTestId(track).locator(`rect[data-pos="${pos}"]`); await m.scrollIntoViewIfNeeded(); await m.click()
  await expect(pinned(page)).toContainText(`Rfam residue ${pos}`)
}

test('family summary: description clamped with Show more / Show less, full text and source kept; short descriptions have no toggle', async ({ page }) => {
  const fam = await J(page, 'RF00167', 'family.json')
  await page.goto('/RF00167'); await ready(page, '4fe5')
  const text = page.getByTestId('family-about-text')
  const clamp = await text.evaluate((e) => ({ h: e.clientHeight, sh: e.scrollHeight, lh: parseFloat(getComputedStyle(e).lineHeight) }))
  expect(clamp.h).toBeLessThanOrEqual(3 * clamp.lh + 1); expect(clamp.sh).toBeGreaterThan(clamp.h)          // about three lines, more available
  await expect(text).toHaveText(fam.comment)                                                                  // the full Rfam text is in the page, unchanged
  await expect(page.getByTestId('rfam-link')).toHaveAttribute('href', fam.rfam_url)
  const toggle = page.getByTestId('family-about-toggle')
  await expect(toggle).toHaveText('Show more'); await expect(toggle).toHaveAttribute('aria-expanded', 'false')
  await toggle.click()
  await expect(toggle).toHaveText('Show less')
  const open = await text.evaluate((e) => ({ h: e.clientHeight, sh: e.scrollHeight }))
  expect(open.h).toBeGreaterThanOrEqual(open.sh - 1)
  await toggle.click(); await expect(toggle).toHaveText('Show more')
  // every statistic is still on the page with its unchanged value
  await expect(page.getByTestId('stat-pdb-entries')).toHaveText(String(fam.counts.pdb_entries))
  await expect(page.getByTestId('stat-rna-instances')).toHaveText(String(fam.counts.instances))
  await expect(page.getByTestId('stat-model-length')).toHaveText(String(fam.clen))
  await expect(page.getByTestId('stat-unique-sequences')).toHaveText(String(fam.counts.unique_mapped_sequences))
  await expect(page.getByTestId('stat-sequence-occupancy')).toHaveText(`${(fam.coverage_summary.sequence_coverage_mean * 100).toFixed(1)} %`)
  await expect(page.getByTestId('stat-structural-occupancy')).toHaveText(`${(fam.coverage_summary.structural_coverage_mean * 100).toFixed(1)} %`)
  // RF00177: one-line description, so no toggle
  await page.goto('/RF00177')
  await expect(page.getByTestId('family-about-text')).toBeVisible()
  await expect(page.getByTestId('family-about-toggle')).toHaveCount(0)
})

test('selection hint: shown only while nothing is selected; nucleotide and loop selection, toggling and Clear still work', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  const hint = page.getByTestId('selection-hint'), status = page.getByTestId('selection-status')
  await expect(hint).toHaveText('Click a nucleotide or loop to highlight it in 2D and 3D.')
  await expect(page.getByTestId('instance-toolbar').locator('input')).toHaveCount(0)                    // no text-entry field
  await selectColumn(page, 50)
  await expect(status).toContainText('Rfam position 50 → 4FE5 B'); await expect(hint).toHaveCount(0)
  await selectColumn(page, 50)                                                                             // click again: deselect
  await expect(status).toHaveCount(0); await expect(hint).toBeVisible()
  await page.getByTestId('loops-toggle').click()
  await page.getByTestId('loop-block-IL_4FE5_001-0').click()
  await expect(status).toContainText('Internal loop · IL_4FE5_001'); await expect(hint).toHaveCount(0)
  await page.getByTestId('clear-selection').first().click()
  await expect(hint).toBeVisible()
})

test('annotation details: measured quantity, numerator / denominator / unit, absence wording and missing data for every category', async ({ page }) => {
  const tr = await J(page, 'RF00167', 'tracks.json'), lig = await J(page, 'RF00167', 'ligands.json'), sasa = await J(page, 'RF00167', 'sasa_tracks.json')
  const pk = (await J(page, 'RF00167', 'pocket_tracks.json')).tracks.all_score_ge_0, md = await J(page, 'RF00167', 'modification_tracks.json')
  await page.goto('/RF00167'); await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-loops')).toBeVisible()
  const tip = page.getByTestId('annotation-tooltip')
  const toTab = async (t: 'ligands' | 'annotations') => { await page.mouse.move(5, 5); await page.getByTestId(`tab-${t}`).click(); await expect(page.getByTestId(t === 'ligands' ? 'category-ligands' : 'category-loops')).toBeVisible() }

  // loops (unit: RNA chain instance): k of n, the rest "not annotated by the source", the 4 instances without an observed nucleotide not counted
  const hl = tr.positions[52].hairpin_loops.denominators
  await pin(page, 'track-HL', 53)
  await expect(row(page, 'observed')).toContainText(`Observed in ${hl.n_positive_observable_and_evaluated} of ${hl.n_observable_and_evaluated} RNA instances (${pct(hl.n_positive_observable_and_evaluated, hl.n_observable_and_evaluated)})`)
  await expect(row(page, 'not-annotated')).toContainText(`${hl.n_not_annotated_by_source} RNA instances (not annotated by the source; not a confirmed absence)`)
  const noObs = tr.n_instances - hl.n_observable_instances
  await expect(row(page, 'not-counted')).toHaveText(noObs ? `Not counted${noObs} RNA instances: no observed nucleotide at this position` : 'Not countedNone')
  await expect(row(page, 'counting-unit')).toContainText('RNA chain instance')

  // ligands (Ligands tab; unit: PDB entry): numerator and denominator are entries with bound HPA that have an observed nucleotide here (tracks.json by_ccd_entries);
  // the numerator equals the contacts recorded in ligands.json (no HPA entry lacks a nucleotide at 55)
  await toTab('ligands')
  const hpa = new Set(lig.bound_molecules.filter((b: { in_default_track: boolean; ccd_id: string; consensus_positions: number[] }) => b.in_default_track && b.ccd_id === 'HPA' && b.consensus_positions.includes(55)).map((b: { pdb_id: string }) => b.pdb_id))
  const e55 = tr.positions[54].ligand_binding_sites.by_ccd_entries.HPA
  expect(e55.n_positive_observable_and_evaluated).toBe(hpa.size); expect(e55.n_observable_and_evaluated).toBe(lig.by_ccd.HPA.n_entries_in_default_track)
  await pin(page, 'track-lig-HPA', 55)
  await expect(row(page, 'recorded-contact')).toHaveText(`Recorded contact${hpa.size} of ${e55.n_observable_and_evaluated} observable PDB entries with bound HPA (${pct(hpa.size, e55.n_observable_and_evaluated)})`)
  await expect(row(page, 'not-counted')).toHaveText('Not countedNone')
  await expect(row(page, 'counting-unit')).toContainText('PDB entry (chains of one entry counted once)')
  await expect(row(page, 'no-contact')).toContainText('not a confirmed absence')
  expect(await page.getByTestId('track-lig-HPA').locator('rect[data-pos="55"]').getAttribute('data-value')).toBe(e55.frequency.toFixed(3))
  // GNG (3 entries): one entry has no observed nucleotide at 25, 81 and 82 and is not counted; recomputed here from the mappings' coordinates
  const gngEntries = [...new Set(lig.bound_molecules.filter((b: { in_default_track: boolean; ccd_id: string }) => b.in_default_track && b.ccd_id === 'GNG').map((b: { pdb_id: string }) => b.pdb_id))] as string[]
  expect(gngEntries.length).toBe(lig.by_ccd.GNG.n_entries_in_default_track)
  const inst = await J(page, 'RF00167', 'instances.json'), instList = (Array.isArray(inst) ? inst : inst.instances) as { pdb_id: string; seq_name: string }[]
  const observedAt = async (sn: string, p: number) => {
    const m = await J(page, 'RF00167', `mappings/${sn}.json`)
    const col = m.consensus_columns.find((c: { rfam_consensus_pos: number }) => c.rfam_consensus_pos === p)
    return col.state === 'match' && !!m.residues.find((r: { polymer_seq_position: number }) => r.polymer_seq_position === col.polymer_seq_position)?.has_coordinates
  }
  for (const [p, k, n] of [[25, 2, 2], [81, 2, 2], [82, 1, 2]]) {
    const contact = new Set(lig.bound_molecules.filter((b: { in_default_track: boolean; ccd_id: string; consensus_positions: number[] }) => b.in_default_track && b.ccd_id === 'GNG' && b.consensus_positions.includes(p)).map((b: { pdb_id: string }) => b.pdb_id))
    const observable: string[] = []
    for (const en of gngEntries) { let ok = false; for (const i of instList.filter((x) => x.pdb_id === en)) ok ||= await observedAt(i.seq_name, p); if (ok) observable.push(en) }
    expect([[...contact].filter((x) => observable.includes(x)).length, observable.length]).toEqual([k, n])
    expect(tr.positions[p - 1].ligand_binding_sites.by_ccd_entries.GNG).toMatchObject({ n_positive_observable_and_evaluated: k, n_observable_and_evaluated: n, n_not_observable: gngEntries.length - n })
    await pin(page, 'track-lig-GNG', p)
    await expect(row(page, 'recorded-contact')).toHaveText(`Recorded contact${k} of ${n} observable PDB entries with bound GNG (${pct(k, n)})`)
    await expect(row(page, 'not-counted')).toHaveText('Not counted1 PDB entry: no observed nucleotide at this position')
    await expect(pinned(page)).toContainText(`of the ${n} (of ${gngEntries.length}) PDB entries containing bound GNG that have an observed nucleotide here`)
    expect(await page.getByTestId('track-lig-GNG').locator(`rect[data-pos="${p}"]`).getAttribute('data-value')).toBe((k / n).toFixed(3))
  }
  await page.getByTestId('track-lig-GNG').locator('rect[data-pos="25"]').hover()
  await expect(tip).toContainText('Recorded contact in 2 of 2 observable PDB entries (100%)')
  // tooltips stay short (title + two lines)
  await page.getByTestId('track-lig-HPA').locator('rect[data-pos="55"]').hover()
  await expect(tip.locator('p')).toHaveCount(3)

  // SASA (continuous; unit Å²; median, IQR and contributing instances; no frequency wording)
  await toTab('annotations')
  const s53 = sasa.positions[52]
  await pin(page, 'track-sasa-median', 53)
  await expect(row(page, 'median')).toHaveText(`Median${s53.median_sasa_A2.toFixed(1)} Å² (interquartile range ${s53.q1_sasa_A2.toFixed(1)}–${s53.q3_sasa_A2.toFixed(1)} Å²) over ${s53.n_contributing_instances} RNA instances`)
  await expect(row(page, 'measurement')).toContainText('(Å²)')
  await expect(row(page, 'observed')).toHaveCount(0)
  const sNo = sasa.summary.n_instances_measured - s53.n_mapped_instances_with_coordinates
  await expect(row(page, 'not-counted')).toContainText(sNo ? `${sNo} RNA instances: no observed nucleotide at this position` : 'None')

  // pockets (prediction: "no pocket" is not a confirmed absence; failed or missing instances not counted, with reasons)
  await toTab('ligands')
  const p53 = pk[52]
  await pin(page, 'track-pocket-occurrence', 53)
  await expect(row(page, 'observed')).toContainText(`In a predicted pocket in ${p53.n_supporting_instances} of ${p53.n_evaluable_instances} RNA instances`)
  await expect(row(page, 'no-pocket')).toContainText(`${p53.n_evaluable_instances - p53.n_supporting_instances} RNA instances (prediction ran, no pocket here; not a confirmed absence of binding)`)
  await expect(row(page, 'not-counted')).toContainText(`${p53.unevaluable_reasons.no_coordinates} RNA instances: no coordinates`)

  // nonstandard nucleotides (deposited sequence: a standard nucleotide is a determined identity; no residue is not counted)
  await toTab('annotations')
  const m87 = md.positions[86]
  const ccc = m87.ccd.find((x: { ccd_id: string }) => x.ccd_id === 'CCC')
  await pin(page, 'track-mod-CCC', 87)
  await expect(row(page, 'observed')).toContainText(`Observed in ${ccc.n_instances} of ${m87.denominator} RNA instances`)
  await expect(row(page, 'other-identities')).toContainText(`${m87.n_negative} RNA instances with a standard nucleotide`)
  if (m87.n_no_residue) await expect(row(page, 'not-counted')).toContainText(`${m87.n_no_residue} RNA instance${m87.n_no_residue === 1 ? '' : 's'}: no residue at this position`)

  // evidence details are still behind the expandable section
  await pinned(page).locator('summary').click()
  await expect(pinned(page)).toContainText('Instances without a residue at this position are not counted.')
  // colour values unchanged
  expect(await page.getByTestId('track-HL').locator('rect[data-pos="53"]').getAttribute('data-value')).toBe(tr.positions[52].hairpin_loops.frequency.toFixed(3))
})

test('View in Family Overview: transfers the Rfam position, keeps family, instance and numbering, reports deletions; not an assertion of the annotation', async ({ page }) => {
  await page.goto('/RF00167'); await ready(page, '4fe5')
  await page.getByTestId('numbering-select').selectOption('1y26_1_X_A')
  await page.getByTestId('instance-select').selectOption('1y26_1_X_A'); await ready(page, '1y26')
  await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-loops')).toBeVisible()
  await expect(page.getByTestId('numbering-select')).toHaveValue('1y26_1_X_A')
  // 1Y26 X has a deletion at Rfam 64 (tested elsewhere): the transfer selects 64 and the Overview says nothing is highlighted
  await pin(page, 'track-sasa-median', 64)
  await page.getByTestId('view-in-overview').click()
  await expect(page.getByTestId('tab-overview')).toHaveAttribute('aria-selected', 'true')
  await expect(page).toHaveURL(/\/RF00167$/)
  await expect(page.getByTestId('instance-select')).toHaveValue('1y26_1_X_A')
  await expect(page.getByTestId('numbering-select')).toHaveValue('1y26_1_X_A')
  await expect(page.getByTestId('selection-status')).toContainText('Rfam position 64 is deleted in 1Y26 X: nothing highlighted')
  await expect(page.getByTestId('selection-origin')).toContainText('Opened from Annotations (Median SASA (Å²) at Rfam position 64). Family-level evidence: not necessarily observed in 1Y26 X.')
  // a mapped position: highlighted, the origin note follows the transferred position only
  await page.getByTestId('tab-ligands').click()
  await pin(page, 'track-lig-HPA', 55)
  await page.getByTestId('view-in-overview').click()
  await expect(page.getByTestId('selection-status')).toContainText('Rfam position 55 → 1Y26 X')
  await expect(page.getByTestId('selection-status')).toContainText('highlighted and centred in 2D/3D')
  await expect(page.getByTestId('selection-origin')).toContainText('HPA at Rfam position 55')
  await ready(page, '1y26')
  await selectColumn(page, 30)
  await expect(page.getByTestId('selection-status')).toContainText('Rfam position 30')
  await expect(page.getByTestId('selection-origin')).toHaveCount(0)
  // the Annotations tab itself never changed the selection: opening it again shows no pinned panel and the Overview keeps position 30
  await page.getByTestId('tab-annotations').click(); await expect(pinned(page)).toHaveCount(0)
  await page.getByTestId('tab-overview').click(); await expect(page.getByTestId('selection-status')).toContainText('Rfam position 30')
})

test('RF00177: View in Family Overview focuses the 2D on the transferred position; Clear returns to the whole 2D view', async ({ page }) => {
  await page.goto('/RF00177')
  await page.getByTestId('tab-annotations').click(); await expect(page.getByTestId('category-loops')).toBeVisible({ timeout: 60_000 })
  await pin(page, 'track-sasa-median', 700)
  await page.getByTestId('view-in-overview').click()
  await expect(page.getByTestId('selection-status')).toContainText('Rfam position 700 → 4V9D BA')
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 120_000 })
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-focused', 'true', { timeout: 30_000 })
  await page.getByTestId('clear-selection').first().click()
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-focused', 'false')
  await expect(page.getByTestId('selection-hint')).toBeVisible()
})

test('compact comparison header: one stored metric line; alignment details in the ⓘ; poor-fit chips kept; absent in single-structure mode', async ({ page }) => {
  const sc = await J(page, 'RF00167', 'structure_clustering.json'), orient = await J(page, 'RF00167', 'orientation.json')
  const ids: string[] = sc.dataset.instances, n = ids.length
  const stored = (metric: string, a: string, b: string) => { let i = ids.indexOf(a), j = ids.indexOf(b); if (i > j) [i, j] = [j, i]; return sc.metrics[metric].values_upper[i * n - (i * (i + 1)) / 2 + (j - i - 1)] as number }
  const rep = (metric: string, g: number) => sc.metrics[metric].clustering.groups.find((x: { group_id: number }) => x.group_id === g).representative as string
  const label = (sn: string) => { const [p, , a] = sn.split('_'); return `${p.toUpperCase()} ${a}` }
  const f3 = (v: number) => (Math.abs(v) >= 1 ? v.toFixed(2) : v.toFixed(3))
  const info = async () => { await page.getByRole('button', { name: 'About the 3D alignment' }).click(); const t = (await page.getByTestId('info-popover').innerText()); await page.keyboard.press('Escape'); return t }
  const fitLine = (sn: string) => { const t = orient.transforms[sn]; return `• ${label(sn)}: ${t.fit_quality} fit, RMSD ${t.rmsd.toFixed(2)} Å` }
  await page.goto('/RF00167'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click()
  await expect(page.getByTestId('alignment-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 60_000 })
  // single group: its caption, no comparison line
  await expect(page.getByTestId('compare-summary')).toHaveCount(0); await expect(page.getByTestId('alignment-caption')).toHaveText('Structural group 1 · Showing representative')
  expect(await info()).toContain(fitLine(rep('discrepancy', 1)))
  // groups 1 and 2 (discrepancy): only the metric line under the heading; representatives stay named by the controls and chips
  await page.getByTestId('compare-on').click()
  const [a, b] = [rep('discrepancy', 1), rep('discrepancy', 2)]
  await expect(page.getByTestId('alignment-caption')).toHaveText('Comparing representatives of structural groups 1 and 2')
  await expect(page.getByTestId('compare-summary')).toHaveText(`Geometric discrepancy ${f3(stored('discrepancy', a, b))} · ${sc.dataset.n_shared_positions} shared Rfam positions`)
  expect(+(await page.getByTestId('compare-summary').getAttribute('data-value'))!).toBeCloseTo(stored('discrepancy', a, b), 6)
  await expect(page.getByTestId('alignment-panel')).not.toContainText(`${label(a)} (Group 1) vs`)
  await expect(page.getByTestId('compare-box')).toContainText(label(a)); await expect(page.getByTestId('compare-select')).toContainText(label(b))
  await expect(page.getByTestId(`viewer-label-${a}`)).toBeVisible(); await expect(page.getByTestId(`viewer-label-${b}`)).toBeVisible()
  // good and poor overlays: the poor fit stays flagged on its chip; the ⓘ gives the reference, every fit and the metric distinction
  expect([orient.transforms[a].fit_quality, orient.transforms[b].fit_quality]).toEqual(['good', 'poor'])
  await expect(page.getByTestId(`fit-warning-${b}`)).toHaveText('poor fit'); await expect(page.getByTestId(`fit-warning-${a}`)).toHaveCount(0)
  const t = await info()
  expect(t).toContain(`fitted onto ${label(orient.reference.seq_name)} using C4′ atoms`)
  expect(t).toContain(`${fitLine(a)}\n${fitLine(b)}`)
  expect(t).toContain('may differ from the alignment used to calculate the geometric discrepancy')
  // the header is one short line: title, caption and metric line sit within 80 px above the comparison controls
  const head = (await page.getByTestId('alignment-panel').boundingBox())!, box = (await page.getByTestId('compare-box').boundingBox())!
  expect(box.y - head.y).toBeLessThan(80)   // title, the caption naming the compared groups, and the metric line
  await page.getByTestId('compare-off').click(); await expect(page.getByTestId('compare-summary')).toHaveCount(0)
  // TM-score: group 1 vs 3, then vs 4; the line follows the chosen pair and metric, the ⓘ names the TM-score
  await page.getByTestId('metric-tm').click()
  for (const g of [3, 4]) {
    if (await page.getByTestId('compare-off').count()) await page.getByTestId('compare-off').click()
    await page.getByTestId('compare-select').selectOption(String(g)); await page.getByTestId('compare-on').click()
    const [x, y] = [rep('tm', 1), rep('tm', g)], v = stored('tm', x, y)
    await expect(page.getByTestId('compare-summary')).toHaveText(`TM-score ${f3(v)} (1 − TM-score ${f3(1 - v)}) · ${sc.dataset.n_shared_positions} shared Rfam positions`)
    const ti = await info(); expect(ti).toContain(`${fitLine(x)}\n${fitLine(y)}`); expect(ti).toContain('calculate the TM-score')
  }
  // narrow window: the line wraps inside the panel, no horizontal page scroll
  await page.setViewportSize({ width: 760, height: 900 })
  const pb = (await page.getByTestId('alignment-panel').boundingBox())!, sb = (await page.getByTestId('compare-summary').boundingBox())!
  expect(sb.x + sb.width).toBeLessThanOrEqual(pb.x + pb.width + 1)
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true)
  // clustering results untouched: groups and representatives as stored
  await expect(page.getByTestId('group-representative')).toContainText(label(rep('tm', 1)))
})
