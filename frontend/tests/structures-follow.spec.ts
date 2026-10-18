import { expect, test, type Page } from '@playwright/test'

/* Structures tab, Clustering mode: the 3D viewer shows one structure (the selected group's representative, or the instance chosen in the
   table, heatmap or dendrogram) or exactly two (the representatives of two groups via Compare representatives). Read-only labels and
   Reset; no member, visibility or remove controls. Groups come from structure_clustering.json (discrepancy: group 1 = 72 instances,
   group 2 = 4 instances). */
type G = { group_id: number; label: string; representative: string; members: string[] }
const viewer = (page: Page) => page.getByTestId('alignment-viewer')
const ready = async (page: Page) => expect(viewer(page)).toHaveAttribute('data-status', 'ready', { timeout: 90_000 })
const label = (sn: string) => { const [p, , a] = sn.split('_'); return `${p.toUpperCase()} ${a}` }
const clickRow = async (page: Page, sn: string) => {
  await page.getByTestId('instance-search').fill(sn.split('_')[0])
  await page.getByTestId(`row-${sn}`).locator('[data-testid="cell-title"]').click()                 // not the PDB / chain links
  await page.getByTestId('instance-search').fill('')
}
let groups: G[] = [], excluded: string[] = []
const open = async (page: Page) => {
  const d = await (await page.request.get('/families/RF00167/structure_clustering.json')).json()
  groups = d.metrics.discrepancy.clustering.groups; excluded = d.dataset.excluded.map((e: { seq_name: string }) => e.seq_name)
  await page.goto('/RF00167'); await page.getByTestId('tab-structures').click(); await page.getByTestId('view-clustering').click(); await ready(page)
}
/** viewer content, title, labels (colours as in the viewer) and caption agree */
async function expectShown(page: Page, sns: string[], caption: string) {
  await expect(viewer(page)).toHaveAttribute('data-structures', sns.join(',')); await ready(page)
  await expect(page.getByTestId('alignment-title')).toHaveText(sns.length === 2 ? '3D superposition' : '3D structure')
  await expect(page.getByTestId('viewer-labels')).toHaveAttribute('data-state', sns.length === 2 ? 'pair' : 'single')
  await expect(page.getByTestId('viewer-labels')).toContainText(sns.length === 2 ? 'Comparing:' : 'Displayed:')
  await expect(page.locator('[data-testid^="viewer-label-"]')).toHaveCount(sns.length)
  const swatches = await Promise.all(sns.map((sn) => page.getByTestId(`swatch-${sn}`).getAttribute('data-color')))
  expect(swatches.join(',')).toBe(await viewer(page).getAttribute('data-colors'))
  await expect(page.getByTestId('alignment-caption')).toHaveText(caption)
}

test('no member controls: default representative; another group shows its representative; a table row shows that RNA alone, its group selected', async ({ page }) => {
  await open(page)
  const [g1, g2] = groups
  for (const id of ['add-member', 'add-selected', 'reset-displayed']) await expect(page.getByTestId(id)).toHaveCount(0)
  await expect(page.locator('[data-testid^="visible-"], [data-testid^="remove-"], [data-testid^="legend-"]')).toHaveCount(0)
  await expectShown(page, [g1.representative], `${g1.label} · Showing representative`)
  await expect(page.getByTestId(`viewer-label-${g1.representative}`)).toContainText(`${label(g1.representative)}· Representative`)
  await expect(page.getByTestId('compare-summary')).toHaveCount(0)
  // another group: its representative alone (the 4-member group is no longer overlaid)
  await page.getByTestId('group-chip-2').click()
  await expectShown(page, [g2.representative], `${g2.label} · Showing representative`)
  // a table row of group 1: that RNA alone; group 1 selected again; caption says member
  const a = g1.members.find((m) => m !== g1.representative)!
  await clickRow(page, a)
  await expectShown(page, [a], `${g1.label} · Showing member ${label(a)}`)
  await expect(page.getByTestId('viewer-rep-label')).toHaveCount(0)
  await expect(page.getByTestId('group-chip-1')).toHaveAttribute('aria-pressed', 'true')
  await expect(page.getByTestId('structures-current')).toContainText(a.split('_')[0].toUpperCase())
  // another member replaces it (never accumulates)
  const b = g1.members.find((m) => ![g1.representative, a].includes(m))!
  await clickRow(page, b); await expectShown(page, [b], `${g1.label} · Showing member ${label(b)}`)
  // a member of group 2: that RNA alone, group 2 selected
  const m2 = g2.members.find((m) => m !== g2.representative)!
  await clickRow(page, m2); await expectShown(page, [m2], `${g2.label} · Showing member ${label(m2)}`)
  await expect(page.getByTestId('group-chip-2')).toHaveAttribute('aria-pressed', 'true')
  // a group's own representative chosen in the table: "Showing representative"
  await clickRow(page, g1.representative); await expectShown(page, [g1.representative], `${g1.label} · Showing representative`)
  // not clustered: alone, explained
  const x = excluded[0]
  await clickRow(page, x); await expectShown(page, [x], `${label(x)} is not clustered; shown on its own`)
  await expect(page.getByTestId('selected-not-clustered')).toBeVisible()
})

test('Compare representatives: exactly the two representatives (not the shown member); another group or Back returns to one; Reset restores the group representative', async ({ page }) => {
  await open(page)
  const [g1, g2] = groups
  const a = g1.members.find((m) => m !== g1.representative)!
  await clickRow(page, a); await expectShown(page, [a], `${g1.label} · Showing member ${label(a)}`)
  await page.getByTestId('compare-on').click()
  await expectShown(page, [g1.representative, g2.representative], `Comparing representatives of structural groups ${g1.group_id} and ${g2.group_id}`)
  await expect(page.getByTestId('compare-summary')).toContainText('shared Rfam positions')
  await expect(page.getByTestId('viewer-labels')).toContainText('↔')
  // selecting another group after a comparison: that group's representative alone
  await page.getByTestId('group-chip-2').click()
  await expectShown(page, [g2.representative], `${g2.label} · Showing representative`); await expect(page.getByTestId('compare-summary')).toHaveCount(0)
  // Back to group
  await page.getByTestId('group-chip-1').click(); await page.getByTestId('compare-on').click(); await expect(viewer(page)).toHaveAttribute('data-structures', `${g1.representative},${g2.representative}`)
  await page.getByTestId('compare-off').click(); await expectShown(page, [g1.representative], `${g1.label} · Showing representative`)
  // Reset: from a member and from a comparison back to the selected group's representative; the group is kept
  await clickRow(page, a); await page.getByTestId('viewer-reset').click()
  await expectShown(page, [g1.representative], `${g1.label} · Showing representative`); await expect(page.getByTestId('group-chip-1')).toHaveAttribute('aria-pressed', 'true')
  await page.getByTestId('group-chip-2').click(); await page.getByTestId('compare-on').click(); await expect(viewer(page)).toHaveAttribute('data-structures', `${g2.representative},${g1.representative}`)
  await page.getByTestId('viewer-reset').click()
  await expectShown(page, [g2.representative], `${g2.label} · Showing representative`); await expect(page.getByTestId('group-chip-2')).toHaveAttribute('aria-pressed', 'true')
})

test('rapid consecutive selections leave only the last structure; mode switches never show more than two', async ({ page }) => {
  await open(page)
  const g1 = groups[0], ms = g1.members.filter((m) => m !== g1.representative).slice(0, 4)
  for (const m of ms) { await page.getByTestId('instance-search').fill(m.split('_')[0]); await page.getByTestId(`row-${m}`).locator('[data-testid="cell-title"]').click() }
  await page.getByTestId('instance-search').fill('')
  await expectShown(page, [ms[3]], `${g1.label} · Showing member ${label(ms[3])}`)
  await expect.poll(async () => page.evaluate(() => document.querySelectorAll('[data-testid^="viewer-label-"]').length)).toBe(1)
  // Compare in Clustering, a pair in Structural Dissimilarity, back and forth: each mode keeps its own valid content (<= 2)
  await page.getByTestId('compare-on').click(); await expect(viewer(page)).toHaveAttribute('data-structures', /^[^,]+,[^,]+$/)
  await page.getByTestId('view-dissimilarity').click(); await expect(viewer(page)).toHaveAttribute('data-structures', /^[^,]+$/)
  await page.getByTestId('dm-cell-30-3').click(); await expect(viewer(page)).toHaveAttribute('data-structures', /^[^,]+,[^,]+$/)
  await page.getByTestId('view-clustering').click(); await expect(viewer(page)).toHaveAttribute('data-structures', `${g1.representative},${groups[1].representative}`)
  for (const v of ['view-dissimilarity', 'view-clustering', 'view-dissimilarity']) {
    await page.getByTestId(v).click(); expect(((await viewer(page).getAttribute('data-structures')) ?? '').split(',').length).toBeLessThanOrEqual(2)
  }
  await ready(page)
})
