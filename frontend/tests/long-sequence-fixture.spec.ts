import { expect, test, type Page, type Route } from '@playwright/test'

/* Long-family rendering check with a SYNTHETIC test fixture (not scientific data, not RF00177): the RF00167 family files are served under a
   test-only accession with every per-position array tiled to 1,533 consensus positions. Only the Family Overview track viewer is exercised. */
const N = 1533, ACC = 'RF99999', OUT = '../docs/screenshots'
const tile = <T,>(a: T[]) => Array.from({ length: N }, (_, i) => a[i % a.length])
const stretch = (o: unknown, n0: number): unknown => Array.isArray(o) ? (o.length === n0 && o.every((v) => v === null || typeof v !== 'object') ? tile(o) : o.map((v) => stretch(v, n0)))
  : o && typeof o === 'object' ? Object.fromEntries(Object.entries(o).map(([k, v]) => [k, stretch(v, n0)])) : o
async function serveFixture(page: Page) {
  await page.route(`**/families/${ACC}/**`, async (route: Route) => {
    const url = route.request().url(), real = url.replace(`/families/${ACC}/`, '/families/RF00167/')
    const r = await page.request.get(real)
    if (!r.ok() || (real.endsWith('.json') && !(r.headers()['content-type'] ?? '').includes('json'))) return route.fulfill({ status: 404, body: 'not found' })   // the dev server answers missing files with HTML
    const name = url.split(`/families/${ACC}/`)[1]
    if (!name.endsWith('.json')) return route.fulfill({ status: 200, body: await r.body() })
    let j = await r.json()
    if (name === 'family.json') j = { ...j, rfam_acc: ACC, clen: N, description: 'Synthetic 1,533-position test fixture (not real data)', consensus_sequence: tile([...j.consensus_sequence]).join(''), ss_cons: tile([...j.ss_cons]).join('') }
    if (name === 'consensus_stats.json') j = { ...(stretch(j, 102) as object), clen: N, rfam_acc: ACC }
    if (name.startsWith('mappings/')) j = { ...j, consensus_columns: tile(j.consensus_columns).map((c: { rfam_consensus_pos: number }, i: number) => ({ ...c, rfam_consensus_pos: i + 1 })) }
    return route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify(j) })
  })
}
const cellPx = async (page: Page) => +(await page.getByTestId('zoom-controls').getAttribute('data-cell-px'))!

test('1,533 synthetic positions: Fit all overview, zoom to nucleotide detail, aligned tracks, selection while scrolled, timing', async ({ page }) => {
  await serveFixture(page)
  const t0 = Date.now()
  await page.goto(`/${ACC}`)
  await expect(page.getByTestId('consensus-tracks')).toBeVisible({ timeout: 60_000 })
  const firstRender = Date.now() - t0
  const n = await page.getByTestId('consensus-tracks').locator('[data-row="seq"]').count(); expect(n).toBeGreaterThan(1400)   // all positions present
  // Fit all: the whole model in the visible width, no horizontal scrolling, no letters, no insertion markers, sparse axis labels
  const sc = page.getByTestId('tracks-scroll')
  expect(await sc.evaluate((e) => e.scrollWidth <= e.clientWidth + 1)).toBe(true)
  expect(await cellPx(page)).toBeLessThan(1)
  await page.getByTestId('reference-toggle').click()
  await expect(page.getByTestId('reference-row-rf')).toHaveAttribute('data-letters', 'false')
  const ticks = await page.getByTestId('consensus-axis').locator('text').allTextContents()
  expect(ticks.length).toBeLessThan(20); expect(ticks.map(Number).every((v, i, a) => i === 0 || v > a[i - 1])).toBe(true)
  await page.getByTestId('loops-toggle').click(); await expect(page.getByTestId('loop-tracks')).toBeVisible()
  await page.getByTestId('consensus-track-area').screenshot({ path: `${OUT}/m17_long_fixture_fit_all.png` })
  // zoom in step by step (timing each re-render) to letter-level detail; the axis numbering stays the consensus numbering
  const times: number[] = []
  while ((await page.getByTestId('reference-row-rf').getAttribute('data-letters')) !== 'true') {
    const before = await cellPx(page), s = Date.now()
    await page.getByTestId('zoom-in').click(); await expect.poll(() => cellPx(page)).toBeGreaterThan(before); times.push(Date.now() - s)
  }
  expect(Math.max(...times)).toBeLessThan(3000)
  expect(await sc.evaluate((e) => e.scrollWidth > 5 * e.clientWidth)).toBe(true)
  // scroll into the middle of the model: every strip moves together (one scroll container) and the columns still coincide
  const mid0 = Math.round(N / 2), cell0 = +(await page.getByTestId('consensus-tracks').getAttribute('data-cell'))!
  // centre position mid0 in the visible data area (right of the sticky label column)
  await sc.evaluate((e, [m, c]) => { e.scrollLeft = (140 + (m - 1) * c) * 1.2 - (140 * 1.2 + e.clientWidth) / 2 }, [mid0, cell0])
  const g = async (id: string) => { const s = page.getByTestId(id), b = (await s.boundingBox())!, vb = +(await s.getAttribute('viewBox'))!.split(' ')[2]; return { x: b.x, k: b.width / vb, cell: +(await s.getAttribute('data-cell'))! } }
  const mid = Math.round(N / 2)
  const xs = await Promise.all(['consensus-tracks', 'loop-tracks', 'reference-tracks', 'consensus-axis'].map(async (id) => { const q = await g(id); return q.x + (140 + (mid - 1) * q.cell + q.cell / 2) * q.k }))
  expect(Math.max(...xs) - Math.min(...xs)).toBeLessThan(0.5)
  // a click in the scrolled view selects the consensus position under the pointer (not a pixel offset)
  const tb = (await page.getByTestId('consensus-tracks').boundingBox())!
  await page.mouse.click(xs[0], tb.y + 3)
  await expect.poll(() => page.evaluate(() => (window as unknown as { __rnaExplorer: { selection: { positions: number[] } } }).__rnaExplorer.selection.positions)).toEqual([mid])
  await page.getByTestId('consensus-track-area').screenshot({ path: `${OUT}/m17_long_fixture_zoomed.png` })
  // back to Fit all keeps the selection
  await page.getByTestId('zoom-fit').click(); await expect.poll(() => cellPx(page)).toBeLessThan(1)
  await expect(page.getByTestId('consensus-tracks').locator(`rect[opacity="0.14"]`)).toHaveCount(1)
  console.log(`fixture: first render ${firstRender} ms; zoom steps ${times.join(', ')} ms`)
})
