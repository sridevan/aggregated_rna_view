import { expect, test } from '@playwright/test'

/* Home page: introduction, one card per available family with numbers read from that family's family.json, links to the family pages,
   and the app title on a family page links back home. */
const ACCS = ['RF00167', 'RF00162', 'RF00059', 'RF00005', 'RF00177', 'RF00011']

test('home lists every family with its real numbers and opens the family page', async ({ page }) => {
  await page.goto('/')
  await expect(page.getByTestId('home-intro')).toContainText('Rfam consensus coordinate system')
  await expect(page.getByTestId('family-cards').locator('[data-testid^="family-card-"]')).toHaveCount(ACCS.length)
  for (const acc of ACCS) {
    const f = await (await page.request.get(`/families/${acc}/family.json`)).json()
    const card = page.getByTestId(`family-card-${acc}`)
    await expect(card).toHaveAttribute('data-state', 'ready')
    await expect(card).toHaveAttribute('href', `/${acc}`)
    await expect(card.getByTestId('card-description')).toHaveText(`Rfam: ${f.description}`)
    await expect(card.getByTestId('card-entries')).toHaveText(f.counts.pdb_entries.toLocaleString('en-US'))
    await expect(card.getByTestId('card-chains')).toHaveText(f.counts.instances.toLocaleString('en-US'))
    await expect(card.getByTestId('card-clen')).toHaveText(f.clen.toLocaleString('en-US'))
    const rep = f.display_representative ?? f.representative
    await expect(card.getByTestId('card-representative')).toHaveText(`Shown by default: ${rep.pdb_id.toUpperCase()} chain ${rep.auth_asym_id}`)
    // curated selection flagged only where the build used one (16S)
    await expect(card.getByTestId('card-curated')).toHaveCount(f.provenance.selection ? 1 : 0)
  }
  // curated subsets: a small label with the explanation in its tooltip; the chain count comes from the data (provenance.selection)
  for (const [acc, n] of [['RF00005', '60'], ['RF00177', '60'], ['RF00011', '50']]) {
    const lab = page.getByTestId(`family-card-${acc}`).getByTestId('card-curated')
    await expect(lab).toHaveText('Curated subset'); await expect(lab).toHaveAttribute('data-n-chains', n)
    await expect(lab).toHaveAttribute('title', 'Selected experimental RNA structures; counts refer to the displayed dataset, not all available PDB structures for this family.')
  }
  await expect(page.getByText(/Curated selection of/)).toHaveCount(0)
  await expect(page.getByTestId('family-card-RF00011').getByTestId('card-chains')).toHaveText('50')
  await expect(page.getByTestId('page-parts')).toContainText('Family Overview'); await expect(page.getByTestId('page-parts')).toContainText('Structures'); await expect(page.getByTestId('page-parts')).toContainText('Annotations')
  await expect(page.getByTestId('home-footer')).toContainText('Rfam 15.1')
  // footer: data sources, then visualisation software, each linked
  for (const [name, href] of [['Rfam 15.1', 'https://rfam.org'], ['PDBe', 'https://www.ebi.ac.uk/pdbe/'], ['BGSU RNA 3D Motif Atlas', 'http://rna.bgsu.edu/rna3dhub/motifs'],
    ['PDBe Mol*', 'https://github.com/molstar/pdbe-molstar'], ['PDBe RNA Viewer', 'https://github.com/PDBeurope/pdb-rna-viewer']] as const)
    await expect(page.getByTestId('home-footer').getByRole('link', { name, exact: true })).toHaveAttribute('href', href)
  await expect(page.getByTestId('home-footer')).toHaveText(/^Data: .*Visualisation: /)
  // tab cards and notes
  await expect(page.getByRole('heading', { name: 'Explore an RNA family' })).toBeVisible()
  await expect(page.getByTestId('home-notes').locator('li')).toHaveText([/^Rfam correspondence: /, /^Residue numbering: /, /^Structural groups: /])
  await expect(page.getByTestId('family-card-RF00011')).toContainText('Catalytic RNA component of bacterial RNase P, which processes the 5′ ends of precursor tRNAs.')
  // only the four small family.json files are loaded on the home page
  // open a family, then come back via the app title
  await page.getByTestId('family-card-RF00059').click()
  await expect(page).toHaveURL(/\/RF00059$/); await expect(page.getByTestId('family-badge')).toHaveText('RF00059')
  await page.getByTestId('app-title').click()
  await expect(page).toHaveURL(/:5173\/$/); await expect(page.getByTestId('home-page')).toBeVisible()
})

test('home page loads only the family summaries and fits a phone screen', async ({ page }) => {
  const urls: string[] = []; page.on('request', (r) => urls.push(r.url()))
  await page.setViewportSize({ width: 390, height: 844 })
  await page.goto('/')
  for (const acc of ACCS) await expect(page.getByTestId(`family-card-${acc}`)).toHaveAttribute('data-state', 'ready')
  const data = urls.filter((u) => u.includes('/families/'))
  expect(data.every((u) => u.endsWith("/family.json"))).toBe(true); expect(new Set(data).size).toBe(ACCS.length)   // dev StrictMode may request twice
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1)).toBe(true)   // no horizontal scroll
})

test('cards are whole-card links: click anywhere, keyboard Enter, visible focus, consistent hover; intro spans the content width', async ({ page }) => {
  await page.goto('/')
  for (const acc of ACCS) await expect(page.getByTestId(`family-card-${acc}`)).toHaveAttribute('data-state', 'ready')
  await expect(page.getByText('Open family page')).toHaveCount(0)
  await expect(page.getByTestId('home-intro')).toHaveText('Explore experimentally determined RNA structures organised by Rfam family. Compare structural variation, nucleotide interactions and residue-level annotations using a shared Rfam consensus coordinate system.')
  await expect(page.getByTestId('page-parts').locator('h3')).toHaveText(['Family Overview', 'Structures', 'Interactions', 'Annotations'])   // in tab order
  await expect(page.getByTestId('home-intro').getByRole('link', { name: 'Rfam' })).toHaveAttribute('href', 'https://rfam.org')
  // intro as wide as the card grid (same container)
  const intro = (await page.getByTestId('home-intro').boundingBox())!, grid = (await page.getByTestId('family-cards').boundingBox())!
  expect(Math.abs(intro.x - grid.x)).toBeLessThan(1); expect(Math.abs(intro.width - grid.width)).toBeLessThan(2)
  // "Good to know" list spans the same width
  const notes = (await page.getByTestId('home-notes').locator('ul').boundingBox())!, sec = (await page.getByTestId('home-notes').boundingBox())!
  expect(Math.abs(sec.width - grid.width)).toBeLessThan(2); expect(notes.x + notes.width).toBeGreaterThan(grid.x + grid.width - 2)
  for (const acc of ACCS) {
    const card = page.getByTestId(`family-card-${acc}`)
    expect(await card.evaluate((e) => e.tagName)).toBe('A')
    expect(await card.locator('a, button, [tabindex]').count()).toBe(0)                    // no nested interactive elements (accession is a badge)
    expect(await card.evaluate((e) => getComputedStyle(e).cursor)).toBe('pointer')
  }
  // hover: same subtle change on every card (border colour + background), no movement
  const look = (acc: string) => page.getByTestId(`family-card-${acc}`).evaluate((e) => { const s = getComputedStyle(e); return { border: s.borderTopColor, bg: s.backgroundColor, transform: s.transform } })
  const before = await look('RF00167')
  const hovered: Record<string, unknown> = {}
  for (const acc of ACCS) { await page.getByTestId(`family-card-${acc}`).hover(); await page.waitForTimeout(250); hovered[acc] = await look(acc) }
  expect(hovered.RF00167).not.toEqual(before); expect((hovered.RF00167 as { transform: string }).transform).toBe('none')
  for (const acc of ACCS) expect(hovered[acc]).toEqual(hovered.RF00167)
  // no leftover gap where the link was: last content line close to the card bottom
  for (const acc of ['RF00167', 'RF00177']) {
    const card = (await page.getByTestId(`family-card-${acc}`).boundingBox())!
    const last = (await page.getByTestId(`family-card-${acc}`).locator('p').last().boundingBox())!
    expect(card.y + card.height - (last.y + last.height)).toBeLessThan(acc === 'RF00177' ? 24 : 60)   // equal-height rows: shorter card may end higher
  }
  // click on a non-text area (the stats block) navigates
  await page.getByTestId('family-card-RF00162').getByTestId('card-stats').click({ position: { x: 5, y: 5 } })
  await expect(page).toHaveURL(/\/RF00162$/); await expect(page.getByTestId('family-badge')).toHaveText('RF00162')
  // keyboard: Tab to a card, visible focus ring, Enter opens it
  await page.goto('/')
  await expect(page.getByTestId('family-card-RF00177')).toHaveAttribute('data-state', 'ready')
  let focused = ''
  for (let i = 0; i < 10 && focused !== 'family-card-RF00177'; i++) { await page.keyboard.press('Tab'); focused = await page.evaluate(() => document.activeElement?.getAttribute('data-testid') ?? '') }
  expect(focused).toBe('family-card-RF00177')
  expect(await page.getByTestId('family-card-RF00177').evaluate((e) => getComputedStyle(e).boxShadow)).toMatch(/rgb/)   // focus ring drawn
  await page.keyboard.press('Enter')
  await expect(page).toHaveURL(/\/RF00177$/)
  for (const acc of ACCS) { await page.goto(`/`); await page.getByTestId(`family-card-${acc}`).click(); await expect(page).toHaveURL(new RegExp(`/${acc}$`)) }
})

test('representative thumbnails: right chain, square and uncropped, part of the card link, placeholder when missing, no viewers', async ({ page }) => {
  const urls: string[] = []; page.on('request', (r) => urls.push(r.url()))
  await page.goto('/')
  const manifest = await (await page.request.get('/thumbnails/manifest.json')).json()
  for (const acc of ACCS) {
    const f = await (await page.request.get(`/families/${acc}/family.json`)).json(), rep = f.display_representative ?? f.representative
    // the image was generated from this family's default representative chain (label and author chain ids)
    expect(manifest[acc]).toMatchObject({ seq_name: rep.seq_name, pdb_id: rep.pdb_id, label_asym_id: rep.label_asym_id, auth_asym_id: rep.auth_asym_id, file: `${acc}.webp` })
    expect(manifest[acc].coordinate_source).toContain(`/${rep.pdb_id}/atoms?label_asym_id=${rep.label_asym_id}&model_nums=1`)
    expect(manifest[acc].bytes).toBeLessThan(40_000)
    const img = page.getByTestId(`card-thumbnail-${acc}`)
    await expect(img).toHaveAttribute('alt', new RegExp(`${rep.pdb_id.toUpperCase()} chain ${rep.auth_asym_id}$`))
    await expect.poll(() => img.evaluate((e: HTMLImageElement) => e.complete && e.naturalWidth)).toBe(240)
    expect(await img.evaluate((e: HTMLImageElement) => e.naturalHeight)).toBe(240)
    const b = (await img.boundingBox())!; expect(Math.abs(b.width - b.height)).toBeLessThan(1); expect(b.width).toBeGreaterThanOrEqual(96); expect(b.width).toBeLessThanOrEqual(130)
    // image inside the card link (one link, no nested interactive element)
    expect(await img.evaluate((e) => e.closest('a')?.getAttribute('data-testid'))).toBe(`family-card-${acc}`)
  }
  // equal-height cards in each row
  const h = await Promise.all(ACCS.map(async (a) => (await page.getByTestId(`family-card-${a}`).boundingBox())!.height))
  expect(Math.abs(h[0] - h[1])).toBeLessThan(1); expect(Math.abs(h[2] - h[3])).toBeLessThan(1)
  // no molecular viewer on the home page: no canvas, no 3Dmol / Mol* script, no coordinate downloads
  expect(await page.locator('canvas').count()).toBe(0)
  // (the app is one bundle, so Mol* code may be downloaded as a module; it is never initialised here)
  const fetched = urls.filter((u) => !/\/src\/|\/node_modules\/|\/assets\/|\/@vite|\/@react-refresh/.test(u))
  expect(fetched.some((u) => /3Dmol|model-server|\.cif|\.bcif|pdbe\/static/i.test(u))).toBe(false)
  expect(await page.locator('.msp-plugin, [data-testid="molstar-viewer"], [data-testid="alignment-viewer"]').count()).toBe(0)
  // clicking the thumbnail opens the family
  await page.getByTestId('card-thumbnail-RF00177').click()
  await expect(page).toHaveURL(/\/RF00177$/)
})

test('a missing thumbnail shows a neutral placeholder and the card still works', async ({ page }) => {
  await page.route('**/thumbnails/RF00059.webp', (r) => r.fulfill({ status: 404, body: '' }))
  await page.goto('/')
  await expect(page.getByTestId('card-thumbnail-placeholder-RF00059')).toBeVisible()
  await expect(page.getByTestId('card-thumbnail-RF00167')).toBeVisible()
  const ph = (await page.getByTestId('card-thumbnail-placeholder-RF00059').boundingBox())!, im = (await page.getByTestId('card-thumbnail-RF00167').boundingBox())!
  expect(Math.abs(ph.width - im.width)).toBeLessThan(1)
  await page.getByTestId('family-card-RF00059').click(); await expect(page).toHaveURL(/\/RF00059$/)
})
