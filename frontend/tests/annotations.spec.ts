import { expect, test, type Page } from '@playwright/test'

/** Browser tests for the synchronised consensus / 2D / 3D selection (Milestone 2). Needs network access to PDBe for coordinates and 2D topology. */
const INSTANCE_1Y26 = '1y26_1_X_A'
const state = (page: Page) => page.evaluate(() => (window as unknown as { __rnaExplorer: Record<string, unknown> }).__rnaExplorer)
const waitViewers = async (page: Page, pdb: string, need2d = true) => {
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-status', need2d ? 'ready' : /ready|unavailable/)
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready')
  await expect(page.getByTestId('rna2d-viewer')).toHaveAttribute('data-pdb', pdb)
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb)
}
const toggle = async (page: Page, pos: number) => { await page.getByTestId('goto-input').fill(String(pos)); await page.getByTestId('goto-button').click() }
const rna2dSelectedIds = (page: Page) => page.evaluate(() => [...document.querySelectorAll('[data-testid="rna2d-viewer"] text.rnaviewEle')]
  .filter((e) => { const f = (e as SVGElement).getAttribute('fill') || ''; return f !== '' && f !== 'none' && f !== 'black' && f !== '#000000' }).map((e) => +(e.className as unknown as SVGAnimatedString).baseVal.split('_').pop()!))

test.describe('family overview and structures', () => {
  test('family statistics match the backend JSON', async ({ page }) => {
    await page.goto('/')
    const fam = await (await page.request.get('/families/RF00167/family.json')).json()
    await expect(page.getByTestId('stat-pdb-entries')).toHaveText(String(fam.counts.pdb_entries))
    await expect(page.getByTestId('stat-rna-instances')).toHaveText(String(fam.counts.instances))
    await expect(page.getByTestId('stat-unique-sequences')).toHaveText(String(fam.counts.unique_mapped_sequences))
    await expect(page.getByTestId('stat-model-length')).toHaveText(String(fam.clen))
    await expect(page.getByTestId('display-representative')).toContainText((fam.display_representative ?? fam.representative).pdb_id, { ignoreCase: true })
    for (const [m, n] of Object.entries(fam.experimental_methods)) await expect(page.getByTestId('methods-list')).toContainText(`${m}${n}`)
    await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready')
  })
  test('structures table rows correspond to instances.json and show pending panels', async ({ page }) => {
    await page.goto('/')
    await page.getByTestId('tab-structures').click()
    const inst = await (await page.request.get('/families/RF00167/instances.json')).json()
    await expect(page.locator('[data-testid="instance-table"] tbody tr')).toHaveCount(10)   // paginated, 10 per page
    await expect(page.getByTestId('table-count')).toContainText(`of ${inst.length} instances`)
    await expect(page.getByTestId('pending-panel')).toHaveCount(3)
    await page.getByTestId('instance-search').fill('1y26')
    await expect(page.locator('[data-testid="instance-table"] tbody tr')).toHaveCount(inst.filter((r: { pdb_id: string }) => r.pdb_id === '1y26').length)
    await page.getByTestId(`row-${INSTANCE_1Y26}`).click()
    await expect(page.getByTestId('structures-current')).toContainText('1Y26 chain X')
    await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', '1y26')
    await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready')
  })
})

test.describe('annotations: synchronised selection', () => {
  test.beforeEach(async ({ page }) => {
    page.on('pageerror', (e) => console.log('PAGEERROR', e.message, e.stack?.split('\n').slice(0, 4).join(' | ')))
    await page.goto('/')
    await page.getByTestId('tab-annotations').click()
    await page.getByTestId('instance-select').selectOption(INSTANCE_1Y26)
  })

  test('1y26 loads in both viewers; consensus 50 -> polymer 34 highlighted in 2D and 3D', async ({ page }) => {
    await waitViewers(page, '1y26')
    expect(await page.locator('[data-testid="rna2d-viewer"] text.rnaviewEle').count()).toBe(71)
    await toggle(page, 50)
    await expect(page.getByTestId('sel-highlightable')).toContainText('50 → G34 (auth 46')
    await expect.poll(() => rna2dSelectedIds(page)).toEqual([34])
    await expect.poll(async () => (await state(page)).lastMolstarSelection).toEqual({ pdbId: '1y26', labelAsymId: 'A', labelSeqIds: [34] })
  })

  test('clicking a nucleotide in 2D updates consensus state and the 3D selection', async ({ page }) => {
    await waitViewers(page, '1y26')
    await page.locator('[data-testid="rna2d-viewer"] text.rnaview_1y26_20').dispatchEvent('click')
    await expect.poll(async () => (await state(page)).selection).toMatchObject({ positions: [36] })
    await expect.poll(async () => (await state(page)).lastMolstarSelection).toEqual({ pdbId: '1y26', labelAsymId: 'A', labelSeqIds: [20] })
    await expect(page.getByTestId('sel-highlightable')).toContainText('36 → G20 (auth 32')
  })

  test('a 3D residue click (PDB.molstar.click event) updates consensus state and the 2D selection', async ({ page }) => {
    await waitViewers(page, '1y26')
    await page.evaluate(() => {   // same event shape as pdbe-molstar lociDetails emits for a residue click
      const ev = new MouseEvent('PDB.molstar.click', { bubbles: true }) as MouseEvent & { eventData?: unknown }
      ev.eventData = { entry_id: '1y26', entity_id: '1', auth_asym_id: 'X', seq_id: 10, residueNumber: 10, auth_seq_id: 22, ins_code: '' }
      document.querySelector('[data-testid="molstar-viewer"]')!.dispatchEvent(ev)
    })
    await expect.poll(async () => (await state(page)).selection).toMatchObject({ positions: [26] })
    await expect.poll(() => rna2dSelectedIds(page)).toEqual([10])
    // a click on a residue of another chain must be ignored
    await page.evaluate(() => {
      const ev = new MouseEvent('PDB.molstar.click', { bubbles: true }) as MouseEvent & { eventData?: unknown }
      ev.eventData = { entry_id: '1y26', entity_id: '2', auth_asym_id: 'B', seq_id: 1, auth_seq_id: 1 }
      document.querySelector('[data-testid="molstar-viewer"]')!.dispatchEvent(ev)
    })
    await page.waitForTimeout(300)
    expect((await state(page)).selection).toMatchObject({ positions: [26] })
  })

  test('switching instances preserves the consensus selection and re-translates it', async ({ page }) => {
    await waitViewers(page, '1y26')
    await toggle(page, 50); await toggle(page, 36)
    await expect(page.getByTestId('sel-highlightable')).toContainText('50 → G34')
    await page.getByTestId('instance-select').selectOption('4tzx_1_X_A')   // 69-nt adenine riboswitch, polymer 2-70 (consensus 50 -> polymer 34 here too)
    await waitViewers(page, '4tzx', false)
    expect((await state(page)).selection).toMatchObject({ positions: [50, 36], instance: '4tzx_1_X_A' })
    await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', '4tzx')
    await expect(page.getByTestId('sel-highlightable')).toContainText('50 → G34 (auth 46')
    await expect.poll(async () => (await state(page)).lastMolstarSelection).toEqual({ pdbId: '4tzx', labelAsymId: 'A', labelSeqIds: [20, 34] })
  })

  test('deleted, uncovered and unobserved positions are reported and never highlighted', async ({ page }) => {
    await page.getByTestId('instance-select').selectOption('9nbh_1_A_A')   // consensus 39,51,64,72 deleted; polymer 31-102 without coordinates
    await waitViewers(page, '9nbh', false)
    await toggle(page, 39); await toggle(page, 50); await toggle(page, 1)
    await expect(page.getByTestId('sel-deleted')).toHaveText('39')
    await expect(page.getByTestId('sel-unobserved')).toContainText('50 → ')
    await expect(page.getByTestId('sel-highlightable')).toContainText('1 → U16')
    await expect.poll(async () => (await state(page)).lastMolstarSelection).toEqual({ pdbId: '9nbh', labelAsymId: 'A', labelSeqIds: [16] })
    await page.getByTestId('instance-select').selectOption('6p2h_1_A_A')   // covers consensus 26-80 only
    await waitViewers(page, '6p2h', false)
    await expect(page.getByTestId('sel-uncovered')).toHaveText('1')
    await expect(page.getByTestId('sel-highlightable')).toContainText('39 → ')
    await expect(page.getByTestId('sel-highlightable')).toContainText('50 → ')
    const s = await state(page)
    expect(s.selection).toMatchObject({ positions: [39, 50, 1] })
    expect((s.lastMolstarSelection as { labelSeqIds: number[] }).labelSeqIds).toHaveLength(2)
  })
})
