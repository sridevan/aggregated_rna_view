import { expect, test, type Page } from '@playwright/test'
import { selectColumn } from './helpers'

/* Family Overview Mol*: "Complete structure" keeps the selected RNA chain prominent (opaque slate cartoon) while other nucleic-acid chains
   and proteins are subdued semi-transparent grey layers; ligands / non-standard residues / ions stay visible, water hidden. A mode switch
   keeps the camera and the selection; Reset fits the selected chain; another instance rebuilds the layers for its own chain. */
const SLATE = 0x64748b
type Layer = { label: string | null; hidden: boolean; alpha: number | null; theme: string | null; value: number | null; elements: number | null }
type H = { __rnaExplorer: { molstarLayers?: Layer[]; molstarComponents?: { label?: string; hidden: boolean }[]; getCamera: () => unknown; molstarModeApplied?: { mode: string; fit: boolean; n: number }
  lastMolstarSelection?: { labelSeqIds: number[] }; focusRepresentation?: number[]; lastView?: { atoms: number } } }
const ready = async (page: Page, pdb: string) => { await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-pdb', pdb); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-status', 'ready', { timeout: 90_000 }) }
const hook = <T,>(page: Page, f: (h: H['__rnaExplorer']) => T) => page.evaluate((src) => new Function('h', `return (${src})(h)`)((window as unknown as H).__rnaExplorer), f.toString()) as Promise<T>
const toComplete = async (page: Page) => {
  const n = await hook(page, (h) => h.molstarModeApplied?.n ?? 0)
  await page.getByTestId('mode-complete').click()
  await expect.poll(() => hook(page, (h) => h.molstarModeApplied?.n ?? 0), { timeout: 60_000 }).toBeGreaterThan(n)
}

test('tRNA with its synthetase (1C0A B): Rfam region by default; Complete structure layers, non-standard residues, water hidden; camera kept; selection kept; reset fits the chain', async ({ page }) => {
  test.setTimeout(240_000)
  await page.goto('/RF00005'); await ready(page, page.url().includes('RF00005') ? (await page.getByTestId('molstar-viewer').getAttribute('data-pdb'))! : '')
  await page.getByTestId('instance-select').selectOption('1c0a_1_B_A'); await ready(page, '1c0a')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-mode', 'region')                          // default
  // a nucleotide selected in region mode
  await selectColumn(page, 30)
  await expect.poll(() => hook(page, (h) => h.focusRepresentation?.length ?? 0)).toBeGreaterThan(0)
  const sel = await hook(page, (h) => h.lastMolstarSelection?.labelSeqIds), focus = await hook(page, (h) => h.focusRepresentation)
  // the camera after the focus animation has settled
  let cam = '', prev = '-'
  await expect.poll(async () => { prev = cam; cam = JSON.stringify(await hook(page, (h) => h.getCamera())); return cam === prev }, { intervals: [400] }).toBe(true)
  await toComplete(page)
  // camera unchanged (no zoom-out to the whole complex); selection and its focus representation unchanged
  expect(JSON.stringify(await hook(page, (h) => h.getCamera()))).toBe(cam)
  expect(await hook(page, (h) => h.molstarModeApplied)).toMatchObject({ mode: 'complete', fit: false })
  expect(await hook(page, (h) => h.lastMolstarSelection?.labelSeqIds)).toEqual(sel); expect(await hook(page, (h) => h.focusRepresentation)).toEqual(focus)
  // layers: the selected chain opaque slate; proteins subdued grey and semi-transparent
  const L = await hook(page, (h) => h.molstarLayers ?? [])
  const byLabel = Object.fromEntries(L.filter((l) => l.label && l.label !== 'Null').map((l) => [l.label!, l]))
  expect(byLabel['Selected RNA chain']).toMatchObject({ hidden: false, alpha: 1, theme: 'uniform', value: SLATE })
  expect(byLabel['Protein chains']).toMatchObject({ hidden: false, theme: 'uniform' }); expect(byLabel['Protein chains'].alpha!).toBeLessThan(0.5)
  expect(byLabel['Protein chains'].value).not.toBe(SLATE)
  // ligands / non-standard residues visible, water and the default polymer component hidden
  const C = await hook(page, (h) => h.molstarComponents ?? [])
  const vis = (label: string) => C.filter((c) => c.label === label).map((c) => !c.hidden)
  // non-standard residues visible; water, the default polymer component and the default ligand component hidden ("Show ligands" is off)
  expect(vis('Non-standard')).toContain(true); expect(vis('Water').every((v) => !v)).toBe(true); expect(vis('Polymer').every((v) => !v)).toBe(true)
  expect(vis('Ligand').every((v) => !v)).toBe(true)
  // Reset: the selected chain is fitted, not the whole deposited model
  await page.getByTestId('reset-orientation').click()
  await expect.poll(() => hook(page, (h) => h.lastView?.atoms ?? 0)).toBe(byLabel['Selected RNA chain'].elements)
  // back to Rfam region: layers hidden, the selection still there
  await page.getByTestId('mode-region').click(); await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-mode', 'region')
  await expect.poll(async () => (await hook(page, (h) => h.molstarLayers ?? [])).every((l) => l.hidden)).toBe(true)
  expect(await hook(page, (h) => h.lastMolstarSelection?.labelSeqIds)).toEqual(sel)
})

test('16S in a 70S ribosome (4V9D Y): other rRNAs and proteins subdued; switching instances in Complete structure rebuilds the layers for the new chain', async ({ page }) => {
  test.setTimeout(300_000)
  await page.goto('/RF00177'); await ready(page, '4v9d')
  await toComplete(page)
  const L = await hook(page, (h) => h.molstarLayers ?? [])
  const byLabel = Object.fromEntries(L.map((l) => [l.label!, l]))
  expect(byLabel['Selected RNA chain']).toMatchObject({ alpha: 1, value: SLATE }); expect(byLabel['Selected RNA chain'].elements!).toBeGreaterThan(30_000)
  for (const k of ['Other nucleic-acid chains', 'Protein chains']) { expect(byLabel[k].alpha!).toBeLessThan(0.5); expect(byLabel[k].value).not.toBe(SLATE); expect(byLabel[k].elements!).toBeGreaterThan(10_000) }
  // another instance while in Complete structure: new layers for the new chain, nothing stale
  await page.getByTestId('instance-select').selectOption('7k00_1_A_A'); await ready(page, '7k00')
  await expect(page.getByTestId('molstar-viewer')).toHaveAttribute('data-mode', 'complete')
  await expect.poll(async () => (await hook(page, (h) => h.molstarLayers ?? [])).find((l) => l.label === 'Selected RNA chain')?.elements ?? 0, { timeout: 90_000 }).toBeGreaterThan(30_000)
  const L2 = await hook(page, (h) => h.molstarLayers ?? [])
  expect(L2.filter((l) => l.label === 'Selected RNA chain')).toHaveLength(1)
  expect(L2.find((l) => l.label === 'Selected RNA chain')!.elements).not.toBe(byLabel['Selected RNA chain'].elements)
})

type Lig = { show: boolean; mode: string; built: boolean; nearVisible: boolean; otherVisible: boolean; near: string[] }
const ligCtl = (page: Page) => page.getByTestId('show-ligands-control')
const candidatesReady = (page: Page) => expect(ligCtl(page)).not.toHaveAttribute('data-candidates', '', { timeout: 60_000 })

test('Show ligands: off by default; chain-specific ligand instances near the selected chain in both modes; kept across modes, reset on instance switch', async ({ page }) => {
  test.setTimeout(240_000)
  await page.goto('/RF00167'); await ready(page, '4fe5'); await candidatesReady(page)
  await expect(page.getByTestId('show-ligands')).not.toBeChecked()
  const cam = JSON.stringify(await hook(page, (h) => h.getCamera()))
  await page.getByTestId('show-ligands').check()
  await expect(ligCtl(page)).toHaveAttribute('data-shown', 'HPA B101')                     // 4FE5 chain B: its own hypoxanthine copy only (not chain A's)
  await expect.poll(() => hook(page, (h) => (h as unknown as { molstarLigands?: Lig }).molstarLigands)).toMatchObject({ show: true, mode: 'region', nearVisible: true, otherVisible: false, near: ['HPA B101'] })
  // sticks, carbons in the chain colour (Family Overview slate 0x64748b), other atoms CPK
  await expect.poll(() => hook(page, (h) => (h as unknown as { molstarLigands?: { style: unknown } }).molstarLigands?.style)).toEqual({ type: 'ball-and-stick', aspect: 1, color: 'element-symbol', carbon: 'uniform', carbonValue: 0x64748b })
  expect(JSON.stringify(await hook(page, (h) => h.getCamera()))).toBe(cam)                   // toggling never moves the camera
  await toComplete(page)
  await expect(page.getByTestId('show-ligands')).toBeChecked()                                   // kept across the mode switch
  await expect.poll(() => hook(page, (h) => (h as unknown as { molstarLigands?: Lig }).molstarLigands)).toMatchObject({ show: true, mode: 'complete', nearVisible: true })
  // PDBe's default ligand component stays hidden (replaced by the toggle's layers)
  expect((await hook(page, (h) => h.molstarComponents ?? [])).filter((c) => c.label === 'Ligand').every((c) => c.hidden)).toBe(true)
  await page.getByTestId('show-ligands').uncheck()
  await expect.poll(() => hook(page, (h) => (h as unknown as { molstarLigands?: Lig }).molstarLigands)).toMatchObject({ show: false, nearVisible: false, otherVisible: false })
  await page.getByTestId('show-ligands').check()
  // another instance: the toggle starts off again and is re-evaluated
  await page.getByTestId('instance-select').selectOption('1y26_1_X_A'); await ready(page, '1y26'); await candidatesReady(page)
  await expect(page.getByTestId('show-ligands')).not.toBeChecked()
  await page.getByTestId('show-ligands').check(); await expect(ligCtl(page)).toHaveAttribute('data-shown', /^ADE X\d+$/)
})

test('Show ligands: TPP bound to 2GDI X only; additives excluded in 16S (7K00: paromomycin shown, polyamines not); disabled when no ligand is near (1EHZ tRNA)', async ({ page }) => {
  test.setTimeout(300_000)
  await page.goto('/RF00059'); await page.getByTestId('instance-select').selectOption('2gdi_1_X_A'); await ready(page, '2gdi'); await candidatesReady(page)
  await page.getByTestId('show-ligands').check(); await expect(ligCtl(page)).toHaveAttribute('data-shown', 'TPP X100')
  await page.goto('/RF00177'); await page.getByTestId('instance-select').selectOption('7k00_1_A_A'); await ready(page, '7k00'); await candidatesReady(page)
  expect(Number(await ligCtl(page).getAttribute('data-candidates'))).toBeGreaterThan(1)
  await page.getByTestId('show-ligands').check(); await expect(ligCtl(page)).toHaveAttribute('data-shown', 'PAR A1601')
  await page.goto('/RF00005'); await ready(page, '1ehz'); await candidatesReady(page)
  await expect(ligCtl(page)).toHaveAttribute('data-candidates', '0'); await expect(page.getByTestId('show-ligands')).toBeDisabled()
  await expect(ligCtl(page)).toHaveText('Show ligands'); await expect(ligCtl(page)).toHaveAttribute('title', /^No small-molecule ligand within 4 Å/)   // the note is in the tooltip (toolbar stays on one line)
})
