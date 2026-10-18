import type { Page } from '@playwright/test'

/** Geometry of a Family Overview track strip: every strip shares x(pos) = left + (pos - 1) * cell (track units, read from the strip itself,
 *  so tests follow the current zoom) and is drawn at `scale` px per unit. */
export async function stripGeom(page: Page, testId: string) {
  const svg = page.getByTestId(testId)
  const b = (await svg.boundingBox())!, vb = +(await svg.getAttribute('viewBox'))!.split(' ')[2]
  const cell = +(await svg.getAttribute('data-cell'))!, left = +(await svg.getAttribute('data-left'))!
  return { svg, b, scale: b.width / vb, cell, left }
}
/** screen x of the centre of consensus column `pos` (or of the boundary after `pos` when `boundary`) in a strip */
export async function colX(page: Page, testId: string, pos: number, boundary = false): Promise<number> {
  const g = await stripGeom(page, testId)
  return g.b.x + (g.left + (pos - 1) * g.cell + (boundary ? g.cell : g.cell / 2)) * g.scale
}

/** Select (or, when already selected, clear) a consensus position the way a user does: click its column in the Family Overview tracks
 *  (top band of the primary strip, away from bars, chain cells and insertion markers). */
export async function selectColumn(page: Page, pos: number): Promise<void> {
  const g = await stripGeom(page, 'consensus-tracks')
  await g.svg.click({ position: { x: (g.left + (pos - 1) * g.cell + g.cell / 2) * g.scale, y: 2 } })
}
