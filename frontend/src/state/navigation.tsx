import { createContext, useContext, type ReactNode } from 'react'

/**
 * Explicit cross-tab navigation (no continuous synchronisation between tabs). Provided by the family page shell, which owns the active tab.
 * `viewInOverview` transfers one Rfam consensus position to the Family Overview: it selects that position for the current RNA instance
 * (the instance and the position numbering are left as they are) and opens the tab. `origin` records where the position came from, so
 * the Overview can say that the annotation is family-level evidence and not necessarily present in the selected chain.
 */
export interface NavigationOrigin { from: 'annotations'; pos: number; label: string }
export interface NavigationApi { viewInOverview: (pos: number, label: string) => void; origin: NavigationOrigin | null }
const Ctx = createContext<NavigationApi | null>(null)
export function NavigationProvider({ value, children }: { value: NavigationApi; children: ReactNode }) {
  return <Ctx.Provider value={value}>{children}</Ctx.Provider>
}
/** null outside a family page shell (e.g. a component rendered on its own in a test): callers then hide navigation actions */
export const useNavigation = (): NavigationApi | null => useContext(Ctx)
