export type TabId = 'overview' | 'structures' | 'ligands' | 'interactions' | 'annotations'
const TABS: { id: TabId; label: string }[] = [{ id: 'overview', label: 'Family Overview' }, { id: 'structures', label: 'Structures' }, { id: 'ligands', label: 'Ligands' }, { id: 'interactions', label: 'Interactions' }, { id: 'annotations', label: 'Annotations' }]
export function Tabs({ value, onChange }: { value: TabId; onChange: (t: TabId) => void }) {
  return (
    // horizontal scrolling only (narrow screens), without a visible scrollbar: overflow-x alone made the 1 px tab underline a vertical overflow
    <nav className="flex gap-1 overflow-x-auto overflow-y-hidden border-b border-stone-200 [scrollbar-width:none] [&::-webkit-scrollbar]:hidden" role="tablist" data-testid="tab-bar">
      {TABS.map((t) => (
        <button key={t.id} role="tab" aria-selected={value === t.id} data-testid={`tab-${t.id}`} onClick={() => onChange(t.id)}
          className={`-mb-px shrink-0 whitespace-nowrap border-b-2 px-3 py-2 text-sm ${value === t.id ? 'border-blue-600 font-medium text-blue-700' : 'border-transparent text-stone-600 hover:text-stone-900'}`}>{t.label}</button>
      ))}
    </nav>
  )
}
