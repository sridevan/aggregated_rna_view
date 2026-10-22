export type TabId = 'overview' | 'structures' | 'interactions' | 'ligands' | 'annotations'
const TABS: { id: TabId; label: string }[] = [{ id: 'overview', label: 'Family Overview' }, { id: 'structures', label: 'Structures' }, { id: 'interactions', label: 'Interactions' }, { id: 'ligands', label: 'Ligands' }, { id: 'annotations', label: 'Annotations' }]
export function Tabs({ value, onChange }: { value: TabId; onChange: (t: TabId) => void }) {
  return (
    <nav className="flex gap-1 overflow-x-auto border-b border-stone-200" role="tablist">
      {TABS.map((t) => (
        <button key={t.id} role="tab" aria-selected={value === t.id} data-testid={`tab-${t.id}`} onClick={() => onChange(t.id)}
          className={`-mb-px shrink-0 whitespace-nowrap border-b-2 px-3 py-2 text-sm ${value === t.id ? 'border-blue-600 font-medium text-blue-700' : 'border-transparent text-stone-600 hover:text-stone-900'}`}>{t.label}</button>
      ))}
    </nav>
  )
}
