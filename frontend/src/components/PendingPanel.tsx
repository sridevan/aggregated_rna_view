export function PendingPanel({ title, detail }: { title: string; detail: string }) {
  return (
    <section className="rounded border border-dashed border-stone-300 bg-stone-100/60 p-4" data-testid="pending-panel">
      <h3 className="text-sm font-medium text-stone-700">{title}</h3>
      <p className="mt-1 text-sm text-stone-500">Unavailable: {detail}</p>
    </section>
  )
}
