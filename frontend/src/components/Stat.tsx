/** Summary statistic card: a fixed label row (two lines max), a value row shared across cards, and a description row.
 *  `grid-template-rows: 2.25rem auto 1fr` keeps every primary value on the same vertical position whatever the label length. */
export function Stat({ label, value, hint, testId }: { label: string; value: string | number; hint?: string; testId?: string }) {
  return (
    <div className="grid h-full grid-rows-[2.25rem_auto_1fr] rounded border border-stone-200 bg-white p-3" data-testid={`${testId ?? `stat-${label.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`}-card`}>
      <div className="self-start text-xs uppercase leading-[1.125rem] tracking-wide text-stone-500">{label}</div>
      <div className="text-xl font-semibold leading-7" data-testid={testId ?? `stat-${label.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`}>{value}</div>
      {hint ? <div className="mt-0.5 self-start text-xs leading-4 text-stone-500">{hint}</div> : <div />}
    </div>
  )
}
export { pct } from '../lib/format'
