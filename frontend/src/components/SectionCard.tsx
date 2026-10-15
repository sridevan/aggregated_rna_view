import { useId, useState, type ReactNode } from 'react'

export interface SectionCardProps {
  /** small category label above the title (optional; omit when it would add clutter) */
  eyebrow?: string
  /** section heading: the primary visual element */
  title: ReactNode
  /** subtitle / explanatory text under the heading */
  description?: ReactNode
  /** controls rendered at the right of the header (toggles, buttons) */
  actions?: ReactNode
  children?: ReactNode
  /** expandable section; `defaultOpen` controls the initial state */
  collapsible?: boolean
  defaultOpen?: boolean
  /** heading level for the semantic hierarchy (tabs use h1 for the family; sections default to h2) */
  headingLevel?: 2 | 3 | 4
  className?: string
  bodyClassName?: string
  testId?: string
}

/**
 * Reusable section container for the scientific-resource style: light border, white background, restrained typography.
 * Header = optional eyebrow, semibold ~18 px title, optional description, optional right-aligned actions.
 * Collapsible sections use a real button with aria-expanded/aria-controls so they are keyboard accessible.
 */
export function SectionCard({ eyebrow, title, description, actions, children, collapsible = false, defaultOpen = true, headingLevel = 2, className = '', bodyClassName = '', testId }: SectionCardProps) {
  const [open, setOpen] = useState(defaultOpen)
  const bodyId = useId()
  const Heading = (`h${headingLevel}`) as 'h2' | 'h3' | 'h4'
  const heading = <Heading className="text-lg font-semibold leading-snug text-stone-900">{title}</Heading>
  return (
    <section className={`flex flex-col rounded border border-stone-200 bg-white ${className}`} data-testid={testId} data-open={collapsible ? String(open) : undefined}>
      <header className="px-4 pt-2.5 pb-1.5">
        {eyebrow && <p className="mb-1"><span className="inline-block rounded border border-stone-200 bg-stone-100 px-1.5 py-px text-[10px] font-medium uppercase tracking-wide text-stone-600" data-testid="section-eyebrow">{eyebrow}</span></p>}
        <div className="flex flex-wrap items-center justify-between gap-x-4 gap-y-1">
          {collapsible ? (
            <button type="button" aria-expanded={open} aria-controls={bodyId} onClick={() => setOpen((o) => !o)} className="flex items-center gap-2 rounded text-left focus:outline-none focus:ring-2 focus:ring-blue-500">
              <span aria-hidden="true" className={`inline-block text-xs text-stone-500 transition-transform ${open ? 'rotate-90' : ''}`}>▶</span>
              {heading}
            </button>
          ) : heading}
          {actions && <div className="flex shrink-0 items-center gap-2">{actions}</div>}
        </div>
        {description && <div className="mt-0.5 text-xs leading-5 text-stone-600">{description}</div>}
      </header>
      {(!collapsible || open) && <div id={bodyId} className={`flex flex-1 flex-col px-4 pb-3 ${bodyClassName}`}>{children}</div>}
    </section>
  )
}
