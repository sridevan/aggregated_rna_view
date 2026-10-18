import { useCallback, useEffect, useId, useLayoutEffect, useRef, useState, type ReactNode } from 'react'
import { createPortal } from 'react-dom'

const PAD = 14           // viewport padding (px)
const GAP = 6            // gap between icon and popover
const MAX_W = 320        // popover max width (px)

type Placement = { left: number; top: number; maxHeight: number; side: 'below' | 'above' }

/**
 * Information icon with a popover. Hover/focus shows it transiently; click pins it (and a second click, Escape or an outside click closes it).
 * The popover is rendered through a portal with fixed positioning, so no scroll container or overflow-hidden parent can clip it; it opens
 * below the icon, flips above when there is not enough room, is shifted horizontally to stay inside the viewport (PAD), is limited to the
 * available height (scrolling internally) and is repositioned on window resize and on any scroll.
 */
export function InfoTip({ text, label = 'Information', title }: { text: string; label?: string; /** optional bold heading above the text */ title?: string }) {
  const id = useId()
  const btn = useRef<HTMLButtonElement>(null)
  const pop = useRef<HTMLSpanElement>(null)
  const [open, setOpen] = useState<false | 'hover' | 'pinned'>(false)
  const [pl, setPl] = useState<Placement | null>(null)

  const place = useCallback(() => {
    const r = btn.current?.getBoundingClientRect(); const el = pop.current
    if (!r) return
    const vw = window.innerWidth, vh = window.innerHeight
    const w = Math.min(MAX_W, vw - 2 * PAD)
    const h = el ? el.offsetHeight : 0
    const below = vh - r.bottom - GAP - PAD, above = r.top - GAP - PAD
    const side: Placement['side'] = h <= below || below >= above ? 'below' : 'above'
    const maxHeight = Math.max(80, side === 'below' ? below : above)
    const top = side === 'below' ? r.bottom + GAP : Math.max(PAD, r.top - GAP - Math.min(h, maxHeight))
    const left = Math.min(Math.max(r.left + r.width / 2 - w / 2, PAD), vw - PAD - w)
    setPl({ left, top, maxHeight, side })
  }, [])
  useLayoutEffect(() => { if (open) place() }, [open, text, title, place])
  useEffect(() => {
    if (!open) return
    const onDoc = (e: MouseEvent) => { if (!btn.current?.contains(e.target as Node) && !pop.current?.contains(e.target as Node)) setOpen(false) }
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') setOpen(false) }
    window.addEventListener('resize', place); window.addEventListener('scroll', place, true); document.addEventListener('mousedown', onDoc); document.addEventListener('keydown', onKey)
    return () => { window.removeEventListener('resize', place); window.removeEventListener('scroll', place, true); document.removeEventListener('mousedown', onDoc); document.removeEventListener('keydown', onKey) }
  }, [open, place])

  return (
    <span className="inline-block align-middle">
      <button ref={btn} type="button" aria-label={label} aria-describedby={open ? id : undefined} aria-expanded={!!open} data-testid="info-tip" data-pinned={open === 'pinned' ? 'true' : undefined}
        className="ml-1 inline-flex h-4 w-4 items-center justify-center rounded-full border border-stone-400 text-[10px] font-semibold normal-case leading-none text-stone-600 hover:bg-stone-200 focus:outline-none focus:ring-2 focus:ring-blue-500"
        onMouseEnter={() => setOpen((o) => o || 'hover')} onMouseLeave={() => setOpen((o) => (o === 'hover' ? false : o))} onFocus={() => setOpen((o) => o || 'hover')} onBlur={() => setOpen((o) => (o === 'hover' ? false : o))}
        onKeyDown={(e) => { if (e.key === 'Escape') setOpen(false) }} onClick={(e) => { e.stopPropagation(); setOpen((o) => (o === 'pinned' ? false : 'pinned')) }}>i</button>
      {open && createPortal(
        <span ref={pop} role="tooltip" id={id} data-testid="info-popover" data-side={pl?.side} style={{ left: pl?.left ?? -9999, top: pl?.top ?? -9999, width: Math.min(MAX_W, (typeof window !== 'undefined' ? window.innerWidth : MAX_W) - 2 * PAD), maxHeight: pl?.maxHeight, visibility: pl ? 'visible' : 'hidden' }}
          className="fixed z-[60] block overflow-y-auto whitespace-pre-line rounded border border-stone-300 bg-white p-2 text-left text-xs font-normal normal-case leading-4 tracking-normal text-stone-800 shadow">{title && <span className="block font-semibold" data-testid="info-popover-title">{title}</span>}{text}</span>,
        document.body)}
    </span>
  )
}

/**
 * Non-interactive hover card anchored to an element's rectangle (same look as the InfoTip popover). Rendered in a portal with fixed
 * positioning: below the anchor, flipped above when there is no room, shifted horizontally to stay PAD px inside the viewport.
 */
export function HoverCard({ anchor, children, testId = 'hover-card', interactive = false, onDismiss, ignoreSelector, label }: {
  anchor: DOMRect | null; children: ReactNode; testId?: string
  /** interactive popover (pointer events, dialog role): Escape or a mousedown outside it (and outside `ignoreSelector`) calls onDismiss */
  interactive?: boolean; onDismiss?: () => void; ignoreSelector?: string; label?: string
}) {
  const el = useRef<HTMLDivElement>(null)
  useEffect(() => {
    if (!interactive || !anchor || !onDismiss) return
    const onDoc = (e: MouseEvent) => { const t = e.target as Element; if (el.current?.contains(t) || (ignoreSelector && t.closest?.(ignoreSelector))) return; onDismiss() }
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') onDismiss() }
    document.addEventListener('mousedown', onDoc); document.addEventListener('keydown', onKey)
    return () => { document.removeEventListener('mousedown', onDoc); document.removeEventListener('keydown', onKey) }
  }, [interactive, anchor, onDismiss, ignoreSelector])
  const [pos, setPos] = useState<{ left: number; top: number; side: 'below' | 'above' } | null>(null)
  useLayoutEffect(() => {
    if (!anchor || !el.current) { setPos(null); return }
    const vw = window.innerWidth, vh = window.innerHeight
    const w = el.current.offsetWidth, h = el.current.offsetHeight
    const below = vh - anchor.bottom - GAP - PAD >= h
    const top = below ? anchor.bottom + GAP : Math.max(PAD, anchor.top - GAP - h)
    const left = Math.min(Math.max(anchor.left + anchor.width / 2 - w / 2, PAD), Math.max(PAD, vw - PAD - w))
    setPos({ left, top, side: below ? 'below' : 'above' })
  }, [anchor, children])
  if (!anchor) return null
  return createPortal(
    <div ref={el} role={interactive ? 'dialog' : 'tooltip'} aria-label={label} data-testid={testId} data-side={pos?.side} style={{ left: pos?.left ?? -9999, top: pos?.top ?? -9999, maxWidth: Math.min(340, window.innerWidth - 2 * PAD), visibility: pos ? 'visible' : 'hidden' }}
      className={`${interactive ? '' : 'pointer-events-none '}fixed z-[60] rounded border border-stone-300 bg-white px-2 py-1.5 text-left text-xs font-normal leading-4 text-stone-800 shadow`}>{children}</div>,
    document.body)
}
