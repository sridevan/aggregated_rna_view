import { useId, useRef, useState } from 'react'

/** Information icon with a tooltip reachable by mouse (hover) and keyboard (focus / Escape). Rendered with fixed positioning so scroll containers never clip it. */
export function InfoTip({ text, label = 'Information' }: { text: string; label?: string }) {
  const id = useId()
  const btn = useRef<HTMLButtonElement>(null)
  const [pos, setPos] = useState<{ left: number; top: number } | null>(null)
  const open = () => { const r = btn.current?.getBoundingClientRect(); if (r) setPos({ left: Math.min(Math.max(r.left + r.width / 2, 140), window.innerWidth - 140), top: r.bottom + 6 }) }
  const close = () => setPos(null)
  return (
    <span className="inline-block align-middle">
      <button ref={btn} type="button" aria-label={label} aria-describedby={pos ? id : undefined} aria-expanded={!!pos} data-testid="info-tip"
        className="ml-1 inline-flex h-4 w-4 items-center justify-center rounded-full border border-stone-400 text-[10px] font-semibold normal-case leading-none text-stone-600 hover:bg-stone-200 focus:outline-none focus:ring-2 focus:ring-blue-500"
        onMouseEnter={open} onMouseLeave={close} onFocus={open} onBlur={close} onKeyDown={(e) => { if (e.key === 'Escape') close() }} onClick={(e) => { e.stopPropagation(); if (pos) close(); else open() }}>i</button>
      {pos && <span role="tooltip" id={id} style={{ left: pos.left, top: pos.top }} className="fixed z-50 w-64 -translate-x-1/2 whitespace-normal rounded border border-stone-300 bg-white p-2 text-left text-xs font-normal normal-case tracking-normal text-stone-800 shadow">{text}</span>}
    </span>
  )
}
