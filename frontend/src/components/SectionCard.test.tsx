import { describe, expect, it } from 'vitest'
import { renderToStaticMarkup } from 'react-dom/server'
import { SectionCard } from './SectionCard'

describe('SectionCard', () => {
  it('renders eyebrow, heading level, description, actions and body', () => {
    const html = renderToStaticMarkup(<SectionCard eyebrow="Diversity" title="Species distribution" description="desc" actions={<button>toggle</button>} headingLevel={3} testId="x">body</SectionCard>)
    expect(html).toContain('Diversity'); expect(html).toContain('<h3'); expect(html).toContain('Species distribution'); expect(html).toContain('desc'); expect(html).toContain('toggle'); expect(html).toContain('body')
  })
  it('collapsible sections render an accessible button and hide the body when closed', () => {
    const closed = renderToStaticMarkup(<SectionCard title="Consensus" collapsible defaultOpen={false}>hidden-body</SectionCard>)
    expect(closed).toContain('aria-expanded="false"'); expect(closed).not.toContain('hidden-body')
    const open = renderToStaticMarkup(<SectionCard title="Consensus" collapsible>shown-body</SectionCard>)
    expect(open).toContain('aria-expanded="true"'); expect(open).toContain('shown-body')
  })
})
