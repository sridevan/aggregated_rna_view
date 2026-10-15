/// <reference types="vite/client" />
/** Vendored PDBe RNA viewer (public/vendor/pdb-rna-viewer, loaded by a <script> tag in index.html). API from src/app/index.ts of the viewer. */
interface PdbRnaViewerPluginInstance {
  render(target: HTMLElement, options: { pdbId: string; entityId: string; chainId: string; subscribeEvents?: boolean; theme?: Record<string, unknown>; apiData?: unknown; FR3DData?: unknown }): Promise<void>
  selectResidue(label_seq_id: number, color?: string, doNotReport?: boolean): void
  selectResidueRange(start: number, end: number, color?: string, doNotReport?: boolean): void
  clearSelection(label_seq_id?: number, doNotReport?: boolean): void
  highlightResidue(label_seq_id: number, color?: string, doNotReport?: boolean): void
  clearHighlight(doNotReport?: boolean): void
  apiData?: { sequence: string; label_seq_ids: (number | null)[]; auth_seq_ids: (number | null)[]; unobserved_label_seq_ids: number[] }
}
interface Window {
  PdbRnaViewerPlugin?: new () => PdbRnaViewerPluginInstance
  __rnaExplorer?: Record<string, unknown>
}
