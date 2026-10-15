# Frontend architecture (Milestone 2)

Stack: React 19, TypeScript, Vite, Tailwind CSS v4 (`@tailwindcss/vite`), Mol* 5 (`createPluginUI` + React 18 renderer, prebuilt `molstar.css`).
Data: static JSON produced by the Python pipeline. Vite's `publicDir` points at `../data`, so `/families/RF00167/*.json` is served in dev and copied into `dist/` on build. No API, no Python, no Infernal at runtime; Mol* fetches coordinates from PDBe (`*_updated.cif`).

```
src/
  main.tsx, App.tsx                 app shell: FamilyProvider + SelectionProvider + three tabs
  types/backend.ts                  TS types mirroring rna_explorer/schemas/*.schema.json
  api/data.ts                       loadFamily / loadInstances / loadStats at start; loadMapping(seq_name) lazily, memoised per chain
  lib/consensus.ts                  pure translation: consensus positions -> structural residues via a chain mapping
                                    (highlightable [has coordinates], unobserved, deleted, uncovered, not-in-mapping)
  state/family.tsx                  FamilyContext: the three start-up JSON files + load status
  state/selection.tsx               SelectionContext: { positions: Set<consensus pos>, instance: seq_name | null }
                                    positions are family coordinates and never depend on the instance (3D and future 2D viewers both consume it)
  components/
    Tabs.tsx                        three tabs
    FamilyOverview.tsx              Tab 1: accession/name/description, counts, methods, CLEN, coverage summary, representative, Mol* viewer
    StructuresTab.tsx               Tab 2: InstanceTable + Mol* viewer + PendingPanel x3 (similarity matrix, clustering/dendrogram, groups)
    InstanceTable.tsx               searchable, sortable table from instances.json; row click selects the instance
    PendingPanel.tsx                explicit "unavailable: not generated yet" state
    AnnotationsTab.tsx              Tab 3: InstanceSelector + ConsensusTracks + SelectionSummary + Mol* viewer (highlight = translation.highlightable)
    ConsensusTracks.tsx             SVG tracks on the consensus axis: axis, sequence coverage, structural coverage, PP mean (display bands), instance state track
                                    (match with coordinates / match without coordinates / deletion / uncovered / insert anchors); click + shift-click selection
    MolstarViewer.tsx               reusable: props { pdbId, labelAsymId, highlight?: { labelAsymId, seqIds } }; loads entry, focuses chain, selects residues by label_asym_id + label_seq_id
```

Selection flow: user clicks consensus column -> SelectionContext.positions -> AnnotationsTab translates through the *selected* chain's mapping (`translateConsensus`) -> MolstarViewer receives label_asym_id + label_seq_id list of residues **with coordinates only**. Switching the instance re-runs the translation; `positions` is untouched. Deleted / uncovered / unobserved positions are listed in the summary and drawn distinctly in the instance track, never sent to Mol*.

Identifiers: Mol* queries use `label_asym_id` and `label_seq_id` (= `polymer_seq_position`). Author chain / number / insertion code are shown as text only.
