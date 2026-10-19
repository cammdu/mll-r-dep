# Methods

The full write-up, including results and discussion against the prior
co-essentiality literature, is in [manuscript.md](manuscript.md). This file is
the procedural summary.

## 1. Inputs (`src/load_seeds.py`)

Seven per-gene co-dependency exports from the DepMap portal are read and
concatenated into a long table of 877 seed–gene records over 625 distinct
genes. Every input must carry the same release label or the pipeline stops.
Profile lengths differ by seed: DOT1L 174, MLLT1 144, MLLT3 129, MLLT10 128,
MEN1 117, MLLT6 103, KMT2A 82.

## 2. Node selection (`src/build_network.py`)

Pivoting to a gene × seed matrix gives 454 genes with at least one positive
correlation. A gene is admitted if positively correlated with ≥ 2 seeds, or
r ≥ 0.35 with a single seed: 127 nodes. This is a reproducibility filter, not a
significance test — the exports carry no p-values. Its false-discovery
behaviour is characterised in `stats/`.

## 3. Topology

Symbols are mapped with STRING `get_string_ids` (one synonym substitution:
CDIN1 → C15orf41, the same gene) and the induced subnetwork retrieved at
combined score ≥ 0.4, giving 505 edges. **The STRING endpoint is pinned to
v12.0.** The unversioned host round-robins across servers running different
builds; an identical query returned 505 edges from one and 490 from another on
the same day.

## 4. Modules

Weighted label propagation with a fixed seed. At score ≥ 0.4 the chromatin
genes collapse into one 41-gene component whose only shared enrichment term is
the uninformative "nucleoplasm", so that component is re-resolved at ≥ 0.7 and
its largest survivor again at ≥ 0.9. That separates the H3K79 writer arm,
MENIN–MLL–LEDGF, PRC2, MYB–p300 and ncBAF.

Modules are named by marker genes (two or more markers present), with the
STRING enrichment term computed alongside and written to
`stats/module_provenance.csv` so every name can be checked against independent
evidence. Naming was previously keyed on the enrichment string itself, which
proved too unstable to drive a layout.

## 5. Statistics

See [../stats/README.md](../stats/README.md).

## 6. Outputs

Node table, typed edge list, GraphML (written without a graph library, see
`src/export_graph.py`), and the two figures. Both figures run a geometric
text-collision check before saving.

## Deliberate choices worth knowing

- **Negative results are censored, not tested.** Each export is truncated, so
  "absent from all seven lists" means the correlation ranked below roughly
  100th everywhere. Per-seed thresholds are in `stats/censoring_thresholds.csv`.
- **No bias correction.** Removing dominant principal components before
  correlating is standard in genome-scale co-essentiality work, and is not
  possible from rank lists.
- **Module membership is not a claim of physical interaction.** STRING edges
  combine experimental, database and co-expression evidence.
