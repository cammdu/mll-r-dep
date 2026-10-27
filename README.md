# DOT1L–H3K79 co-dependency network

A map of the genetic-dependency neighbourhood around DOT1L in cancer cell lines,
built from DepMap CRISPR co-dependency profiles and STRING protein-association
topology, and the selection step is calibrated against an explicit null.

![Co-dependency network](figures/dot1l_h3k79_dependency_network.png)

## The question

DOT1L is the only enzyme that methylates histone H3 at lysine 79, and
KMT2A-rearranged (MLL-r) leukaemia depends on that activity. Catalytic DOT1L
inhibitors nonetheless produced only modest single-agent clinical activity,
which raises a practical question: **what else does a cell need in order to
need DOT1L?** Genes whose CRISPR knockout fitness profiles track DOT1L's
across many cell lines are candidates for acting in the same process, and
therefore candidates for combination targeting.

Concretely I ask three things:

1. Which genes recur across the co-dependency neighbourhoods of DOT1L and of
   six genes tied to it physically or functionally?
2. Do those genes form coherent complexes, or an unstructured list?
3. How much of the resulting structure is real, and how much is chance?

## How

Seven per-gene co-dependency profiles were exported from the DepMap portal
(DOT1L, MLLT10/AF10, MLLT6/AF17, MLLT1/ENL, MLLT3/AF9, KMT2A, MEN1), giving
877 seed–gene correlation records over 625 distinct genes. A gene was admitted
to the network if it was positively correlated with at least two of the seven
seeds, or reached r ≥ 0.35 with one, yielding 127 nodes. Protein-association
edges among those nodes came from STRING v12.0 at combined score ≥ 0.4 (505
edges), and modules were found by weighted label propagation, with the densely
connected chromatin component re-resolved at stricter confidence. Module names
are anchored to marker genes and cross-checked against STRING functional
enrichment.

Because the portal exports carry no p-values and are truncated at roughly the
hundredth-ranked partner, the selection step was calibrated separately — see
[stats/README.md](stats/README.md). 

## Main results

**The MENIN–MLL axis, not the rest of the writer complex, is closest to DOT1L.**
The strongest relationships in the whole network are KMT2A–MEN1 (r = 0.58),
DOT1L–MEN1 (0.57) and DOT1L–KMT2A (0.52), with LEDGF/PSIP1 at 0.53. Partial
correlations show this is a genuine three-way relationship rather than a chain
through any one member: controlling for KMT2A, DOT1L–MEN1 only falls to 0.38.
The cells that need DOT1L are largely the cells that need MENIN, which is an
argument for combination and a warning about shared resistance.

**Only two of four DotCom paralogues travel with DOT1L.** AF10 (0.45) and ENL
(0.38) do; AF17 reaches just 0.17 and AF9 returns no correlation with DOT1L at
all. Neither AF17 nor AF9 contributed a node of its own. This is fitness-level
redundancy, not evidence against a physical interaction.

**KAT6A is the most reproducible partner in the analysis**, appearing in six of
seven seed profiles — the only gene at that level (q = 3.0 × 10⁻⁴) — inside an
11-gene MOZ/MORF–HBO1 acetyltransferase module. Its correlations with KMT2A and
MLLT1 exceed what transitivity through DOT1L would produce, so the association
is not an artefact of DOT1L's own links. **Important counterpoint:** a direct
test found the KAT6A/B inhibitor WM-1119 effective in KAT6A-rearranged but
*not* KMT2A-rearranged AML, so a shared dependency profile should not be read
as catalytic tractability in MLL-r disease.

**A biochemically required cofactor arm is missing.** The H2B-ubiquitination
ligases RNF20 and RNF40, which stimulate DOT1L, appear in none of the seven
profiles, while the SAGA H2B-deubiquitinase module appears in full (12 of 12
members). Required for the chemistry, apparently not for survival — though see
the censoring caveat below before reading too much into any absence.

**Most of the permissive network is chance, and the calibration says which
part.** At the recurrence ≥ 2 threshold we started with, the false-discovery
rate is 85%. At ≥ 3 it is 15%, at ≥ 4 it is 1.5%. Applying the defensible
threshold, the chromatin modules retain 71% of their members and the metabolic
modules 23% — confirming on statistical rather than rhetorical grounds that the
mitochondrial, nucleotide and NAD blocks reflect lineage and proliferation
covariance among haematopoietic lines rather than H3K79 biology.

![Statistical calibration](figures/dot1l_statistical_calibration.png)

## Main limitations

- **These are pan-cancer correlations.** The therapeutically decisive contrast,
  co-dependency structure *within* KMT2A-rearranged lines versus everything
  else, is not modelled here. This is the single biggest caveat.
- **Correlations are uncorrected for systematic covariance.** Best practice
  removes dominant principal components from the gene-effect matrix before
  correlating. That requires the matrix, not rank lists, so any inflation
  shared across the axis remains in every value reported.
- **Absence is censoring, not a null result.** Each export is truncated at
  roughly the hundredth partner. A gene missing from a list only means its
  correlation fell below that seed's threshold (0.142–0.193 depending on seed;
  see `stats/censoring_thresholds.csv`). No absence here is a tested negative.
- **Module names inherit STRING's literature bias.** 42 of 127 nodes fall
  outside any named module, disproportionately poorly characterised genes.
- **The gene universe is estimated, not known.** The capture–recapture figure
  of ~2,380 assumes homogeneous capture probability, which is false; it is a
  lower bound, making the reported false-discovery rates conservative upper
  bounds. Two internally unrelated seeds (AF17, AF9) were used to check the
  null is not anti-conservative.

## Next steps

1. **Restrict to KMT2A-rearranged lines.** With `CRISPRGeneEffect.csv`,
   `Model.csv` and `OmicsFusionFiltered.csv` from the same release, compute
   MLL-r versus non-MLL-r differential dependency with effect sizes and formal
   FDR, and recompute co-dependency within KMT2A-r lines only. This separates
   genuine module structure from the lineage signal that dominates the
   metabolic blocks.
2. **Apply principal-component bias correction** once the matrix is available,
   and re-derive every correlation reported here.
3. **Test KAT6A by degradation rather than catalytic inhibition**, since the
   co-dependency may report a requirement for complex integrity that an
   active-site inhibitor does not reproduce.
4. **Add an expression layer** to confirm the module genes are expressed in the
   dependent lines and to relate co-dependency to HOXA9/MEIS1 activation.

## Repository layout

```
data/      node table, edge list and GraphML network (main outputs)
           raw/ is where you place the DepMap exports; it is not committed
methods/   the manuscript, and how each step was done
stats/     calibration outputs: censoring, null, overlap, partial correlations
figures/   the two manuscript figures
src/       the pipeline that produces everything above
```

## Reproducing

```bash
pip install -r requirements.txt
# place the seven DepMap exports in data/raw/ - see data/README.md
cd src && python run_pipeline.py
```

Runs in about a minute. Reaches the STRING API; everything else is local. The
STRING endpoint is pinned to v12.0 deliberately: the unversioned host
round-robins across servers running different builds, and returned 505 edges
from one and 490 from another on the same day for an identical query.

## Data source and credit

Dependency data are from the **Cancer Dependency Map (DepMap)**, DepMap Public
26Q1+Score, Chronos gene effect, generated by the Broad Institute and the
Wellcome Sanger Institute. DepMap data are released under CC BY 4.0. This
repository does not redistribute the portal export files; see
[data/README.md](data/README.md) for attribution and how to obtain them.

Protein-association data are from **STRING v12.0**, also CC BY 4.0.

If you use this work, please cite the DepMap and STRING sources listed in
`data/README.md` alongside this repository.

## Licence

Code is MIT licensed (see `LICENSE`). Derived data files inherit CC BY 4.0 from
DepMap and STRING; attribution requirements are in `data/README.md`.
