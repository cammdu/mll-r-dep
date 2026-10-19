# Statistical calibration

The DepMap portal exports carry no p-values, no confidence intervals and no
multiple-testing correction, and they are truncated at roughly the
hundredth-ranked partner per seed. The selection step therefore needs its own
null. This directory holds that calibration. It changed a conclusion, so it is
worth reading before using the network.

## The problem

The network admits a gene if it recurs across ≥ 2 of 7 seed lists. How often
would that happen by chance? The answer depends almost entirely on how many
genes are *eligible* to appear in a top-100 correlation list, which is far
fewer than the ~18,000 genes screened: a gene needs appreciable variance in its
gene-effect profile before it can correlate strongly with anything.

## What was done

**Effective universe (capture–recapture).** Treating the seven lists as
repeated captures: 625 distinct genes observed, 489 seen exactly once, 67
exactly twice, giving a bias-corrected Chao estimate of **≈ 2,380 eligible
genes**. Heterogeneous capture probability makes this a lower bound on
richness, so the false-discovery rates below are conservative upper bounds.

**Recurrence null.** Under independent sampling of each list, a gene's
recurrence count is Poisson-binomial. Convolving the exact distribution gives
the expected number of genes reaching each threshold by chance, and a p-value
per gene. These are converted to q-values by Benjamini–Hochberg over the full
estimated universe, with unobserved genes entering as recurrence-0 tests.
Because the p-values are heavily tied, all genes at a given recurrence level
share their tie group's q — a naive cumulative minimum over an arbitrary
ordering would assign different q to identical p.

| recurrence | observed genes | expected by chance | q (BH) |
|---|---|---|---|
| ≥ 2 | 136 | 115.5 | 0.849 |
| ≥ 3 | 69 | 10.2 | 0.147 |
| ≥ 4 | 35 | 0.5 | 0.015 |
| ≥ 5 | 11 | 0.04 | 0.0016 |
| ≥ 6 | 1 | 0.001 | 0.0003 |

**Is the null calibrated?** MLLT3 (AF9) and MLLT6 (AF17) share no appreciable
dependency with DOT1L, so they act as internal negative controls. Every pair
involving them overlapped at or *below* the independence expectation
(observed/expected 0.13–0.86), while all ten pairs among the five related seeds
were enriched 3.4- to 8.7-fold (q < 10⁻⁵). The null is approximately
calibrated, or mildly conservative, for unrelated genes.

**Censoring.** The smallest reported positive correlation per seed is the
truncation threshold: 0.142 (MLLT10) to 0.193 (DOT1L). Absence from a list is
an interval-censored observation, never a tested null.

**Module stability.** Label propagation repeated with 400 random orderings;
mean within-module co-assignment ranges from 0.87 (H3K79 writer, whose members
bridge to adjacent blocks) to 1.00 (SAGA, Mediator, PRC2, MENIN–MLL–LEDGF,
ncBAF).

**Partial correlations.** Each seed-pair correlation recomputed controlling for
the third member, to test whether the MENIN–MLL–DOT1L triad is a genuine
three-way relationship or a chain.

## What this changed

The ≥ 2 threshold carries an 85% false-discovery rate and is **not**
defensible. Every gene is therefore tiered — A (≥ 4), B (= 3), C (everything
else) — and a 66-gene high-confidence subset defined at recurrence ≥ 3.

The effect is strongly non-uniform: chromatin modules retain 32 of 45 members
at the defensible threshold, metabolic modules only 7 of 31. The suspicion that
the mitochondrial, nucleotide and NAD blocks reflect lineage covariance rather
than H3K79 biology is confirmed, rather than merely asserted.

## A note on nominal significance

It does not discriminate here. The smallest correlation in any export is 0.142;
for a panel of 800+ cell lines that is p < 10⁻⁴, and the Bonferroni threshold
for 18,000 tests falls at r ≈ 0.135–0.165 over the same range. The portal's
truncation already sits at roughly the genome-wide significance boundary, so
essentially every listed correlation is nominally significant and the
discrimination must come from effect size and recurrence.

## Files

| File | Contents |
|---|---|
| `censoring_thresholds.csv` | Per-seed truncation bounds |
| `recurrence_null.csv` | Observed vs expected counts, p and q per recurrence level |
| `universe_sensitivity.csv` | How the FDR moves with the assumed universe size |
| `pairwise_overlap.csv` | Hypergeometric test per seed pair, BH corrected |
| `partial_correlations.csv` | Seed-pair correlations controlling for a third seed |
| `module_provenance.csv` | Each module's members and its STRING enrichment term |

## Caveats

Both the discreteness of the p-values and the positive dependence among seeds
make the BH correction conservative rather than exact. And none of this
addresses bias in the correlations themselves — principal-component correction
requires the gene-effect matrix, not rank lists.
