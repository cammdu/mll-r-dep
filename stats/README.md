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


---

# MLL-r-restricted differential dependency

A separate analysis, added after the calibration above, using the full Chronos
gene-effect matrix rather than the portal's rank lists. It tests the premise the
network rests on: that DOT1L is *selectively* required in KMT2A-rearranged
disease.

## Design

Lines are classified by two independent sources that must agree — a canonical
leukaemogenic KMT2A fusion partner, and DepMap's own subtype annotation. They
agreed on all 13 calls. Fourteen lines with other KMT2A evidence (passenger
fusions in solid tumours; the partial-tandem-duplication lines EOL-1 and
Kasumi-6) are marked **ambiguous and excluded**, not pooled into either arm.
Sister lines from one patient (MOLM-13/MOLM-14, SEM/SEMK2) are averaged, giving
11 independent KMT2A-r patients against 38 other acute leukaemias. Each gene is
fitted by OLS with lineage indicators (AML reference, B-ALL, T-ALL), so a gene
cannot score by being a myeloid or lymphoid dependency; p-values are BH
corrected across 17,874 genes.

## Results

| gene | β (KMT2A-r) | 95% CI | q | rank |
|---|---|---|---|---|
| ZFP64 | −0.41 | — | 0.002 | 1 |
| FLT3 | −0.66 | — | 0.005 | 2 |
| JMJD1C | −0.44 | — | 0.039 | 4 |
| MEN1 | −0.72 | — | 0.133 | 6 |
| KMT2A | −0.49 | — | 0.250 | 13 |
| **DOT1L** | **−0.18** | **[−0.44, +0.08]** | **0.947** | **3,532** |

Permutation tests (10,000 label shuffles) confirm the parametric p-values:
DOT1L p = 0.19; MEN1, KMT2A, ZFP64 and JMJD1C all p < 10⁻⁴.

## Reading the DOT1L null correctly

DOT1L *is* essential in KMT2A-r lines (mean effect −0.55, p = 4 × 10⁻⁴ against
zero). What it is not is *selectively* essential relative to other acute
leukaemias. And the minimum detectable effect at conventional power is ~0.375
Chronos units, so the design resolves a MENIN-sized window (−0.72) but not a
modest one. The interval excludes effects more negative than −0.44 while
remaining compatible with zero.

**This bounds DOT1L selectivity; it does not refute it.** Eleven patients is the
entire public KMT2A-r panel, so no reframing of the analysis escapes that cap.

## Known confounds

- **FLT3** ranks second, but 2 of 11 KMT2A-r patients carry FLT3-ITD and no
  comparator does. Those two drive the group mean (−1.54 versus −0.59 in the
  ITD-negative KMT2A-r patients). Treat the FLT3 rank as partly genotype
  confounding, not pure MLL-r biology.
- **ZFP64 and JMJD1C** are both established MLL-r dependencies, which is
  reassurance that the design detects real signal rather than noise.
- Correlations and effects remain uncorrected for dominant principal
  components, as in the pan-cancer analysis.

## Files

| File | Contents |
|---|---|
| `mllr_differential_dependency.csv` | All 17,874 genes: β, SE, p, q, rank, group means, Hedges' g |
| `mllr_module_selectivity.csv` | Mean β per module with a one-sample test against zero |
| `mllr_within_group_seed_correlations.csv` | Seed-seed co-dependency computed inside KMT2A-r only (n = 11; indicative) |
| `../data/mllr_cellline_classification.csv` | Every screened line's class and fusion partner |
