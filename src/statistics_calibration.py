"""Step 3 - calibrate the recurrence filter.

The portal exports carry no p-values and are truncated at roughly the hundredth
ranked partner per seed, so the selection step needs its own null. Four analyses:

  censoring_thresholds   what "absent from a list" actually bounds
  effective_universe     capture-recapture estimate of genes eligible to appear
  recurrence_null        Poisson-binomial null -> per-gene BH q-values
  pairwise_overlap       hypergeometric test per seed pair, BH corrected

Outputs: stats/*.csv
"""
import itertools

import numpy as np
import pandas as pd
from scipy import stats

from config import SEEDS


def censoring_thresholds(long):
    """Smallest reported |r| per seed. Absence from a list implies |r| below this."""
    rows = []
    for seed, df in long.groupby("seed"):
        pos = df[df.correlation > 0].correlation
        neg = df[df.correlation < 0].correlation
        rows.append({"seed": seed, "n_listed": len(df), "n_pos": len(pos), "n_neg": len(neg),
                     "min_pos_r": pos.min() if len(pos) else np.nan,
                     "max_neg_r": neg.max() if len(neg) else np.nan,
                     "max_r": pos.max() if len(pos) else np.nan})
    out = pd.DataFrame(rows).set_index("seed").loc[SEEDS]
    out["note"] = "absence from this seed's list implies min_pos_r > r > max_neg_r"
    return out


def effective_universe(long):
    """Chao capture-recapture estimate of the number of genes eligible to be listed.

    The seven lists are treated as repeated captures. This is far smaller than the
    screening library because a gene needs appreciable variance in its gene-effect
    profile before it can correlate strongly with anything. Heterogeneous capture
    probability makes this a lower bound on richness.
    """
    captures = long.groupby("gene").seed.nunique()
    s_obs = len(captures)
    f1 = int((captures == 1).sum())
    f2 = int((captures == 2).sum())
    chao = s_obs + f1 ** 2 / (2 * f2) if f2 else np.nan
    chao_bc = s_obs + f1 * (f1 - 1) / (2 * (f2 + 1))      # bias-corrected
    return {"observed": s_obs, "f1": f1, "f2": f2,
            "chao1": chao, "chao1_bias_corrected": chao_bc,
            "spectrum": captures.value_counts().sort_index().to_dict()}


def poisson_binomial(probs):
    """Exact distribution of a sum of independent non-identical Bernoullis."""
    dist = np.array([1.0])
    for p in probs:
        dist = np.convolve(dist, [1 - p, p])
    return dist


def recurrence_null(long, universe_n):
    """Per-gene q-value for recurrence under independent sampling of each list.

    Returns (per_gene_q: dict {recurrence: q}, table of expected vs observed).
    BH is applied over the full estimated universe, so genes never listed enter
    as recurrence-0 tests. p-values are heavily tied, therefore all genes at a
    given recurrence level share their tie group's q-value - the BH step-up is
    applied to the sorted vector and the per-level minimum taken, which is the
    correct treatment of ties (a naive cumulative minimum over an arbitrary
    ordering assigns different q to identical p).
    """
    sizes = long.groupby("seed").size().reindex(SEEDS).values
    captures = long.groupby("gene").seed.nunique()
    m = int(round(universe_n))
    dist = poisson_binomial(sizes / universe_n)
    p_by_r = {k: float(dist[k:].sum()) for k in range(len(dist))}

    r_vec = np.concatenate([captures.values, np.zeros(max(m - len(captures), 0), dtype=int)])
    p_vec = np.array([p_by_r.get(int(k), 0.0) for k in r_vec])
    q_vec = stats.false_discovery_control(p_vec, method="bh")
    q_by_r = pd.Series(q_vec).groupby(pd.Series(r_vec)).min().to_dict()

    rows = []
    for k in range(1, 7):
        rows.append({"recurrence": k,
                     "observed_genes": int((captures >= k).sum()),
                     "expected_by_chance": float(m * dist[k:].sum()),
                     "p_value": p_by_r.get(k, 0.0),
                     "q_BH": float(q_by_r.get(k, 1.0))})
    return q_by_r, pd.DataFrame(rows)


def universe_sensitivity(long, universe_sizes):
    """How the filter's false-discovery rate moves with the assumed universe size."""
    sizes = long.groupby("seed").size().reindex(SEEDS).values
    captures = long.groupby("gene").seed.nunique()
    rows = []
    for n in universe_sizes:
        dist = poisson_binomial(sizes / n)
        row = {"universe_N": int(n)}
        for k in (2, 3, 4, 5):
            obs = int((captures >= k).sum())
            exp = n * dist[k:].sum()
            row[f"expected_ge{k}"] = round(float(exp), 2)
            row[f"observed_ge{k}"] = obs
            row[f"fdr_ge{k}"] = round(float(min(1.0, exp / max(obs, 1))), 4)
        rows.append(row)
    return pd.DataFrame(rows)


def pairwise_overlap(long, universe_n, matrix=None):
    """Hypergeometric test of list overlap for every seed pair, BH corrected.

    Seeds with no real relationship to the axis (MLLT3, MLLT6) act as internal
    negative controls: if the null is calibrated their overlaps sit at or below
    the independence expectation.
    """
    sets = {s: set(df.gene) for s, df in long.groupby("seed")}
    rows = []
    for a, b in itertools.combinations(SEEDS, 2):
        obs = len(sets[a] & sets[b])
        exp = len(sets[a]) * len(sets[b]) / universe_n
        p = stats.hypergeom.sf(obs - 1, int(universe_n), len(sets[a]), len(sets[b]))
        seed_r = np.nan
        if matrix is not None and a in matrix.index and b in matrix.columns:
            seed_r = matrix.loc[a, b]
        rows.append({"pair": f"{a}-{b}", "observed": obs, "expected": round(exp, 2),
                     "ratio": round(obs / exp, 2), "p_hypergeom": p, "seed_seed_r": seed_r})
    out = pd.DataFrame(rows).sort_values("ratio", ascending=False)
    out["q_BH"] = stats.false_discovery_control(out.p_hypergeom.values, method="bh")
    return out


def partial_correlation(r_ab, r_ac, r_bc):
    """Correlation of a and b with c partialled out."""
    return (r_ab - r_ac * r_bc) / np.sqrt((1 - r_ac ** 2) * (1 - r_bc ** 2))


def module_stability(nodes, edges, modules, runs, label_propagation_fn):
    """Mean within-module co-assignment frequency over `runs` label-propagation restarts."""
    index = {n: i for i, n in enumerate(nodes)}
    co = np.zeros((len(nodes), len(nodes)))
    for seed in range(runs):
        labels = label_propagation_fn(nodes, edges, seed=seed)
        vec = np.array([labels[n] for n in nodes])
        co += (vec[:, None] == vec[None, :])
    co /= runs
    per_gene, per_module = {}, {}
    for module in set(modules.values()):
        members = [g for g, m in modules.items() if m == module and g in index]
        if len(members) < 2:
            continue
        ii = [index[g] for g in members]
        block = co[np.ix_(ii, ii)]
        per_module[module] = float(block[np.triu_indices(len(ii), 1)].mean())
        for g in members:
            peers = [index[p] for p in members if p != g]
            per_gene[g] = float(np.mean([co[index[g], p] for p in peers]))
    return per_gene, per_module
