"""MLL-r-restricted differential dependency.

Tests whether each gene is MORE essential in KMT2A-rearranged lines than in
other acute leukaemias, using the full Chronos gene-effect matrix. This is the
analysis the per-gene co-dependency exports cannot do, and the one that
separates MLL-r-specific biology from pan-cancer lineage covariance.

Inputs (place in data/raw/, not redistributed - see data/README.md):
    CRISPRGeneEffect.csv          Chronos gene effect, lines x genes
    Model.csv                     cell-line annotation
    OmicsFusionFiltered.csv       filtered fusion calls

Outputs:
    stats/mllr_differential_dependency.csv
    stats/mllr_module_selectivity.csv
    stats/mllr_within_group_seed_correlations.csv
    data/mllr_cellline_classification.csv
    figures/mllr_selectivity.png

    cd src && python mllr_selectivity.py
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats

from config import RAW, DATA, STATS, FIGS, SEEDS, RANDOM_SEED

# Fusion partners that define a canonical leukaemogenic KMT2A rearrangement.
# A KMT2A fusion to an arbitrary passenger gene in a solid tumour is not MLL-r.
CANONICAL_PARTNERS = {"AFF1", "MLLT3", "MLLT1", "MLLT10", "AFDN", "ELL", "EPS15",
                      "MLLT6", "SEPTIN6", "SEPTIN9", "MLLT11", "FOXO4", "CLTC",
                      "AFF3", "AFF4"}
ACUTE_LEUKAEMIA = {"Acute Myeloid Leukemia", "B-Cell Acute Lymphoblastic Leukemia",
                   "T-Lymphoblastic Leukemia/Lymphoma", "Acute Leukemias of Ambiguous Lineage"}
LINEAGE = {"Acute Myeloid Leukemia": "AML", "B-Cell Acute Lymphoblastic Leukemia": "BALL",
           "T-Lymphoblastic Leukemia/Lymphoma": "TALL"}
N_PERMUTATIONS = 10000

# Modules carried over from the co-dependency network, for module-level testing.
MODULES = {
    "MENIN-MLL-LEDGF": ["KMT2A", "MEN1", "PSIP1"],
    "H3K79 writer": ["DOT1L", "MLLT10", "MLLT1", "AFF1"],
    "MOZ/MORF-HBO1": ["KAT6A", "KAT6B", "BRPF1", "KAT7", "JADE1", "JADE2", "ING5",
                      "MEAF6", "MORF4L1", "MRGBP", "MSL3"],
    "SAGA": ["SGF29", "TADA1", "TADA2B", "TADA3", "TAF5L", "TAF6L", "TAF10",
             "SUPT7L", "SUPT20H", "USP22", "ATXN7", "ATXN7L3"],
    "PRC2": ["EZH2", "EED", "SUZ12"],
    "Mediator": ["MED12", "MED13", "MED13L", "MED15", "MED16", "MED19", "MED23",
                 "MED24", "CCNC"],
    "Fusion output": ["HOXA9", "HOXA10", "MEIS1", "PBX3", "MEF2C", "JMJD1C", "ZFP64"],
    "H2Bub ligase": ["RNF20", "RNF40", "UBE2A", "UBE2B"],
    "OXPHOS": ["SDHA", "SDHB", "UQCRC1", "UQCRC2", "COX5B", "COX6B1", "CYC1", "SURF1"],
}


def load_inputs():
    need = {"CRISPRGeneEffect.csv", "Model.csv", "OmicsFusionFiltered.csv"}
    missing = [f for f in need if not (RAW / f).exists()]
    if missing:
        raise FileNotFoundError(f"Missing from {RAW}: {missing}. See data/README.md")
    ge = pd.read_csv(RAW / "CRISPRGeneEffect.csv", index_col=0)
    ge.columns = [c.split(" (")[0] for c in ge.columns]
    return ge, pd.read_csv(RAW / "Model.csv"), pd.read_csv(RAW / "OmicsFusionFiltered.csv")


def classify_lines(ge, model, fusion):
    """Label each screened line KMT2A_r / non_r / ambiguous.

    Two independent sources must agree: a canonical fusion call, and DepMap's
    own subtype annotation. Lines with any other KMT2A evidence (passenger
    fusions, partial tandem duplications) are marked ambiguous and EXCLUDED
    rather than pooled into either group.
    """
    m = model.set_index("ModelID")
    f = fusion[fusion.IsDefaultEntryForModel == "Yes"].copy()
    f["g1"] = f.Gene1.str.split(" ").str[0]
    f["g2"] = f.Gene2.str.split(" ").str[0]
    k = f[(f.g1 == "KMT2A") | (f.g2 == "KMT2A")].copy()
    k["partner"] = np.where(k.g1 == "KMT2A", k.g2, k.g1)
    k = k[k.partner != "KMT2A"]                       # drop PTD-like self-fusions
    partners = k.groupby("ModelID").partner.apply(lambda s: ",".join(sorted(set(s))))

    by_fusion = set(k[k.partner.isin(CANONICAL_PARTNERS)].ModelID)
    annot_cols = ["ModelSubtypeFeatures", "PatientSubtypeFeatures"]
    by_annot = set(m.index[m[annot_cols].astype(str).apply(
        lambda s: s.str.contains(r"KMT2A(?:::|-)(?!PTD)[A-Z]", regex=True)).any(axis=1)])

    ids = ge.index
    rearranged = (by_fusion | by_annot) & set(ids)
    other_kmt2a = ((set(k.ModelID) | by_annot) & set(ids)) - rearranged
    cls = pd.Series("non_r", index=ids)
    cls[list(rearranged)] = "KMT2A_r"
    cls[list(other_kmt2a)] = "ambiguous"

    info = m.loc[ids, ["CellLineName", "PatientID", "OncotreeLineage", "OncotreePrimaryDisease"]].copy()
    info["cls"] = cls
    info["acute_leukaemia"] = info.OncotreePrimaryDisease.isin(ACUTE_LEUKAEMIA)
    info["KMT2A_partner"] = info.index.map(partners)
    disagree = (by_fusion ^ by_annot) & set(ids)
    return info, sorted(disagree)


def collapse_to_patients(ge, info):
    """Average sister lines from the same patient - they are not independent."""
    keep = info.index[info.cls != "ambiguous"]
    sub = info.loc[keep]
    if (sub.groupby("PatientID").cls.nunique() > 1).any():
        raise ValueError("A patient has lines with conflicting KMT2A status.")
    profiles = ge.loc[keep].groupby(sub.PatientID).mean()
    meta = (sub.groupby("PatientID")
               .agg(cls=("cls", "first"), disease=("OncotreePrimaryDisease", "first"),
                    acute=("acute_leukaemia", "first"), n_lines=("CellLineName", "size"),
                    lines=("CellLineName", lambda s: "/".join(s))))
    meta["lineage"] = meta.disease.map(LINEAGE).fillna("other")
    return profiles, meta


def _design(meta):
    """Intercept + KMT2A-r indicator + lineage dummies (AML as reference)."""
    dummies = pd.get_dummies(meta.lineage, drop_first=False).astype(float)
    dummies = dummies.drop(columns=[c for c in ["AML"] if c in dummies])
    x = (meta.cls == "KMT2A_r").astype(float).values
    return np.column_stack([np.ones(len(meta)), x, dummies.values]), list(dummies.columns)


def differential_dependency(profiles, meta, min_group=5):
    """OLS per gene: effect of KMT2A-r status on gene effect, adjusted for lineage."""
    acute = meta[meta.acute]
    Y = profiles.loc[acute.index]
    X, lineage_cols = _design(acute)
    is_r = (acute.cls == "KMT2A_r").values

    def fit(y):
        ok = ~np.isnan(y)
        Xo, yo = X[ok], y[ok]
        if ok.sum() < X.shape[1] + min_group or Xo[:, 1].sum() < min_group:
            return (np.nan,) * 4
        beta, *_ = np.linalg.lstsq(Xo, yo, rcond=None)
        resid = yo - Xo @ beta
        dof = len(yo) - Xo.shape[1]
        cov = (resid @ resid / dof) * np.linalg.pinv(Xo.T @ Xo)
        se = np.sqrt(cov[1, 1])
        return beta[1], se, 2 * stats.t.sf(abs(beta[1] / se), dof), dof

    res = pd.DataFrame([fit(Y[g].values) for g in Y.columns], index=Y.columns,
                       columns=["beta_KMT2Ar", "se", "p", "dof"])
    res["mean_KMT2Ar"] = Y[is_r].mean()
    res["mean_nonr_AL"] = Y[~is_r].mean()
    pooled = np.sqrt((Y[is_r].var() * (is_r.sum() - 1) +
                      Y[~is_r].var() * ((~is_r).sum() - 1)) / (len(is_r) - 2))
    res["g_raw"] = (res.mean_KMT2Ar - res.mean_nonr_AL) / pooled
    res = res.dropna(subset=["p"])
    res["q"] = stats.false_discovery_control(res.p.values, method="bh")
    res = res.sort_values("p")
    res["rank"] = np.arange(1, len(res) + 1)
    res.index.name = "gene"
    return res, X, Y, lineage_cols


def permutation_test(Y, X, gene, n=N_PERMUTATIONS, seed=RANDOM_SEED):
    """Exact-ish null for one gene by permuting the KMT2A-r label."""
    y = Y[gene].values
    ok = ~np.isnan(y)
    Xo, yo = X[ok], y[ok]
    observed = np.linalg.lstsq(Xo, yo, rcond=None)[0][1]
    rng = np.random.default_rng(seed)
    labels = Xo[:, 1].copy()
    draws = np.empty(n)
    for i in range(n):
        Xp = Xo.copy()
        Xp[:, 1] = rng.permutation(labels)
        draws[i] = np.linalg.lstsq(Xp, yo, rcond=None)[0][1]
    return observed, float((np.abs(draws) >= abs(observed)).mean())


def minimum_detectable_effect(se, power_z=2.8):
    """Effect size detectable at ~80% power, nominal two-sided alpha=0.05."""
    return power_z * se


def module_selectivity(res):
    rows = []
    for name, genes in MODULES.items():
        betas = res.reindex(genes).beta_KMT2Ar.dropna()
        rows.append({"module": name, "n": len(betas), "mean_beta": betas.mean(),
                     "median_beta": betas.median(),
                     "p": stats.ttest_1samp(betas, 0).pvalue,
                     "n_q_lt_0.25": int((res.reindex(genes).q < 0.25).sum())})
    return pd.DataFrame(rows).sort_values("mean_beta")


def within_group_correlations(profiles, meta, min_sd=0.05):
    """Co-dependency computed INSIDE the KMT2A-r group only."""
    Y = profiles.loc[meta.index[meta.cls == "KMT2A_r"]]
    Y = Y.loc[:, (Y.std() > min_sd) & Y.notna().all()]
    present = [s for s in SEEDS if s in Y.columns]
    return pd.DataFrame({s: Y.corrwith(Y[s]) for s in present}).reindex(present), len(Y)
