"""Step 1 - read the per-gene DepMap co-dependency exports into one long table.

Inputs  : data/raw/<SEED>_top100_codependencies_26Q1.csv  (not redistributed; see data/README.md)
Outputs : returns a long DataFrame; nothing written, the merged raw table is
          deliberately not committed because it is a verbatim concatenation of
          the DepMap exports.
"""
import pandas as pd
from config import RAW, SEEDS, DATASET_LABEL


def load_seed(seed):
    """Read one portal export. Accepts the portal's original filename or the renamed form."""
    matches = sorted(RAW.glob(f"{seed}[_']*.csv")) + sorted(RAW.glob(f"{seed}.csv"))
    if not matches:
        raise FileNotFoundError(
            f"No export found for {seed} in {RAW}. "
            f"Download it from the DepMap portal gene page (Perturbation Effects -> "
            f"Co-dependencies -> download) and save as {seed}_top100_codependencies_26Q1.csv"
        )
    df = pd.read_csv(matches[0])
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={"Gene": "gene", "Entrez Id": "entrez",
                            "Correlation": "correlation", "Dataset": "dataset"})
    missing = {"gene", "correlation", "dataset"} - set(df.columns)
    if missing:
        raise ValueError(f"{matches[0].name} is missing columns: {sorted(missing)}")
    df["correlation"] = df["correlation"].astype(float)
    df["seed"] = seed
    return df[["seed", "gene", "entrez", "correlation", "dataset"]]


def load_all():
    long = pd.concat([load_seed(s) for s in SEEDS], ignore_index=True)
    labels = set(long.dataset.unique())
    if labels != {DATASET_LABEL}:
        raise ValueError(
            f"Inputs mix DepMap releases: {sorted(labels)}. "
            f"All seeds must come from '{DATASET_LABEL}' or correlations are not comparable."
        )
    if long.duplicated(["seed", "gene"]).any():
        raise ValueError("Duplicate seed/gene rows in the exports.")
    return long


def correlation_matrix(long):
    """gene x seed matrix of reported correlations; NaN means the gene was not listed."""
    return long.pivot_table(index="gene", columns="seed", values="correlation").reindex(columns=SEEDS)


if __name__ == "__main__":
    long = load_all()
    print(f"{len(long)} seed-gene records over {long.gene.nunique()} distinct genes")
    print(long.groupby("seed").size().reindex(SEEDS).to_string())
