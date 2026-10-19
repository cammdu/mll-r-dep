"""Run the whole analysis end to end.

    cd src && python run_pipeline.py

Requires the seven DepMap co-dependency exports in data/raw/ (see data/README.md).
Reaches the STRING API; everything else is local.
"""
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

from config import (SEEDS, DATA, STATS, FIGS, STRING_INTER_MODULE, CONSENSUS_RUNS,
                    HIGH_CONFIDENCE_RECURRENCE, RANDOM_SEED)
import load_seeds
import build_network as bn
import statistics_calibration as sc
import export_graph
import figures

# Modules are detected and named automatically, but the manuscript uses short
# display names and groups them for layout. Mapping from the STRING enrichment
# term to the name used in the paper; unmatched modules keep their STRING term.
SUPERGROUP = {
    "H3K79 writer (DotCom)": "axis", "MENIN-MLL-LEDGF": "axis",
    "PRC2": "chrom", "MOZ/MORF-HBO1 HAT": "chrom", "SAGA complex": "chrom",
    "Mediator complex": "chrom", "ncBAF": "chrom",
    "Lymphoid/myeloid TFs": "tf", "MYB-p300 coactivator": "tf",
    "Mitochondrial OXPHOS": "met", "Nucleotide biosynthesis": "met",
    "NAD synthesis": "met", "Heme / sugar metabolism": "met",
}
SHORT = {"H3K79 writer (DotCom)": "H3K79 writer", "MOZ/MORF-HBO1 HAT": "MOZ/MORF-HBO1", "MYB-p300 coactivator": "MYB-p300",
         "Lymphoid/myeloid TFs": "Haematopoietic TFs", "Nucleotide biosynthesis": "Nucleotide synthesis",
         "Heme / sugar metabolism": "Heme/sugar metab."}
CHROMATIN = ["SAGA complex", "MOZ/MORF-HBO1 HAT", "Mediator complex",
             "MENIN-MLL-LEDGF", "PRC2", "H3K79 writer (DotCom)"]
METABOLIC = ["Mitochondrial OXPHOS", "Nucleotide biosynthesis",
             "NAD synthesis", "Heme / sugar metabolism"]
LAYOUT = [(0, ["H3K79 writer (DotCom)", "MENIN-MLL-LEDGF", "PRC2", "ncBAF"]),
          (1, ["MOZ/MORF-HBO1 HAT", "MYB-p300 coactivator", "Lymphoid/myeloid TFs"]),
          (2, ["SAGA complex", "Nucleotide biosynthesis", "NAD synthesis"]),
          (3, ["Mediator complex", "Heme / sugar metabolism"]),
          (4, ["Mitochondrial OXPHOS"])]


def main():
    for d in (DATA, STATS, FIGS):
        d.mkdir(parents=True, exist_ok=True)

    # ---- 1. inputs
    long = load_seeds.load_all()
    matrix = load_seeds.correlation_matrix(long)
    print(f"[1] {len(long)} records, {long.gene.nunique()} distinct genes")

    # ---- 2. node selection and topology
    nodes = bn.select_nodes(matrix)
    kept = nodes[nodes.in_network]
    print(f"[2] {len(kept)} nodes admitted of {len(nodes)} candidates")
    idmap, edges = bn.string_graph(list(kept.index))
    renamed = idmap[idmap.queryItem != idmap.preferredName]
    if len(renamed):
        print(f"    STRING resolved {len(renamed)} symbol(s) to a synonym: "
              + ", ".join(f"{a}->{b}" for a, b in zip(renamed.queryItem, renamed.preferredName)))
    print(f"    {len(edges)} STRING edges at score >= 0.4")

    labels = bn.hierarchical_modules(idmap, edges)
    groups = {}
    for gene, lab in labels.items():
        groups.setdefault(lab, []).append(gene)
    module_of, provenance = {}, []
    for lab, members in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        name = bn.name_module(members)
        if name is None:
            continue
        for gene in members:
            module_of[gene] = name
        term = bn.enrich(members) if len(members) >= 4 else {"term": None, "fdr": float("nan"), "n": 0}
        provenance.append({"module": name, "n_genes": len(members),
                           "string_enrichment_term": term["term"], "enrichment_fdr": term["fdr"],
                           "genes_in_term": term["n"],
                           "members": " ".join(sorted(members))})
        time.sleep(0.6)
    pd.DataFrame(provenance).to_csv(STATS / "module_provenance.csv", index=False)
    nodes["module"] = pd.Series(module_of).reindex(nodes.index)
    print(f"    {len(provenance)} named modules covering {len(module_of)} genes")

    # ---- 3. statistics
    cens = sc.censoring_thresholds(long)
    cens.to_csv(STATS / "censoring_thresholds.csv")
    uni = sc.effective_universe(long)
    universe_n = uni["chao1_bias_corrected"]
    print(f"[3] effective universe (Chao, bias-corrected): {universe_n:.0f} genes")

    q_by_r, null_tbl = sc.recurrence_null(long, universe_n)
    null_tbl.to_csv(STATS / "recurrence_null.csv", index=False)
    sens = sc.universe_sensitivity(long, [universe_n, 4000, 6000, 10000, 14000, 18000])
    sens.to_csv(STATS / "universe_sensitivity.csv", index=False)
    overlap = sc.pairwise_overlap(long, universe_n, matrix)
    overlap.to_csv(STATS / "pairwise_overlap.csv", index=False)
    for k in (2, 3, 4):
        print(f"    recurrence >= {k}: q = {q_by_r.get(k, float('nan')):.4f}")

    string_nodes = list(idmap.preferredName)
    per_gene_stab, per_module_stab = sc.module_stability(
        string_nodes, edges, {g: m for g, m in module_of.items() if g in string_nodes},
        CONSENSUS_RUNS, bn.label_propagation)
    print(f"    module stability over {CONSENSUS_RUNS} restarts: "
          f"min {min(per_module_stab.values()):.2f}, max {max(per_module_stab.values()):.2f}")

    # partial correlations within the seed set
    triads = [("DOT1L", "MEN1", "KMT2A"), ("DOT1L", "KMT2A", "MEN1"), ("KMT2A", "MEN1", "DOT1L"),
              ("DOT1L", "MLLT10", "MEN1"), ("DOT1L", "MLLT10", "KMT2A"), ("DOT1L", "MLLT1", "KMT2A")]
    rows = []
    for a, b, c in triads:
        r_ab, r_ac, r_bc = matrix.loc[a, b], matrix.loc[a, c], matrix.loc[b, c]
        if pd.notna([r_ab, r_ac, r_bc]).all():
            rows.append({"pair": f"{a}-{b}", "controlling_for": c, "r": round(float(r_ab), 4),
                         "partial_r": round(float(sc.partial_correlation(r_ab, r_ac, r_bc)), 4)})
    pd.DataFrame(rows).to_csv(STATS / "partial_correlations.csv", index=False)

    # ---- 4. node table
    degree = pd.Series(0, index=nodes.index, dtype=int)
    for a, b in zip(edges.A, edges.B):
        for g in (a, b):
            if g in degree.index:
                degree[g] += 1
    nodes["string_degree"] = degree
    nodes["recurrence_q_BH"] = nodes.recurrence.map(lambda k: round(float(q_by_r.get(int(k), 1.0)), 6))
    nodes["tier"] = np.where(nodes.recurrence >= 4, "A", np.where(nodes.recurrence >= 3, "B", "C"))
    nodes["high_confidence"] = nodes.recurrence >= HIGH_CONFIDENCE_RECURRENCE
    nodes["module_stability"] = pd.Series(per_gene_stab).reindex(nodes.index).round(3)
    nodes["is_seed"] = nodes.index.isin(SEEDS)
    nodes = nodes.join(matrix.add_prefix("r_"))
    nodes.index.name = "gene"
    front = ["tier", "high_confidence", "recurrence", "recurrence_q_BH", "max_r", "mean_r",
             "module", "module_stability", "in_network", "is_seed", "string_degree"]
    nodes = nodes[front + [c for c in nodes.columns if c not in front]]
    nodes.sort_values(["tier", "max_r"], ascending=[True, False]).round(4).to_csv(
        DATA / "dot1l_h3k79_node_table.csv")

    # ---- 5. edges and graph
    edge_rows = [(a, b, "string_association", round(float(s), 3))
                 for a, b, s in zip(edges.A, edges.B, edges.score)]
    for seed in SEEDS:
        for gene, r in matrix[seed].dropna().items():
            if gene in kept.index and gene != seed:
                edge_rows.append((seed, gene, "codependency", round(float(r), 4)))
    edge_df = pd.DataFrame(edge_rows, columns=["source", "target", "edge_type", "weight"])
    edge_df.to_csv(DATA / "dot1l_h3k79_network_edges.csv", index=False)
    net = nodes[nodes.in_network]
    n_nodes, n_edges = export_graph.write_graphml(
        DATA / "dot1l_h3k79_network.graphml", net, edge_df)
    print(f"[5] graphml: {n_nodes} nodes, {n_edges} edges")

    # ---- 6. figures
    drawn = {m for _, ms in LAYOUT for m in ms}
    shown = net[net.module.isin(drawn)]
    inter = [(a, b, s) for a, b, s in zip(edges.A, edges.B, edges.score)
             if s >= STRING_INTER_MODULE and a in shown.index and b in shown.index
             and module_of.get(a) != module_of.get(b)]
    fig1 = figures.network_figure(
        shown, inter, LAYOUT, SUPERGROUP,
        {m: SHORT.get(m, m) for m in drawn},
        FIGS / "dot1l_h3k79_dependency_network.png",
        n_omitted=int(len(net) - len(shown)))
    counts = {k: int((nodes.recurrence >= k).sum()) for k in (2, 3, 4, 5)}
    before = {m: int((nodes.module == m).sum()) for m in CHROMATIN + METABOLIC}
    after = {m: int(((nodes.module == m) & (nodes.recurrence >= HIGH_CONFIDENCE_RECURRENCE)).sum())
             for m in CHROMATIN + METABOLIC}
    fig2 = figures.calibration_figure(
        counts, long.groupby("seed").size().reindex(SEEDS).values,
        [universe_n, 4000, 6000, 10000, 18000], universe_n,
        before, after, CHROMATIN, METABOLIC, sc.poisson_binomial,
        FIGS / "dot1l_statistical_calibration.png")
    for name, fig in (("network", fig1), ("calibration", fig2)):
        bad = figures.check_overlaps(fig)
        print(f"[6] {name} figure: {len(bad)} text collisions")

    print("\ndone. outputs in data/, stats/ and figures/")
    return nodes


if __name__ == "__main__":
    main()
