"""Step 2 - select nodes, fetch STRING topology, resolve modules.

Outputs: data/dot1l_h3k79_network_edges.csv
         data/dot1l_h3k79_network.graphml   (written by export_graph.py)
         module assignment returned to the caller
"""
import io
import random
import time
from collections import defaultdict

import pandas as pd
import requests

from config import (SEEDS, STRING_API, STRING_SPECIES, STRING_MIN_SCORE, CALLER_IDENTITY,
                    STRING_RESOLVE_1, STRING_RESOLVE_2, MIN_RECURRENCE, SINGLE_SEED_R,
                    RANDOM_SEED)


def select_nodes(matrix):
    """Admit a gene if positively correlated with >=MIN_RECURRENCE seeds, or >=SINGLE_SEED_R with one."""
    pos = matrix.where(matrix > 0)
    summary = pd.DataFrame({
        "recurrence": pos.notna().sum(axis=1),
        "max_r": pos.max(axis=1),
        "mean_r": pos.mean(axis=1),
        "mean_r_DOT1L_complex": pos[["DOT1L", "MLLT10", "MLLT6", "MLLT1", "MLLT3"]].mean(axis=1),
        "mean_r_MENIN_MLL": pos[["KMT2A", "MEN1"]].mean(axis=1),
        "n_neg": (matrix < 0).sum(axis=1),
    })
    summary = summary[summary.recurrence > 0]
    keep = (summary.recurrence >= MIN_RECURRENCE) | (summary.max_r >= SINGLE_SEED_R)
    summary["in_network"] = keep
    return summary


def _post(endpoint, data, timeout):
    r = requests.post(f"{STRING_API}/tsv/{endpoint}", data=data, timeout=timeout)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text), sep="\t")


def string_graph(genes, min_score=STRING_MIN_SCORE):
    """Map symbols to STRING ids and fetch the induced subnetwork. Returns (idmap, edges)."""
    idmap = _post("get_string_ids",
                  {"identifiers": "\r".join(sorted(genes)), "species": STRING_SPECIES,
                   "limit": 1, "echo_query": 1, "caller_identity": CALLER_IDENTITY},
                  timeout=300).drop_duplicates("queryItem")
    edges = _post("network",
                  {"identifiers": "\r".join(idmap.stringId), "species": STRING_SPECIES,
                   "add_nodes": 0, "caller_identity": CALLER_IDENTITY},
                  timeout=600)
    edges = edges.rename(columns={"preferredName_A": "A", "preferredName_B": "B"})
    edges = edges[(edges.score >= min_score) & (edges.A != edges.B)]
    return idmap, edges.drop_duplicates(subset=["A", "B"])


def label_propagation(nodes, edges, seed=RANDOM_SEED, max_iter=100):
    """Weighted label propagation. Deterministic given `seed`."""
    adj = defaultdict(dict)
    for a, b, w in zip(edges.A, edges.B, edges.score):
        adj[a][b] = w
        adj[b][a] = w
    labels = {n: i for i, n in enumerate(nodes)}
    rng = random.Random(seed)
    order = list(nodes)
    for _ in range(max_iter):
        rng.shuffle(order)
        changed = 0
        for node in order:
            neighbours = adj.get(node)
            if not neighbours:
                continue
            tally = {}
            for nb, w in neighbours.items():
                if nb in labels:
                    tally[labels[nb]] = tally.get(labels[nb], 0) + w
            # sort before max so ties resolve deterministically
            best = max(sorted(tally.items()), key=lambda kv: kv[1])[0]
            if best != labels[node]:
                labels[node] = best
                changed += 1
        if changed == 0:
            break
    return labels


def _groups(labels):
    out = defaultdict(list)
    for node, lab in labels.items():
        out[lab].append(node)
    return out


def hierarchical_modules(idmap, edges, min_component=10):
    """Partition the STRING graph, then re-resolve any component too large and
    generically connected to be one complex.

    At score >= 0.4 the chromatin genes collapse into a single component whose
    only shared enrichment term is "nucleoplasm". That component is re-resolved
    at >= STRING_RESOLVE_1 and its largest surviving component again at
    >= STRING_RESOLVE_2, which separates the writer, MENIN-MLL, PRC2, MYB-p300
    and ncBAF groups. Components not implicated keep their original label.

    Returns {gene: label}.
    """
    nodes = list(idmap.preferredName)
    final = dict(label_propagation(nodes, edges))
    for cut in (STRING_RESOLVE_1, STRING_RESOLVE_2):
        groups = _groups(final)
        biggest = max(groups.values(), key=len)
        if len(biggest) < min_component:
            break
        sub_edges = edges[edges.A.isin(biggest) & edges.B.isin(biggest) & (edges.score >= cut)]
        sub = label_propagation(biggest, sub_edges)
        offset = max(final.values()) + 1
        for gene, lab in sub.items():
            final[gene] = offset + lab
    return final


# Display names are anchored to marker genes rather than to the STRING
# enrichment string, which is not stable enough to key a layout on. The
# enrichment term is still computed and reported alongside, so every name
# can be checked against independent evidence.
MODULE_MARKERS = {
    "H3K79 writer (DotCom)": {"DOT1L", "MLLT10"},
    "MENIN-MLL-LEDGF": {"KMT2A", "MEN1"},
    "PRC2": {"EZH2", "SUZ12"},
    "ncBAF": {"BRD9", "SMARCD1"},
    "MYB-p300 coactivator": {"MYB", "EP300"},
    "Lymphoid/myeloid TFs": {"IKZF1", "TCF3"},
    "MOZ/MORF-HBO1 HAT": {"KAT6A", "BRPF1"},
    "SAGA complex": {"TADA1", "SUPT7L"},
    "Mediator complex": {"MED13", "MED23"},
    "Mitochondrial OXPHOS": {"SDHA", "UQCRC1"},
    "Nucleotide biosynthesis": {"ATIC", "PAICS"},
    "NAD synthesis": {"NAMPT", "NMNAT1"},
    "Heme / sugar metabolism": {"ALAD", "ALAS1"},
}


def name_module(members):
    """Marker-anchored name for a module, or None if no marker set matches."""
    members = set(members)
    best, best_hits = None, 0
    for name, markers in MODULE_MARKERS.items():
        hits = len(markers & members)
        if hits > best_hits or (hits == best_hits and hits and name < (best or "~")):
            best, best_hits = name, hits
    return best if best_hits >= 2 else None


def enrich(genes, categories=("Component", "Process", "RCTM", "KEGG", "COMPARTMENTS")):
    """Lowest-FDR STRING functional enrichment term for a gene set."""
    tbl = _post("enrichment",
                {"identifiers": "\r".join(sorted(genes)), "species": STRING_SPECIES,
                 "caller_identity": CALLER_IDENTITY}, timeout=300)
    tbl = tbl[tbl.category.isin(categories)].sort_values("fdr")
    if tbl.empty:
        return {"term": "(no enrichment)", "fdr": float("nan"), "n": 0}
    top = tbl.iloc[0]
    return {"term": top.description, "fdr": float(top.fdr), "n": int(top.number_of_genes)}


def name_modules(groups, min_size=4, pause=1.0):
    """Name each module by its lowest-FDR STRING enrichment term."""
    named = {}
    for lab, genes in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        if len(genes) >= min_size:
            named[lab] = enrich(genes)
            time.sleep(pause)   # courtesy to the shared STRING endpoint
        else:
            named[lab] = {"term": None, "fdr": float("nan"), "n": 0}
    return named
