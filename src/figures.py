"""Step 4 - the two manuscript figures."""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from config import FIGS

SUPERGROUP_COLOUR = {"axis": "#14477a", "chrom": "#3c7fb1", "tf": "#4f8a47", "met": "#b5762f"}
GREYS = ["#c9c9c9", "#ababab", "#8d8d8d", "#6f6f6f"]


def _style():
    mpl.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "font.size": 8,
        "axes.titlesize": 8, "axes.labelsize": 8,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "savefig.bbox": "tight",
    })


def check_overlaps(fig):
    """Geometric legibility check: returns colliding visible text pairs."""
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    return [(a.get_text()[:24], b.get_text()[:24])
            for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]


NETWORK_TITLE = ("DOT1L co-dependencies resolve into the MENIN-MLL axis, four chromatin "
                 "co-activator\ncomplexes, and lineage-independent metabolic modules")


def network_figure(nodes, inter_edges, columns, supergroup, short_names,
                   path=None, title=NETWORK_TITLE, caption=None, n_omitted=0):
    """Module-block layout: one column per group of modules, nodes ordered by max_r."""
    _style()
    x_step, dy, title_gap, module_gap = 1.95, 0.118, 0.20, 0.30
    pos, blocks = {}, {}
    for col, (x, modules) in enumerate(columns):
        y = 0.0
        for module in modules:
            members = nodes[nodes.module == module].sort_values("max_r", ascending=False).index.tolist()
            blocks[module] = (col * x_step, y, len(members))
            y -= title_gap
            for gene in members:
                pos[gene] = (col * x_step, y)
                y -= dy
            y -= module_gap

    smax, smin = nodes.loc[list(pos)].max_r.max(), nodes.loc[list(pos)].max_r.min()
    fig, ax = plt.subplots(figsize=(9.0, 6.6))
    for a, b, score in inter_edges:
        if a in pos and b in pos:
            (x1, y1), (x2, y2) = pos[a], pos[b]
            ax.annotate("", xy=(x2, y2), xytext=(x1, y1), zorder=1,
                        arrowprops=dict(arrowstyle="-", color="#9fb6c9", lw=0.55 * score,
                                        alpha=0.55, connectionstyle="arc3,rad=0.16",
                                        shrinkA=3, shrinkB=3))
    for module, (x, y, n) in blocks.items():
        ax.text(x - 0.045, y + 0.05, f"{short_names.get(module, module)}  ({n})", fontsize=7,
                fontweight="bold", color=SUPERGROUP_COLOUR[supergroup[module]], ha="left", va="bottom")
    for gene, (x, y) in pos.items():
        row = nodes.loc[gene]
        ax.scatter([x], [y], s=16 + 180 * (row.max_r - smin) / (smax - smin),
                   c=SUPERGROUP_COLOUR[supergroup[row.module]], edgecolors="white",
                   linewidths=1.5 if row.recurrence >= 4 else 0.4, zorder=4)
        ax.text(x + 0.085, y, gene, fontsize=6.2, style="italic", va="center", ha="left",
                zorder=5, fontweight="bold" if row.get("is_seed", False) else "normal")
    # key
    kx, ky = 4 * x_step, -2.80
    ax.text(kx - 0.045, ky + 0.33, "Key", fontsize=7, fontweight="bold", color="#3a3a3a")
    for i, (rv, lbl) in enumerate([(smax, f"strongest correlation with a seed, r = {smax:.2f}"),
                                   (0.35, "r = 0.35"), (0.20, "r = 0.20")]):
        ax.scatter([kx], [ky - i * 0.19], s=16 + 180 * (rv - smin) / (smax - smin),
                   c="#8a8a8a", edgecolors="white", linewidths=0.4)
        ax.text(kx + 0.11, ky - i * 0.19, lbl, fontsize=6.2, va="center")
    ax.scatter([kx], [ky - 0.66], s=95, c="#8a8a8a", edgecolors="white", linewidths=1.5)
    ax.text(kx + 0.11, ky - 0.66, "thick ring: co-dependent with $\\geq$4 of 7 seeds",
            fontsize=6.2, va="center")
    ax.plot([kx - 0.04, kx + 0.06], [ky - 0.85, ky - 0.85], color="#9fb6c9", lw=0.55)
    ax.text(kx + 0.11, ky - 0.85, "STRING association $\\geq$0.7, between modules",
            fontsize=6.2, va="center")
    ax.text(kx - 0.045, ky - 1.04, "bold gene = CRISPR seed profile", fontsize=6.2, va="center")

    ax.set_xlim(-0.22, 4 * x_step + 2.15)
    ax.set_ylim(-4.75, 0.62)
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=8.5, loc="left", x=0.0, y=1.005)
    if caption is None:
        caption = (
            "Nodes: genes positively correlated with $\\geq$2 of 7 seed co-dependency profiles "
            "(DOT1L, MLLT10, MLLT6, MLLT1, MLLT3, KMT2A, MEN1), or r $\\geq$0.35 with a single "
            "seed; dot area scales with the strongest of those\ncorrelations. Modules from "
            "hierarchical weighted label propagation on the STRING v12.0 graph (score "
            "$\\geq$0.4, chromatin core re-resolved at $\\geq$0.7 then $\\geq$0.9), named by "
            "marker genes and verified\nagainst STRING enrichment (stats/module_provenance.csv). "
            f"{n_omitted} qualifying nodes fell outside a named module and are omitted here; all "
            "appear in the node table. Seed profiles are DepMap Public\n26Q1+Score (Chronos) "
            "pan-cancer correlations; MLLT6 (AF17) and MLLT3 (AF9) returned no co-dependency "
            "with DOT1L and contribute no nodes of their own.")
        ax.text(-0.22, -4.70, caption, fontsize=5.8, color="#5a5a5a", va="bottom", ha="left")
    if path:
        fig.savefig(path)
    return fig


def calibration_figure(recurrence_counts, seed_sizes, universe_sizes, focal_n,
                       module_before, module_after, chromatin, metabolic,
                       poisson_binomial_fn, path=None):
    """Panel a: FDR vs recurrence threshold. Panel b: module survival."""
    _style()
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.3), gridspec_kw={"width_ratios": [1, 1.2]})
    ths = [2, 3, 4, 5]

    ax = axes[0]
    for i, n in enumerate(universe_sizes):
        dist = poisson_binomial_fn(seed_sizes / n)
        fdr = [min(1, n * dist[k:].sum() / recurrence_counts[k]) for k in ths]
        focal = (n == focal_n)
        ax.plot(ths, fdr, "o-", lw=2.0 if focal else 1.0, ms=5 if focal else 3.2,
                color="#14477a" if focal else GREYS[min(i - 1, len(GREYS) - 1)],
                zorder=4 if focal else 2,
                label=(f"{round(n):,} (estimated from the data)" if focal else f"{round(n):,}"))
    ax.axhline(0.10, ls=(0, (3, 3)), lw=0.7, color="#b5762f")
    ax.text(5.06, 0.105, "10%", fontsize=6, color="#b5762f", va="bottom", ha="right")
    ax.set_yscale("log")
    ax.set_ylim(4e-4, 2.2)
    ax.set_xticks(ths)
    ax.margins(x=0.09)
    ax.set_xlabel("recurrence threshold (seed lists with the gene)")
    ax.set_ylabel("false-discovery rate of the filter")
    ax.set_yticks([1e-3, 1e-2, 1e-1, 1])
    ax.set_yticklabels(["0.1%", "1%", "10%", "100%"])
    leg = ax.legend(title="assumed gene universe", fontsize=5.6, title_fontsize=6,
                    loc="lower left", handlelength=1.5, borderaxespad=0.1)
    leg._legend_box.align = "left"
    ax.set_title("A gene in $\\geq$2 lists is mostly chance;\n$\\geq$3 is defensible", loc="left")

    ax = axes[1]
    order = chromatin + metabolic
    before = [module_before[m] for m in order]
    after = [module_after[m] for m in order]
    y = np.arange(len(order))[::-1]
    ax.barh(y, before, height=0.64, color="#e3e8ed", zorder=2)
    ax.barh(y, after, height=0.64, zorder=3,
            color=["#14477a"] * len(chromatin) + ["#b5762f"] * len(metabolic))
    for yy, a, b in zip(y, before, after):
        ax.text(a + 0.4, yy, f"{b}/{a}", va="center", fontsize=6, color="#4a4a4a", zorder=4)
        if b == 0:
            ax.plot([0.08], [yy], marker="|", ms=6, color="#b5762f", zorder=5)
    ax.set_yticks(y)
    ax.set_yticklabels([m.replace(" (DotCom)", "") for m in order], fontsize=6.4)
    ax.set_xlabel("genes in module")
    ax.set_xlim(0, max(before) + 3.5)
    ax.margins(y=0.05)
    ax.legend(handles=[Patch(fc="#e3e8ed", label="all nodes in the module"),
                       Patch(fc="#14477a", label="surviving at $\\geq$3: chromatin"),
                       Patch(fc="#b5762f", label="surviving at $\\geq$3: metabolic")],
              fontsize=6, loc="lower right", borderaxespad=0.2, handlelength=1.1)
    ax.set_title("Chromatin modules survive the filter,\nmetabolic modules largely do not", loc="left")
    for ax_, letter in zip(axes, "ab"):
        ax_.text(-0.14, 1.06, letter, transform=ax_.transAxes, fontsize=10, fontweight="bold")
    fig.subplots_adjust(wspace=0.50)
    if path:
        fig.savefig(path)
    return fig
