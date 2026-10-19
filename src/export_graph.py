"""GraphML writer - no third-party graph library required."""
import pandas as pd

NODE_ATTRS = [("module", "string"), ("tier", "string"), ("recurrence", "int"),
              ("recurrence_q_BH", "double"), ("max_r", "double"), ("mean_r", "double"),
              ("mean_r_DOT1L_complex", "double"), ("mean_r_MENIN_MLL", "double"),
              ("string_degree", "int"), ("module_stability", "double"),
              ("is_seed", "boolean"), ("high_confidence", "boolean")]


def _esc(v):
    return str(v).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_graphml(path, nodes, edges, graph_id="DOT1L_H3K79_dependency_network"):
    """nodes: DataFrame indexed by gene. edges: columns source,target,edge_type,weight."""
    ids = set(nodes.index)
    with open(path, "w") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n'
                 '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n')
        for key, typ in NODE_ATTRS:
            if key in nodes.columns:
                fh.write(f'  <key id="{key}" for="node" attr.name="{key}" attr.type="{typ}"/>\n')
        fh.write('  <key id="edge_type" for="edge" attr.name="edge_type" attr.type="string"/>\n')
        fh.write('  <key id="weight" for="edge" attr.name="weight" attr.type="double"/>\n')
        fh.write(f'  <graph id="{graph_id}" edgedefault="undirected">\n')
        for gene, row in nodes.iterrows():
            fh.write(f'    <node id="{_esc(gene)}">\n')
            for key, _ in NODE_ATTRS:
                if key not in nodes.columns:
                    continue
                val = row[key]
                if pd.isna(val):
                    continue
                if isinstance(val, (bool,)):
                    val = str(val).lower()
                fh.write(f'      <data key="{key}">{_esc(val)}</data>\n')
            fh.write("    </node>\n")
        n = 0
        for i, (src, tgt, etype, w) in enumerate(edges.itertuples(index=False)):
            if src in ids and tgt in ids:
                fh.write(f'    <edge id="e{i}" source="{_esc(src)}" target="{_esc(tgt)}">'
                         f'<data key="edge_type">{etype}</data>'
                         f'<data key="weight">{w}</data></edge>\n')
                n += 1
        fh.write("  </graph>\n</graphml>\n")
    return len(ids), n
