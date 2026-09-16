import networkx as nx
import json
from pathlib import Path
from datetime import date

PROC_DIR = Path("data/processed")
TODAY = date.today().strftime("%Y%m%d")
JSON_PATH = PROC_DIR / f"cvd_pgx_network_{TODAY}.json"

def verify_network():
    if not JSON_PATH.exists():
        print(f"Error: Network file not found at {JSON_PATH}. Run build_network.py first.")
        return

    with open(JSON_PATH, "r") as fp:
        data = json.load(fp)

    G = nx.node_link_graph(data)

    print("=== Network Topological Metrics ===")
    print(f"Total Nodes: {G.number_of_nodes()}")
    print(f"Total Edges: {G.number_of_edges()}")

    node_types = {}
    for _, attrs in G.nodes(data=True):
        ntype = attrs.get("node_type", "Unknown")
        node_types[ntype] = node_types.get(ntype, 0) + 1

    print("\nNode Breakdown by Category:")
    for ntype, count in node_types.items():
        print(f"  - {ntype}: {count}")

    print("\nDegree Centrality (Top 5 Hub Genes/Drugs):")
    centrality = nx.degree_centrality(G)
    sorted_centrality = sorted(centrality.items(), key=lambda x: x[1], reverse=True)[:5]
    for node, score in sorted_centrality:
        print(f"  - {node}: {score:.4f}")

if __name__ == "__main__":
    verify_network()
