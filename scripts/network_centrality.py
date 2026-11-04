#!/usr/bin/env python3
"""
scripts/network_centrality.py
Week 8: Network Centrality Computation on Week 7 STRING Ego-Network
Conforms strictly to W8.1 & W8.2 SOP guidelines.
"""

import os
import sys
import hashlib
import yaml
import networkx as nx
import pandas as pd
from scipy.stats import spearmanr

CONFIG_FILE = "config/config.yaml"
GRAPHML_FILE = "results/network/network.graphml"
OUT_CSV = "results/tables/centrality.csv"
EXPECTED_SHA = "33bcd28abc5583d264ebda841b09cdba2364ceeccd5d4d63b5a6501b068c8b23"

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def main():
    print("=== WEEK 8: NETWORK CENTRALITY PIPELINE ===")

    # Verify input graph checksum
    actual_sha = compute_sha256(GRAPHML_FILE)
    print(f"[*] Input GraphML SHA-256: {actual_sha}")
    if actual_sha != EXPECTED_SHA:
        print(f"[!] STOP: GraphML checksum mismatch! Expected {EXPECTED_SHA}, got {actual_sha}")
        sys.exit(1)

    with open(CONFIG_FILE, "r") as f:
        cfg = yaml.safe_load(f)
    panel = cfg.get("gene_panel", [])

    G = nx.read_graphml(GRAPHML_FILE)
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    print(f"[*] Read graph: {n_nodes} nodes, {n_edges} edges.")

    # 1. Degree & Betweenness (unweighted, exact)
    deg_dict = dict(G.degree())
    dc_dict = nx.degree_centrality(G)
    bc_dict = nx.betweenness_centrality(G, normalized=True, weight=None)

    # 2. Build rows
    rows = []
    for nid, data in G.nodes(data=True):
        sym = data.get("gene_symbol", "")
        role = data.get("node_role", "")
        is_panel = (sym in panel)
        is_iso = bool(deg_dict[nid] == 0)

        rows.append({
            "node_id": nid,
            "gene_symbol": sym,
            "node_role": role,
            "degree": int(deg_dict[nid]),
            "degree_centrality": float(dc_dict[nid]),
            "betweenness_centrality": float(bc_dict[nid]),
            "is_panel_gene": is_panel,
            "is_isolated": is_iso,
            "network_sha256": actual_sha,
            "network_approval_status": "provisional_pending_supervisor",
            "networkx_version": nx.__version__
        })

    df = pd.DataFrame(rows)

    # QC checks
    assert len(df) == n_nodes, f"Row count mismatch: expected {n_nodes}, got {len(df)}"
    assert df["degree_centrality"].isna().sum() == 0, "Missing degree centrality values!"
    assert df["betweenness_centrality"].isna().sum() == 0, "Missing betweenness centrality values!"
    assert df["node_id"].duplicated().sum() == 0, "Duplicate node IDs detected!"

    os.makedirs("results/tables", exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"[+] Saved centrality metrics to: {OUT_CSV} ({len(df)} rows)")

    # Diagnostic: Spearman correlation between degree and betweenness
    corr, pval = spearmanr(df["degree"], df["betweenness_centrality"])
    print(f"\n[*] Diagnostic: Spearman correlation (Degree vs. Betweenness): rho = {corr:.4f} (p = {pval:.4e})")

    # Panel gene ranking table
    print("\n=== PANEL PHARMACOGENES CENTRALITY SUMMARY ===")
    panel_df = df[df["is_panel_gene"]].sort_values(by="degree", ascending=False)
    cols_show = ["gene_symbol", "degree", "degree_centrality", "betweenness_centrality", "is_isolated"]
    print(panel_df[cols_show].to_string(index=False))

    # W8.2 CYP2C19 Informal Check
    cyp2c19_row = df[df["gene_symbol"] == "CYP2C19"].iloc[0]
    cyp2c19_id = cyp2c19_row["node_id"]
    neighbors = [G.nodes[nbr].get("gene_symbol", nbr) for nbr in G.neighbors(cyp2c19_id)]
    print("\n=== CYP2C19 INFORMAL CONNECTIVITY CHECK ===")
    print(f"Node ID:               {cyp2c19_id}")
    print(f"Degree:                {cyp2c19_row['degree']}")
    print(f"Degree Centrality:     {cyp2c19_row['degree_centrality']:.4f}")
    print(f"Betweenness Centrality:{cyp2c19_row['betweenness_centrality']:.4f}")
    print(f"First 10 Neighbors:    {neighbors[:10]}")

if __name__ == "__main__":
    main()
