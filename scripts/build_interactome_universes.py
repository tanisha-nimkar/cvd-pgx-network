#!/usr/bin/env python3
"""
scripts/build_interactome_universes.py
Task 4: Construct Functional, TM-Excluded, and Physical STRING Interactomes
Conforms strictly to W9-3 and Part IV deliverables.
"""

import os
import sys
import gzip
import hashlib
import pandas as pd
import numpy as np
import networkx as nx

CONFIG_FILE = "config/config.yaml"
LINKS_FULL = "data/raw/9606.protein.links.full.v12.5.txt.gz"
LINKS_PHYS = "data/raw/9606.protein.physical.links.v12.5.txt.gz"
INFO_FILE = "data/raw/9606.protein.info.v12.5.txt.gz"
OUT_DIR = "data/processed"
AUDIT_OUT = "results/tables/network_construction_manifest.csv"

PANEL_GENES = ["CYP2C19", "CYP2C9", "VKORC1", "CYP4F2", "SLCO1B1", "ABCG2"]

def sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def main():
    print("=== TASK 4: BUILD INTERACTOME UNIVERSES ===")
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs("results/tables", exist_ok=True)

    # 1. Load Symbol Mapping from protein.info
    print("[*] Loading protein info identifier map...")
    info_map = {}
    with gzip.open(INFO_FILE, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                info_map[parts[0]] = parts[1]
    reverse_map = {v: k for k, v in info_map.items()}

    manifest_rows = []

    # =========================================================================
    # 2. Network A: Primary Functional Interactome (All channels >= 700)
    # =========================================================================
    print("\n[*] Building Network A: Primary Functional (All channels >= 700)...")
    edges_func_all = []
    chunk_size = 1_000_000
    with gzip.open(LINKS_FULL, "rt") as f:
        for chunk in pd.read_csv(f, sep=" ", chunksize=chunk_size):
            hi = chunk[chunk["combined_score"] >= 700][["protein1", "protein2", "combined_score"]]
            edges_func_all.append(hi)

    df_func_all = pd.concat(edges_func_all, ignore_index=True)
    # Remove self-loops
    df_func_all = df_func_all[df_func_all["protein1"] != df_func_all["protein2"]]
    # Deduplicate undirected edges
    df_func_all["edge_key"] = [tuple(sorted([p1, p2])) for p1, p2 in zip(df_func_all["protein1"], df_func_all["protein2"])]
    df_func_all = df_func_all.drop_duplicates(subset=["edge_key"]).copy()

    G_func = nx.Graph()
    for _, r in df_func_all.iterrows():
        G_func.add_edge(r["protein1"], r["protein2"], weight=int(r["combined_score"]))

    lcc_nodes_func = max(nx.connected_components(G_func), key=len)
    G_func_lcc = G_func.subgraph(lcc_nodes_func).copy()

    # Save edges
    out_func_edges = os.path.join(OUT_DIR, "string_functional_universe_edges.tsv")
    rows_func = [{
        "protein1": u, "protein2": v,
        "symbol1": info_map.get(u, u), "symbol2": info_map.get(v, v),
        "combined_score": d["weight"]
    } for u, v, d in G_func_lcc.edges(data=True)]
    pd.DataFrame(rows_func).to_csv(out_func_edges, sep="\t", index=False)

    panel_in_func = [g for g in PANEL_GENES if reverse_map.get(g) in G_func_lcc]
    manifest_rows.append({
        "network_name": "Primary Functional Network",
        "description": "STRING v12.5 functional associations, all channels >= 700",
        "total_nodes_raw": G_func.number_of_nodes(),
        "total_edges_raw": G_func.number_of_edges(),
        "lcc_nodes": G_func_lcc.number_of_nodes(),
        "lcc_edges": G_func_lcc.number_of_edges(),
        "lcc_percentage": round(G_func_lcc.number_of_nodes() / G_func.number_of_nodes() * 100, 2),
        "panel_genes_in_lcc": ";".join(panel_in_func),
        "output_edge_file": out_func_edges,
        "output_sha256": sha256(out_func_edges)
    })
    print(f"[+] Network A (Functional LCC): {G_func_lcc.number_of_nodes():,} nodes, {G_func_lcc.number_of_edges():,} edges.")

    # =========================================================================
    # 3. Network B: Sensitivity Functional Interactome (TM-excluded >= 700)
    # =========================================================================
    print("\n[*] Building Network B: Sensitivity Functional (Text-Mining Excluded >= 700)...")
    non_text = ["neighborhood", "neighborhood_transferred", "fusion", "cooccurence",
                "homology", "coexpression", "coexpression_transferred",
                "experiments", "experiments_transferred", "database", "database_transferred"]
    text_cols = ["textmining", "textmining_transferred"]

    edges_func_tm_ex = []
    with gzip.open(LINKS_FULL, "rt") as f:
        for chunk in pd.read_csv(f, sep=" ", chunksize=chunk_size):
            hi = chunk[chunk["combined_score"] >= 700].copy()
            if len(hi) == 0:
                continue
            is_tm_only = (hi[text_cols].max(axis=1) > 0) & (hi[non_text].max(axis=1) == 0)
            hi_retained = hi[~is_tm_only][["protein1", "protein2", "combined_score"]]
            edges_func_tm_ex.append(hi_retained)

    df_func_tm = pd.concat(edges_func_tm_ex, ignore_index=True)
    df_func_tm = df_func_tm[df_func_tm["protein1"] != df_func_tm["protein2"]]
    df_func_tm["edge_key"] = [tuple(sorted([p1, p2])) for p1, p2 in zip(df_func_tm["protein1"], df_func_tm["protein2"])]
    df_func_tm = df_func_tm.drop_duplicates(subset=["edge_key"]).copy()

    G_tm = nx.Graph()
    for _, r in df_func_tm.iterrows():
        G_tm.add_edge(r["protein1"], r["protein2"], weight=int(r["combined_score"]))

    lcc_nodes_tm = max(nx.connected_components(G_tm), key=len)
    G_tm_lcc = G_tm.subgraph(lcc_nodes_tm).copy()

    out_tm_edges = os.path.join(OUT_DIR, "string_functional_tm_excluded_edges.tsv")
    rows_tm = [{
        "protein1": u, "protein2": v,
        "symbol1": info_map.get(u, u), "symbol2": info_map.get(v, v),
        "combined_score": d["weight"]
    } for u, v, d in G_tm_lcc.edges(data=True)]
    pd.DataFrame(rows_tm).to_csv(out_tm_edges, sep="\t", index=False)

    panel_in_tm = [g for g in PANEL_GENES if reverse_map.get(g) in G_tm_lcc]
    manifest_rows.append({
        "network_name": "Sensitivity Functional Network (TM-Excluded)",
        "description": "STRING v12.5 functional, combined >= 700, text-mining-only edges excluded",
        "total_nodes_raw": G_tm.number_of_nodes(),
        "total_edges_raw": G_tm.number_of_edges(),
        "lcc_nodes": G_tm_lcc.number_of_nodes(),
        "lcc_edges": G_tm_lcc.number_of_edges(),
        "lcc_percentage": round(G_tm_lcc.number_of_nodes() / G_tm.number_of_nodes() * 100, 2),
        "panel_genes_in_lcc": ";".join(panel_in_tm),
        "output_edge_file": out_tm_edges,
        "output_sha256": sha256(out_tm_edges)
    })
    print(f"[+] Network B (TM-Excluded LCC): {G_tm_lcc.number_of_nodes():,} nodes, {G_tm_lcc.number_of_edges():,} edges.")

    # =========================================================================
    # 4. Network C: Physical Interactome (Physical Subnetwork >= 700)
    # =========================================================================
    if os.path.exists(LINKS_PHYS):
        print("\n[*] Building Network C: Physical Subnetwork (Combined >= 700)...")
        edges_phys = []
        with gzip.open(LINKS_PHYS, "rt") as f:
            for chunk in pd.read_csv(f, sep=" ", chunksize=chunk_size):
                hi = chunk[chunk["combined_score"] >= 700][["protein1", "protein2", "combined_score"]]
                edges_phys.append(hi)

        df_phys = pd.concat(edges_phys, ignore_index=True)
        df_phys = df_phys[df_phys["protein1"] != df_phys["protein2"]]
        df_phys["edge_key"] = [tuple(sorted([p1, p2])) for p1, p2 in zip(df_phys["protein1"], df_phys["protein2"])]
        df_phys = df_phys.drop_duplicates(subset=["edge_key"]).copy()

        G_phys = nx.Graph()
        for _, r in df_phys.iterrows():
            G_phys.add_edge(r["protein1"], r["protein2"], weight=int(r["combined_score"]))

        lcc_nodes_phys = max(nx.connected_components(G_phys), key=len)
        G_phys_lcc = G_phys.subgraph(lcc_nodes_phys).copy()

        out_phys_edges = os.path.join(OUT_DIR, "string_physical_universe_edges.tsv")
        rows_phys = [{
            "protein1": u, "protein2": v,
            "symbol1": info_map.get(u, u), "symbol2": info_map.get(v, v),
            "combined_score": d["weight"]
        } for u, v, d in G_phys_lcc.edges(data=True)]
        pd.DataFrame(rows_phys).to_csv(out_phys_edges, sep="\t", index=False)

        panel_in_phys = [g for g in PANEL_GENES if reverse_map.get(g) in G_phys_lcc]
        manifest_rows.append({
            "network_name": "Physical Interaction Network",
            "description": "STRING v12.5 physical binding subnetwork, combined >= 700",
            "total_nodes_raw": G_phys.number_of_nodes(),
            "total_edges_raw": G_phys.number_of_edges(),
            "lcc_nodes": G_phys_lcc.number_of_nodes(),
            "lcc_edges": G_phys_lcc.number_of_edges(),
            "lcc_percentage": round(G_phys_lcc.number_of_nodes() / G_phys.number_of_nodes() * 100, 2),
            "panel_genes_in_lcc": ";".join(panel_in_phys),
            "output_edge_file": out_phys_edges,
            "output_sha256": sha256(out_phys_edges)
        })
        print(f"[+] Network C (Physical LCC): {G_phys_lcc.number_of_nodes():,} nodes, {G_phys_lcc.number_of_edges():,} edges.")
    else:
        print(f"[!] Physical links file not found at {LINKS_PHYS}. Skipping Network C build.")

    # Save Construction Manifest
    df_man = pd.DataFrame(manifest_rows)
    df_man.to_csv(AUDIT_OUT, index=False)
    print(f"\n[+] Saved network construction manifest to {AUDIT_OUT}")
    print(df_man[["network_name", "lcc_nodes", "lcc_edges", "lcc_percentage", "panel_genes_in_lcc"]].to_string(index=False))

if __name__ == "__main__":
    main()
