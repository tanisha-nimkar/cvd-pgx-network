#!/usr/bin/env python3
"""
scripts/build_string_universe.py
Week 8: Build STRING Background Interactome (Universe LCC)
Conforms strictly to W8.3 SOP guidelines.
"""

import os
import sys
import hashlib
import gzip
import yaml
import pandas as pd
import numpy as np
import networkx as nx

CONFIG_FILE = "config/config.yaml"
AUDIT_FILE = "results/network/string_edges_audit.tsv"
LINKS_FILE = "data/raw/9606.protein.links.full.v12.5.txt.gz"
INFO_FILE = "data/raw/9606.protein.info.v12.5.txt.gz"
OUT_EDGES = "data/processed/string_universe_edges.tsv"

EXPECTED_LINKS_SHA = "c645be3e636a896f82b0acf63d181aa58a1320b8df9b4106ae9348e68e086ad8"
EXPECTED_INFO_SHA = "144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c"

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def main():
    print("=== WEEK 8: STRING UNIVERSE LCC CONSTRUCTION ===")

    # 1. Verify Checksums
    print("[*] Verifying raw file SHA-256 checksums...")
    sha_links = compute_sha256(LINKS_FILE)
    sha_info = compute_sha256(INFO_FILE)
    if sha_links != EXPECTED_LINKS_SHA or sha_info != EXPECTED_INFO_SHA:
        print("[!] STOP: Checksum mismatch on STRING raw bulk files!")
        sys.exit(1)
    print("  [Pass] Both STRING v12.5 bulk files match verified hashes.")

    with open(CONFIG_FILE, "r") as f:
        cfg = yaml.safe_load(f)
    panel = cfg.get("gene_panel", [])

    # 2. Test File-Derived Classification vs. Week 7 Audit File
    print("\n[*] Validating bulk file channel mapping against Week 7 audit file...")
    df_audit = pd.read_csv(AUDIT_FILE, sep="\t")
    # In audit: stringId_A, stringId_B, edge_status, score
    audit_retained_pairs = set(
        tuple(sorted([r["stringId_A"], r["stringId_B"]]))
        for _, r in df_audit[df_audit["edge_status"] == "retained"].iterrows()
    )
    audit_tm_pairs = set(
        tuple(sorted([r["stringId_A"], r["stringId_B"]]))
        for _, r in df_audit[df_audit["edge_status"] == "excluded_text_mining_only"].iterrows()
    )

    print(f"    Week 7 audit contains {len(audit_retained_pairs)} retained seed edges, {len(audit_tm_pairs)} text-mining-only edges.")

    # 3. Stream and filter links.full in chunks
    print("\n[*] Streaming and filtering protein.links.full (combined_score >= 700)...")
    non_text_cols = [
        "neighborhood", "neighborhood_transferred", "fusion", "cooccurence",
        "homology", "coexpression", "coexpression_transferred",
        "experiments", "experiments_transferred", "database", "database_transferred"
    ]
    text_cols = ["textmining", "textmining_transferred"]

    retained_edges = []
    chunk_size = 1_000_000
    total_processed = 0

    with gzip.open(LINKS_FILE, "rt") as f:
        for chunk in pd.read_csv(f, sep=" ", chunksize=chunk_size):
            total_processed += len(chunk)
            # Threshold >= 700
            hi_conf = chunk[chunk["combined_score"] >= 700].copy()
            if len(hi_conf) == 0:
                continue

            # Identify text-mining only
            max_non_text = hi_conf[non_text_cols].max(axis=1)
            max_text = hi_conf[text_cols].max(axis=1)
            is_tm_only = (max_text > 0) & (max_non_text == 0)

            # Filter out text-mining-only
            retained = hi_conf[~is_tm_only][["protein1", "protein2", "combined_score"]]
            retained_edges.append(retained)
            print(f"    Processed {total_processed:,} lines | Retained high-conf candidates: {sum(len(x) for x in retained_edges):,}", end="\r")

    print(f"\n[+] Total rows processed: {total_processed:,}")
    df_all_retained = pd.concat(retained_edges, ignore_index=True)

    # 4. Deduplicate undirected pairs
    print("[*] Deduplicating undirected edges...")
    df_all_retained["edge_key"] = [
        tuple(sorted([p1, p2])) for p1, p2 in zip(df_all_retained["protein1"], df_all_retained["protein2"])
    ]
    df_dedup = df_all_retained.drop_duplicates(subset=["edge_key"]).copy()
    print(f"[+] Unique undirected high-confidence edges: {len(df_dedup):,}")

    # Check consistency against Week 7 seed edges
    universe_pairs = set(df_dedup["edge_key"])
    missing_from_universe = audit_retained_pairs - universe_pairs
    print(f"[*] Audit concordance check: {len(missing_from_universe)} retained audit seed edges missing from universe.")
    if len(missing_from_universe) > 0:
        print(f"    Missing sample: {list(missing_from_universe)[:3]}")

    tm_in_universe = audit_tm_pairs.intersection(universe_pairs)
    print(f"[*] Zero text-mining audit edges retained in universe: {len(tm_in_universe) == 0}")

    # 5. Build Graph and extract Largest Connected Component (LCC)
    print("\n[*] Constructing NetworkX universe graph...")
    G = nx.Graph()
    for _, r in df_dedup.iterrows():
        G.add_edge(r["protein1"], r["protein2"], combined_score=int(r["combined_score"]))

    n_nodes_total = G.number_of_nodes()
    n_edges_total = G.number_of_edges()
    print(f"    Raw Universe Graph: {n_nodes_total:,} nodes, {n_edges_total:,} edges.")

    components = sorted(nx.connected_components(G), key=len, reverse=True)
    lcc_nodes = components[0]
    G_lcc = G.subgraph(lcc_nodes).copy()

    n_lcc_nodes = G_lcc.number_of_nodes()
    n_lcc_edges = G_lcc.number_of_edges()
    lcc_share = (n_lcc_nodes / n_nodes_total) * 100.0
    print(f"    Largest Connected Component (LCC): {n_lcc_nodes:,} nodes ({lcc_share:.2f}%), {n_lcc_edges:,} edges.")

    # 6. Load protein.info to map ENSP -> Gene Symbol
    print("\n[*] Loading protein.info to map IDs...")
    info_map = {}
    with gzip.open(INFO_FILE, "rt") as f:
        df_info = pd.read_csv(f, sep="\t")
        for _, r in df_info.iterrows():
            info_map[r["#string_protein_id"]] = r["preferred_name"]

    # Check panel genes in LCC
    print("\n=== PANEL GENES IN STRING UNIVERSE LCC ===")
    reverse_info = {v: k for k, v in info_map.items()}
    panel_lcc_status = {}
    for g in panel:
        ensp = reverse_info.get(g)
        in_lcc = (ensp in G_lcc) if ensp else False
        deg = G_lcc.degree(ensp) if in_lcc else 0
        panel_lcc_status[g] = in_lcc
        print(f"  - {g:<8} ({ensp}): in_LCC = {in_lcc} | LCC Degree = {deg}")

    # 7. Save edge list to data/processed/string_universe_edges.tsv
    print(f"\n[*] Saving universe LCC edges to {OUT_EDGES}...")
    os.makedirs("data/processed", exist_ok=True)
    out_rows = []
    for u, v, d in G_lcc.edges(data=True):
        out_rows.append({
            "protein1": u,
            "protein2": v,
            "gene_symbol1": info_map.get(u, u),
            "gene_symbol2": info_map.get(v, v),
            "combined_score": d["combined_score"]
        })
    df_out = pd.DataFrame(out_rows)
    df_out.to_csv(OUT_EDGES, sep="\t", index=False)
    print(f"[+] Saved {len(df_out):,} LCC edges to {OUT_EDGES}")

    # Summary statistics
    degrees = [d for _, d in G_lcc.degree()]
    print("\n=== UNIVERSE LCC TOPOLOGY METRICS ===")
    print(f"Total Nodes:     {n_lcc_nodes:,}")
    print(f"Total Edges:     {n_lcc_edges:,}")
    print(f"Mean Degree:     {np.mean(degrees):.2f}")
    print(f"Median Degree:   {np.median(degrees):.1f}")
    print(f"Min / Max Deg:   {np.min(degrees)} / {np.max(degrees)}")
    print(f"Nodes with k>=100: {sum(1 for d in degrees if d >= 100):,} (Degree binning requirement)")

if __name__ == "__main__":
    main()
