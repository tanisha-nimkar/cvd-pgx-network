#!/usr/bin/env python3
"""
scripts/build_network.py
Week 7: First-Degree STRING PPI Network Construction & Audit Pipeline
Conforms strictly to W7.3 - W7.13 SOP guidelines.
"""

import os
import sys
import time
import hashlib
from datetime import datetime
import yaml
import requests
import pandas as pd
import numpy as np
import networkx as nx

CONFIG_FILE = "config/config.yaml"
OUT_DIR = "results/network"
MAP_FILE = os.path.join(OUT_DIR, "string_id_mapping.tsv")
AUDIT_FILE = os.path.join(OUT_DIR, "string_edges_audit.tsv")
GRAPHML_FILE = os.path.join(OUT_DIR, "network.graphml")

BASE_URL = "https://string-db.org/api"
CALLER = "cvd-pgx-network"
LIMIT = 1000
STRING_VERSION = "12.5"
NON_TEXT_CHANNELS = ["nscore", "fscore", "pscore", "ascore", "escore", "dscore"]

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def main():
    print("=== WEEK 7: STRING PPI NETWORK PIPELINE ===")
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Load config and panel
    with open(CONFIG_FILE, "r") as f:
        cfg = yaml.safe_load(f)

    panel_genes = cfg.get("gene_panel", [])
    if not panel_genes:
        print("[!] STOP: gene_panel in config.yaml is empty.")
        sys.exit(1)

    species = int(cfg.get("string_species", 9606))
    network_type = str(cfg.get("string_network_type", "functional"))
    conf_thresh_raw = float(cfg.get("string_confidence_threshold", 700))
    thresh = conf_thresh_raw / 1000.0 if conf_thresh_raw > 1.0 else conf_thresh_raw

    print(f"[*] Panel Seeds ({len(panel_genes)}): {panel_genes}")
    print(f"[*] Species: {species} | Network Type: {network_type} | Confidence Threshold: {thresh} ({conf_thresh_raw})")

    # 2. Gene Symbol -> STRING ID Mapping (W7.3)
    print("\n[*] Step 1: Resolving STRING protein identifiers...")
    mapping_rows = []
    seed_id_map = {} # gene_symbol -> string_id

    for sym in panel_genes:
        url = f"{BASE_URL}/tsv/get_string_ids"
        params = {
            "identifiers": sym,
            "species": species,
            "limit": 5,
            "caller_identity": CALLER
        }
        res = requests.get(url, params=params, timeout=15)
        if not res.ok:
            print(f"[!] Error querying get_string_ids for {sym}: HTTP {res.status_code}")
            mapping_rows.append({
                "input_symbol": sym, "string_id": "", "preferred_name": "",
                "n_candidates": 0, "mapping_status": "api_query_failure"
            })
            continue

        lines = res.text.strip().splitlines()
        if len(lines) <= 1:
            mapping_rows.append({
                "input_symbol": sym, "string_id": "", "preferred_name": "",
                "n_candidates": 0, "mapping_status": "unmapped"
            })
            continue

        header = lines[0].split("\t")
        candidates = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
        n_cand = len(candidates)

        # Exact match rule
        top_cand = candidates[0]
        pref_name = top_cand.get("preferredName", "")
        str_id = top_cand.get("stringId", "")
        taxon = str(top_cand.get("ncbiTaxonId", ""))

        if pref_name.upper() == sym.upper() and taxon == str(species):
            status = "exact_match"
            seed_id_map[sym] = str_id
        else:
            status = "ambiguous_mapping"

        mapping_rows.append({
            "input_symbol": sym,
            "string_id": str_id,
            "preferred_name": pref_name,
            "n_candidates": n_cand,
            "mapping_status": status
        })
        time.sleep(1.0) # rate limit guard

    df_map = pd.DataFrame(mapping_rows)
    df_map.to_csv(MAP_FILE, sep="\t", index=False)
    print(f"[+] Saved mapping audit to {MAP_FILE}")
    print(df_map.to_string(index=False))

    if len(seed_id_map) != len(panel_genes):
        print(f"[!] STOP: Only {len(seed_id_map)}/{len(panel_genes)} seeds resolved unambiguously.")
        sys.exit(1)

    # 3. Query First-Degree Interactors for Seeds (W7.4)
    print("\n[*] Step 2: Querying interaction partners (seeds only, no recursion)...")
    all_edges_raw = []

    for sym, str_id in seed_id_map.items():
        print(f"  -> Fetching partners for {sym} ({str_id})...")
        url = f"{BASE_URL}/tsv/interaction_partners"
        params = {
            "identifiers": str_id,
            "species": species,
            "limit": LIMIT,
            "required_score": 0,
            "network_type": network_type,
            "caller_identity": CALLER
        }
        res = requests.get(url, params=params, timeout=25)
        if not res.ok:
            print(f"[!] API error for {sym}: HTTP {res.status_code}")
            continue

        lines = res.text.strip().splitlines()
        if len(lines) <= 1:
            print(f"     No partners returned for {sym}.")
            continue

        header = lines[0].split("\t")
        rows = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
        print(f"     Retrieved {len(rows)} raw edge records.")

        # Truncation check
        if len(rows) == LIMIT:
            scores = [float(r["score"]) for r in rows if "score" in r]
            if scores and min(scores) >= thresh:
                print(f"[!] STOP: Truncation detected for {sym}! Hit LIMIT={LIMIT} with min score >= {thresh}")
                sys.exit(1)

        all_edges_raw.extend(rows)
        time.sleep(1.0) # rate limit guard

    df_edges = pd.DataFrame(all_edges_raw)
    print(f"[+] Total raw edge records retrieved across seeds: {len(df_edges)}")

    # 4. Classification & Audit (W7.5)
    print("\n[*] Step 3: Classifying edges against threshold and text-mining filters...")
    score_cols = ["score", "nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore"]
    for sc in score_cols:
        df_edges[sc] = pd.to_numeric(df_edges[sc], errors="coerce").fillna(0.0)

    df_edges["above_threshold"] = df_edges["score"].round(3) >= thresh
    df_edges["max_non_text_score"] = df_edges[NON_TEXT_CHANNELS].max(axis=1)
    df_edges["text_mining_only"] = (df_edges["tscore"] > 0) & (df_edges["max_non_text_score"] == 0.0)

    # Edge status assignment
    df_edges["edge_status"] = "retained"
    df_edges.loc[df_edges["text_mining_only"], "edge_status"] = "excluded_text_mining_only"
    df_edges.loc[~df_edges["above_threshold"], "edge_status"] = "excluded_below_threshold"

    df_edges.to_csv(AUDIT_FILE, sep="\t", index=False)
    print(f"[+] Saved complete edge audit table to {AUDIT_FILE}")

    print("\n--- Edge Classification Summary ---")
    print(df_edges["edge_status"].value_counts())

    # 5. Build NetworkX Graph (W7.6)
    print("\n[*] Step 4: Constructing NetworkX Graph...")
    G = nx.Graph()

    # Add all panel seeds as baseline nodes
    for sym, str_id in seed_id_map.items():
        G.add_node(
            str_id,
            gene_symbol=sym,
            string_id=str_id,
            node_role="panel_gene",
            seed_genes=sym,
            isolated=True
        )

    # Process retained edges
    retained_df = df_edges[df_edges["edge_status"] == "retained"].copy()
    seen_edges = set()

    # Map string_id -> preferred name
    name_map = {}
    for _, r in df_edges.iterrows():
        name_map[r["stringId_A"]] = r["preferredName_A"]
        name_map[r["stringId_B"]] = r["preferredName_B"]

    for _, row in retained_df.iterrows():
        id_a = row["stringId_A"]
        id_b = row["stringId_B"]
        if id_a == id_b:
            continue # no self-loops

        edge_key = tuple(sorted([id_a, id_b]))
        if edge_key in seen_edges:
            continue
        seen_edges.add(edge_key)

        # Add partner node if not in G
        for nid in [id_a, id_b]:
            if not G.has_node(nid):
                sym_partner = name_map.get(nid, nid)
                G.add_node(
                    nid,
                    gene_symbol=sym_partner,
                    string_id=nid,
                    node_role="first_degree_interactor",
                    seed_genes="",
                    isolated=False
                )

        # Add edge with scalar GraphML attributes
        G.add_edge(
            id_a, id_b,
            combined_score=float(row["score"]),
            nscore=float(row["nscore"]),
            fscore=float(row["fscore"]),
            pscore=float(row["pscore"]),
            ascore=float(row["ascore"]),
            escore=float(row["escore"]),
            dscore=float(row["dscore"]),
            tscore=float(row["tscore"]),
            text_mining_contribution=bool(row["tscore"] > 0),
            max_non_text_score=float(row["max_non_text_score"]),
            source="STRING",
            string_version=STRING_VERSION
        )

    # Update node attributes: isolated & seed_genes tracking
    seed_str_ids = set(seed_id_map.values())
    for n in G.nodes():
        deg = G.degree(n)
        G.nodes[n]["isolated"] = bool(deg == 0)
        if n not in seed_str_ids:
            # Record which seeds this interactor is adjacent to
            neighbors = set(G.neighbors(n))
            adj_seeds = [seed_id_map[sym] for sym in panel_genes if seed_id_map[sym] in neighbors]
            adj_symbols = [name_map.get(sid, sid) for sid in adj_seeds]
            G.nodes[n]["seed_genes"] = ";".join(sorted(adj_symbols))

    # Graph-level attributes
    G.graph["string_version"] = STRING_VERSION
    G.graph["access_date"] = datetime.now().strftime("%Y-%m-%d")
    G.graph["threshold"] = float(thresh)
    G.graph["species"] = str(species)
    G.graph["network_type"] = network_type

    # Write and read back GraphML to verify integrity
    nx.write_graphml(G, GRAPHML_FILE)
    G_read = nx.read_graphml(GRAPHML_FILE)

    assert G.number_of_nodes() == G_read.number_of_nodes(), "GraphML node count mismatch!"
    assert G.number_of_edges() == G_read.number_of_edges(), "GraphML edge count mismatch!"

    print(f"[+] Saved GraphML to {GRAPHML_FILE}")
    print(f"    Graph contains {G.number_of_nodes()} nodes and {G.number_of_edges()} retained edges.")

    # 6. Network QC Audit (W7.8)
    print("\n=== NETWORK QUALITY CONTROL AUDIT ===")
    all_seeds_in = all(G.has_node(sid) for sid in seed_str_ids)
    print(f"1. All panel seeds present in graph:        {'PASS' if all_seeds_in else 'FAIL'}")

    first_degree_valid = True
    for u, v in G.edges():
        if (u not in seed_str_ids) and (v not in seed_str_ids):
            first_degree_valid = False
            break
    print(f"2. Strict first-degree restriction verified:  {'PASS' if first_degree_valid else 'FAIL'}")

    min_retained_score = min(d["combined_score"] for _, _, d in G.edges(data=True)) if G.edges else 0.0
    print(f"3. Combined score >= {thresh}:                 {'PASS' if min_retained_score >= thresh else 'FAIL'} (min: {min_retained_score:.4f})")

    tm_only_in_graph = any(d["max_non_text_score"] == 0.0 and d["tscore"] > 0 for _, _, d in G.edges(data=True))
    print(f"4. Zero text-mining-only edges in graph:     {'PASS' if not tm_only_in_graph else 'FAIL'}")

    # Isolated genes report (W7.9)
    print("\n=== PANEL SEED ISOLATION & DEGREE REPORT ===")
    for sym in panel_genes:
        sid = seed_id_map[sym]
        deg = G.degree(sid)
        raw_edges_sym = df_edges[df_edges["stringId_A"] == sid]
        below_count = (raw_edges_sym["edge_status"] == "excluded_below_threshold").sum()
        tm_count = (raw_edges_sym["edge_status"] == "excluded_text_mining_only").sum()
        ret_count = (raw_edges_sym["edge_status"] == "retained").sum()
        print(f"  - {sym:<8} (Degree: {deg:<3}) | Retained: {ret_count:<3} | Below Thresh: {below_count:<3} | Text-mining-only: {tm_count:<3} | Isolated: {G.nodes[sid]['isolated']}")

    # 7. Reproducible Validation Edges (W7.10)
    print("\n=== MANUAL REPRODUCIBLE VALIDATION (W7.10) ===")
    retained_for_sampling = retained_df.rename(columns={"preferredName_A": "seed", "preferredName_B": "partner"})
    if len(retained_for_sampling) >= 2:
        chosen = retained_for_sampling.sort_values(["seed", "partner"]).sample(n=2, random_state=7)
        cols_show = ["seed", "partner", "score", "escore", "dscore", "tscore", "max_non_text_score"]
        print(chosen[cols_show].to_string(index=False))
    else:
        print("Fewer than 2 retained edges available for sampling.")

    # Checksum of GraphML
    graphml_sha = compute_sha256(GRAPHML_FILE)
    print(f"\n[+] network.graphml SHA-256: {graphml_sha}")

if __name__ == "__main__":
    main()
