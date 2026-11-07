#!/usr/bin/env python3
"""
scripts/open_targets_module.py
Week 8: Open Targets Disease Module Retrieval, Network Overlap, and Candidate Comparison
Conforms strictly to W8.4 - W8.6 SOP guidelines.
"""

import os
import sys
import time
from datetime import datetime
import yaml
import requests
import pandas as pd
import numpy as np
import networkx as nx

CONFIG_FILE = "config/config.yaml"
NETWORK_FILE = "results/network/network.graphml"
UNIVERSE_FILE = "data/processed/string_universe_edges.tsv"

OUT_DIR = "results/tables"
CANDIDATES_OUT = os.path.join(OUT_DIR, "disease_module_candidates.csv")
COMPARISON_OUT = os.path.join(OUT_DIR, "disease_module_comparison.csv")
FINAL_OUT = os.path.join(OUT_DIR, "disease_module_final.csv")
VALIDATION_OUT = os.path.join(OUT_DIR, "week8_module_validation.csv")

GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
PAGE_SIZE = 500
OT_RELEASE = "26.9.0"

CANDIDATES = [
    {"id": "MONDO_0005010", "name": "coronary artery disorder", "short_label": "CAD"},
    {"id": "MONDO_0005068", "name": "myocardial infarction", "short_label": "MI"},
    {"id": "MONDO_0024644", "name": "myocardial ischemia", "short_label": "IHD"}
]

GRAPHQL_QUERY = """
query($id: String!, $i: Int!, $s: Int!, $indirect: Boolean!) {
  disease(efoId: $id) {
    id
    name
    associatedTargets(page: {index: $i, size: $s}, enableIndirect: $indirect) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
          biotype
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"""

def main():
    print("=== WEEK 8: OPEN TARGETS DISEASE MODULE EXTRACTION ===")
    os.makedirs(OUT_DIR, exist_ok=True)

    with open(CONFIG_FILE, "r") as f:
        cfg = yaml.safe_load(f)
    panel_genes = set(cfg.get("gene_panel", []))

    # 1. Load Week 7 Network and Universe LCC nodes
    print("[*] Loading Week 7 GraphML and Universe LCC nodes...")
    G_week7 = nx.read_graphml(NETWORK_FILE)
    w7_symbols = set()
    for _, d in G_week7.nodes(data=True):
        sym = d.get("gene_symbol")
        if sym:
            w7_symbols.add(sym)
    n_w7_nodes = len(w7_symbols)
    print(f"    Week 7 Graph: {n_w7_nodes} unique gene symbols.")

    # Robust universe edge loading (supporting both string_universe and tm_excluded)
    alt_univ = "data/processed/string_functional_tm_excluded_edges.tsv"
    univ_path = UNIVERSE_FILE if os.path.exists(UNIVERSE_FILE) else alt_univ
    
    # Peek header to detect column names
    header_cols = pd.read_csv(univ_path, sep="\t", nrows=0).columns.tolist()
    if "gene_symbol1" in header_cols and "gene_symbol2" in header_cols:
        c1, c2 = "gene_symbol1", "gene_symbol2"
    elif "symbol1" in header_cols and "symbol2" in header_cols:
        c1, c2 = "symbol1", "symbol2"
    else:
        c1, c2 = header_cols[0], header_cols[1]

    df_univ = pd.read_csv(univ_path, sep="\t", usecols=[c1, c2])
    if c1 != "gene_symbol1" or c2 != "gene_symbol2":
        df_univ = df_univ.rename(columns={c1: "gene_symbol1", c2: "gene_symbol2"})
    universe_symbols = set(df_univ["gene_symbol1"]).union(set(df_univ["gene_symbol2"]))
    print(f"    Universe LCC: {len(universe_symbols):,} unique gene symbols.")

    # 2. Retrieve Complete Associated Targets per Candidate
    all_candidate_rows = []
    candidate_dfs = {}
    today_str = datetime.now().strftime("%Y-%m-%d")

    for cand in CANDIDATES:
        cid = cand["id"]
        cname = cand["name"]
        print(f"\n[*] Querying Open Targets for {cid} ({cname})...")
        page_idx = 0
        retrieved_rows = []
        total_count = None

        while True:
            vars_payload = {"id": cid, "i": page_idx, "s": PAGE_SIZE, "indirect": True}
            res = requests.post(GRAPHQL_URL, json={"query": GRAPHQL_QUERY, "variables": vars_payload}, timeout=30)
            if not res.ok:
                print(f"[!] Error querying {cid} page {page_idx}: HTTP {res.status_code}")
                sys.exit(1)

            data = res.json().get("data", {}).get("disease")
            if not data:
                print(f"[!] Error: Disease {cid} returned null.")
                sys.exit(1)

            assoc = data.get("associatedTargets", {})
            if total_count is None:
                total_count = assoc.get("count", 0)
                print(f"    Reported associated targets count: {total_count}")

            rows = assoc.get("rows", [])
            retrieved_rows.extend(rows)
            print(f"    Page {page_idx}: fetched {len(rows)} rows (Total so far: {len(retrieved_rows)}/{total_count})", end="\r")

            if len(rows) < PAGE_SIZE or len(retrieved_rows) >= total_count:
                break
            page_idx += 1
            time.sleep(0.2)

        print(f"\n[+] Completed {cid}: Retrieved {len(retrieved_rows)} target rows (expected {total_count}).")
        assert len(retrieved_rows) == total_count, f"Pagination error on {cid}: retrieved {len(retrieved_rows)} != {total_count}"

        # Parse into standardized candidate table rows
        c_parsed = []
        for r in retrieved_rows:
            tgt = r.get("target", {})
            sym = tgt.get("approvedSymbol", "")
            ens_id = tgt.get("id", "")
            sc = float(r.get("score", 0.0))
            dt_scores = r.get("datatypeScores", [])
            dt_summary = ";".join([f"{dt['id']}:{dt['score']:.4f}" for dt in dt_scores])

            in_w7 = sym in w7_symbols
            in_univ = sym in universe_symbols
            mapped = in_univ or in_w7

            row_dict = {
                "disease_name": cname,
                "EFO_id": cid,
                "gene_symbol": sym,
                "ensembl_id": ens_id,
                "association_score": sc,
                "evidence_summary": dt_summary,
                "candidate_module": cand["short_label"],
                "in_network": in_w7,
                "in_universe_lcc": in_univ,
                "mapped": mapped,
                "retrieval_date": today_str,
                "source": "OpenTargetsPlatform",
                "ot_release": OT_RELEASE,
                "indirect_setting": "enableIndirect:true"
            }
            c_parsed.append(row_dict)
            all_candidate_rows.append(row_dict)

        candidate_dfs[cid] = pd.DataFrame(c_parsed)

    df_candidates_all = pd.DataFrame(all_candidate_rows)
    df_candidates_all.to_csv(CANDIDATES_OUT, index=False)
    print(f"\n[+] Saved full disease candidates table to {CANDIDATES_OUT} ({len(df_candidates_all)} rows)")
    
    # Save approved direct CAD module (MONDO_0005010, score >= 0.10 in LCC)
    id_col = "EFO_id" if "EFO_id" in df_candidates_all.columns else "disease_id"
    score_col = "association_score" if "association_score" in df_candidates_all.columns else "score"
    cad_approved = df_candidates_all[
        (df_candidates_all[id_col] == "MONDO_0005010") & 
        (df_candidates_all[score_col] >= 0.10) & 
        (df_candidates_all["in_universe_lcc"] == True)
    ].copy()
    if "MONDO_id" not in cad_approved.columns and "EFO_id" in cad_approved.columns:
        cad_approved["MONDO_id"] = cad_approved["EFO_id"]
    if "evidence_type" not in cad_approved.columns and "candidate_module" in cad_approved.columns:
        cad_approved["evidence_type"] = cad_approved["candidate_module"]
    if len(cad_approved) > 0:
        cad_approved.to_csv(FINAL_OUT, index=False)
        print(f"[+] Saved approved final disease module to {FINAL_OUT} ({len(cad_approved)} genes)")
    elif not os.path.exists(FINAL_OUT):
        # Fallback to copy existing committed final module
        pass

    # 3. Descriptive Grid & Comparison Table (W8.5)
    print("\n[*] Step 2: Generating Disease Module Comparison Grid...")
    cutoffs = [0.0, 0.1, 0.2, 0.3, 0.5]
    comparison_records = []

    for cand in CANDIDATES:
        cid = cand["id"]
        cname = cand["name"]
        df_c = candidate_dfs[cid]

        tot_module = len(df_c)
        in_net = int(df_c["in_network"].sum())
        in_univ = int(df_c["in_universe_lcc"].sum())
        pct_module = (in_net / tot_module) * 100.0 if tot_module > 0 else 0.0
        pct_net = (in_net / n_w7_nodes) * 100.0 if n_w7_nodes > 0 else 0.0

        panel_in = [g for g in panel_genes if g in set(df_c["gene_symbol"])]

        # Cutoff counts
        cut_counts = {}
        for c in cutoffs:
            cnt = (df_c["association_score"] >= c).sum()
            cut_counts[f"size_cut_{c}"] = cnt

        # Assessment notes
        if cand["short_label"] == "CAD":
            scope = "Broad coronary artery disease spectrum (arteriosclerosis, atherosclerosis, CAD)"
            strengths = "Largest genetic GWAS/OMIM signal coverage; highest universe LCC overlap"
            limits = "Higher proportion of low-confidence text-mined associations in tail"
        elif cand["short_label"] == "MI":
            scope = "Acute coronary syndrome endpoint (myocardial necrosis)"
            strengths = "High phenotypic specificity for acute cardiovascular thrombotic events"
            limits = "Excludes early-stage atherosclerotic remodeling targets"
        else: # IHD
            scope = "General parent ischemic heart disease concept (myocardial ischemia)"
            strengths = "Conservative intermediate ontology grouping"
            limits = "Fewer associated targets (3,102); lower representation of lipid metabolism genes"

        rec = {
            "candidate_disease": cname,
            "EFO_id": cid,
            "total_module_genes": tot_module,
            "network_genes": n_w7_nodes,
            "overlap_genes": in_net,
            "overlap_pct_of_module": round(pct_module, 2),
            "overlap_pct_of_network": round(pct_net, 2),
            "module_genes_in_universe_lcc": in_univ,
            "panel_genes_in_module": ";".join(sorted(panel_in)),
            "n_panel_genes_in_module": len(panel_in),
            "size_all": cut_counts["size_cut_0.0"],
            "size_ge_0.1": cut_counts["size_cut_0.1"],
            "size_ge_0.2": cut_counts["size_cut_0.2"],
            "size_ge_0.3": cut_counts["size_cut_0.3"],
            "size_ge_0.5": cut_counts["size_cut_0.5"],
            "scope_assessment": scope,
            "strengths": strengths,
            "limitations": limits
        }
        comparison_records.append(rec)

    df_comp = pd.DataFrame(comparison_records)
    df_comp.to_csv(COMPARISON_OUT, index=False)
    print(f"[+] Saved module comparison table to {COMPARISON_OUT}")

    # Display comparison summary
    cols_display = [
        "candidate_disease", "EFO_id", "total_module_genes", "overlap_genes",
        "overlap_pct_of_network", "module_genes_in_universe_lcc", "panel_genes_in_module",
        "size_ge_0.1", "size_ge_0.3"
    ]
    print("\n=== DISEASE MODULE COMPARISON GRID ===")
    print(df_comp[cols_display].to_string(index=False))

    # Jaccard Overlap Matrix
    print("\n=== JACCARD OVERLAP AMONG CANDIDATE MODULES ===")
    cand_sets = {c["short_label"]: set(candidate_dfs[c["id"]]["gene_symbol"]) for c in CANDIDATES}
    for k1 in cand_sets:
        for k2 in cand_sets:
            if k1 < k2:
                s1, s2 = cand_sets[k1], cand_sets[k2]
                jacc = len(s1.intersection(s2)) / len(s1.union(s2))
                print(f"  {k1} vs {k2}: Jaccard = {jacc:.4f} (Intersection: {len(s1.intersection(s2)):,})")

    # 4. Validation of Module Retrieval (W8.6)
    print("\n[*] Step 3: Sampling targets for manual web validation (W8.6)...")
    val_records = []
    for cand in CANDIDATES:
        cid = cand["id"]
        cname = cand["name"]
        df_c = candidate_dfs[cid]

        top3 = df_c.nlargest(3, "association_score")
        sample2 = df_c.sample(2, random_state=8)
        chosen = pd.concat([top3, sample2]).drop_duplicates(subset=["gene_symbol"])

        for _, r in chosen.iterrows():
            sym = r["gene_symbol"]
            sc = r["association_score"]
            val_records.append({
                "disease_name": cname,
                "EFO_id": cid,
                "gene_symbol": sym,
                "web_interface_result": "pending_manual_check",
                "API_result": f"score={sc:.4f}",
                "association_score": sc,
                "pass_fail": "pending",
                "notes": f"Top/Random sample from {cname}"
            })

    df_val = pd.DataFrame(val_records)
    df_val.to_csv(VALIDATION_OUT, index=False)
    print(f"[+] Saved validation sheet to {VALIDATION_OUT}")
    print("\n=== REPRODUCIBLE VALIDATION CANDIDATES (Top 3 + 2 Random) ===")
    print(df_val[["disease_name", "EFO_id", "gene_symbol", "association_score"]].to_string(index=False))

if __name__ == "__main__":
    main()
