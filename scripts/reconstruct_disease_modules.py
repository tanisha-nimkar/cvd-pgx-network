#!/usr/bin/env python3
"""
scripts/reconstruct_disease_modules.py
Task 3: Open Targets Disease Module Reconstruction (Direct vs. Indirect & Threshold Matrix)
Conforms strictly to W9-2 SOP guidelines.
"""

import os
import sys
import time
import requests
import pandas as pd
import numpy as np

OUT_DIR = "results/tables"
UNIVERSE_FILE = "data/processed/string_universe_edges.tsv"
GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
PAGE_SIZE = 500
OT_RELEASE = "26.9.0"

PANEL_GENES = {"CYP2C19", "CYP2C9", "VKORC1", "CYP4F2", "SLCO1B1", "ABCG2"}

DISEASES = [
    {"id": "MONDO_0005010", "name": "coronary artery disorder", "short": "CAD"},
    {"id": "MONDO_0005068", "name": "myocardial infarction", "short": "MI"},
    {"id": "MONDO_0024644", "name": "myocardial ischemia", "short": "IHD"}
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

def fetch_targets(disease_id, indirect_setting):
    page_idx = 0
    retrieved = []
    total_count = None
    while True:
        payload = {"id": disease_id, "i": page_idx, "s": PAGE_SIZE, "indirect": indirect_setting}
        res = requests.post(GRAPHQL_URL, json={"query": GRAPHQL_QUERY, "variables": payload}, timeout=30)
        if not res.ok:
            print(f"[!] Error querying {disease_id} (indirect={indirect_setting}) page {page_idx}: HTTP {res.status_code}")
            sys.exit(1)
        assoc = res.json().get("data", {}).get("disease", {}).get("associatedTargets", {})
        if total_count is None:
            total_count = assoc.get("count", 0)
        rows = assoc.get("rows", [])
        retrieved.extend(rows)
        if len(rows) < PAGE_SIZE or len(retrieved) >= total_count:
            break
        page_idx += 1
        time.sleep(0.15)
    return retrieved

def main():
    print("=== TASK 3: OPEN TARGETS DISEASE MODULE RECONSTRUCTION ===")
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Load universe LCC symbols for mapping coverage
    df_univ = pd.read_csv(UNIVERSE_FILE, sep="\t", usecols=["gene_symbol1", "gene_symbol2"])
    universe_symbols = set(df_univ["gene_symbol1"]).union(set(df_univ["gene_symbol2"]))
    n_universe = len(universe_symbols)
    print(f"[*] Universe LCC contains {n_universe:,} unique gene symbols.")

    all_disease_records = []

    # 2. Fetch Direct and Indirect modules for CAD (MONDO_0005010)
    print("\n[*] Fetching CAD (MONDO_0005010) DIRECT associations (Primary Module)...")
    raw_cad_direct = fetch_targets("MONDO_0005010", indirect_setting=False)
    print(f"[+] Retrieved {len(raw_cad_direct):,} DIRECT CAD target associations.")

    print("\n[*] Fetching CAD (MONDO_0005010) INDIRECT associations (Sensitivity Module)...")
    raw_cad_indirect = fetch_targets("MONDO_0005010", indirect_setting=True)
    print(f"[+] Retrieved {len(raw_cad_indirect):,} INDIRECT CAD target associations.")

    # Parse CAD direct into DataFrame
    direct_rows = []
    for r in raw_cad_direct:
        tgt = r["target"]
        dt_scores = ";".join([f"{dt['id']}:{dt['score']:.4f}" for dt in r.get("datatypeScores", [])])
        sym = tgt["approvedSymbol"]
        direct_rows.append({
            "disease_name": "coronary artery disorder",
            "MONDO_id": "MONDO:0005010",
            "gene_symbol": sym,
            "ensembl_id": tgt["id"],
            "association_score": float(r["score"]),
            "evidence_summary": dt_scores,
            "evidence_type": "direct",
            "in_universe_lcc": sym in universe_symbols,
            "is_panel_gene": sym in PANEL_GENES,
            "ot_release": OT_RELEASE
        })
    df_cad_direct = pd.DataFrame(direct_rows)
    df_cad_direct.to_csv(os.path.join(OUT_DIR, "disease_module_direct_mondo0005010.csv"), index=False)
    print(f"[+] Saved {len(df_cad_direct)} direct CAD targets to {OUT_DIR}/disease_module_direct_mondo0005010.csv")

    # Parse CAD indirect into DataFrame
    indirect_rows = []
    for r in raw_cad_indirect:
        tgt = r["target"]
        dt_scores = ";".join([f"{dt['id']}:{dt['score']:.4f}" for dt in r.get("datatypeScores", [])])
        sym = tgt["approvedSymbol"]
        indirect_rows.append({
            "disease_name": "coronary artery disorder",
            "MONDO_id": "MONDO:0005010",
            "gene_symbol": sym,
            "ensembl_id": tgt["id"],
            "association_score": float(r["score"]),
            "evidence_summary": dt_scores,
            "evidence_type": "indirect_inclusive",
            "in_universe_lcc": sym in universe_symbols,
            "is_panel_gene": sym in PANEL_GENES,
            "ot_release": OT_RELEASE
        })
    df_cad_indirect = pd.DataFrame(indirect_rows)
    df_cad_indirect.to_csv(os.path.join(OUT_DIR, "disease_module_indirect_mondo0005010.csv"), index=False)
    print(f"[+] Saved {len(df_cad_indirect)} indirect CAD targets to {OUT_DIR}/disease_module_indirect_mondo0005010.csv")

    # 3. Construct Threshold Sensitivity Matrix across Cutoffs
    print("\n[*] Constructing threshold sensitivity modules ({unfiltered, >=0.05, >=0.10, >=0.20, >=0.30})...")
    cutoffs = [0.0, 0.05, 0.10, 0.20, 0.30]
    sensitivity_summary = []

    # Evaluate Direct CAD (Primary)
    for cut in cutoffs:
        subset = df_cad_direct[df_cad_direct["association_score"] >= cut]
        in_lcc = subset[subset["in_universe_lcc"]]
        cov_pct = (len(in_lcc) / len(subset) * 100.0) if len(subset) > 0 else 0.0
        panel_overlap = subset[subset["is_panel_gene"]]["gene_symbol"].tolist()

        sensitivity_summary.append({
            "disease": "CAD (MONDO:0005010)",
            "evidence_mode": "DIRECT (Primary)",
            "score_threshold": cut,
            "total_genes_raw": len(subset),
            "genes_in_universe_lcc": len(in_lcc),
            "interactome_coverage_pct": round(cov_pct, 2),
            "panel_genes_in_module": ";".join(sorted(panel_overlap)),
            "n_panel_genes": len(panel_overlap)
        })

    # Evaluate Indirect CAD (Sensitivity)
    for cut in cutoffs:
        subset = df_cad_indirect[df_cad_indirect["association_score"] >= cut]
        in_lcc = subset[subset["in_universe_lcc"]]
        cov_pct = (len(in_lcc) / len(subset) * 100.0) if len(subset) > 0 else 0.0
        panel_overlap = subset[subset["is_panel_gene"]]["gene_symbol"].tolist()

        sensitivity_summary.append({
            "disease": "CAD (MONDO:0005010)",
            "evidence_mode": "INDIRECT (Descendant-Inclusive)",
            "score_threshold": cut,
            "total_genes_raw": len(subset),
            "genes_in_universe_lcc": len(in_lcc),
            "interactome_coverage_pct": round(cov_pct, 2),
            "panel_genes_in_module": ";".join(sorted(panel_overlap)),
            "n_panel_genes": len(panel_overlap)
        })

    df_sens = pd.DataFrame(sensitivity_summary)
    df_sens.to_csv(os.path.join(OUT_DIR, "disease_modules_threshold_sensitivity.csv"), index=False)
    print(f"[+] Saved threshold sensitivity matrix to {OUT_DIR}/disease_modules_threshold_sensitivity.csv")

    print("\n=== CAD DISEASE MODULE THRESHOLD SENSITIVITY GRID ===")
    print(df_sens[["evidence_mode", "score_threshold", "total_genes_raw", "genes_in_universe_lcc", "interactome_coverage_pct", "panel_genes_in_module"]].to_string(index=False))

if __name__ == "__main__":
    main()
