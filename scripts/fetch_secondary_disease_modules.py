#!/usr/bin/env python3
"""
scripts/fetch_secondary_disease_modules.py
Task 3 Deliverable: Fetch MI (MONDO:0005068) and IHD (MONDO:0024644) modules
Conforms strictly to W9-2 and Part IV deliverables.
"""

import os
import sys
import time
import requests
import pandas as pd

OUT_DIR = "results/tables"
UNIVERSE_FILE = "data/processed/string_universe_edges.tsv"
GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
PAGE_SIZE = 500
OT_RELEASE = "26.9.0"
PANEL_GENES = {"CYP2C19", "CYP2C9", "VKORC1", "CYP4F2", "SLCO1B1", "ABCG2"}

SECONDARY_DISEASES = [
    {"id": "MONDO_0005068", "name": "myocardial infarction", "label": "mi"},
    {"id": "MONDO_0024644", "name": "myocardial ischemia", "label": "ihd"}
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

def fetch_targets(disease_id, indirect):
    page_idx = 0
    retrieved = []
    total_count = None
    while True:
        payload = {"id": disease_id, "i": page_idx, "s": PAGE_SIZE, "indirect": indirect}
        res = requests.post(GRAPHQL_URL, json={"query": GRAPHQL_QUERY, "variables": payload}, timeout=30)
        if not res.ok:
            print(f"[!] HTTP error {res.status_code} on {disease_id} (indirect={indirect})")
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
    print("=== TASK 3: FETCHING SECONDARY DISEASE MODULES (MI & IHD) ===")
    df_univ = pd.read_csv(UNIVERSE_FILE, sep="\t", usecols=["gene_symbol1", "gene_symbol2"])
    universe_symbols = set(df_univ["gene_symbol1"]).union(set(df_univ["gene_symbol2"]))

    for dis in SECONDARY_DISEASES:
        did = dis["id"]
        dname = dis["name"]
        dlabel = dis["label"]

        for indirect, mode in [(False, "direct"), (True, "indirect")]:
            print(f"[*] Fetching {dname} ({did}) {mode} targets...")
            raw_rows = fetch_targets(did, indirect=indirect)
            parsed = []
            for r in raw_rows:
                tgt = r["target"]
                sym = tgt["approvedSymbol"]
                dt_scores = ";".join([f"{dt['id']}:{dt['score']:.4f}" for dt in r.get("datatypeScores", [])])
                parsed.append({
                    "disease_name": dname,
                    "MONDO_id": did.replace("_", ":"),
                    "gene_symbol": sym,
                    "ensembl_id": tgt["id"],
                    "association_score": float(r["score"]),
                    "evidence_summary": dt_scores,
                    "evidence_type": mode,
                    "in_universe_lcc": sym in universe_symbols,
                    "is_panel_gene": sym in PANEL_GENES,
                    "ot_release": OT_RELEASE
                })
            df_out = pd.DataFrame(parsed)
            out_file = os.path.join(OUT_DIR, f"disease_module_{dlabel}_{mode}.csv")
            df_out.to_csv(out_file, index=False)
            print(f"[+] Saved {len(df_out)} {mode} targets to {out_file}")

if __name__ == "__main__":
    main()
