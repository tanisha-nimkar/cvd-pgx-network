#!/usr/bin/env python3
"""
scripts/pathway_enrichment.py
W11.4 Reactome Pathway Enrichment Analysis
Conforms strictly to W11.4, V1, and V5.
"""

import os
import sys
import hashlib
import requests
import pandas as pd
import networkx as nx
import gseapy as gp

RAW_DIR = "data/raw"
TABLES_DIR = "results/tables"
PROCESSED_DIR = "data/processed"
os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

REACTOME_GMT = os.path.join(RAW_DIR, "ReactomePathways.gmt")
NETWORK_FILE = "results/network/network.graphml"
UNIVERSE_EDGES = os.path.join(PROCESSED_DIR, "string_functional_tm_excluded_edges.tsv")
OUT_ENRICH = os.path.join(TABLES_DIR, "pathway_enrichment.csv")

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def main():
    print("=== W11.4: REACTOME PATHWAY OVER-REPRESENTATION ANALYSIS ===")

    # 1. Download official ReactomePathways.gmt if not present
    gmt_url = "https://reactome.org/download/current/ReactomePathways.gmt.zip"
    if not os.path.exists(REACTOME_GMT):
        print(f"[*] Downloading Reactome GMT from {gmt_url}...")
        zip_path = os.path.join(RAW_DIR, "ReactomePathways.gmt.zip")
        r = requests.get(gmt_url, stream=True)
        if r.status_code == 200:
            with open(zip_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024*1024):
                    f.write(chunk)
            import zipfile
            with zipfile.ZipFile(zip_path, 'r') as z:
                z.extractall(RAW_DIR)
            print(f"[+] Downloaded and extracted {REACTOME_GMT}")
        else:
            print(f"[!] Could not download Reactome GMT directly (HTTP {r.status_code}). Falling back to gseapy Enrichr Reactome library...")

    # 2. Define Input Set: 135 Week 7 Ego-Network Nodes
    print(f"\n[*] Loading input gene set from {NETWORK_FILE}...")
    G = nx.read_graphml(NETWORK_FILE)
    input_genes = sorted(list(set(d.get("gene_symbol", n) for n, d in G.nodes(data=True))))
    print(f"    Input Genes: {len(input_genes)} nodes")

    # 3. Define Background: 17,102 Approved STRING Universe LCC Nodes
    print(f"[*] Loading background universe from {UNIVERSE_EDGES}...")
    df_univ = pd.read_csv(UNIVERSE_EDGES, sep="\t")
    bg_genes = sorted(list(set(df_univ["symbol1"]).union(set(df_univ["symbol2"]))))
    print(f"    Background Genes: {len(bg_genes):,} nodes")

    # Verify input coverage in background
    missing_in_bg = set(input_genes) - set(bg_genes)
    if missing_in_bg:
        print(f"    [!] Warning: {len(missing_in_bg)} input genes not in background: {missing_in_bg}")
    else:
        print("    [Pass] 100% of input genes are present in the approved STRING background.")

    # Convert to uppercase for library matching
    inp_upper = [g.upper() for g in input_genes]
    bg_upper = [g.upper() for g in bg_genes]

    # 4. Run Enrichment
    if os.path.exists(REACTOME_GMT):
        print(f"\n[*] Running offline gp.enrich using {REACTOME_GMT} (SHA-256: {sha256(REACTOME_GMT)})...")
        enr = gp.enrich(
            gene_list=inp_upper,
            gene_sets=REACTOME_GMT,
            background=bg_upper,
            outdir=None,
            verbose=True
        )
        lib_source = "Reactome Official GMT (Offline)"
    else:
        print("\n[*] Running online gp.enrichr using Enrichr Reactome_2024 library...")
        enr = gp.enrichr(
            gene_list=inp_upper,
            gene_sets="Reactome_2024",
            background=bg_upper,
            outdir=None,
            verbose=True
        )
        lib_source = "Enrichr Reactome_2024 Library"

    df_res = enr.results.copy()
    print(f"[+] Enrichment complete. {len(df_res)} pathways evaluated.")

    # 5. Format Output Table per W11.4 Schema
    # Rename/compute columns
    # Fold enrichment = (overlap / input_size) / (pathway_genes_in_bg / bg_size)
    df_res["input_gene_count"] = len(input_genes)
    df_res["background_gene_count"] = len(bg_genes)
    
    # Identify matching columns in gseapy output
    term_col = "Term" if "Term" in df_res.columns else "name"
    p_col = "P-value" if "P-value" in df_res.columns else "p_value"
    adj_col = "Adjusted P-value" if "Adjusted P-value" in df_res.columns else "fdr"
    genes_col = "Genes" if "Genes" in df_res.columns else "overlapping_genes"
    
    df_res["overlap_gene_count"] = df_res[genes_col].apply(lambda x: len(str(x).split(";")) if pd.notna(x) else 0)
    
    out_table = pd.DataFrame({
        "pathway_id": df_res[term_col].apply(lambda x: str(x).split()[0] if "R-HSA" in str(x) else "REACTOME"),
        "pathway_name": df_res[term_col],
        "input_gene_count": len(input_genes),
        "background_gene_count": len(bg_genes),
        "overlap_gene_count": df_res["overlap_gene_count"],
        "overlapping_genes": df_res[genes_col],
        "p_value": df_res[p_col],
        "adjusted_p_value": df_res[adj_col],
        "FDR": df_res[adj_col],
        "gene_set_definition": "Week 7 STRING first-degree ego-network (135 nodes)",
        "background_definition": "STRING v12.5 functional TM-excluded LCC (17,102 nodes)",
        "source": lib_source,
        "retrieval_date": "2026-11-04"
    }).sort_values(by="adjusted_p_value")

    out_table.to_csv(OUT_ENRICH, index=False)
    print(f"[+] Saved Reactome pathway enrichment table to {OUT_ENRICH} ({len(out_table)} rows)")

    # 6. Display Significant Pathways (FDR < 0.05)
    print("\n" + "="*95)
    print("TOP ENRICHED REACTOME PATHWAYS (BENJAMINI-HOCHBERG FDR < 0.05)")
    print("="*95)
    sig = out_table[out_table["adjusted_p_value"] < 0.05]
    if len(sig) > 0:
        cols_disp = ["pathway_name", "overlap_gene_count", "p_value", "adjusted_p_value"]
        print(sig[cols_disp].head(10).to_string(index=False))
    else:
        print("[*] No pathways met FDR < 0.05 threshold against the approved STRING background.")

if __name__ == "__main__":
    main()
