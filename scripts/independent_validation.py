#!/usr/bin/env python3
"""
scripts/independent_validation.py
W11.5 Integrated Multi-Stream Validation Pipeline
Combines BioGRID, gnomAD, and PubMed streams into a unified long-format table.
Conforms strictly to W11.5 and V3 vocabulary rules.
"""

import os
import pandas as pd

TABLES_DIR = "results/tables"
BG_FILE = os.path.join(TABLES_DIR, "validation_detail_biogrid.csv")
GN_FILE = os.path.join(TABLES_DIR, "validation_detail_gnomad.csv")
PM_FILE = os.path.join(TABLES_DIR, "validation_detail_pubmed.csv")
OUT_VAL = os.path.join(TABLES_DIR, "validation_results.csv")

def main():
    print("=== W11.5: GENERATING INTEGRATED VALIDATION RESULTS TABLE ===")
    rows = []

    # 1. Stream: BioGRID (Physical PPI Centrality & Proximity)
    if os.path.exists(BG_FILE):
        df_bg = pd.read_csv(BG_FILE)
        for _, r in df_bg.iterrows():
            gene = r["gene"]
            z_bg = r["biogrid_z"]
            z_str = r["string_z"]
            res_str = f"STRING z={z_str:.2f}, BioGRID z={z_bg:.2f}" if pd.notna(z_bg) else "No BioGRID metrics"
            rows.append({
                "evidence_stream": "biogrid",
                "candidate_id": gene,
                "gene_symbol": gene,
                "variant_id": "gene_level",
                "input_source": "BioGRID Physical Interactions",
                "source_version": "5.0.263",
                "retrieval_date": "2026-11-04",
                "comparison_method": "Degree-matched closest proximity (B=1000) & ego-betweenness",
                "result": res_str,
                "limitations": "Physical binding lacks dynamic metabolic links; enzyme interactors depleted.",
                "interpretation": r["notes"],
                "status": r["status"]
            })

    # 2. Stream: gnomAD v4 SAS (Independent Allele Frequency Cross-Check)
    if os.path.exists(GN_FILE):
        df_gn = pd.read_csv(GN_FILE)
        for _, r in df_gn.iterrows():
            v_id = r["variant_id"]
            gene = r["gene_symbol"]
            dtype = r["gnomad_data_type"]
            af_gn = r["gnomad_SAS_af"]
            p_raw = r["primary_cohen_h"]
            res_str = f"Primary Cohen's h={p_raw:.4f}; gnomAD {dtype} SAS AF={af_gn:.5f}" if pd.notna(af_gn) else f"Absent from gnomAD {dtype}"
            rows.append({
                "evidence_stream": "gnomad",
                "candidate_id": v_id,
                "gene_symbol": gene,
                "variant_id": v_id,
                "input_source": f"gnomAD v4.1 ({dtype})",
                "source_version": "v4.1",
                "retrieval_date": "2026-11-04",
                "comparison_method": f"Direction of allele frequency enrichment (SAS vs NFE) via {dtype}",
                "result": res_str,
                "limitations": r["independence_note"],
                "interpretation": r["notes"],
                "status": r["cross_check_status"]
            })

    # 3. Stream: PubMed (Gene-Drug Relevance Literature Assessment)
    if os.path.exists(PM_FILE):
        df_pm = pd.read_csv(PM_FILE)
        for _, r in df_pm.iterrows():
            gene = r["gene_symbol"]
            drug = r["drug"]
            pmid = r["PMID"]
            rows.append({
                "evidence_stream": "pubmed",
                "candidate_id": f"{gene}_{drug}",
                "gene_symbol": gene,
                "variant_id": r["variant_or_mechanism"],
                "input_source": "NCBI PubMed E-Utilities",
                "source_version": "PubMed API 2026",
                "retrieval_date": "2026-11-04",
                "comparison_method": "Systematic query screening (Q1-Q4) across reviews and regional studies",
                "result": f"PMID:{pmid} ({r['publication_year']}): {r['publication_title'][:70]}...",
                "limitations": r["limitations"],
                "interpretation": f"{r['relationship_to_project']}; {r['interpretation']}",
                "status": "consistent" if "supported" in r["interpretation"] else "insufficient_data"
            })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(OUT_VAL, index=False)
    print(f"[+] Saved integrated validation results to {OUT_VAL} ({len(df_out)} rows)")
    
    print("\nSummary by evidence stream and status:")
    print(pd.crosstab(df_out["evidence_stream"], df_out["status"]))

if __name__ == "__main__":
    main()
