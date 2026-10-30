#!/usr/bin/env python3
"""
Step W4.8, W4.9 & W4.10: Week 4 Quality Control & Validation Spot-Check.
Evaluates gene symbol mismatches, extracts a targeted 10-variant audit panel,
and generates summary statistics for documentation.
"""

import pandas as pd

def main():
    path = "results/tables/variants_annotated.csv"
    df = pd.read_csv(path, low_memory=False)
    print(f"[*] Loaded {len(df)} variants from {path}")

    # 1. Audit gene mismatches
    if "vep_gene_mismatch_flag" in df.columns:
        mismatches = df[df["vep_gene_mismatch_flag"] == True]
        print(f"[*] VEP Gene Mismatches flagged: {len(mismatches)}")
        if len(mismatches) > 0:
            print(mismatches[["variant_key", "gene_symbol", "vep_gene_symbol"]].head(10))

    # 2. Select targeted audit variants (Key star alleles + high/moderate consequences)
    star_targets = ["CYP2C19*2", "CYP2C9*2", "SLCO1B1*5", "CYP2C19*3", "CYP4F2*2"]
    
    # Identify star allele rows
    star_rows = pd.DataFrame()
    for col in ["star_allele", "star_id", "haplotype", "name"]:
        if col in df.columns:
            star_rows = df[df[col].astype(str).str.contains("|".join([f"\\{s}" for s in star_targets]), regex=True, na=False)]
            if len(star_rows) > 0:
                break
    
    # High impact rows
    high_rows = df[df["vep_impact"] == "HIGH"]
    mod_rows = df[df["vep_impact"] == "MODERATE"]

    spot_candidates = pd.concat([star_rows, high_rows, mod_rows]).drop_duplicates(subset=["variant_key"])
    spot_10 = spot_candidates.head(10).copy()

    spot_cols = [c for c in ["variant_key", "gene_symbol", "vep_consequence", "vep_impact", 
                             "vep_feature_transcript", "clinvar_sig_category"] if c in spot_10.columns]
    
    spot_out_path = "results/tables/week4_spotcheck.csv"
    spot_10[spot_cols].to_csv(spot_out_path, index=False)
    print(f"[+] Saved 10-variant audit spot-check to {spot_out_path}:")
    print(spot_10[spot_cols].to_string(index=False))

    # 3. Generate QC Summary by Gene
    summary = df.groupby("gene_symbol").agg(
        total_variants=("variant_key", "count"),
        vep_annotated=("vep_impact", lambda x: x.notna().sum()),
        high_impact=("vep_impact", lambda x: (x == "HIGH").sum()),
        moderate_impact=("vep_impact", lambda x: (x == "MODERATE").sum()),
        vus_count=("clinvar_sig_category", lambda x: (x == "VUS").sum()),
        pathogenic_count=("clinvar_sig_category", lambda x: (x == "Pathogenic_Likely_Pathogenic").sum())
    ).reset_index()

    qc_out_path = "results/tables/week4_qc_summary.csv"
    summary.to_csv(qc_out_path, index=False)
    print(f"\n[+] Saved QC summary table to {qc_out_path}:")
    print(summary.to_string(index=False))

if __name__ == "__main__":
    main()
