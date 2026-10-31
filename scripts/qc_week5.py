#!/usr/bin/env python3
"""
scripts/qc_week5.py
Automated Quality Control & Integrity Audit for Week 5 Deliverables:
- Invariant maintenance (N=465)
- Allele frequency probability boundary checks [0.0, 1.0]
- Delta AF boundary checks [-1.0, 1.0]
- Cross-gene population summary matrix
"""

import sys
import pandas as pd
import numpy as np

def main():
    file_path = "results/tables/variants_annotated_w5.csv"
    qc_out = "results/tables/week5_qc_summary.csv"

    print(f"[*] Loading Week 5 annotated catalog: {file_path}")
    df = pd.read_csv(file_path, low_memory=False)

    print("\n=== 5-POINT QUALITY CONTROL AUDIT: WEEK 5 ===")

    # 1. Exact Row Invariant
    n_rows = len(df)
    c1 = (n_rows == 465)
    print(f"1. Row Count Invariant (N=465):        {'PASS' if c1 else 'FAIL'} ({n_rows}/465)")

    # 2. Probability Bounds on Allele Frequencies [0, 1]
    af_cols = [c for c in df.columns if c.startswith("af_")]
    invalid_af = 0
    for col in af_cols:
        vals = df[col].dropna()
        out_of_bounds = ((vals < 0.0) | (vals > 1.0)).sum()
        invalid_af += out_of_bounds

    c2 = (invalid_af == 0)
    print(f"2. Allele Frequency Bounds [0.0, 1.0]:  {'PASS' if c2 else 'FAIL'} (Invalid: {invalid_af})")

    # 3. Delta AF Bounds [-1.0, 1.0]
    valid_deltas = df["delta_af_sas"].dropna()
    invalid_deltas = ((valid_deltas < -1.0) | (valid_deltas > 1.0)).sum()
    c3 = (invalid_deltas == 0)
    print(f"3. Delta AF Bounds [-1.0, 1.0]:        {'PASS' if c3 else 'FAIL'} (Invalid: {invalid_deltas})")

    # 4. Canonical Pharmacogenomic Markers Check
    test_markers = ["rs4149056", "rs4244285", "rs1057910", "rs9923231"]
    matched_markers = df[df["rsid"].isin(test_markers) & df["af_sas"].notna()]["rsid"].tolist()
    c4 = len(matched_markers) > 0
    print(f"4. Core Pharmacogene Loci Verified:    {'PASS' if c4 else 'FAIL'} ({len(matched_markers)}/{len(test_markers)} anchors present)")

    # 5. Non-null Enrichment Ratio when SAS AF is present
    has_sas_null_ratio = df[df["af_sas"].notna()]["sas_enrichment_ratio"].isna().sum()
    c5 = (has_sas_null_ratio == 0)
    print(f"5. Enrichment Ratio Consistency:       {'PASS' if c5 else 'FAIL'} (Missing: {has_sas_null_ratio})")

    # Per-gene Summary Table
    summary = df.groupby("gene_symbol").agg(
        total_variants=("variant_key", "count"),
        variants_with_global_af=("af_global", lambda x: x.notna().sum()),
        variants_with_sas_af=("af_sas", lambda x: x.notna().sum()),
        common_in_sas=("rarity_tier_sas", lambda x: (x == "Common").sum()),
        rare_in_sas=("rarity_tier_sas", lambda x: (x.isin(["Rare", "Ultra_Rare"])).sum()),
        mean_delta_af_sas=("delta_af_sas", lambda x: x.dropna().mean())
    ).reset_index()

    summary.to_csv(qc_out, index=False)
    print(f"\n[+] QC Summary table saved to: {qc_out}\n")
    print(summary.to_string(index=False))

    all_passed = all([c1, c2, c3, c4, c5])
    print(f"\nOVERALL AUDIT VERDICT: {'PASSED (Ready for Week 5 Commit)' if all_passed else 'FAILED'}")

if __name__ == "__main__":
    main()
