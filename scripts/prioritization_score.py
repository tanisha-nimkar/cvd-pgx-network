#!/usr/bin/env python3
"""
scripts/prioritization_score.py
Week 10 Composite Variant Prioritization Score + Sensitivity Analysis
Conforms strictly to W10.1 - W10.10.
"""

import os
import sys
import json
import itertools
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, kendalltau
from sklearn.preprocessing import MinMaxScaler

TABLES_DIR = "results/tables"
CONFIG_FILE = "config/config.yaml"

VARIANTS_FILE = os.path.join(TABLES_DIR, "variants_annotated.csv")
POP_FREQ_FILE = os.path.join(TABLES_DIR, "population_frequencies.csv")
METRICS_FILE = os.path.join(TABLES_DIR, "network_metrics.csv")
CPIC_FILE = os.path.join(TABLES_DIR, "annotated_cpic_pairs.tsv")

OUT_RANKED = os.path.join(TABLES_DIR, "prioritization_ranked.csv")
OUT_SENS_LONG = os.path.join(TABLES_DIR, "sensitivity_analysis.csv")
OUT_SENS_SUMM = os.path.join(TABLES_DIR, "sensitivity_summary.csv")

ORDER = ["population_enrichment", "functional_consequence", "network_centrality", "disease_proximity"]
INPUTS = ["P_n", "F_n", "C_n", "D_n"]
W0 = {
    "population_enrichment": 0.25,
    "functional_consequence": 0.25,
    "network_centrality": 0.25,
    "disease_proximity": 0.25
}

def cohens_h(p_target, p_nfe):
    """Signed Cohen's h = 2*asin(sqrt(p_target)) - 2*asin(sqrt(p_nfe))."""
    if pd.isna(p_target) or pd.isna(p_nfe):
        return np.nan
    p_t = np.clip(p_target, 0.0, 1.0)
    p_n = np.clip(p_nfe, 0.0, 1.0)
    return float(2.0 * np.arcsin(np.sqrt(p_t)) - 2.0 * np.arcsin(np.sqrt(p_n)))

def map_vep_impact(impact_str):
    """Ordinal VEP IMPACT: MODIFIER=0, LOW=1, MODERATE=2, HIGH=3."""
    if pd.isna(impact_str):
        return np.nan
    mapping = {"MODIFIER": 0.0, "LOW": 1.0, "MODERATE": 2.0, "HIGH": 3.0}
    return mapping.get(str(impact_str).upper(), np.nan)

def composite(df_in, w):
    """The ONLY scoring function: weighted dot product."""
    return df_in[INPUTS].to_numpy() @ np.array([w[k] for k in ORDER])

def main():
    print("=== WEEK 10: COMPOSITE VARIANT PRIORITIZATION PIPELINE ===")

    # 1. Load inputs
    df_var = pd.read_csv(VARIANTS_FILE)
    df_pop = pd.read_csv(POP_FREQ_FILE)
    df_met = pd.read_csv(METRICS_FILE)

    # Filter df_pop to primary India-source cohort (GIH)
    target_pops = df_pop["target_population"].unique()
    primary_pop = "GIH" if "GIH" in target_pops else target_pops[0]
    print(f"[*] Using primary India-source population cohort: {primary_pop}")
    df_pop_target = df_pop[df_pop["target_population"] == primary_pop].drop_duplicates(subset=["variant_key"]).copy()

    # Merge Population Frequencies into Variant Candidates
    merged = df_var.merge(df_pop_target[["variant_key", "target_af", "nfe_af"]], on="variant_key", how="left")

    # Merge Network Metrics (Gene-Level)
    df_met_panel = df_met[df_met["is_panel_gene"]].copy()
    gene_col = "gene_symbol"
    merged = merged.merge(df_met_panel[[gene_col, "degree", "betweenness_centrality", "proximity_z_score", "proximity_status"]],
                          left_on="gene_symbol", right_on=gene_col, how="left")

    # Load Guideline annotations for descriptive baseline & labeling
    if os.path.exists(CPIC_FILE):
        df_cpic = pd.read_csv(CPIC_FILE, sep="\t")
        gene_to_drugs = df_cpic.groupby("genesymbol")["drugname"].apply(lambda x: ";".join(sorted(set(x)))).to_dict()
        merged["gene_drug_pair"] = merged["gene_symbol"].map(gene_to_drugs).fillna("None")
        gene_to_cpic = df_cpic.groupby("genesymbol")["cpiclevel"].apply(lambda x: ";".join(sorted(set(x)))).to_dict()
        merged["guideline_covered"] = merged["gene_symbol"].map(lambda g: g in gene_to_cpic)
    else:
        merged["gene_drug_pair"] = "None"
        merged["guideline_covered"] = False

    # 2. Compute Raw Components
    print("\n[*] Calculating raw score components (P_raw, F_raw, C_raw, D_raw)...")
    merged["P_raw"] = [cohens_h(t, n) for t, n in zip(merged["target_af"], merged["nfe_af"])]
    merged["F_raw"] = merged["vep_impact"].map(map_vep_impact)
    
    # Gene-level raw centrality C_raw: unweighted degree & betweenness on 135-node graph
    deg_norm = merged["degree"] / 134.0
    merged["C_raw"] = (deg_norm + merged["betweenness_centrality"]) / 2.0
    
    # Proximity D_raw = -z (closer = higher score)
    merged["D_raw"] = -merged["proximity_z_score"]

    # 3. Quality Control on Components (W10.4)
    print("\nPre-normalization checks:")
    for c in ["P_raw", "F_raw", "C_raw", "D_raw"]:
        desc = merged[c].describe().to_dict()
        print(f"  {c:<8}: min={desc.get('min', np.nan):.4f}, max={desc.get('max', np.nan):.4f} | NaN: {merged[c].isna().sum()} | nunique: {merged[c].nunique()}")
        if merged[c].nunique() <= 1:
            raise ValueError(f"Constant column detected for {c}! STOP.")

    # 4. Normalization with MinMaxScaler (W10.4)
    ok_var = merged["P_raw"].notna() & merged["F_raw"].notna()
    sc_var = MinMaxScaler().fit(merged.loc[ok_var, ["P_raw", "F_raw"]])
    
    genes_unique = merged[["gene_symbol", "C_raw", "D_raw"]].drop_duplicates(subset=["gene_symbol"]).dropna()
    sc_gene = MinMaxScaler().fit(genes_unique[["C_raw", "D_raw"]])

    fit_log = {
        "P_min_max": (float(sc_var.data_min_[0]), float(sc_var.data_max_[0])),
        "F_min_max": (float(sc_var.data_min_[1]), float(sc_var.data_max_[1])),
        "C_min_max": (float(sc_gene.data_min_[0]), float(sc_gene.data_max_[0])),
        "D_min_max": (float(sc_gene.data_min_[1]), float(sc_gene.data_max_[1]))
    }
    print(f"\n[+] Scaler fit log: {fit_log}")

    # Transform
    merged[["P_n", "F_n"]] = sc_var.transform(merged[["P_raw", "F_raw"]])
    merged[["C_n", "D_n"]] = sc_gene.transform(merged[["C_raw", "D_raw"]])

    # 5. Missing-Value Policy & Candidate Status (S6)
    complete_case = merged[INPUTS].notna().all(axis=1)
    merged["candidate_status"] = "ranked"
    merged.loc[~complete_case, "candidate_status"] = "unranked_missing_components"
    print(f"\n[*] Candidate completeness: {complete_case.sum()} ranked, {(~complete_case).sum()} unranked")

    # 6. Composite Score & Baseline Rank (W10.5)
    print("\n[*] Calculating composite scores with equal weights (0.25 x 4)...")
    assert abs(sum(W0.values()) - 1.0) < 1e-12 and all(abs(W0[k] - 0.25) < 1e-12 for k in ORDER)
    
    merged["composite_score"] = np.nan
    merged["rank"] = np.nan
    
    scores = np.round(composite(merged.loc[complete_case], W0), 12)
    merged.loc[complete_case, "composite_score"] = scores
    merged.loc[complete_case, "rank"] = rankdata(-scores, method="min")

    # 7. Mandatory Guideline Exclusion Test (W10.6)
    print("\n[*] Executing W10.6 Guideline Exclusion Assertions...")
    assert not any("guideline" in c.lower() for c in INPUTS)
    sub_ranked = merged.loc[complete_case].copy()
    base_scores = composite(sub_ranked, W0)
    
    # Test flipping guideline_covered to all False, all True, and inverted
    for flip_mask in [~sub_ranked["guideline_covered"], pd.Series(True, index=sub_ranked.index), pd.Series(False, index=sub_ranked.index)]:
        mut = sub_ranked.copy()
        mut["guideline_covered"] = flip_mask
        assert np.array_equal(composite(mut, W0), base_scores), "Guideline leakage into score calculation!"
        assert np.array_equal(rankdata(-composite(mut, W0), method="min"), sub_ranked["rank"].to_numpy()), "Guideline leakage into rank!"
    print("    [PASS] All W10.6 guideline mutation assertions passed successfully.")

    # 8. Output results/tables/prioritization_ranked.csv (W10.2, W10.10)
    merged["candidate_id"] = merged["variant_key"]
    merged["gene_level_component"] = True
    for k in ORDER:
        merged[f"weight_{k}"] = W0[k]

    out_cols = [
        "candidate_id", "gene_symbol", "variant_key", "rsid", "gene_drug_pair",
        "guideline_covered", "P_raw", "P_n", "F_raw", "F_n", "C_raw", "C_n",
        "D_raw", "D_n", "weight_population_enrichment", "weight_functional_consequence",
        "weight_network_centrality", "weight_disease_proximity", "composite_score",
        "rank", "candidate_status", "annotation_conflict_flag", "vep_consequence", "clinvar_clnsig", "gene_level_component"
    ]
    merged[out_cols].to_csv(OUT_RANKED, index=False)
    print(f"[+] Saved ranked candidates to {OUT_RANKED} ({len(merged)} rows)")

    # 9. Systematic Sensitivity Analysis (+-50% Relative Grid; W10.7)
    print("\n[*] Running full factorial sensitivity grid (+-50% relative; 5 levels)...")
    levels = [-0.50, -0.25, 0.0, 0.25, 0.50]
    vecs = {}
    for rel in itertools.product(levels, repeat=4):
        w = np.array([0.25 * (1.0 + r) for r in rel])
        w = w / w.sum()
        vecs[tuple(np.round(w, 9))] = rel

    scen = [(i, np.array(k), rel) for i, (k, rel) in enumerate(vecs.items())]
    print(f"    Constructed {len(scen)} unique normalized weight scenarios.")

    sens_rows = []
    base_ranks = sub_ranked["rank"].to_numpy()

    for sid, w, rel in scen:
        s = np.round(sub_ranked[INPUTS].to_numpy() @ w, 12)
        r = rankdata(-s, method="min")
        sens_rows.append(pd.DataFrame({
            "scenario_id": sid,
            "population_enrichment_weight": w[0],
            "functional_consequence_weight": w[1],
            "network_centrality_weight": w[2],
            "disease_proximity_weight": w[3],
            "variant_key": sub_ranked["variant_key"],
            "gene_symbol": sub_ranked["gene_symbol"],
            "composite_score": s,
            "rank": r,
            "rank_change_from_baseline": r - base_ranks
        }))

    df_sens_long = pd.concat(sens_rows, ignore_index=True)
    df_sens_long.to_csv(OUT_SENS_LONG, index=False)
    print(f"[+] Saved long sensitivity table to {OUT_SENS_LONG} ({len(df_sens_long)} rows)")

    # Sensitivity Summary per Candidate (W10.7)
    summary_list = []
    for v_key, grp in df_sens_long.groupby("variant_key"):
        ranks = grp["rank"]
        top10_freq = (ranks <= 10).mean()
        q75, q25 = np.percentile(ranks, [75, 25])
        summary_list.append({
            "variant_key": v_key,
            "gene_symbol": grp["gene_symbol"].iloc[0],
            "top10_frequency": round(top10_freq, 4),
            "median_rank": float(np.median(ranks)),
            "min_rank": int(np.min(ranks)),
            "max_rank": int(np.max(ranks)),
            "rank_iqr": float(q75 - q25)
        })

    df_sens_summ = pd.DataFrame(summary_list).sort_values(by="median_rank")
    df_sens_summ.to_csv(OUT_SENS_SUMM, index=False)
    print(f"[+] Saved sensitivity summary table to {OUT_SENS_SUMM} ({len(df_sens_summ)} rows)")

    # 10. Baselines and Rank Correlations (W10.8)
    print("\n[*] Calculating baseline rank correlations (Spearman rho & Kendall tau)...")
    sub_ranked["evidence_ord"] = sub_ranked["star_allele"].map(lambda x: 4.0 if "*2" in str(x) or "*3" in str(x) else 2.0)
    
    corr_results = []
    for b_name, b_col in [("frequency_only", "P_raw"), ("evidence_only", "evidence_ord")]:
        rho, p_rho = spearmanr(sub_ranked["composite_score"], sub_ranked[b_col])
        tau, p_tau = kendalltau(sub_ranked["composite_score"], sub_ranked[b_col])
        corr_results.append({
            "baseline": b_name,
            "comparison_column": b_col,
            "spearman_rho": round(rho, 4),
            "spearman_p": p_rho,
            "kendall_tau": round(tau, 4),
            "kendall_p": p_tau
        })
    df_corr = pd.DataFrame(corr_results)
    print(df_corr.to_string(index=False))

    # Display Top 10 Ranked Candidates
    print("\n" + "="*95)
    print("TOP 10 PRIORITIZED CARDIOVASCULAR PHARMACOGENOMIC VARIANTS (BASELINE)")
    print("="*95)
    top10 = merged.sort_values(by="rank").head(10)
    cols_show = ["rank", "variant_key", "rsid", "gene_symbol", "star_allele", "P_n", "F_n", "C_n", "D_n", "composite_score"]
    print(top10[cols_show].to_string(index=False))

if __name__ == "__main__":
    main()
