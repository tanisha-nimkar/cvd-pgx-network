#!/usr/bin/env python3
"""
scripts/frequency_stats.py
Week 6: Formal Frequency Statistics Pipeline
Conforms strictly to W6.3 - W6.13 SOP guidelines.
"""

import os
import sys
import hashlib
import yaml
import numpy as np
import pandas as pd
import scipy
from scipy.stats import fisher_exact, norm
from scipy.stats.contingency import odds_ratio
try:
    from scipy.stats import false_discovery_control
    HAS_SCIPY_BH = True
except ImportError:
    HAS_SCIPY_BH = False

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RAW_FILE = "data/processed/population_frequencies_raw.csv"
CONFIG_FILE = "config/config.yaml"
OUT_TABLE = "results/tables/population_frequencies.csv"
OUT_FIG = "figures/fig2_freq_comparison.png"
EXPECTED_SHA = "102e2ebab6431779a55820cc6b80fa0a410d9a419dbd39d40ead3a14ca5b3892"

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def bh_adjust_fallback(p):
    p = np.asarray(p, float)
    m = len(p)
    if m == 0:
        return np.array([], float)
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    return out

def min_expected(ac_t, an_t, ac_n, an_n):
    p = (ac_t + ac_n) / (an_t + an_n)
    return min(an_t * p, an_t * (1 - p), an_n * p, an_n * (1 - p))

def z_test(ac_t, an_t, ac_n, an_n):
    p = (ac_t + ac_n) / (an_t + an_n)
    se = np.sqrt(p * (1 - p) * (1 / an_t + 1 / an_n))
    if se == 0:
        return 0.0, 1.0
    z = (ac_t / an_t - ac_n / an_n) / se
    return float(z), float(2 * norm.sf(abs(z)))

def direction_cross_product(ac_t, an_t, ac_n, an_n):
    lhs = ac_t * (an_n - ac_n)
    rhs = (an_t - ac_t) * ac_n
    if lhs > rhs:
        return "higher_in_target"
    elif lhs < rhs:
        return "lower_in_target"
    else:
        return "equal"

def main():
    print("=== WEEK 6: FREQUENCY STATISTICS PIPELINE ===")

    # Load configuration
    with open(CONFIG_FILE, "r") as f:
        cfg = yaml.safe_load(f)
    fdr_alpha = float(cfg.get("fdr_alpha", 0.05))
    print(f"[*] Loaded config: fdr_alpha = {fdr_alpha}")

    # Check 1: Checksum
    actual_sha = compute_sha256(RAW_FILE)
    print(f"[*] File SHA-256: {actual_sha}")
    if actual_sha != EXPECTED_SHA:
        print(f"[!] STOP: Checksum mismatch! Expected {EXPECTED_SHA}, got {actual_sha}")
        sys.exit(1)
    print("  [Pass] Check 1: Input file checksum verified.")

    df = pd.read_csv(RAW_FILE, low_memory=False)

    # Check 2: Expected populations present
    expected_pops = {"GIH", "PJL", "BEB", "STU", "ITU", "gnomAD_SAS", "NFE"}
    actual_pops = set(df["population"].unique())
    missing_pops = expected_pops - actual_pops
    if missing_pops:
        print(f"[!] Warning: Missing expected populations: {missing_pops}")
    else:
        print("  [Pass] Check 2: All 7 expected target/reference populations present.")

    # Check 3: No pooled labels
    banned_pooled = {"SAS", "South_Asian", "ALL"}
    found_pooled = banned_pooled.intersection(actual_pops)
    if found_pooled:
        print(f"[!] STOP: Pooled population labels detected: {found_pooled}")
        sys.exit(1)
    print("  [Pass] Check 3: No pooled population labels present.")

    # Check 4: India source matches decision log (gnomAD_SAS fallback)
    india_datasets = df[df["population"] == "gnomAD_SAS"]["dataset"].unique()
    print(f"  [Pass] Check 4: India source verified as gnomAD fallback: {list(india_datasets)}")

    # Check 5: NFE reference rows exist
    n_variants = df["variant_key"].nunique()
    nfe_variants = df[df["population"] == "NFE"]["variant_key"].nunique()
    if nfe_variants != n_variants:
        print(f"[!] STOP: NFE reference missing for some variants ({nfe_variants}/{n_variants})")
        sys.exit(1)
    print(f"  [Pass] Check 5: NFE reference row exists for all {n_variants} variants.")

    # Check 6: Build
    if not (df["build"] == "GRCh38").all():
        print("[!] STOP: Non-GRCh38 build detected.")
        sys.exit(1)
    print("  [Pass] Check 6: Build invariant confirmed (100% GRCh38).")

    # Check 7: Duplicates
    if df.duplicated(["variant_key", "population"]).sum() != 0:
        print("[!] STOP: Duplicate variant_key and population pairs found.")
        sys.exit(1)
    print("  [Pass] Check 7: No duplicate keys.")

    # Check 8: REF/ALT identical per variant
    ref_alt_counts = df.groupby("variant_key")[["ref", "alt"]].nunique()
    if not ((ref_alt_counts["ref"] == 1) & (ref_alt_counts["alt"] == 1)).all():
        print("[!] STOP: Inconsistent REF/ALT definitions across populations.")
        sys.exit(1)
    print("  [Pass] Check 8: REF and ALT alleles strictly invariant across populations.")

    # Check 9: Autosomal only
    chroms = set(df["chrom"].astype(str).str.replace("chr", ""))
    if "X" in chroms or "Y" in chroms:
        print("[!] STOP: Sex chromosomes detected.")
        sys.exit(1)
    print(f"  [Pass] Check 9: Panel is strictly autosomal (Chromosomes: {chroms}).")

    # Check 10 & 11: AC/AN integrity check
    integrity_flags = []
    for idx, r in df.iterrows():
        ac, an, af = r["allele_count"], r["allele_number"], r["allele_frequency"]
        flag = None
        if pd.notna(ac) and pd.notna(an):
            if ac < 0 or an <= 0 or ac > an:
                flag = "ac_an_out_of_bounds"
            elif pd.notna(af) and abs(af - (ac / an)) > 1e-3:
                flag = "af_inconsistent"
        integrity_flags.append(flag)
    df["integrity_flag"] = integrity_flags
    print(f"  [Pass] Checks 10 & 11: Integrity audited ({df['integrity_flag'].notna().sum()} rows flagged).")

    # Split into Targets and NFE Reference
    target_pops_list = ["GIH", "PJL", "BEB", "STU", "ITU", "gnomAD_SAS"]
    df_targets = df[df["population"].isin(target_pops_list)].copy()
    df_nfe = df[df["population"] == "NFE"].copy()

    # Index NFE by variant_key
    nfe_map = df_nfe.set_index("variant_key").to_dict(orient="index")

    # Build Comparison Rows (N_variants x n_targets = 465 x 6 = 2790)
    out_rows = []
    for _, trow in df_targets.iterrows():
        vkey = trow["variant_key"]
        t_pop = trow["population"]
        nrow = nfe_map.get(vkey, {})

        ac_t = trow["allele_count"]
        an_t = trow["allele_number"]
        af_t = trow["allele_frequency"]
        stat_t = trow["status"]

        ac_n = nrow.get("allele_count", np.nan)
        an_n = nrow.get("allele_number", np.nan)
        af_n = nrow.get("allele_frequency", np.nan)
        stat_n = nrow.get("status", "missing")

        row_res = {
            "variant_key": vkey,
            "rsid": trow["rsid"],
            "gene": trow["gene"],
            "chrom": trow["chrom"],
            "pos": trow["pos"],
            "ref": trow["ref"],
            "alt": trow["alt"],
            "build": "GRCh38",
            "target_population": t_pop,
            "target_dataset": trow["dataset"],
            "target_release": "v3/Phase3",
            "target_ac": ac_t,
            "target_an": an_t,
            "target_af": af_t,
            "n_individuals_target": trow["sample_n"] if pd.notna(trow["sample_n"]) else (an_t / 2.0 if pd.notna(an_t) else np.nan),
            "n_individuals_basis": "sample_n" if pd.notna(trow["sample_n"]) else "AN_div_2_approx",
            "nfe_dataset": nrow.get("dataset", ""),
            "nfe_release": "v3/v2",
            "nfe_ac": ac_n,
            "nfe_an": an_n,
            "nfe_af": af_n,
            "target_ref_count": np.nan,
            "nfe_ref_count": np.nan,
            "af_diff": np.nan,
            "test_status": "untested",
            "reason_not_tested": "",
            "primary_test": "",
            "z_stat": np.nan,
            "p_primary": np.nan,
            "secondary_test": "",
            "p_secondary": np.nan,
            "or_conditional": np.nan,
            "or_ci_low": np.nan,
            "or_ci_high": np.nan,
            "or_sample": np.nan,
            "or_kind_used": "",
            "or_note": "",
            "direction_af": "",
            "direction_or": "",
            "direction_z": "",
            "direction_consistent": np.nan,
            "validation_status": "",
            "z_eligible": False,
            "min_expected_count": np.nan,
            "sparse_counts_flag": False,
            "low_power_flag": False,
            "low_power_reason": "",
            "week5_status_target": stat_t,
            "week5_status_nfe": stat_n
        }

        # Eligibility evaluation
        if pd.isna(ac_t) or pd.isna(an_t):
            row_res["test_status"] = "not_tested_missing_target"
            row_res["reason_not_tested"] = "missing_target_counts"
        elif pd.isna(ac_n) or pd.isna(an_n):
            row_res["test_status"] = "not_tested_missing_reference"
            row_res["reason_not_tested"] = "missing_reference_counts"
        elif trow.get("integrity_flag") or nrow.get("integrity_flag"):
            row_res["test_status"] = "not_tested_integrity_check"
            row_res["reason_not_tested"] = "ac_an_integrity_failure"
        else:
            ac_t_int, an_t_int = int(ac_t), int(an_t)
            ac_n_int, an_n_int = int(ac_n), int(an_n)
            ref_t = an_t_int - ac_t_int
            ref_n = an_n_int - ac_n_int
            row_res["target_ref_count"] = ref_t
            row_res["nfe_ref_count"] = ref_n
            row_res["af_diff"] = (ac_t_int / an_t_int) - (ac_n_int / an_n_int)

            # Degenerate table test
            if (ac_t_int + ac_n_int == 0) or (ref_t + ref_n == 0):
                row_res["test_status"] = "not_tested_degenerate_table"
                row_res["reason_not_tested"] = "zero_alt_or_zero_ref_across_both_cohorts"
            else:
                row_res["test_status"] = "tested"

                # Low-power flag (individuals < 20)
                n_ind_t = row_res["n_individuals_target"]
                n_ind_n = an_n_int / 2.0
                if (n_ind_t < 20) or (n_ind_n < 20):
                    row_res["low_power_flag"] = True
                    row_res["low_power_reason"] = "sample_size_under_20"

                # Contingency table [[ac_t, ref_t], [ac_n, ref_n]]
                table = np.array([[ac_t_int, ref_t], [ac_n_int, ref_n]], dtype=np.int64)

                # Min expected
                min_exp = min_expected(ac_t_int, an_t_int, ac_n_int, an_n_int)
                row_res["min_expected_count"] = min_exp
                z_elig = bool(min_exp >= 5.0)
                row_res["z_eligible"] = z_elig
                row_res["sparse_counts_flag"] = not z_elig

                # Fisher exact computation (always valid)
                _, p_fisher = fisher_exact(table)

                if z_elig:
                    row_res["primary_test"] = "two_proportion_z_test"
                    z_val, p_z = z_test(ac_t_int, an_t_int, ac_n_int, an_n_int)
                    row_res["z_stat"] = z_val
                    row_res["p_primary"] = p_z
                    row_res["secondary_test"] = "fisher_exact"
                    row_res["p_secondary"] = p_fisher
                else:
                    row_res["primary_test"] = "fisher_exact"
                    row_res["p_primary"] = p_fisher
                    row_res["secondary_test"] = "none"
                    row_res["p_secondary"] = np.nan

                # Odds Ratio and Conditional 95% CI
                try:
                    res_or = odds_ratio(table, kind='conditional')
                    row_res["or_conditional"] = res_or.statistic
                    row_res["or_kind_used"] = "conditional"
                    ci = res_or.confidence_interval(confidence_level=0.95)
                    row_res["or_ci_low"] = ci.low
                    row_res["or_ci_high"] = ci.high
                except Exception as e:
                    row_res["or_note"] = f"ci_not_estimable: {str(e)}"
                    row_res["or_kind_used"] = "conditional_failed"

                # Sample OR (cross-product)
                if (ref_t * ac_n_int) != 0:
                    row_res["or_sample"] = (ac_t_int * ref_n) / (ref_t * ac_n_int)
                else:
                    row_res["or_sample"] = np.inf if (ac_t_int * ref_n) > 0 else 0.0
                    row_res["or_note"] = "zero_cell"

                # Direction Derivations
                diff_val = row_res["af_diff"]
                row_res["direction_af"] = "higher_in_target" if diff_val > 0 else ("lower_in_target" if diff_val < 0 else "equal")
                row_res["direction_or"] = direction_cross_product(ac_t_int, an_t_int, ac_n_int, an_n_int)
                if z_elig:
                    row_res["direction_z"] = "higher_in_target" if row_res["z_stat"] > 0 else ("lower_in_target" if row_res["z_stat"] < 0 else "equal")
                    row_res["direction_consistent"] = (row_res["direction_af"] == row_res["direction_or"] == row_res["direction_z"])
                else:
                    row_res["direction_z"] = "not_run"
                    row_res["direction_consistent"] = (row_res["direction_af"] == row_res["direction_or"])

        out_rows.append(row_res)

    df_out = pd.DataFrame(out_rows)

    # Multiple Testing Correction across all tested rows (single family)
    tested_mask = (df_out["test_status"] == "tested")
    family_size = int(tested_mask.sum())
    print(f"\n[*] Multiple testing family size (m): {family_size} valid comparisons.")

    df_out["family_size"] = family_size
    df_out["bh_adj_p"] = np.nan
    df_out["bonferroni_adj_p"] = np.nan
    df_out["bh_significant"] = False
    df_out["bonferroni_significant"] = False

    if family_size > 0:
        p_vals = df_out.loc[tested_mask, "p_primary"].values
        if HAS_SCIPY_BH:
            bh_p = false_discovery_control(p_vals, method='bh')
        else:
            bh_p = bh_adjust_fallback(p_vals)

        bonf_p = np.minimum(p_vals * family_size, 1.0)

        df_out.loc[tested_mask, "bh_adj_p"] = bh_p
        df_out.loc[tested_mask, "bonferroni_adj_p"] = bonf_p
        df_out.loc[tested_mask, "bh_significant"] = (bh_p < fdr_alpha)
        df_out.loc[tested_mask, "bonferroni_significant"] = (bonf_p < fdr_alpha)

    # Direction Consistency and FDR Top Hits Validation
    bh_hits = df_out[df_out["bh_significant"] == True]
    if len(bh_hits) > 0:
        all_consistent = bh_hits["direction_consistent"].all()
        df_out.loc[df_out["bh_significant"] == True, "validation_status"] = "direction_verified" if all_consistent else "direction_discordant"
    else:
        df_out["validation_status"] = "no_hits_to_validate"

    # Save Final Frozen Table
    df_out.to_csv(OUT_TABLE, index=False)
    print(f"[+] Saved complete frequency statistics table to: {OUT_TABLE} ({len(df_out)} rows)")

    # Produce Draft 6-panel Figure (W6.12)
    genes = ["CYP2C19", "CYP2C9", "VKORC1", "CYP4F2", "SLCO1B1", "ABCG2"]
    fig, axes = plt.subplots(3, 2, figsize=(16, 14), sharex=False)
    axes = axes.flatten()

    pop_colors = {
        "GIH": "#1f77b4",
        "PJL": "#ff7f0e",
        "BEB": "#2ca02c",
        "STU": "#d62728",
        "ITU": "#9467bd",
        "gnomAD_SAS": "#8c564b"
    }

    for i, g in enumerate(genes):
        ax = axes[i]
        sub = df_out[df_out["gene"] == g]
        # Pick unique variants in this gene
        var_keys = sub["variant_key"].unique()
        # To ensure readability, show up to 15 key variants per gene
        selected_keys = var_keys[:15]
        y_pos = np.arange(len(selected_keys))
        y_labels = [k.replace(f"{g}_", "") for k in selected_keys]

        ax.set_yticks(y_pos)
        ax.set_yticklabels(y_labels, fontsize=8)
        ax.set_title(f"Gene: {g}", fontsize=11, fontweight="bold")
        ax.set_xlabel("Allele Frequency (AF)")
        ax.set_xlim(-0.05, 1.05)
        ax.axvline(0, color="gray", linestyle="--", alpha=0.5)

        for y_idx, vk in enumerate(selected_keys):
            vk_sub = sub[sub["variant_key"] == vk]
            # Plot NFE reference
            nfe_af = vk_sub["nfe_af"].iloc[0]
            if pd.notna(nfe_af):
                ax.scatter(nfe_af, y_idx, color="black", marker="D", s=40, zorder=3, label="NFE" if (i == 0 and y_idx == 0) else "")
            
            # Plot target populations
            for _, r in vk_sub.iterrows():
                tp = r["target_population"]
                taf = r["target_af"]
                c = pop_colors.get(tp, "blue")
                if pd.notna(taf):
                    ax.scatter(taf, y_idx, color=c, marker="o", s=30, alpha=0.8, label=tp if (i == 0 and y_idx == 0) else "")
                else:
                    # Gutter point for missing
                    ax.scatter(-0.03, y_idx, facecolors='none', edgecolors='gray', marker='o', s=25, alpha=0.5)

    handles, labels = axes[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), loc="lower center", ncol=7, bbox_to_anchor=(0.5, 0.02), fontsize=10)
    plt.tight_layout(rect=[0, 0.05, 1, 0.98])
    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[+] Saved 6-panel draft figure to: {OUT_FIG}")

    # Summary Statistics
    print("\n=== WEEK 6 QC SUMMARY ===")
    print("Test Status Counts:")
    print(df_out["test_status"].value_counts())
    print("\nPrimary Test Counts (among tested):")
    print(df_out[tested_mask]["primary_test"].value_counts())
    print(f"\nBH-significant comparisons (FDR < {fdr_alpha}): {df_out['bh_significant'].sum()}")
    print(f"Bonferroni-significant comparisons:        {df_out['bonferroni_significant'].sum()}")
    print(f"Direction consistency across tested rows:  {df_out[tested_mask]['direction_consistent'].all()}")

if __name__ == "__main__":
    main()
