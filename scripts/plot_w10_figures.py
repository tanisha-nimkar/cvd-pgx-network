#!/usr/bin/env python3
"""
scripts/plot_w10_figures.py
Week 10 Deliverables: fig6_ranking.png and fig7_sensitivity.png
Conforms strictly to W10.10.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TABLES_DIR = "results/tables"
FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial"]
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

# Load tables
df_ranked = pd.read_csv(os.path.join(TABLES_DIR, "prioritization_ranked.csv"))
df_sens = pd.read_csv(os.path.join(TABLES_DIR, "sensitivity_analysis.csv"))

# If star_allele is not in df_ranked, merge it from variants_annotated.csv
if "star_allele" not in df_ranked.columns:
    df_var = pd.read_csv(os.path.join(TABLES_DIR, "variants_annotated.csv"))
    df_ranked = df_ranked.merge(df_var[["variant_key", "star_allele"]], on="variant_key", how="left")

# Filter to complete-case ranked candidates
df_ok = df_ranked[df_ranked["candidate_status"] == "ranked"].sort_values(by="rank").copy()

# =============================================================================
# 1. FIGURE 6: Top 30 Horizontal Stacked Bars (fig6_ranking.png)
# =============================================================================
print("[*] Generating figures/fig6_ranking.png (Top 30 Stacked Decompositions)...")
top30 = df_ok.head(30).iloc[::-1].copy() # Reverse so Rank 1 is at top

fig, ax = plt.subplots(figsize=(10, 11), dpi=300)

y_pos = np.arange(len(top30))

# Weighted components (each weight is 0.25 at baseline)
p_bar = 0.25 * top30["P_n"].to_numpy()
f_bar = 0.25 * top30["F_n"].to_numpy()
c_bar = 0.25 * top30["C_n"].to_numpy()
d_bar = 0.25 * top30["D_n"].to_numpy()

# Stack bars horizontally
ax.barh(y_pos, p_bar, color="#2b5c8f", height=0.68, edgecolor="none", label="Population Enrichment (0.25 × P)")
ax.barh(y_pos, f_bar, left=p_bar, color="#d95f02", height=0.68, edgecolor="none", label="Functional Consequence (0.25 × F)")
ax.barh(y_pos, c_bar, left=p_bar + f_bar, color="#7570b3", height=0.68, edgecolor="none", label="Network Centrality (0.25 × C)")
ax.barh(y_pos, d_bar, left=p_bar + f_bar + c_bar, color="#1b9e77", height=0.68, edgecolor="none", label="Disease Proximity (0.25 × D)")

# Construct descriptive y-labels
labels = []
for _, r in top30.iterrows():
    raw_star = r.get("star_allele", "")
    star = str(raw_star).split(";")[0] if pd.notna(raw_star) and str(raw_star) != "nan" else ""
    star_str = f" {star}" if star else ""
    cpic_tag = " [CPIC]" if r.get("guideline_covered", False) else ""
    labels.append(f"{r['gene_symbol']}{star_str} ({r['rsid']}){cpic_tag}")

ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel("Composite Prioritization Score (Baseline Equal Weights = 0.25 × 4)", fontsize=9, fontweight="bold")
ax.set_xlim(0, 1.0)

ax.set_title("Top 30 Prioritized Cardiovascular Pharmacogenomic Variants (India vs. NFE)\nprioritization index, not a clinical risk estimate", fontsize=11, pad=12)

ax.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=8)
ax.grid(axis="x", linestyle=":", alpha=0.5)

plt.tight_layout()
fig6_path = os.path.join(FIG_DIR, "fig6_ranking.png")
plt.savefig(fig6_path)
plt.close()
print(f"[+] Saved {fig6_path}")

# =============================================================================
# 2. FIGURE 7: Sensitivity Analysis Rank Distribution (fig7_sensitivity.png)
# =============================================================================
print("\n[*] Generating figures/fig7_sensitivity.png (Scenario Rank Stability)...")
top_candidates = df_ok.head(40)["variant_key"].tolist()
df_sens_top = df_sens[df_sens["variant_key"].isin(top_candidates)].copy()

order_map = {vk: idx for idx, vk in enumerate(top_candidates)}
df_sens_top["order"] = df_sens_top["variant_key"].map(order_map)
df_sens_top = df_sens_top.sort_values(by="order")

plot_keys = top_candidates[::-1]
data_to_plot = [df_sens_top[df_sens_top["variant_key"] == vk]["rank"].to_numpy() for vk in plot_keys]
baseline_ranks = [df_ok[df_ok["variant_key"] == vk]["rank"].iloc[0] for vk in plot_keys]

fig, ax = plt.subplots(figsize=(9, 12), dpi=300)

bp = ax.boxplot(data_to_plot, orientation='horizontal', patch_artist=True, widths=0.55,
                boxprops=dict(facecolor="#c6dbef", color="#2b5c8f", linewidth=0.8),
                whiskerprops=dict(color="#2b5c8f", linewidth=0.8),
                capprops=dict(color="#2b5c8f", linewidth=0.8),
                medianprops=dict(color="#d95f02", linewidth=1.2),
                flierprops=dict(marker='o', markersize=2, alpha=0.3, markeredgecolor='none', markerfacecolor='#7f7f7f'))

y_coords = np.arange(1, len(plot_keys) + 1)
ax.scatter(baseline_ranks, y_coords, color="black", s=25, zorder=5, marker="D", label="Baseline Rank (Equal Weights)")

labels_sens = []
for vk in plot_keys:
    row = df_ok[df_ok["variant_key"] == vk].iloc[0]
    raw_star = row.get("star_allele", "")
    star = str(raw_star).split(";")[0] if pd.notna(raw_star) and str(raw_star) != "nan" else ""
    star_str = f" {star}" if star else ""
    labels_sens.append(f"{row['gene_symbol']}{star_str} ({row['rsid']})")

ax.set_yticks(y_coords)
ax.set_yticklabels(labels_sens, fontsize=8)
ax.set_xlabel("Rank Across 593 Sensitivity Scenarios (±50% Relative Weight Grid)", fontsize=9, fontweight="bold")
ax.set_title("Rank Stability Across Multi-Criteria Weight Sensitivity Grid\n(Top 40 Ranked Pharmacogenomic Candidates)", fontsize=11, pad=12)

ax.grid(axis="x", linestyle=":", alpha=0.5)
ax.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=8)

plt.tight_layout()
fig7_path = os.path.join(FIG_DIR, "fig7_sensitivity.png")
plt.savefig(fig7_path)
plt.close()
print(f"[+] Saved {fig7_path}")
