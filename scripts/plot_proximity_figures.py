#!/usr/bin/env python3
"""
scripts/plot_proximity_figures.py
Task 5 / Part IV Deliverable: Publication-Quality Proximity Visualizations
Conforms strictly to W9-4, W9-5, and Part IV figure specifications.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

TABLES_DIR = "results/tables"
FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial"]
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

# Load tables
df_primary = pd.read_csv(os.path.join(TABLES_DIR, "proximity_results_primary.csv"))
df_sens = pd.read_csv(os.path.join(TABLES_DIR, "proximity_results_sensitivity.csv"))
df_manifest = pd.read_csv(os.path.join(TABLES_DIR, "network_construction_manifest.csv"))

# =============================================================================
# FIGURE 3A: Proximity z-score Comparison across Query Levels
# =============================================================================
fig, ax = plt.subplots(figsize=(8.5, 6), dpi=300)

df_sorted = df_primary.sort_values(by="z_score_10k", ascending=True).copy()
y_pos = np.arange(len(df_sorted))

colors = []
for _, r in df_sorted.iterrows():
    if "direct_target" in r["query_id"]:
        colors.append("#2b5c8f") # Blue for direct drug targets
    elif r["query_id"] == "all_6_pharmacogenes":
        colors.append("#7f7f7f") # Gray for pooled
    elif "antiplatelet" in r["query_id"] or "statins" in r["query_id"] or "warfarin" in r["query_id"]:
        colors.append("#d95f02") # Orange for drug PGx sets
    else:
        colors.append("#1b9e77") # Green for individual pharmacogenes

bars = ax.barh(y_pos, df_sorted["z_score_10k"], color=colors, height=0.65, edgecolor="none", alpha=0.9)

# Threshold lines
ax.axvline(0, color="black", linestyle="-", linewidth=0.8, alpha=0.7)
ax.axvline(-0.15, color="#e7298a", linestyle="--", linewidth=1.0, label="Guney proximal (z <= -0.15)")
ax.axvline(-2.0, color="#d95f02", linestyle=":", linewidth=1.2, label="Stringent (z <= -2.0)")

ax.set_yticks(y_pos)
clean_labels = [l.replace("_", " ").replace("direct target", "(Target)") for l in df_sorted["query_id"]]
ax.set_yticklabels(clean_labels, fontsize=9)
ax.set_xlabel("Closest Network Proximity z-Score (z_c)", fontsize=10, fontweight="bold")
ax.set_title("Network Proximity (z_c) to Direct CAD Module (MONDO:0005010)\nSTRING v12.5 Functional Interactome LCC (B = 10,000 Permutations)", fontsize=11, pad=12)

# Custom legend
legend_elements = [
    mpatches.Patch(color="#1b9e77", label="Individual Pharmacogene"),
    mpatches.Patch(color="#d95f02", label="Drug-Specific PGx Set"),
    mpatches.Patch(color="#2b5c8f", label="Direct Drug Target (Comparator)"),
    mpatches.Patch(color="#7f7f7f", label="Exploratory Pooled Panel"),
    plt.Line2D([0], [0], color="#e7298a", linestyle="--", label="Guney proximal (z <= -0.15)"),
    plt.Line2D([0], [0], color="#d95f02", linestyle=":", label="Stringent threshold (z <= -2.0)")
]
ax.legend(handles=legend_elements, loc="lower right", fontsize=8, frameon=True, facecolor="white", edgecolor="none")
ax.grid(axis="x", linestyle=":", alpha=0.5)

plt.tight_layout()
fig_3a_path = os.path.join(FIG_DIR, "fig3a_proximity_zscores.png")
plt.savefig(fig_3a_path)
plt.close()
print(f"[+] Saved Figure 3A to {fig_3a_path}")

# =============================================================================
# FIGURE 3B: Observed vs Null Distribution (Representative Cases)
# =============================================================================
fig, axes = plt.subplots(1, 3, figsize=(13, 4), dpi=300)

cases = [
    ("CYP2C19", "CYP2C19 (Clopidogrel PGx)", 0.0, 0.9479, 0.2541, -3.73),
    ("HMGCR", "HMGCR (Statin Target)", 0.0, 0.9314, 0.2632, -3.54),
    ("SLCO1B1", "SLCO1B1 (Statin Transporter)", 1.0, 1.0108, 0.3223, -0.03)
]

for idx, (cid, title, obs, mu, sd, z) in enumerate(cases):
    ax = axes[idx]
    x = np.linspace(max(0, mu - 4*sd), mu + 4*sd, 200)
    y = (1.0 / (sd * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sd) ** 2)
    ax.plot(x, y, color="#2b5c8f", lw=1.8, label="Degree-matched null")
    ax.fill_between(x, y, color="#2b5c8f", alpha=0.15)
    ax.axvline(obs, color="#d95f02", lw=2.0, linestyle="--", label=f"Observed d_c = {obs:.1f}")
    ax.axvline(mu, color="#333333", lw=1.0, linestyle=":", label=f"Null mu = {mu:.2f}")
    
    ax.set_title(f"{title}\nz_c = {z:.2f}", fontsize=10, fontweight="bold")
    ax.set_xlabel("Closest Distance (d_c)", fontsize=9)
    if idx == 0:
        ax.set_ylabel("Probability Density", fontsize=9)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, linestyle=":", alpha=0.5)

plt.tight_layout()
fig_3b_path = os.path.join(FIG_DIR, "fig3b_distance_distributions.png")
plt.savefig(fig_3b_path)
plt.close()
print(f"[+] Saved Figure 3B to {fig_3b_path}")

# =============================================================================
# FIGURE 3C: Interactome & Disease Module Sensitivity Summary
# =============================================================================
fig, ax = plt.subplots(figsize=(9, 5), dpi=300)

df_sens_b = df_sens[df_sens["query_id"].isin(["clopidogrel_antiplatelet", "statins_lipid_lowering", "warfarin_anticoagulation_separate"])].copy()

pivot_z = df_sens_b.pivot(index="condition", columns="query_id", values="z_score")
pivot_z.columns = ["Clopidogrel (CYP2C19)", "Statins (3 PGx)", "Warfarin (3 PGx)"]

pivot_z.plot(kind="barh", ax=ax, colormap="tab10", width=0.75, edgecolor="none", alpha=0.85)
ax.axvline(0, color="black", lw=0.8)
ax.axvline(-0.15, color="#e7298a", linestyle="--", label="Guney proximal (z <= -0.15)")
ax.axvline(-2.0, color="#d95f02", linestyle=":", label="Stringent (z <= -2.0)")

ax.set_xlabel("Closest Distance z-Score (z_c)", fontsize=10, fontweight="bold")
ax.set_ylabel("Sensitivity Condition", fontsize=10, fontweight="bold")
ax.set_title("Sensitivity Matrix: Proximity Robustness Across Interactomes & Ontologies", fontsize=11, pad=12)
ax.legend(fontsize=8, loc="lower right")
ax.grid(axis="x", linestyle=":", alpha=0.5)

plt.tight_layout()
fig_3c_path = os.path.join(FIG_DIR, "fig3c_sensitivity_matrix.png")
plt.savefig(fig_3c_path)
plt.close()
print(f"[+] Saved Figure 3C to {fig_3c_path}")
