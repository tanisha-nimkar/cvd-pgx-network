#!/usr/bin/env python3
"""
scripts/plot_w9_figures.py
Task 5 / W9.7 Deliverables: fig3_network.png and fig5_proximity.png
Conforms strictly to W9.7 visual specifications.
"""

import os
import json
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter

FIG_DIR = "figures"
TABLES_DIR = "results/tables"
NETWORK_FILE = "results/network/network.graphml"
METRICS_FILE = os.path.join(TABLES_DIR, "network_metrics.csv")
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial"]
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

# Load data
G = nx.read_graphml(NETWORK_FILE)
df_metrics = pd.read_csv(METRICS_FILE)
panel_genes = ["CYP2C19", "CYP2C9", "VKORC1", "CYP4F2", "SLCO1B1", "ABCG2"]
query_genes = panel_genes + ["ADGRG6"]

# =============================================================================
# 1. FIGURE 3: Ego-Network Visualization (fig3_network.png)
# =============================================================================
print("[*] Generating figures/fig3_network.png (W9.7 ego-network layout)...")
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

pos = nx.spring_layout(G, seed=42, k=0.25, iterations=50)

# Extract degree sizes from network_metrics.csv
deg_map = dict(zip(df_metrics["gene_symbol"], df_metrics["degree"]))

node_colors = []
node_sizes = []
node_shapes = []

for n, d in G.nodes(data=True):
    sym = d.get("gene_symbol", n)
    k = deg_map.get(sym, 1)
    # Size scaled by degree
    sz = max(40, k * 14)
    node_sizes.append(sz)
    if sym in panel_genes:
        node_colors.append("#d95f02") # Distinct orange for panel
    else:
        node_colors.append("#7570b3") # Purple for interactors

# Draw edges
nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.25, edge_color="#888888", width=0.7)

# Draw interactors
interactor_nodes = [n for n, d in G.nodes(data=True) if d.get("gene_symbol") not in panel_genes]
interactor_sizes = [deg_map.get(G.nodes[n].get("gene_symbol"), 1) * 14 for n in interactor_nodes]
nx.draw_networkx_nodes(G, pos, nodelist=interactor_nodes, node_color="#7570b3",
                       node_size=interactor_sizes, alpha=0.7, ax=ax, label="Interactors")

# Draw panel genes
panel_node_ids = [n for n, d in G.nodes(data=True) if d.get("gene_symbol") in panel_genes]
panel_sizes = [deg_map.get(G.nodes[n].get("gene_symbol"), 1) * 14 for n in panel_node_ids]
nx.draw_networkx_nodes(G, pos, nodelist=panel_node_ids, node_color="#d95f02",
                       node_size=panel_sizes, alpha=0.95, ax=ax, label="Panel Pharmacogenes")

# Label ONLY panel genes
labels = {n: G.nodes[n].get("gene_symbol") for n in panel_node_ids}
nx.draw_networkx_labels(G, pos, labels=labels, font_size=10, font_weight="bold",
                        font_color="#111111", ax=ax)

ax.set_title("Week 7 Cardiovascular Pharmacogenomics Ego-Network (135 Nodes, 225 Edges)\nNode size encodes degree; panel seeds labeled", fontsize=11, pad=12)
ax.axis("off")

# Legend
leg = ax.legend(scatterpoints=1, loc="lower right", frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9)
leg.get_frame().set_alpha(0.9)

plt.tight_layout()
fig3_path = os.path.join(FIG_DIR, "fig3_network.png")
plt.savefig(fig3_path)
plt.close()
print(f"[+] Saved {fig3_path}")

# =============================================================================
# 2. FIGURE 5: Discrete Integer Proximity Null Distributions (fig5_proximity.png)
# =============================================================================
print("\n[*] Generating figures/fig5_proximity.png (W9.7 discrete null distributions)...")
fig, axes = plt.subplots(2, 4, figsize=(14, 7), dpi=300)
axes_flat = axes.flatten()

metrics_map = df_metrics.set_index("gene_symbol").to_dict(orient="index")

for idx, g in enumerate(query_genes):
    ax = axes_flat[idx]
    null_file = os.path.join(TABLES_DIR, f"proximity_null_{g}.csv")
    if not os.path.exists(null_file):
        ax.text(0.5, 0.5, f"{g}\nMissing data", ha="center", va="center")
        ax.set_title(g)
        continue
    
    df_null = pd.read_csv(null_file)
    null_vals = df_null["null_distance"].astype(int)
    counts = Counter(null_vals)
    total = len(null_vals)
    
    # Range of discrete distances
    x_vals = sorted(counts.keys())
    freqs = [counts[x] / total for x in x_vals]
    
    meta = metrics_map.get(g, {})
    obs_d = meta.get("observed_proximity", np.nan)
    z_val = meta.get("proximity_z_score", np.nan)
    p_val = meta.get("empirical_p_value", np.nan)
    status = meta.get("proximity_status", "ok")
    
    # Discrete bar plot
    ax.bar(x_vals, freqs, width=0.45, color="#2b5c8f", alpha=0.75, edgecolor="#1a3b5c", label="Null frequency")
    
    # Mark observed distance
    if not np.isnan(obs_d):
        ax.axvline(obs_d, color="#d95f02", linestyle="--", linewidth=1.8, label=f"Observed d={int(obs_d)}")
        
    ax.set_xticks(range(min(x_vals + [int(obs_d) if not np.isnan(obs_d) else 0]), max(x_vals + [int(obs_d) if not np.isnan(obs_d) else 1]) + 1))
    ax.set_xlabel("Shortest Path Distance", fontsize=8)
    ax.set_ylabel("Relative Frequency", fontsize=8)
    
    # Subplot title & stats as plain text (no stars, no color coding)
    title_label = f"{g} (Negative Control)" if g == "ADGRG6" else f"{g} (Panel)"
    ax.set_title(title_label, fontsize=10, fontweight="bold")
    
    z_str = f"z = {z_val:.2f}" if not np.isnan(z_val) else "z = NaN"
    p_str = f"p_emp = {p_val:.4f}" if not np.isnan(p_val) else "p_emp = NA"
    ax.text(0.95, 0.85, f"{z_str}\n{p_str}\nstatus: {status}", transform=ax.transAxes,
            fontsize=8, ha="right", va="top", bbox=dict(boxstyle="round,pad=0.3", facecolor="#f8f9fa", edgecolor="#cccccc"))
    
    ax.grid(axis="y", linestyle=":", alpha=0.5)
    ax.legend(fontsize=7, loc="upper left")

# Empty slot 8
axes_flat[7].axis("off")

plt.tight_layout()
fig5_path = os.path.join(FIG_DIR, "fig5_proximity.png")
plt.savefig(fig5_path)
plt.close()
print(f"[+] Saved {fig5_path}")
