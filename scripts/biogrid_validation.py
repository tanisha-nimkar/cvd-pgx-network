#!/usr/bin/env python3
"""
scripts/biogrid_validation.py
W11.1 BioGRID Physical PPI Replication Pipeline
Conforms strictly to W11.1, V1, and W9.4 permutation specifications.
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import networkx as nx
from collections import deque
from scipy.stats import rankdata, spearmanr
from sklearn.preprocessing import MinMaxScaler

TABLES_DIR = "results/tables"
PROCESSED_DIR = "data/processed"
BIOGRID_EDGES = os.path.join(PROCESSED_DIR, "biogrid_edges.tsv")
DISEASE_MODULE_FILE = os.path.join(TABLES_DIR, "disease_module_final.csv")
STRING_METRICS_FILE = os.path.join(TABLES_DIR, "network_metrics.csv")
PRIMARY_RANKED_FILE = os.path.join(TABLES_DIR, "prioritization_ranked.csv")
VARIANTS_FILE = os.path.join(TABLES_DIR, "variants_annotated.csv")

OUT_BIOGRID_METRICS = os.path.join(TABLES_DIR, "biogrid_network_metrics.csv")
OUT_BIOGRID_DETAIL = os.path.join(TABLES_DIR, "validation_detail_biogrid.csv")

PANEL_GENES = ["CYP2C19", "CYP2C9", "VKORC1", "CYP4F2", "SLCO1B1", "ABCG2"]
NEG_CONTROL_GENE = "ADGRG6"
QUERY_LIST = PANEL_GENES + [NEG_CONTROL_GENE]

ORDER = ["population_enrichment", "functional_consequence", "network_centrality", "disease_proximity"]
INPUTS_BG = ["P_n", "F_n", "C_bg_n", "D_bg_n"]
W0 = {"population_enrichment": 0.25, "functional_consequence": 0.25, "network_centrality": 0.25, "disease_proximity": 0.25}

# -----------------------------------------------------------------------------
# W9.4 Core Graph & Proximity Functions
# -----------------------------------------------------------------------------
def dist_to_set(adj, sources):
    dist = {s: 0 for s in sources if s in adj}
    dq = deque(dist.keys())
    while dq:
        u = dq.popleft()
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                dq.append(v)
    return dist

def make_degree_bins(deg, min_bin=100):
    by_deg = {}
    for n, k in deg.items():
        by_deg.setdefault(k, []).append(n)
    bins, cur = [], []
    for k in sorted(by_deg):
        cur += by_deg[k]
        if len(cur) >= min_bin:
            bins.append(cur)
            cur = []
    if cur:
        if bins:
            bins[-1] += cur
        else:
            bins.append(cur)
    return bins, {n: i for i, b in enumerate(bins) for n in b}

def draw_matched(rng, bins, node2bin, observed, exclude=frozenset()):
    need = {}
    for n in observed:
        need[node2bin[n]] = need.get(node2bin[n], 0) + 1
    out = []
    for b, k in need.items():
        pool = [x for x in bins[b] if x not in exclude]
        if k > len(pool):
            raise ValueError(f"bin capacity: bin {b} needs {k}, pool has {len(pool)}")
        out += [pool[i] for i in rng.choice(len(pool), size=k, replace=False)]
    return out

def main():
    print("=== W11.1: BIOGRID REPLICATION PIPELINE ===")

    # 1. Load BioGRID Edges & Build Graph
    print(f"[*] Loading BioGRID physical edges from {BIOGRID_EDGES}...")
    df_edges = pd.read_csv(BIOGRID_EDGES, sep="\t")
    G_raw = nx.from_pandas_edgelist(df_edges, source="symbol1", target="symbol2")
    print(f"    Raw BioGRID Physical Graph: {G_raw.number_of_nodes():,} nodes, {G_raw.number_of_edges():,} edges")

    # Extract LCC
    lcc_nodes = max(nx.connected_components(G_raw), key=len)
    G_lcc = G_raw.subgraph(lcc_nodes).copy()
    print(f"    BioGRID LCC: {G_lcc.number_of_nodes():,} nodes, {G_lcc.number_of_edges():,} edges")

    # Check mapping of panel and module genes
    panel_in_lcc = [g for g in PANEL_GENES if g in G_lcc]
    print(f"    Panel genes in BioGRID LCC: {len(panel_in_lcc)}/{len(PANEL_GENES)} {panel_in_lcc}")
    missing_panel = set(PANEL_GENES) - set(panel_in_lcc)
    if missing_panel:
        print(f"    [!] Panel genes missing from BioGRID LCC: {missing_panel}")

    # Build BioGRID Adjacency
    adj_bg = {n: list(G_lcc.neighbors(n)) for n in G_lcc.nodes()}

    # 2. Centrality on BioGRID Ego-Network
    print("\n[*] Constructing BioGRID first-degree ego-network for panel genes...")
    ego_nodes = set(panel_in_lcc)
    for g in panel_in_lcc:
        ego_nodes.update(G_lcc.neighbors(g))
    G_ego = G_lcc.subgraph(ego_nodes).copy()
    print(f"    BioGRID Ego-Network: {G_ego.number_of_nodes():,} nodes, {G_ego.number_of_edges():,} edges")

    deg_ego = dict(G_ego.degree())
    bet_ego = nx.betweenness_centrality(G_ego, normalized=True)

    # 3. Proximity to Disease Module in BioGRID
    print(f"\n[*] Mapping direct CAD disease module to BioGRID LCC...")
    df_mod = pd.read_csv(DISEASE_MODULE_FILE)
    S_obs = [g for g in df_mod["gene_symbol"] if g in G_lcc]
    print(f"    CAD Module Genes in BioGRID LCC: {len(S_obs):,}/{len(df_mod):,}")

    deg_univ = {n: len(adj_bg[n]) for n in G_lcc.nodes()}
    bins, node2bin = make_degree_bins(deg_univ, min_bin=100)
    print(f"    Degree bins constructed: {len(bins)} bins.")

    # Calculate observed distance
    dist_obs = dist_to_set(adj_bg, S_obs)
    d_obs = {g: dist_obs.get(g, np.inf) for g in QUERY_LIST}

    # 1,000 degree-matched permutations
    print("[*] Running 1,000 degree-matched permutations on BioGRID LCC...")
    rng = np.random.default_rng(42)
    null = {g: [] for g in QUERY_LIST}
    panel_set = frozenset(PANEL_GENES)
    done, failed = 0, 0

    for iteration in range(1000):
        try:
            S_rand = draw_matched(rng, bins, node2bin, S_obs)
            dist_r = dist_to_set(adj_bg, S_rand)
            for g in QUERY_LIST:
                if g in G_lcc:
                    g_rand = draw_matched(rng, bins, node2bin, [g], exclude=panel_set)[0]
                    null[g].append(dist_r.get(g_rand, np.inf))
            done += 1
        except Exception as e:
            failed += 1

    print(f"    Permutations completed: {done} | Failed: {failed}")

    # Compute BioGRID statistics
    bg_stats = {}
    for g in QUERY_LIST:
        if g in G_lcc and len(null[g]) > 0:
            mu = float(np.mean(null[g]))
            sd = float(np.std(null[g], ddof=0))
            z = (d_obs[g] - mu) / sd if sd > 0 else np.nan
            p_emp = float((1 + np.sum(np.array(null[g]) <= d_obs[g])) / (1 + len(null[g])))
            bg_stats[g] = {
                "observed_proximity": d_obs[g],
                "permutation_mean": mu,
                "permutation_sd": sd,
                "proximity_z_score": z,
                "empirical_p_value": p_emp,
                "status": "ok" if sd > 0 else "zero_variance_null",
                "coverage": "in_lcc"
            }
        else:
            bg_stats[g] = {
                "observed_proximity": np.nan,
                "permutation_mean": np.nan,
                "permutation_sd": np.nan,
                "proximity_z_score": np.nan,
                "empirical_p_value": np.nan,
                "status": "not_in_biogrid_lcc",
                "coverage": "missing_from_lcc"
            }

    # 4. Save results/tables/biogrid_network_metrics.csv
    bg_metrics_rows = []
    for g in G_ego.nodes():
        row = {
            "gene_symbol": g,
            "degree": deg_ego.get(g, 0),
            "betweenness_centrality": bet_ego.get(g, 0.0),
            "is_panel_gene": g in PANEL_GENES,
            "observed_proximity": bg_stats.get(g, {}).get("observed_proximity", np.nan),
            "permutation_mean": bg_stats.get(g, {}).get("permutation_mean", np.nan),
            "permutation_sd": bg_stats.get(g, {}).get("permutation_sd", np.nan),
            "proximity_z_score": bg_stats.get(g, {}).get("proximity_z_score", np.nan),
            "empirical_p_value": bg_stats.get(g, {}).get("empirical_p_value", np.nan),
            "proximity_status": bg_stats.get(g, {}).get("status", "not_a_query_gene"),
            "biogrid_coverage": bg_stats.get(g, {}).get("coverage", "in_ego")
        }
        bg_metrics_rows.append(row)
    df_bg_metrics = pd.DataFrame(bg_metrics_rows)
    df_bg_metrics.to_csv(OUT_BIOGRID_METRICS, index=False)
    print(f"[+] Saved BioGRID network metrics to {OUT_BIOGRID_METRICS} ({len(df_bg_metrics)} rows)")

    # 5. Load Primary STRING Metrics for Comparison
    df_str = pd.read_csv(STRING_METRICS_FILE)
    df_str_panel = df_str[df_str["is_panel_gene"]].set_index("gene_symbol").to_dict(orient="index")

    # Load Primary Ranked Candidates
    df_ranked = pd.read_csv(PRIMARY_RANKED_FILE)
    top10_variants = df_ranked[df_ranked["candidate_status"] == "ranked"].sort_values(by="rank").head(10)["variant_key"].tolist()
    top10_genes = set(df_ranked[df_ranked["variant_key"].isin(top10_variants)]["gene_symbol"])

    # 6. Secondary Composite Replication on BioGRID
    print("\n[*] Building secondary biogrid_composite_replication on candidate pool...")
    df_rep = df_ranked[df_ranked["candidate_status"] == "ranked"].copy()

    # Attach BioGRID C and D
    max_deg_bg = max(deg_ego.values()) if deg_ego else 1
    c_raw_bg = {}
    d_raw_bg = {}
    for g in PANEL_GENES:
        if g in G_ego:
            c_raw_bg[g] = ((deg_ego.get(g, 0) / max_deg_bg) + bet_ego.get(g, 0.0)) / 2.0
            d_raw_bg[g] = -bg_stats[g]["proximity_z_score"] if not np.isnan(bg_stats[g]["proximity_z_score"]) else np.nan
        else:
            c_raw_bg[g] = np.nan
            d_raw_bg[g] = np.nan

    df_rep["C_bg_raw"] = df_rep["gene_symbol"].map(c_raw_bg)
    df_rep["D_bg_raw"] = df_rep["gene_symbol"].map(d_raw_bg)

    # Scale BioGRID components on unique panel genes
    bg_gene_df = pd.DataFrame([{"gene": g, "C": c_raw_bg[g], "D": d_raw_bg[g]} for g in PANEL_GENES]).dropna()
    if len(bg_gene_df) >= 2:
        sc_bg_gene = MinMaxScaler().fit(bg_gene_df[["C", "D"]])
        # Apply transform
        for idx, r in df_rep.iterrows():
            g = r["gene_symbol"]
            if pd.notna(c_raw_bg.get(g)) and pd.notna(d_raw_bg.get(g)):
                sample_df = pd.DataFrame([[c_raw_bg[g], d_raw_bg[g]]], columns=["C", "D"])
                transformed = sc_bg_gene.transform(sample_df)[0]
                df_rep.loc[idx, "C_bg_n"] = transformed[0]
                df_rep.loc[idx, "D_bg_n"] = transformed[1]
            else:
                df_rep.loc[idx, "C_bg_n"] = np.nan
                df_rep.loc[idx, "D_bg_n"] = np.nan

        # Calculate secondary composite
        complete_bg = df_rep[INPUTS_BG].notna().all(axis=1)
        df_rep["biogrid_composite_replication"] = np.nan
        df_rep["rank_biogrid"] = np.nan

        bg_scores = np.round(df_rep.loc[complete_bg, INPUTS_BG].to_numpy() @ np.array([W0[k] for k in ORDER]), 12)
        df_rep.loc[complete_bg, "biogrid_composite_replication"] = bg_scores
        df_rep.loc[complete_bg, "rank_biogrid"] = rankdata(-bg_scores, method="min")
    else:
        print("[!] Too few panel genes in BioGRID LCC to fit gene scaler.")
        df_rep["biogrid_composite_replication"] = np.nan
        df_rep["rank_biogrid"] = np.nan

    # 7. Build validation_detail_biogrid.csv
    detail_rows = []
    for g in PANEL_GENES + [NEG_CONTROL_GENE]:
        str_m = df_str_panel.get(g, {})
        bg_m = bg_stats.get(g, {})
        in_t10 = g in top10_genes

        str_z = str_m.get("proximity_z_score", np.nan)
        bg_z = bg_m.get("proximity_z_score", np.nan)
        str_p = str_m.get("empirical_p_value", np.nan)
        bg_p = bg_m.get("empirical_p_value", np.nan)

        # Apply V3 status definitions
        if bg_m.get("status") == "not_in_biogrid_lcc":
            status = "unavailable"
            notes = "Gene absent from BioGRID physical LCC; cannot evaluate physical network proximity."
        elif np.isnan(str_z) or np.isnan(bg_z):
            status = "insufficient_data"
            notes = "Missing proximity metrics in one or both networks."
        elif (str_p < 0.05 and bg_p < 0.05) or (str_p >= 0.05 and bg_p >= 0.05):
            status = "consistent"
            notes = f"Both networks yield consistent empirical inference (p_str={str_p:.4f}, p_bg={bg_p:.4f})."
        else:
            status = "discordant"
            notes = f"Inconsistent empirical inference across networks (p_str={str_p:.4f}, p_bg={bg_p:.4f})."

        # Compare ranks in primary vs replication
        g_ranked = df_rep[df_rep["gene_symbol"] == g]
        min_rank_str = g_ranked["rank"].min() if len(g_ranked) > 0 else np.nan
        min_rank_bg = g_ranked["rank_biogrid"].min() if len(g_ranked) > 0 else np.nan

        detail_rows.append({
            "gene": g,
            "in_top10": in_t10,
            "string_degree": str_m.get("degree", np.nan),
            "biogrid_degree": deg_ego.get(g, np.nan),
            "string_betweenness": str_m.get("betweenness_centrality", np.nan),
            "biogrid_betweenness": bet_ego.get(g, np.nan),
            "string_z": str_z,
            "biogrid_z": bg_z,
            "string_p": str_p,
            "biogrid_p": bg_p,
            "string_status": str_m.get("proximity_status", "ok"),
            "biogrid_status": bg_m.get("status", "ok"),
            "rank_string": min_rank_str,
            "rank_biogrid": min_rank_bg,
            "coverage": bg_m.get("coverage", "in_lcc"),
            "status": status,
            "notes": notes
        })

    df_detail = pd.DataFrame(detail_rows)
    df_detail.to_csv(OUT_BIOGRID_DETAIL, index=False)
    print(f"[+] Saved BioGRID validation detail to {OUT_BIOGRID_DETAIL}")

    # 8. Report Summary Table
    print("\n" + "="*100)
    print("W11.1 BIOGRID PHYSICAL PPI REPLICATION SUMMARY")
    print("="*100)
    cols_display = ["gene", "in_top10", "string_degree", "biogrid_degree", "string_z", "biogrid_z", "status"]
    print(df_detail[cols_display].to_string(index=False))

    # Calculate Top-10 Overlap
    if "rank_biogrid" in df_rep.columns and df_rep["rank_biogrid"].notna().sum() > 0:
        top10_primary = set(df_rep.sort_values(by="rank").head(10)["variant_key"])
        top10_rep = set(df_rep.sort_values(by="rank_biogrid").head(10)["variant_key"])
        overlap = len(top10_primary.intersection(top10_rep))
        print(f"\n[*] Top-10 Candidate Variant Overlap (STRING Primary vs. BioGRID Secondary): {overlap}/10 ({overlap * 10}%)")

if __name__ == "__main__":
    main()
