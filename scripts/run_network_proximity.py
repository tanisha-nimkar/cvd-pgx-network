#!/usr/bin/env python3
"""
scripts/run_network_proximity.py
Tasks 5, 6, 7: Optimized Network Proximity Analysis & Degree-Preserved Permutation Testing
Conforms strictly to W9-4, W9-5, W9-6, and Part IV deliverables.
"""

import os
import sys
import json
import time
import random
import numpy as np
import pandas as pd
import networkx as nx
from collections import defaultdict, Counter

TABLES_DIR = "results/tables"
PROCESSED_DIR = "data/processed"
PRIMARY_OUT = os.path.join(TABLES_DIR, "proximity_results_primary.csv")
SENSITIVITY_OUT = os.path.join(TABLES_DIR, "proximity_results_sensitivity.csv")
DIAGNOSTICS_OUT = os.path.join(TABLES_DIR, "randomisation_diagnostics.csv")

SEED = 42
N_PERM_PRIMARY = 10000
N_PERM_GUNEY = 1000
MIN_BIN_SIZE = 100

def load_graph_from_tsv(filepath):
    df = pd.read_csv(filepath, sep="\t")
    G = nx.Graph()
    for _, r in df.iterrows():
        G.add_edge(r["symbol1"], r["symbol2"])
    return G

def get_degree_bins(G, min_bin_size=100):
    deg_dict = dict(G.degree())
    degree_to_nodes = defaultdict(list)
    for node, deg in deg_dict.items():
        degree_to_nodes[deg].append(node)
    
    sorted_degs = sorted(degree_to_nodes.keys())
    bins = []
    curr_nodes = []
    curr_degs = []
    
    for d in sorted_degs:
        curr_nodes.extend(degree_to_nodes[d])
        curr_degs.append(d)
        if len(curr_nodes) >= min_bin_size:
            bins.append((min(curr_degs), max(curr_degs), list(curr_nodes)))
            curr_nodes = []
            curr_degs = []
    
    if curr_nodes:
        if bins:
            prev_min, _, prev_nodes = bins[-1]
            bins[-1] = (prev_min, max(curr_degs), prev_nodes + curr_nodes)
        else:
            bins.append((min(curr_degs), max(curr_degs), list(curr_nodes)))
            
    node_to_bin = {}
    for b_idx, (_, _, b_nodes) in enumerate(bins):
        for n in b_nodes:
            node_to_bin[n] = b_idx
            
    return bins, node_to_bin

def get_bin_counts(nodes, node_to_bin):
    counts = Counter()
    for n in nodes:
        counts[node_to_bin[n]] += 1
    return counts

def degree_matched_sample_fast(bin_counts, bin_pools):
    sampled = []
    for b_idx, k in bin_counts.items():
        pool = bin_pools[b_idx]
        if len(pool) >= k:
            sampled.extend(random.sample(pool, k))
        else:
            sampled.extend(random.choices(pool, k=k))
    return set(sampled)

def closest_distance_fast(adj_list, T_set, S_set):
    valid_T = [t for t in T_set if t in adj_list]
    if not valid_T or not S_set:
        return np.nan
    
    min_distances = []
    for t in valid_T:
        if t in S_set:
            min_distances.append(0.0)
            continue
        
        visited = {t}
        queue = [t]
        found = False
        dist = 0
        idx = 0
        
        while idx < len(queue) and not found:
            next_queue = []
            dist += 1
            for curr in queue[idx:]:
                for nbr in adj_list[curr]:
                    if nbr not in visited:
                        if nbr in S_set:
                            min_distances.append(float(dist))
                            found = True
                            break
                        visited.add(nbr)
                        next_queue.append(nbr)
                if found:
                    break
            queue = next_queue
            idx = 0
            
        if not found:
            min_distances.append(999999.0)
            
    return float(np.mean(min_distances))

def run_permutations_fast(adj_list, T_bin_counts, S_bin_counts, bin_pools, n_perm):
    null_dists = np.empty(n_perm, dtype=np.float64)
    for i in range(n_perm):
        T_rand = degree_matched_sample_fast(T_bin_counts, bin_pools)
        S_rand = degree_matched_sample_fast(S_bin_counts, bin_pools)
        null_dists[i] = closest_distance_fast(adj_list, T_rand, S_rand)
    return null_dists

def benjamini_hochberg(p_values):
    p_vals = np.array(p_values)
    n = len(p_vals)
    sorted_indices = np.argsort(p_vals)
    sorted_p = p_vals[sorted_indices]
    q_vals = np.zeros(n)
    q_vals[-1] = sorted_p[-1]
    for i in range(n - 2, -1, -1):
        q_vals[i] = min(sorted_p[i] * n / (i + 1), q_vals[i + 1])
    q_vals = np.clip(q_vals, 0.0, 1.0)
    adjusted = np.zeros(n)
    adjusted[sorted_indices] = q_vals
    return adjusted.tolist()

def main():
    print("=== WEEK 9: PROXIMITY ANALYSIS & STATISTICAL PERMUTATION TESTING ===")
    random.seed(SEED)
    np.random.seed(SEED)

    # 1. Load Query Sets
    with open(os.path.join(TABLES_DIR, "query_sets_locked.json")) as f:
        query_cfg = json.load(f)
    panel_set = set(query_cfg["level_C_exploratory_whole_panel"]["all_6_pharmacogenes"])

    # 2. Load Networks
    print("\n[*] Loading Interactome Networks...")
    G_primary = load_graph_from_tsv(os.path.join(PROCESSED_DIR, "string_functional_universe_edges.tsv"))
    G_tm = load_graph_from_tsv(os.path.join(PROCESSED_DIR, "string_functional_tm_excluded_edges.tsv"))
    G_phys = load_graph_from_tsv(os.path.join(PROCESSED_DIR, "string_physical_universe_edges.tsv"))

    adj_primary = {n: list(G_primary.neighbors(n)) for n in G_primary.nodes()}
    adj_tm = {n: list(G_tm.neighbors(n)) for n in G_tm.nodes()}
    adj_phys = {n: list(G_phys.neighbors(n)) for n in G_phys.nodes()}

    print(f"    Primary Functional LCC: {G_primary.number_of_nodes():,} nodes, {G_primary.number_of_edges():,} edges.")
    print(f"    Sensitivity TM-Ex LCC:  {G_tm.number_of_nodes():,} nodes, {G_tm.number_of_edges():,} edges.")
    print(f"    Physical LCC:           {G_phys.number_of_nodes():,} nodes, {G_phys.number_of_edges():,} edges.")

    # 3. Load Disease Modules
    print("\n[*] Loading Disease Modules...")
    df_cad_dir = pd.read_csv(os.path.join(TABLES_DIR, "disease_module_direct_mondo0005010.csv"))
    df_cad_indir = pd.read_csv(os.path.join(TABLES_DIR, "disease_module_indirect_mondo0005010.csv"))
    df_mi_dir = pd.read_csv(os.path.join(TABLES_DIR, "disease_module_mi_direct.csv"))
    df_ihd_dir = pd.read_csv(os.path.join(TABLES_DIR, "disease_module_ihd_direct.csv"))

    cad_primary_genes = set(df_cad_dir[df_cad_dir["association_score"] >= 0.10]["gene_symbol"])
    print(f"    Primary CAD Module (Direct >= 0.10): {len(cad_primary_genes)} genes.")

    # Degree Bins for Primary
    print("\n[*] Precomputing degree bins for Primary Functional LCC...")
    bins_prim, n2b_prim = get_degree_bins(G_primary, min_bin_size=MIN_BIN_SIZE)
    print(f"    Constructed {len(bins_prim)} degree bins.")

    diag_rows = [{
        "bin_id": idx, "min_degree": b[0], "max_degree": b[1], "n_nodes": len(b[2])
    } for idx, b in enumerate(bins_prim)]
    pd.DataFrame(diag_rows).to_csv(DIAGNOSTICS_OUT, index=False)
    print(f"[+] Saved degree bin diagnostics to {DIAGNOSTICS_OUT}")

    # Pre-filter candidate pools for Null 1 (panel excluded) and Null 2 (panel included)
    bin_pools_ex = [[n for n in b[2] if n not in panel_set] for b in bins_prim]
    bin_pools_inc = [list(b[2]) for b in bins_prim]

    # =========================================================================
    # PRIMARY RUN (B=10,000 Permutations)
    # =========================================================================
    print("\n" + "="*80)
    print("RUNNING PRIMARY PROXIMITY CALCULATIONS (B=10,000 Permutations)")
    print("="*80)

    S_mapped_prim = set(cad_primary_genes).intersection(set(G_primary.nodes()))
    S_bin_counts_prim = get_bin_counts(S_mapped_prim, n2b_prim)
    print(f"[*] Mapped CAD disease genes in Primary LCC: {len(S_mapped_prim)}/{len(cad_primary_genes)}")

    query_definitions = []
    for sym in sorted(panel_set):
        query_definitions.append(("Level A: Individual", sym, [sym]))
    for k, v in query_cfg["level_B_drug_specific"].items():
        query_definitions.append((f"Level B: Drug-specific ({k})", k, v["genes"]))
        dt_comp = v.get("direct_targets_comparator", [])
        if dt_comp:
            query_definitions.append((f"Level B: Direct Target Comparator ({k})", f"{k}_direct_target", dt_comp))
    query_definitions.append(("Level C: Exploratory Pooled Panel", "all_6_pharmacogenes", list(panel_set)))

    primary_records = []
    for level_label, q_id, q_genes in query_definitions:
        T_mapped = [g for g in q_genes if g in G_primary]
        if not T_mapped:
            continue
        overlap = set(T_mapped).intersection(S_mapped_prim)
        d_obs = closest_distance_fast(adj_primary, T_mapped, S_mapped_prim)

        T_bin_counts = get_bin_counts(T_mapped, n2b_prim)

        print(f"\n[*] Evaluating {q_id} (Size={len(T_mapped)})...")
        print(f"    Observed d_c = {d_obs:.4f} | Overlap: {sorted(list(overlap))}")

        # 10,000 permutations (Panel-excluded null)
        t0 = time.time()
        null_10k = run_permutations_fast(adj_primary, T_bin_counts, S_bin_counts_prim, bin_pools_ex, N_PERM_PRIMARY)
        elapsed_10k = time.time() - t0

        mu_10k = float(np.mean(null_10k))
        sd_10k = float(np.std(null_10k))
        zc_10k = float((d_obs - mu_10k) / sd_10k) if sd_10k > 0 else 0.0
        p_emp_10k = float((1.0 + np.sum(null_10k <= d_obs)) / (N_PERM_PRIMARY + 1.0))
        print(f"    [10k Permutations] mu={mu_10k:.4f}, sd={sd_10k:.4f}, z={zc_10k:.4f}, p_emp={p_emp_10k:.6f} ({elapsed_10k:.1f}s)")

        # 1,000 permutations (Guney-style baseline, without panel exclusion)
        null_1k_guney = run_permutations_fast(adj_primary, T_bin_counts, S_bin_counts_prim, bin_pools_inc, N_PERM_GUNEY)
        mu_1k = float(np.mean(null_1k_guney))
        sd_1k = float(np.std(null_1k_guney))
        zc_1k = float((d_obs - mu_1k) / sd_1k) if sd_1k > 0 else 0.0
        p_emp_1k = float((1.0 + np.sum(null_1k_guney <= d_obs)) / (N_PERM_GUNEY + 1.0))

        # Overlap-excluded sensitivity
        S_no_overlap = S_mapped_prim - set(T_mapped)
        d_no_overlap = closest_distance_fast(adj_primary, T_mapped, S_no_overlap) if S_no_overlap else np.nan

        primary_records.append({
            "query_level": level_label,
            "query_id": q_id,
            "query_genes": ";".join(q_genes),
            "mapped_query_genes": ";".join(T_mapped),
            "n_query_genes": len(T_mapped),
            "n_disease_genes": len(S_mapped_prim),
            "direct_overlap_genes": ";".join(sorted(overlap)),
            "n_direct_overlap": len(overlap),
            "observed_distance": round(d_obs, 4),
            "overlap_excluded_distance": round(d_no_overlap, 4) if not np.isnan(d_no_overlap) else "NA",
            "null_mean_10k": round(mu_10k, 4),
            "null_sd_10k": round(sd_10k, 4),
            "z_score_10k": round(zc_10k, 4),
            "p_empirical_10k": round(p_emp_10k, 6),
            "null_mean_1k_guney": round(mu_1k, 4),
            "null_sd_1k_guney": round(sd_1k, 4),
            "z_score_1k_guney": round(zc_1k, 4),
            "p_empirical_1k_guney": round(p_emp_1k, 6),
            "permutations_run": N_PERM_PRIMARY,
            "random_seed": SEED,
            "benchmark_runtime_sec": round(elapsed_10k, 2)
        })

    p_vals_prim = [r["p_empirical_10k"] for r in primary_records]
    q_vals_prim = benjamini_hochberg(p_vals_prim)
    for r, q in zip(primary_records, q_vals_prim):
        r["fdr_bh_qval_10k"] = round(q, 6)

    df_primary = pd.DataFrame(primary_records)
    df_primary.to_csv(PRIMARY_OUT, index=False)
    print(f"\n[+] Saved primary proximity results to {PRIMARY_OUT}")

    # =========================================================================
    # SENSITIVITY ANALYSIS MATRIX (Tasks 5 & 7, 1k Permutations)
    # =========================================================================
    print("\n" + "="*80)
    print("RUNNING MULTI-DIMENSIONAL SENSITIVITY MATRIX (1,000 Permutations per Cell)")
    print("="*80)

    sens_records = []

    # 1. Network Sensitivity: TM-excluded vs. Physical
    net_configs = [
        ("TM_Excluded_Functional", G_tm, adj_tm),
        ("Physical_Subnetwork", G_phys, adj_phys)
    ]
    for net_name, G_t, adj_t in net_configs:
        bins_t, n2b_t = get_degree_bins(G_t, min_bin_size=MIN_BIN_SIZE)
        pools_t = [[n for n in b[2] if n not in panel_set] for b in bins_t]
        S_t = set(cad_primary_genes).intersection(set(G_t.nodes()))
        if not S_t:
            continue
        S_cnt_t = get_bin_counts(S_t, n2b_t)

        for level_label, q_id, q_genes in query_definitions:
            T_t = [g for g in q_genes if g in G_t]
            if not T_t:
                continue
            d_o = closest_distance_fast(adj_t, T_t, S_t)
            T_cnt_t = get_bin_counts(T_t, n2b_t)
            null_t = run_permutations_fast(adj_t, T_cnt_t, S_cnt_t, pools_t, N_PERM_GUNEY)
            mu_t, sd_t = float(np.mean(null_t)), float(np.std(null_t))
            zc_t = float((d_o - mu_t) / sd_t) if sd_t > 0 else 0.0
            p_t = float((1.0 + np.sum(null_t <= d_o)) / (N_PERM_GUNEY + 1.0))

            sens_records.append({
                "sensitivity_dimension": "Interactome Representation",
                "condition": net_name,
                "query_id": q_id,
                "mapped_genes": ";".join(T_t),
                "n_query": len(T_t),
                "n_disease": len(S_t),
                "d_obs": round(d_o, 4),
                "null_mean": round(mu_t, 4),
                "z_score": round(zc_t, 4),
                "p_empirical": round(p_t, 6)
            })

    # 2. Disease Concept Sensitivity: MI vs. IHD
    secondary_dis = [
        ("Myocardial_Infarction_Direct", df_mi_dir),
        ("Myocardial_Ischemia_Direct", df_ihd_dir)
    ]
    for d_label, df_sec in secondary_dis:
        s_sec = set(df_sec[df_sec["association_score"] >= 0.10]["gene_symbol"]).intersection(set(G_primary.nodes()))
        if not s_sec:
            continue
        s_sec_cnt = get_bin_counts(s_sec, n2b_prim)
        for level_label, q_id, q_genes in query_definitions:
            T_m = [g for g in q_genes if g in G_primary]
            if not T_m:
                continue
            d_o = closest_distance_fast(adj_primary, T_m, s_sec)
            T_cnt = get_bin_counts(T_m, n2b_prim)
            null_t = run_permutations_fast(adj_primary, T_cnt, s_sec_cnt, bin_pools_ex, N_PERM_GUNEY)
            mu_t, sd_t = float(np.mean(null_t)), float(np.std(null_t))
            zc_t = float((d_o - mu_t) / sd_t) if sd_t > 0 else 0.0
            p_t = float((1.0 + np.sum(null_t <= d_o)) / (N_PERM_GUNEY + 1.0))

            sens_records.append({
                "sensitivity_dimension": "Disease Concept",
                "condition": d_label,
                "query_id": q_id,
                "mapped_genes": ";".join(T_m),
                "n_query": len(T_m),
                "n_disease": len(s_sec),
                "d_obs": round(d_o, 4),
                "null_mean": round(mu_t, 4),
                "z_score": round(zc_t, 4),
                "p_empirical": round(p_t, 6)
            })

    # 3. Direct vs. Indirect CAD Sensitivity
    s_indir = set(df_cad_indir[df_cad_indir["association_score"] >= 0.10]["gene_symbol"]).intersection(set(G_primary.nodes()))
    s_indir_cnt = get_bin_counts(s_indir, n2b_prim)
    for level_label, q_id, q_genes in query_definitions:
        T_m = [g for g in q_genes if g in G_primary]
        if not T_m:
            continue
        d_o = closest_distance_fast(adj_primary, T_m, s_indir)
        T_cnt = get_bin_counts(T_m, n2b_prim)
        null_t = run_permutations_fast(adj_primary, T_cnt, s_indir_cnt, bin_pools_ex, N_PERM_GUNEY)
        mu_t, sd_t = float(np.mean(null_t)), float(np.std(null_t))
        zc_t = float((d_o - mu_t) / sd_t) if sd_t > 0 else 0.0
        p_t = float((1.0 + np.sum(null_t <= d_o)) / (N_PERM_GUNEY + 1.0))

        sens_records.append({
            "sensitivity_dimension": "Ontology Evidence (Direct vs Indirect)",
            "condition": "CAD_Indirect_Descendant_Inclusive",
            "query_id": q_id,
            "mapped_genes": ";".join(T_m),
            "n_query": len(T_m),
            "n_disease": len(s_indir),
            "d_obs": round(d_o, 4),
            "null_mean": round(mu_t, 4),
            "z_score": round(zc_t, 4),
            "p_empirical": round(p_t, 6)
        })

    # 4. Association Score Cutoff Sensitivity on Direct CAD (0.05, 0.20, 0.30)
    for cut in [0.05, 0.20, 0.30]:
        s_cut = set(df_cad_dir[df_cad_dir["association_score"] >= cut]["gene_symbol"]).intersection(set(G_primary.nodes()))
        s_cut_cnt = get_bin_counts(s_cut, n2b_prim)
        for level_label, q_id, q_genes in query_definitions:
            T_m = [g for g in q_genes if g in G_primary]
            if not T_m:
                continue
            d_o = closest_distance_fast(adj_primary, T_m, s_cut)
            T_cnt = get_bin_counts(T_m, n2b_prim)
            null_t = run_permutations_fast(adj_primary, T_cnt, s_cut_cnt, bin_pools_ex, N_PERM_GUNEY)
            mu_t, sd_t = float(np.mean(null_t)), float(np.std(null_t))
            zc_t = float((d_o - mu_t) / sd_t) if sd_t > 0 else 0.0
            p_t = float((1.0 + np.sum(null_t <= d_o)) / (N_PERM_GUNEY + 1.0))

            sens_records.append({
                "sensitivity_dimension": "Score Cutoff",
                "condition": f"CAD_Direct_Cutoff_{cut}",
                "query_id": q_id,
                "mapped_genes": ";".join(T_m),
                "n_query": len(T_m),
                "n_disease": len(s_cut),
                "d_obs": round(d_o, 4),
                "null_mean": round(mu_t, 4),
                "z_score": round(zc_t, 4),
                "p_empirical": round(p_t, 6)
            })

    df_sens = pd.DataFrame(sens_records)
    for dim, group in df_sens.groupby("sensitivity_dimension"):
        q_sub = benjamini_hochberg(group["p_empirical"].tolist())
        df_sens.loc[group.index, "fdr_bh_qval"] = [round(x, 6) for x in q_sub]

    df_sens.to_csv(SENSITIVITY_OUT, index=False)
    print(f"\n[+] Saved complete sensitivity matrix to {SENSITIVITY_OUT} ({len(df_sens)} rows)")

    print("\n" + "="*80)
    print("SUMMARY OF PRIMARY PROXIMITY RESULTS (B=10,000 PERMUTATIONS)")
    print("="*80)
    cols_show = ["query_id", "n_query_genes", "direct_overlap_genes", "observed_distance", "null_mean_10k", "z_score_10k", "p_empirical_10k", "fdr_bh_qval_10k"]
    print(df_primary[cols_show].to_string(index=False))

if __name__ == "__main__":
    main()
