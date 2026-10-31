#!/usr/bin/env python3
"""
scripts/extract_population_af.py
Week 5: Systematic Extraction of 1000G and gnomAD Population Allele Frequencies
Uses persistent HTTP session single GET requests with caching and rate limiting.
"""

import os
import sys
import time
import json
import requests
import pandas as pd
import numpy as np

CACHE_FILE = "data/processed/ensembl_pop_cache_single.json"
ENSEMBL_REST_URL = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

def parse_population_frequencies(pop_list, alt_allele):
    af_data = {
        "af_1000g_global": np.nan,
        "af_1000g_sas": np.nan,
        "af_1000g_eur": np.nan,
        "af_1000g_afr": np.nan,
        "af_1000g_eas": np.nan,
        "af_1000g_amr": np.nan,
        "af_gnomad_global": np.nan,
        "af_gnomad_sas": np.nan,
        "af_gnomad_nfe": np.nan,
        "af_gnomad_afr": np.nan,
        "af_gnomad_eas": np.nan,
        "af_gnomad_amr": np.nan
    }
    
    if not pop_list:
        return af_data

    alt_upper = str(alt_allele).strip().upper()

    for p in pop_list:
        p_allele = str(p.get("allele", "")).strip().upper()
        if p_allele != alt_upper:
            continue
            
        pop = str(p.get("population", ""))
        try:
            freq = float(p.get("frequency", np.nan))
        except (ValueError, TypeError):
            continue

        # 1000 Genomes Phase 3
        if pop == "1000GENOMES:phase_3:ALL":
            af_data["af_1000g_global"] = freq
        elif pop == "1000GENOMES:phase_3:SAS":
            af_data["af_1000g_sas"] = freq
        elif pop == "1000GENOMES:phase_3:EUR":
            af_data["af_1000g_eur"] = freq
        elif pop == "1000GENOMES:phase_3:AFR":
            af_data["af_1000g_afr"] = freq
        elif pop == "1000GENOMES:phase_3:EAS":
            af_data["af_1000g_eas"] = freq
        elif pop == "1000GENOMES:phase_3:AMR":
            af_data["af_1000g_amr"] = freq
            
        # gnomAD genomes / exomes
        elif pop in ["gnomADg:ALL", "gnomADe:ALL"] and np.isnan(af_data["af_gnomad_global"]):
            af_data["af_gnomad_global"] = freq
        elif pop in ["gnomADg:sas", "gnomADe:sas"] and np.isnan(af_data["af_gnomad_sas"]):
            af_data["af_gnomad_sas"] = freq
        elif pop in ["gnomADg:nfe", "gnomADe:nfe"] and np.isnan(af_data["af_gnomad_nfe"]):
            af_data["af_gnomad_nfe"] = freq
        elif pop in ["gnomADg:afr", "gnomADe:afr"] and np.isnan(af_data["af_gnomad_afr"]):
            af_data["af_gnomad_afr"] = freq
        elif pop in ["gnomADg:eas", "gnomADe:eas"] and np.isnan(af_data["af_gnomad_eas"]):
            af_data["af_gnomad_eas"] = freq
        elif pop in ["gnomADg:amr", "gnomADe:amr"] and np.isnan(af_data["af_gnomad_amr"]):
            af_data["af_gnomad_amr"] = freq

    return af_data

def categorize_rarity(af):
    if pd.isna(af):
        return "Unknown"
    if af >= 0.05:
        return "Common"
    elif af >= 0.01:
        return "Low_Frequency"
    elif af >= 0.001:
        return "Rare"
    else:
        return "Ultra_Rare"

def main():
    input_path = "results/tables/variants_annotated.csv"
    output_path = "results/tables/variants_annotated_w5.csv"
    intermediate_path = "results/tables/variants_population_af.csv"
    os.makedirs("data/processed", exist_ok=True)

    print(f"[*] Reading annotated catalog: {input_path}")
    df = pd.read_csv(input_path, low_memory=False)
    n_variants = len(df)

    valid_rsids = [str(rs).strip() for rs in df["rsid"].dropna().unique() if str(rs).strip().startswith("rs")]
    print(f"[*] Total variants: {n_variants}, Unique rsIDs to query: {len(valid_rsids)}")

    # Load existing cache if available
    cache = {}
    if os.path.exists(CACHE_FILE):
        try:
            with open(CACHE_FILE, "r") as f:
                cache = json.load(f)
            print(f"[*] Loaded {len(cache)} existing cached payloads from {CACHE_FILE}")
        except Exception:
            cache = {}

    to_query = [rs for rs in valid_rsids if rs not in cache]
    print(f"[*] Remaining rsIDs to fetch: {len(to_query)}")

    session = requests.Session()
    session.headers.update(HEADERS)

    for idx, rs in enumerate(to_query):
        url = f"{ENSEMBL_REST_URL}/variation/homo_sapiens/{rs}?pops=1"
        success = False
        for attempt in range(3):
            try:
                r = session.get(url, timeout=12)
                if r.status_code == 429:
                    wait_time = int(r.headers.get("Retry-After", 2))
                    time.sleep(wait_time)
                    continue
                if r.ok:
                    cache[rs] = r.json().get("populations", [])
                    success = True
                    break
                elif r.status_code == 404:
                    cache[rs] = []
                    success = True
                    break
            except Exception:
                time.sleep(0.5)

        if not success:
            cache[rs] = []

        # Flush cache to disk every 25 records
        if (idx + 1) % 25 == 0 or idx == len(to_query) - 1:
            with open(CACHE_FILE, "w") as f:
                json.dump(cache, f)
            print(f"[*] Progress: {idx + 1}/{len(to_query)} rsIDs processed ({len(cache)} cached)...")

        time.sleep(0.08) # ~12 requests/sec (safely under Ensembl 15 req/sec rate limit)

    with open(CACHE_FILE, "w") as f:
        json.dump(cache, f)

    print(f"\n[+] Fetch complete. Total cached rsID populations: {len(cache)}")

    # Parse frequencies
    parsed_rows = []
    for _, row in df.iterrows():
        rs = str(row.get("rsid", "")).strip()
        alt = str(row.get("alt", "")).strip()
        var_key = row["variant_key"]
        
        pop_list = cache.get(rs, [])
        af_dict = parse_population_frequencies(pop_list, alt)
        af_dict["variant_key"] = var_key
        parsed_rows.append(af_dict)

    df_af = pd.DataFrame(parsed_rows)
    df_af.to_csv(intermediate_path, index=False)
    print(f"[+] Saved population frequency matrix to: {intermediate_path}")

    # Merge maintaining exact row invariant
    df_merged = pd.merge(df, df_af, on="variant_key", how="left")

    # Composite frequencies
    df_merged["af_global"] = df_merged["af_gnomad_global"].combine_first(df_merged["af_1000g_global"])
    df_merged["af_sas"] = df_merged["af_gnomad_sas"].combine_first(df_merged["af_1000g_sas"])

    # Population delta and enrichment metrics
    df_merged["delta_af_sas"] = df_merged["af_sas"] - df_merged["af_global"]
    eps = 1e-5
    df_merged["sas_enrichment_ratio"] = (df_merged["af_sas"].fillna(0) + eps) / (df_merged["af_global"].fillna(0) + eps)

    df_merged["rarity_tier_global"] = df_merged["af_global"].apply(categorize_rarity)
    df_merged["rarity_tier_sas"] = df_merged["af_sas"].apply(categorize_rarity)

    assert len(df_merged) == n_variants, f"Row count drift: expected {n_variants}, got {len(df_merged)}"

    df_merged.to_csv(output_path, index=False)
    print(f"\n[+] Successfully saved Week 5 dataset ({len(df_merged)} variants) to {output_path}")

    print("\n--- Population Frequency Summary ---")
    print(f"Variants with Global AF: {df_merged['af_global'].notna().sum()} / {len(df_merged)}")
    print(f"Variants with SAS AF:    {df_merged['af_sas'].notna().sum()} / {len(df_merged)}")
    print("\nSAS Rarity Breakdown:")
    print(df_merged["rarity_tier_sas"].value_counts())

if __name__ == "__main__":
    main()
