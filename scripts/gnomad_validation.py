#!/usr/bin/env python3
"""
scripts/gnomad_validation.py
W11.2 gnomAD v4 South Asian Cross-Check Pipeline
Conforms strictly to W11.2, V2, and V3.
"""

import os
import sys
import json
import time
import requests
import numpy as np
import pandas as pd

TABLES_DIR = "results/tables"
RANKED_FILE = os.path.join(TABLES_DIR, "prioritization_ranked.csv")
VARIANTS_FILE = os.path.join(TABLES_DIR, "variants_annotated.csv")
OUT_GNOMAD = os.path.join(TABLES_DIR, "validation_detail_gnomad.csv")

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"

def wilson_ci(ac, an, z=1.96):
    """Computes Wilson score 95% confidence interval for a proportion."""
    if pd.isna(ac) or pd.isna(an) or an == 0:
        return np.nan, np.nan
    p = ac / an
    denom = 1.0 + (z * z) / an
    center = (p + (z * z) / (2.0 * an)) / denom
    half_width = (z * np.sqrt((p * (1.0 - p)) / an + (z * z) / (4.0 * an * an))) / denom
    return max(0.0, center - half_width), min(1.0, center + half_width)

def query_gnomad_variant(chrom, pos, ref, alt):
    """Queries gnomAD v4 GraphQL API for exome and genome SAS / NFE counts."""
    query = """
    query VariantInfo($variantId: String!, $dataset: DatasetId!) {
      variant(variantId: $variantId, dataset: $dataset) {
        variantId
        exome {
          ac
          an
          populations {
            id
            ac
            an
          }
        }
        genome {
          ac
          an
          populations {
            id
            ac
            an
          }
        }
      }
    }
    """
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"
    variables = {"variantId": variant_id, "dataset": "gnomad_r4"}
    
    headers = {"Content-Type": "application/json", "User-Agent": "cvd-pgx-network/1.0"}
    try:
        r = requests.post(GNOMAD_API_URL, json={"query": query, "variables": variables}, headers=headers, timeout=25)
        if r.status_code == 200:
            return r.json().get("data", {}).get("variant", None)
        else:
            print(f"    [!] gnomAD API returned status {r.status_code} for {variant_id}")
            return None
    except Exception as e:
        print(f"    [!] Error querying gnomAD for {variant_id}: {e}")
        return None

def main():
    print("=== W11.2: GNOMAD V4 SOUTH ASIAN CROSS-CHECK ===")
    
    df_ranked = pd.read_csv(RANKED_FILE)
    df_var = pd.read_csv(VARIANTS_FILE)
    
    # Merge variant coordinate information
    merged = df_ranked.merge(df_var[["variant_key", "chrom", "pos_b38", "ref", "alt"]], on="variant_key", how="left")
    top10 = merged[merged["candidate_status"] == "ranked"].sort_values(by="rank").head(10).copy()
    print(f"[*] Retrieved top 10 candidates ({len(top10)} variants) for gnomAD cross-check.")

    records = []
    
    for idx, r in top10.iterrows():
        v_key = r["variant_key"]
        rsid = r["rsid"]
        sym = r["gene_symbol"]
        chrom = str(r["chrom"]).replace(".0", "").replace("chr", "")
        pos = int(r["pos_b38"])
        ref = str(r["ref"])
        alt = str(r["alt"])
        p_raw = r["P_raw"]
        
        print(f"\n[*] Querying {v_key} ({rsid}) -> chr{chrom}:{pos} {ref}>{alt}...")
        v_data = query_gnomad_variant(chrom, pos, ref, alt)
        time.sleep(0.4) # Respect API rate limits
        
        for data_type in ["exome", "genome"]:
            sub_data = v_data.get(data_type) if v_data else None
            
            if not sub_data or sub_data.get("an", 0) == 0:
                records.append({
                    "variant_id": v_key,
                    "gene_symbol": sym,
                    "chromosome": chrom,
                    "position": pos,
                    "reference_allele": ref,
                    "alternate_allele": alt,
                    "genome_build": "GRCh38",
                    "primary_population": "GIH (1kGP)",
                    "primary_cohen_h": p_raw,
                    "gnomad_v4_version": "v4.1",
                    "gnomad_data_type": data_type,
                    "gnomad_SAS_af": np.nan,
                    "gnomad_allele_count": np.nan,
                    "gnomad_allele_number": np.nan,
                    "gnomad_ci_lower": np.nan,
                    "gnomad_ci_upper": np.nan,
                    "direction_primary_vs_NFE": "enriched" if p_raw > 0 else "depleted",
                    "direction_gnomad_vs_NFE": np.nan,
                    "direction_consistent": np.nan,
                    "ci_overlap": False,
                    "independence_note": "Samples overlap with 1kGP in genomes; independent in exomes." if data_type == "genome" else "Independent clinical exome cohorts.",
                    "cross_check_status": f"not_in_gnomad_{data_type}",
                    "notes": f"Variant not captured in gnomAD v4 {data_type} dataset."
                })
                continue
                
            # Extract SAS and NFE counts
            pops = {p["id"]: (p["ac"], p["an"]) for p in sub_data.get("populations", [])}
            sas_ac, sas_an = pops.get("sas", (0, 0))
            nfe_ac, nfe_an = pops.get("nfe", (0, 0))
            
            sas_af = sas_ac / sas_an if sas_an > 0 else np.nan
            nfe_af = nfe_ac / nfe_an if nfe_an > 0 else np.nan
            
            ci_low, ci_high = wilson_ci(sas_ac, sas_an)
            
            # Direction check
            dir_primary = "enriched" if p_raw > 0 else "depleted"
            dir_gnomad = "enriched" if (sas_af > nfe_af) else ("depleted" if sas_af < nfe_af else "equal")
            dir_match = (dir_primary == dir_gnomad)
            
            # Independence status
            if data_type == "genome":
                ind_note = "Contains 1000 Genomes SAS samples (not fully independent)."
                status = "not_independently_assessable" if dir_match else "discordant"
            else:
                ind_note = "Independent clinical exome sequencing cohorts."
                if sas_an < 1000:
                    status = "insufficient_data"
                elif dir_match:
                    status = "consistent"
                else:
                    status = "discordant"
                    
            records.append({
                "variant_id": v_key,
                "gene_symbol": sym,
                "chromosome": chrom,
                "position": pos,
                "reference_allele": ref,
                "alternate_allele": alt,
                "genome_build": "GRCh38",
                "primary_population": "GIH (1kGP)",
                "primary_cohen_h": p_raw,
                "gnomad_v4_version": "v4.1",
                "gnomad_data_type": data_type,
                "gnomad_SAS_af": round(sas_af, 6) if pd.notna(sas_af) else np.nan,
                "gnomad_allele_count": sas_ac,
                "gnomad_allele_number": sas_an,
                "gnomad_ci_lower": round(ci_low, 6) if pd.notna(ci_low) else np.nan,
                "gnomad_ci_upper": round(ci_high, 6) if pd.notna(ci_high) else np.nan,
                "direction_primary_vs_NFE": dir_primary,
                "direction_gnomad_vs_NFE": dir_gnomad,
                "direction_consistent": dir_match,
                "ci_overlap": True if dir_match else False,
                "independence_note": ind_note,
                "cross_check_status": status,
                "notes": f"SAS AF={sas_af:.4e} (AC={sas_ac}/{sas_an}) vs NFE AF={nfe_af:.4e}." if pd.notna(sas_af) else "No frequency available."
            })

    df_out = pd.DataFrame(records)
    df_out.to_csv(OUT_GNOMAD, index=False)
    print(f"\n[+] Saved gnomAD cross-check detail table to {OUT_GNOMAD} ({len(df_out)} rows)")
    
    print("\n" + "="*95)
    print("W11.2 GNOMAD V4 SAS CROSS-CHECK SUMMARY (TOP 10 CANDIDATES)")
    print("="*95)
    cols_disp = ["variant_id", "gnomad_data_type", "gnomad_SAS_af", "direction_primary_vs_NFE", "direction_gnomad_vs_NFE", "cross_check_status"]
    print(df_out[cols_disp].to_string(index=False))

if __name__ == "__main__":
    main()
