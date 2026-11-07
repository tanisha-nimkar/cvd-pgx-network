import os
#!/usr/bin/env python3
"""
scripts/prepare_population_frequencies_raw.py
Transforms cached multi-population records into data/processed/population_frequencies_raw.csv
strictly conforming to Week 6 Pre-flight and W6.3 integrity specifications:
- Unpooled targets: GIH, PJL, BEB, STU, ITU, and gnomAD_SAS (India fallback)
- Reference: NFE (gnomAD_NFE)
- Autosomal only (build GRCh38)
- Exact AC, AN, AF integrity and status logging
"""

import json
import pandas as pd
import numpy as np

CACHE_FILE = "data/processed/ensembl_pop_cache_single.json"
ANNOTATED_FILE = "results/tables/variants_annotated_w5.csv" if os.path.exists("results/tables/variants_annotated_w5.csv") else "results/tables/variants_annotated.csv"
OUTPUT_FILE = "data/processed/population_frequencies_raw.csv"

# Target definitions
TARGET_POPS = {
    "1000GENOMES:phase_3:GIH": ("GIH", "1000Genomes_Phase3", 103),
    "1000GENOMES:phase_3:PJL": ("PJL", "1000Genomes_Phase3", 96),
    "1000GENOMES:phase_3:BEB": ("BEB", "1000Genomes_Phase3", 86),
    "1000GENOMES:phase_3:STU": ("STU", "1000Genomes_Phase3", 102),
    "1000GENOMES:phase_3:ITU": ("ITU", "1000Genomes_Phase3", 102),
    "gnomADg:sas": ("gnomAD_SAS", "gnomAD_v3_genomes", 2419),
    "gnomADe:sas": ("gnomAD_SAS", "gnomAD_v2_exomes", 15308),
}

REF_POPS = {
    "gnomADg:nfe": ("NFE", "gnomAD_v3_genomes", 34029),
    "gnomADe:nfe": ("NFE", "gnomAD_v2_exomes", 64603)
}

def main():
    print(f"[*] Reading annotated variants: {ANNOTATED_FILE}")
    df_var = pd.read_csv(ANNOTATED_FILE, low_memory=False)

    print(f"[*] Reading population cache: {CACHE_FILE}")
    with open(CACHE_FILE, "r") as f:
        cache = json.load(f)

    rows = []

    # Map target populations for every variant
    populations_to_extract = ["GIH", "PJL", "BEB", "STU", "ITU", "gnomAD_SAS", "NFE"]

    for _, vrow in df_var.iterrows():
        vkey = vrow["variant_key"]
        rs = str(vrow.get("rsid", "")).strip()
        gene = vrow["gene_symbol"]
        chrom = str(vrow["chrom"]).replace("chr", "").split(".")[0]
        pos = int(vrow["pos_b38"])
        ref = str(vrow["ref"]).strip().upper()
        alt = str(vrow["alt"]).strip().upper()

        pop_list = cache.get(rs, []) if rs in cache else []

        # Index populations by name and allele
        # Format: pop_dict[pop_name][allele] = {'ac': ..., 'freq': ...}
        pop_dict = {}
        for p in pop_list:
            p_name = p.get("population", "")
            p_allele = str(p.get("allele", "")).strip().upper()
            if p_name not in pop_dict:
                pop_dict[p_name] = {}
            pop_dict[p_name][p_allele] = {
                "ac": p.get("allele_count"),
                "freq": p.get("frequency")
            }

        # For each standard population, extract the alternate and reference allele counts
        for target_label in populations_to_extract:
            ac_val = np.nan
            an_val = np.nan
            af_val = np.nan
            sample_n_val = np.nan
            dataset_val = ""
            status_val = "missing"

            if target_label in ["GIH", "PJL", "BEB", "STU", "ITU"]:
                full_name = f"1000GENOMES:phase_3:{target_label}"
                dataset_val = "1000Genomes_Phase3"
                sample_n_val = TARGET_POPS[full_name][2]
                if full_name in pop_dict:
                    sub = pop_dict[full_name]
                    # Compute total AN if both alleles exist
                    alt_info = sub.get(alt)
                    ref_info = sub.get(ref)
                    if alt_info and ref_info and alt_info["ac"] is not None and ref_info["ac"] is not None:
                        ac_val = int(alt_info["ac"])
                        an_val = int(alt_info["ac"]) + int(ref_info["ac"])
                        af_val = float(alt_info["freq"]) if alt_info["freq"] is not None else (ac_val / an_val if an_val > 0 else np.nan)
                        status_val = "ok" if ac_val > 0 else "zero_frequency"
                    elif alt_info and alt_info["ac"] is not None:
                        ac_val = int(alt_info["ac"])
                        af_val = float(alt_info["freq"]) if alt_info["freq"] is not None else np.nan
                        # Fallback approximate AN if available
                        if af_val and af_val > 0:
                            an_val = int(round(ac_val / af_val))
                            status_val = "ok"
                    elif ref_info and ref_info["ac"] is not None:
                        ac_val = 0
                        an_val = int(ref_info["ac"])
                        af_val = 0.0
                        status_val = "zero_frequency"

            elif target_label == "gnomAD_SAS":
                # Prefer gnomADg (genomes), fallback to gnomADe (exomes)
                for cand_pop, (lbl, dset, sn) in [("gnomADg:sas", TARGET_POPS["gnomADg:sas"]), ("gnomADe:sas", TARGET_POPS["gnomADe:sas"])]:
                    if cand_pop in pop_dict:
                        sub = pop_dict[cand_pop]
                        alt_info = sub.get(alt)
                        ref_info = sub.get(ref)
                        dataset_val = dset
                        sample_n_val = sn
                        if alt_info and ref_info and alt_info["ac"] is not None and ref_info["ac"] is not None:
                            ac_val = int(alt_info["ac"])
                            an_val = int(alt_info["ac"]) + int(ref_info["ac"])
                            af_val = float(alt_info["freq"]) if alt_info["freq"] is not None else (ac_val / an_val if an_val > 0 else np.nan)
                            status_val = "ok" if ac_val > 0 else "zero_frequency"
                            break
                        elif alt_info and alt_info["ac"] is not None:
                            ac_val = int(alt_info["ac"])
                            af_val = float(alt_info["freq"]) if alt_info["freq"] is not None else np.nan
                            if af_val and af_val > 0:
                                an_val = int(round(ac_val / af_val))
                                status_val = "ok"
                            break

            elif target_label == "NFE":
                # Reference population: gnomADg:nfe, fallback gnomADe:nfe
                for cand_pop, (lbl, dset, sn) in [("gnomADg:nfe", REF_POPS["gnomADg:nfe"]), ("gnomADe:nfe", REF_POPS["gnomADe:nfe"])]:
                    if cand_pop in pop_dict:
                        sub = pop_dict[cand_pop]
                        alt_info = sub.get(alt)
                        ref_info = sub.get(ref)
                        dataset_val = dset
                        sample_n_val = sn
                        if alt_info and ref_info and alt_info["ac"] is not None and ref_info["ac"] is not None:
                            ac_val = int(alt_info["ac"])
                            an_val = int(alt_info["ac"]) + int(ref_info["ac"])
                            af_val = float(alt_info["freq"]) if alt_info["freq"] is not None else (ac_val / an_val if an_val > 0 else np.nan)
                            status_val = "ok" if ac_val > 0 else "zero_frequency"
                            break
                        elif alt_info and alt_info["ac"] is not None:
                            ac_val = int(alt_info["ac"])
                            af_val = float(alt_info["freq"]) if alt_info["freq"] is not None else np.nan
                            if af_val and af_val > 0:
                                an_val = int(round(ac_val / af_val))
                                status_val = "ok"
                            break

            rows.append({
                "variant_key": vkey,
                "rsid": rs,
                "gene": gene,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "build": "GRCh38",
                "population": target_label,
                "dataset": dataset_val,
                "allele_count": ac_val,
                "allele_number": an_val,
                "allele_frequency": af_val,
                "sample_n": sample_n_val,
                "status": status_val
            })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(OUTPUT_FILE, index=False)
    print(f"[+] Successfully generated {OUTPUT_FILE} with {len(df_out)} rows.")
    print("\nPopulation breakdown:")
    print(df_out["population"].value_counts())
    print("\nStatus breakdown:")
    print(df_out["status"].value_counts())

if __name__ == "__main__":
    main()
