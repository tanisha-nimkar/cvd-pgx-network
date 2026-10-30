#!/usr/bin/env python3
"""
Patch unmapped GRCh38 coordinates for PharmVar-derived variants in CYP2C19, CYP2C9, and SLCO1B1
using pyliftover from hg19 (GRCh37) to hg38 (GRCh38).
Generates the fully mapped N=465 vep_input_full.vcf.
"""

import pandas as pd
from pyliftover import LiftOver

def main():
    csv_path = "results/tables/variants_annotated_week3_backup.csv"
    df = pd.read_csv(csv_path, low_memory=False)
    print(f"[*] Loaded baseline variants: {len(df)}")

    # Ensure columns can accept string and numeric mixed values without lossy upcast errors
    for col in ["chrom", "pos_b38", "pos"]:
        if col in df.columns:
            df[col] = df[col].astype("object")
        else:
            df[col] = None
            df[col] = df[col].astype("object")

    lo = LiftOver('hg19', 'hg38')

    chrom_map = {
        "CYP2C19": "10",
        "CYP2C9": "10",
        "SLCO1B1": "12",
        "ABCG2": "4",
        "CYP4F2": "19",
        "VKORC1": "16"
    }

    patched_count = 0
    unresolved_count = 0
    vcf_records = []

    for idx, row in df.iterrows():
        gene = str(row["gene_symbol"]).strip()
        rsid = str(row["rsid"]).strip() if pd.notna(row["rsid"]) else "-"
        ref = str(row["ref"]).strip()
        alt = str(row["alt"]).strip()
        
        var_key = f"{gene}_{rsid}"
        if rsid in ["-", "", "nan", "None"]:
            var_key = f"{gene}_{idx}"

        pos_b38 = row.get("pos_b38")
        chrom = str(row.get("chrom", "")).strip()

        # Check if coordinate needs liftOver
        needs_lift = pd.isna(pos_b38) or str(pos_b38).strip() in ["unmapped", "nan", "", "."]

        if needs_lift:
            pos_b37 = row.get("pos_b37")
            chr_num = chrom_map.get(gene, "10")
            chr_str = f"chr{chr_num}"
            
            if pd.notna(pos_b37):
                pos_b37_int = int(pos_b37)
                res = lo.convert_coordinate(chr_str, pos_b37_int - 1)
                if res and len(res) > 0:
                    new_chr = str(res[0][0]).replace("chr", "")
                    new_pos = int(res[0][1]) + 1
                    df.at[idx, "chrom"] = str(new_chr)
                    df.at[idx, "pos_b38"] = int(new_pos)
                    df.at[idx, "pos"] = int(new_pos)
                    patched_count += 1
                    chrom = str(new_chr)
                    pos_b38 = new_pos
                else:
                    unresolved_count += 1
            else:
                unresolved_count += 1
        else:
            if not chrom or chrom in ["nan", "None", ""]:
                chrom = chrom_map.get(gene, "")
                df.at[idx, "chrom"] = str(chrom)
            pos_b38 = int(float(pos_b38))
            df.at[idx, "pos_b38"] = pos_b38
            df.at[idx, "pos"] = pos_b38

        if pd.notna(pos_b38) and str(pos_b38) not in ["unmapped", "nan", ""]:
            vcf_records.append({
                "CHROM": str(chrom),
                "POS": int(pos_b38),
                "ID": var_key,
                "REF": ref,
                "ALT": alt,
                "QUAL": ".",
                "FILTER": ".",
                "INFO": f"GENE={gene}"
            })

    print(f"[+] LiftOver complete: Patched {patched_count} coordinates, Unresolved: {unresolved_count}")
    
    # Save patched catalog maintaining exact 465 rows
    df.to_csv("results/tables/master_variant_annotations_b38_patched.tsv", sep="\t", index=False)
    print(f"[+] Wrote patched catalog ({len(df)} rows) to results/tables/master_variant_annotations_b38_patched.tsv")

    # Generate full VCF for VEP
    vcf_df = pd.DataFrame(vcf_records)
    vcf_df["POS"] = vcf_df["POS"].astype(int)
    vcf_df = vcf_df.sort_values(by=["CHROM", "POS"])
    
    vcf_path = "data/processed/vep_input_full.vcf"
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, r in vcf_df.iterrows():
            f.write(f"{r['CHROM']}\t{r['POS']}\t{r['ID']}\t{r['REF']}\t{r['ALT']}\t{r['QUAL']}\t{r['FILTER']}\t{r['INFO']}\n")

    print(f"[+] Wrote {len(vcf_df)} records to {vcf_path}")

if __name__ == "__main__":
    main()
