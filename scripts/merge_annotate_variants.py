#!/usr/bin/env python3
"""
scripts/merge_annotate_variants.py
Master SOP Week 3 compliant variant merge & annotation pipeline.
Builds master_variant_annotations.tsv uniting PharmVar star-alleles,
ClinVar clinical significance, and CPIC drug-gene recommendations.
"""

import os
import re
from pathlib import Path
from io import StringIO
import pandas as pd
import yaml

def load_config():
    repo_root = Path(__file__).resolve().parent.parent
    possible_paths = [
        repo_root / "config" / "config.yaml",
        repo_root / "config.yaml"
    ]
    for path in possible_paths:
        if path.exists():
            with open(path, "r") as f:
                return yaml.safe_load(f), repo_root
    raise FileNotFoundError("config.yaml not found in ./config/ or root directory.")

def clean_rsid(val):
    if pd.isna(val) or val is None:
        return ""
    s = str(val).strip().lower()
    if s in [".", "none", "nan", ""]:
        return ""
    m = re.search(r"rs(\d+)", s)
    if m:
        return f"rs{m.group(1)}"
    if s.isdigit():
        return f"rs{s}"
    return s

def parse_clinvar_vcf(vcf_path):
    """Extract variants and clinical significance from ClinVar VCF."""
    records = []
    if not os.path.exists(vcf_path):
        return pd.DataFrame()

    with open(vcf_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, rsid, ref, alt, qual, filt, info = parts[:8]
            
            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v
                else:
                    info_dict[item] = True

            final_rs = rsid if rsid != "." else info_dict.get("RS", ".")
            
            records.append({
                "chrom": str(chrom).replace("chr", ""),
                "pos_b38": str(pos),
                "rsid": clean_rsid(final_rs),
                "ref": ref,
                "alt": alt,
                "clinical_significance": info_dict.get("CLNSIG", "not_specified"),
                "review_status": info_dict.get("CLNREVSTAT", "not_specified"),
                "gene_symbol": info_dict.get("GENEINFO", "").split(":")[0] if "GENEINFO" in info_dict else ""
            })
    return pd.DataFrame(records)

def parse_pharmvar_tsv(tsv_path, gene_name):
    """Extract star alleles and coordinates from PharmVar TSV."""
    if not os.path.exists(tsv_path):
        return pd.DataFrame()

    lines = []
    with open(tsv_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("#"):
                lines.append(line)
    if not lines:
        return pd.DataFrame()

    df = pd.read_csv(StringIO("".join(lines)), sep="\t")
    df["gene_symbol"] = gene_name
    df["rsid"] = df["rsID"].apply(clean_rsid) if "rsID" in df.columns else ""
    df["star_allele"] = df["Haplotype Name"] if "Haplotype Name" in df.columns else "none"
    df["variant_type"] = df["Type"] if "Type" in df.columns else "not_specified"
    df["ref"] = df["Reference Allele"] if "Reference Allele" in df.columns else ""
    df["alt"] = df["Variant Allele"] if "Variant Allele" in df.columns else ""
    df["pos_b37"] = df["Variant Start"].dropna().astype(str).str.split(".").str[0] if "Variant Start" in df.columns else ""
    
    # Return only defined variant rows
    return df[df["rsid"] != ""][["gene_symbol", "rsid", "star_allele", "variant_type", "ref", "alt", "pos_b37"]]

def parse_cpic_json(json_path):
    """Load CPIC gene-drug pair records."""
    if not os.path.exists(json_path):
        return pd.DataFrame()
    import json
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, list):
        return pd.DataFrame(data)
    elif isinstance(data, dict):
        return pd.DataFrame([data])
    return pd.DataFrame()

def main():
    config, repo_root = load_config()
    genes = config.get("gene_panel", ["CYP2C19", "CYP2C9", "VKORC1", "CYP4F2", "SLCO1B1", "ABCG2"])
    
    raw_dir = repo_root / "data" / "raw"
    out_dir = repo_root / "results" / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[*] Processing panel genes: {genes}")

    cpic_frames = []
    clinvar_frames = []
    pharmvar_frames = []

    for gene in genes:
        # CPIC
        cpic_files = list(raw_dir.glob(f"cpic_{gene}_*.json"))
        if cpic_files:
            df_c = parse_cpic_json(sorted(cpic_files)[-1])
            if not df_c.empty:
                cpic_frames.append(df_c)

        # ClinVar
        cv_files = list(raw_dir.glob(f"clinvar_{gene}_*.vcf"))
        if cv_files:
            df_v = parse_clinvar_vcf(sorted(cv_files)[-1])
            if not df_v.empty:
                if "gene_symbol" in df_v.columns:
                    df_v["gene_symbol"] = df_v["gene_symbol"].replace("", gene)
                clinvar_frames.append(df_v)

        # PharmVar
        pv_files = list(raw_dir.glob(f"pharmvar_{gene}_*.tsv"))
        if pv_files:
            df_p = parse_pharmvar_tsv(sorted(pv_files)[-1], gene)
            if not df_p.empty:
                pharmvar_frames.append(df_p)

    # 1. Output CPIC Pairs Table
    if cpic_frames:
        master_cpic = pd.concat(cpic_frames, ignore_index=True)
        if "cpiclevel" in master_cpic.columns:
            master_cpic = master_cpic[master_cpic["cpiclevel"].astype(str).str.upper().isin(["A", "B", "A/B"])]
        cpic_out = out_dir / "annotated_cpic_pairs.tsv"
        master_cpic.to_csv(cpic_out, sep="\t", index=False)
        print(f"[✓] Generated CPIC table: {cpic_out} ({len(master_cpic)} pairs)")

    # 2. Output ClinVar Variants Table
    df_clinvar = pd.concat(clinvar_frames, ignore_index=True) if clinvar_frames else pd.DataFrame()
    if not df_clinvar.empty:
        cv_out = out_dir / "annotated_clinvar_variants.tsv"
        df_clinvar.to_csv(cv_out, sep="\t", index=False)
        print(f"[✓] Generated ClinVar table: {cv_out} ({len(df_clinvar)} variants)")

    # 3. Build Unified Master Variant Annotations Table
    # Collapse PharmVar records by rsid and gene
    if pharmvar_frames:
        df_pv_all = pd.concat(pharmvar_frames, ignore_index=True)
        pv_agg = df_pv_all.groupby(["gene_symbol", "rsid"]).agg({
            "star_allele": lambda s: ";".join(sorted(set(str(x) for x in s if x != "none"))),
            "variant_type": lambda s: ";".join(sorted(set(str(x) for x in s if str(x) != "nan"))),
            "ref": "first",
            "alt": "first",
            "pos_b37": "first"
        }).reset_index()
    else:
        pv_agg = pd.DataFrame(columns=["gene_symbol", "rsid", "star_allele", "variant_type", "ref", "alt", "pos_b37"])

    # Create mapping lookups from ClinVar
    cv_sig_map = {}
    cv_rev_map = {}
    cv_pos_map = {}
    if not df_clinvar.empty:
        for _, row in df_clinvar[df_clinvar["rsid"] != ""].iterrows():
            key = (row["gene_symbol"], row["rsid"])
            cv_sig_map[key] = row["clinical_significance"]
            cv_rev_map[key] = row["review_status"]
            cv_pos_map[key] = row["pos_b38"]

    # Annotate PharmVar core alleles with ClinVar data
    pv_records = []
    for _, row in pv_agg.iterrows():
        key = (row["gene_symbol"], row["rsid"])
        pv_records.append({
            "gene_symbol": row["gene_symbol"],
            "rsid": row["rsid"],
            "star_allele": row["star_allele"] if row["star_allele"] else "none",
            "variant_type": row["variant_type"],
            "ref": row["ref"],
            "alt": row["alt"],
            "pos_b37": row["pos_b37"],
            "pos_b38": cv_pos_map.get(key, "unmapped"),
            "clinical_significance": cv_sig_map.get(key, "unclassified_in_clinvar"),
            "review_status": cv_rev_map.get(key, "none")
        })

    master_df = pd.DataFrame(pv_records)

    # Append non-PharmVar ClinVar variants for panel genes without PharmVar (e.g., VKORC1, CYP4F2, ABCG2)
    non_pv_genes = [g for g in genes if g not in ["CYP2C19", "CYP2C9", "SLCO1B1"]]
    if not df_clinvar.empty:
        extra_clinvar = df_clinvar[df_clinvar["gene_symbol"].isin(non_pv_genes)].copy()
        if not extra_clinvar.empty:
            extra_clinvar["star_allele"] = "none"
            extra_clinvar["pos_b37"] = "unmapped"
            rename_dict = {"chrom": "chrom", "pos_b38": "pos_b38"}
            extra_records = extra_clinvar[[
                "gene_symbol", "rsid", "star_allele", "ref", "alt", "pos_b37", "pos_b38", "clinical_significance", "review_status"
            ]].to_dict("records")
            master_df = pd.concat([master_df, pd.DataFrame(extra_records)], ignore_index=True)

    master_out = out_dir / "master_variant_annotations.tsv"
    master_df.to_csv(master_out, sep="\t", index=False)
    print(f"[✓] Generated Master Variant Annotations: {master_out} ({len(master_df)} unified variants)")
    print("[✓] Week 3 merge & annotation complete.")

if __name__ == "__main__":
    main()
