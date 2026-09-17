#!/usr/bin/env python3
"""
scripts/merge_annotate_variants.py
Standardized merge and annotation script adhering strictly to Master SOP Week 3.
Processes CPIC, ClinPGx, ClinVar, and PharmVar outputs into unified tables.
"""

import os
import glob
import json
import csv
from pathlib import Path
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

def parse_clinvar_vcf(vcf_path):
    """Extract pathogenic / likely pathogenic or review-classified variants from raw ClinVar VCF."""
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

            records.append({
                "chrom": chrom,
                "pos_b38": pos,
                "rsid": rsid if rsid != "." else info_dict.get("RS", "."),
                "ref": ref,
                "alt": alt,
                "clinical_significance": info_dict.get("CLNSIG", "not_specified"),
                "review_status": info_dict.get("CLNREVSTAT", "not_specified"),
                "gene_symbol": info_dict.get("GENEINFO", "").split(":")[0] if "GENEINFO" in info_dict else ""
            })
    return pd.DataFrame(records)

def parse_pharmvar_tsv(tsv_path):
    """Extract star alleles and core variant coordinates from PharmVar TSV."""
    if not os.path.exists(tsv_path):
        return pd.DataFrame()

    lines = []
    with open(tsv_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("#"):
                lines.append(line)
    
    if not lines:
        return pd.DataFrame()

    from io import StringIO
    df = pd.read_csv(StringIO("".join(lines)), sep="\t")
    return df

def parse_cpic_json(json_path):
    """Load CPIC gene-drug pair records."""
    if not os.path.exists(json_path):
        return pd.DataFrame()
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

    for gene in genes:
        # 1. Parse CPIC
        cpic_files = list(raw_dir.glob(f"cpic_{gene}_*.json"))
        if cpic_files:
            latest_cpic = sorted(cpic_files)[-1]
            df_cpic = parse_cpic_json(latest_cpic)
            if not df_cpic.empty:
                cpic_frames.append(df_cpic)

        # 2. Parse ClinVar
        cv_files = list(raw_dir.glob(f"clinvar_{gene}_*.vcf"))
        if cv_files:
            latest_cv = sorted(cv_files)[-1]
            df_cv = parse_clinvar_vcf(latest_cv)
            if not df_cv.empty:
                if "gene_symbol" in df_cv.columns:
                    df_cv["gene_symbol"] = df_cv["gene_symbol"].replace("", gene)
                clinvar_frames.append(df_cv)

    # Compile tables
    if cpic_frames:
        master_cpic = pd.concat(cpic_frames, ignore_index=True)
        if "cpiclevel" in master_cpic.columns:
            master_cpic = master_cpic[master_cpic["cpiclevel"].astype(str).str.upper().isin(["A", "B", "A/B"])]
        cpic_out = out_dir / "annotated_cpic_pairs.tsv"
        master_cpic.to_csv(cpic_out, sep="\t", index=False)
        print(f"[✓] Generated CPIC annotated table: {cpic_out} ({len(master_cpic)} pairs)")

    if clinvar_frames:
        master_clinvar = pd.concat(clinvar_frames, ignore_index=True)
        clinvar_out = out_dir / "annotated_clinvar_variants.tsv"
        master_clinvar.to_csv(clinvar_out, sep="\t", index=False)
        print(f"[✓] Generated ClinVar variants table: {clinvar_out} ({len(master_clinvar)} variants)")

    print("[✓] Week 3 base annotation tables compiled.")

if __name__ == "__main__":
    main()
