import json
import re
import pandas as pd
from pathlib import Path
from datetime import date
import logging

RAW_DIR = Path("data/raw")
PROC_DIR = Path("data/processed")
QC_DIR = Path("data/qc_reports")

PROC_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

TODAY = date.today().strftime("%Y%m%d")

# 1. Normalize CPIC Guidelines
def normalize_cpic():
    records = []
    for f in RAW_DIR.glob("cpic_*.json"):
        with open(f, "r") as fp:
            data = json.load(fp)
            for item in data:
                # Handle dictionary response structure
                if isinstance(item, dict):
                    records.append({
                        "source": "CPIC",
                        "gene_symbol": str(item.get("genesymbol", "NA")).upper(),
                        "drug_name": str(item.get("drugname", "NA")).lower(),
                        "guideline_name": item.get("guidelinename", "NA"),
                        "lookup_method": item.get("lookupmethod", "NA")
                    })
    df = pd.DataFrame(records).drop_duplicates().fillna("NA")
    out_path = PROC_DIR / f"cpic_normalized_{TODAY}.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    logging.info(f"Normalized CPIC: {len(df)} records -> {out_path}")
    return len(df)

# 2. Normalize ClinPGx Annotations
def normalize_clinpgx():
    records = []
    for f in RAW_DIR.glob("clinpgx_*.json"):
        with open(f, "r") as fp:
            data = json.load(fp)
            # Standardize list vs single object
            if isinstance(data, dict):
                data = [data]
            for item in data:
                records.append({
                    "source": "ClinPGx",
                    "gene_id": item.get("id", "NA"),
                    "gene_symbol": str(item.get("symbol", "NA")).upper(),
                    "gene_name": item.get("name", "NA"),
                    "allele_file": item.get("alleleFile", "NA"),
                    "allele_function_source": item.get("alleleFunctionSource", "NA")
                })
    df = pd.DataFrame(records).drop_duplicates().fillna("NA")
    out_path = PROC_DIR / f"clinpgx_normalized_{TODAY}.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    logging.info(f"Normalized ClinPGx: {len(df)} records -> {out_path}")
    return len(df)

# 3. Normalize PharmVar Star-Allele Tables
def normalize_pharmvar():
    records = []
    for f in RAW_DIR.glob("pharmvar_*.tsv"):
        gene_match = re.search(r"pharmvar_([A-Z0-9]+)_", f.name)
        gene = gene_match.group(1) if gene_match else "NA"
        
        try:
            df_pv = pd.read_csv(f, sep="\t")
            for _, row in df_pv.iterrows():
                star_allele = str(row.get("Gene", row.get("gene", "NA"))).strip()
                # Ensure star allele notation (e.g. *2)
                if not star_allele.startswith("*") and star_allele != "NA":
                    # Parse if embedded in symbol string
                    m = re.search(r"(\*[0-9]+(\.[0-9]+)?)", str(row))
                    star_allele = m.group(1) if m else star_allele

                records.append({
                    "source": "PharmVar",
                    "gene_symbol": gene.upper(),
                    "star_allele": star_allele,
                    "rsid": str(row.get("RSID", row.get("rsid", "NA"))),
                    "ref_seq": str(row.get("Reference", row.get("ref", "NA"))),
                    "pos_grch38": str(row.get("Position", row.get("pos", "NA"))),
                    "variant_type": str(row.get("Type", row.get("type", "NA")))
                })
        except Exception as e:
            logging.error(f"Error parsing PharmVar file {f}: {e}")

    df = pd.DataFrame(records).drop_duplicates().fillna("NA")
    out_path = PROC_DIR / f"pharmvar_normalized_{TODAY}.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    logging.info(f"Normalized PharmVar: {len(df)} records -> {out_path}")
    return len(df)

# 4. Normalize ClinVar VCF Slices
def normalize_clinvar():
    records = []
    for f in RAW_DIR.glob("clinvar_*.vcf"):
        gene_match = re.search(r"clinvar_([A-Z0-9]+)_", f.name)
        gene = gene_match.group(1) if gene_match else "NA"
        
        with open(f, "r") as fp:
            for line in fp:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue
                
                chrom, pos, id_, ref, alt, qual, filter_, info = parts[:8]
                
                # Extract RSID from INFO if ID is missing '.'
                rsid = id_
                if rsid == "." and "rsID=" in info:
                    m = re.search(r"rsID=([0-9]+)", info)
                    if m:
                        rsid = f"rs{m.group(1)}"
                
                # Extract Clinical Significance (CLNSIG)
                clnsig = "NA"
                if "CLNSIG=" in info:
                    m = re.search(r"CLNSIG=([^;]+)", info)
                    if m:
                        clnsig = m.group(1)

                records.append({
                    "source": "ClinVar",
                    "gene_symbol": gene.upper(),
                    "chrom": chrom,
                    "pos_grch38": pos,
                    "rsid": rsid,
                    "ref": ref,
                    "alt": alt,
                    "clinical_significance": clnsig
                })

    df = pd.DataFrame(records).drop_duplicates().fillna("NA")
    out_path = PROC_DIR / f"clinvar_normalized_{TODAY}.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    logging.info(f"Normalized ClinVar: {len(df)} records -> {out_path}")
    return len(df)

def main():
    logging.info("Starting Phase 2 Data Normalization Pipeline...")
    cpic_cnt = normalize_cpic()
    clinpgx_cnt = normalize_clinpgx()
    pharmvar_cnt = normalize_pharmvar()
    clinvar_cnt = normalize_clinvar()
    
    summary = {
        "date": TODAY,
        "normalized_counts": {
            "cpic": cpic_cnt,
            "clinpgx": clinpgx_cnt,
            "pharmvar": pharmvar_cnt,
            "clinvar": clinvar_cnt
        }
    }
    
    qc_summary_path = QC_DIR / f"qc_summary_{TODAY}.json"
    with open(qc_summary_path, "w") as fp:
        json.dump(summary, fp, indent=2)
    logging.info(f"Phase 2 Normalization Complete. QC Summary saved to {qc_summary_path}")

if __name__ == "__main__":
    main()
