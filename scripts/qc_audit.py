import json
import pandas as pd
from pathlib import Path
from datetime import date

PROC_DIR = Path("data/processed")
QC_DIR = Path("data/qc_reports")
TODAY = date.today().strftime("%Y%m%d")

def run_qc_audit():
    audit_results = {
        "date": TODAY,
        "null_value_check": {},
        "format_checks": {},
        "benchmark_validation_CYP2C19": {}
    }
    
    # Load processed TSVs
    cpic = pd.read_csv(PROC_DIR / f"cpic_normalized_{TODAY}.tsv", sep="\t")
    clinpgx = pd.read_csv(PROC_DIR / f"clinpgx_normalized_{TODAY}.tsv", sep="\t")
    pharmvar = pd.read_csv(PROC_DIR / f"pharmvar_normalized_{TODAY}.tsv", sep="\t")
    clinvar = pd.read_csv(PROC_DIR / f"clinvar_normalized_{TODAY}.tsv", sep="\t")

    # 1. Null/NA Check
    for name, df in [("cpic", cpic), ("clinpgx", clinpgx), ("pharmvar", pharmvar), ("clinvar", clinvar)]:
        na_count = (df == "NA").sum().sum() + df.isna().sum().sum()
        audit_results["null_value_check"][name] = {
            "total_rows": len(df),
            "unhandled_na_cells": int(na_count)
        }

    # 2. Format Checks
    audit_results["format_checks"]["cpic_lowercase_drugs"] = cpic["drug_name"].str.islower().all()
    audit_results["format_checks"]["gene_symbols_uppercase"] = all(
        df["gene_symbol"].str.isupper().all() for df in [cpic, clinpgx, pharmvar, clinvar]
    )
    audit_results["format_checks"]["clinvar_chrom_valid"] = clinvar["chrom"].astype(str).str.contains(r"^(chr)?[0-9XY]+$").all()

    # 3. Benchmark Validation (CYP2C19)
    audit_results["benchmark_validation_CYP2C19"] = {
        "cpic_cyp2c19_records": int((cpic["gene_symbol"] == "CYP2C19").sum()),
        "clinpgx_cyp2c19_records": int((clinpgx["gene_symbol"] == "CYP2C19").sum()),
        "pharmvar_cyp2c19_records": int((pharmvar["gene_symbol"] == "CYP2C19").sum()),
        "clinvar_cyp2c19_variants": int((clinvar["gene_symbol"] == "CYP2C19").sum()),
        "clopidogrel_guideline_present": "clopidogrel" in cpic["drug_name"].values
    }

    # Save detailed audit report
    audit_path = QC_DIR / f"qc_audit_report_{TODAY}.json"
    with open(audit_path, "w") as fp:
        json.dump(audit_results, fp, indent=2)

    print(json.dumps(audit_results, indent=2))

if __name__ == "__main__":
    run_qc_audit()
