import pandas as pd

df = pd.read_csv("results/tables/variants_annotated_week3_backup.csv")
cv = pd.read_csv("results/tables/annotated_clinvar_variants.tsv", sep="\t")

# Build map from ClinVar for any missing coordinates by rsid
cv_valid = cv.dropna(subset=["chrom", "pos_b38", "ref", "alt"]).copy()
cv_map = cv_valid.set_index("rsid")[["chrom", "pos_b38", "ref", "alt"]].to_dict("index")

for idx, r in df.iterrows():
    rs = str(r.get("rsid", ""))
    if pd.isna(r.get("pos")) or str(r.get("pos")) in ["nan", "unmapped", "."]:
        if rs in cv_map:
            df.at[idx, "chrom"] = cv_map[rs]["chrom"]
            df.at[idx, "pos"] = cv_map[rs]["pos_b38"]
            df.at[idx, "ref"] = cv_map[rs]["ref"]
            df.at[idx, "alt"] = cv_map[rs]["alt"]

# Also include ClinVar variants directly to ensure all panel genes are represented
cv_extra = cv_valid[~cv_valid["rsid"].isin(df["rsid"])].copy()
if not cv_extra.empty:
    cv_extra["gene"] = cv_extra["gene_symbol"]
    cv_extra["pos"] = cv_extra["pos_b38"]
    cv_extra["variant_key"] = cv_extra["gene"] + "_" + cv_extra["rsid"]
    df = pd.concat([df, cv_extra], ignore_index=True)

# Drop any without full coordinates
vep = df.dropna(subset=["chrom", "pos", "ref", "alt"]).copy()
vep = vep[~vep["pos"].astype(str).isin([".", "unmapped", "nan"])].copy()

n_dropped = len(df) - len(vep)
print(f"Total candidate variants evaluated: {len(df)}")
print(f"Rows lacking GRCh38 coordinates (cannot go to VEP): {n_dropped}")
print(f"Variants ready for VEP: {len(vep)}")

# Format chrom and pos as pure integers
vep["chrom"] = vep["chrom"].astype(str).str.replace(r"^chr", "", regex=True).str.split(".").str[0]
vep["pos"] = vep["pos"].astype(str).str.split(".").str[0]

# Standard 8-column VCF structure
out = vep[["chrom", "pos", "variant_key", "ref", "alt"]].drop_duplicates(subset=["variant_key"]).copy()
out["QUAL"] = "."
out["FILTER"] = "."
out["INFO"] = "."

out.to_csv("data/processed/vep_input_raw.tsv", sep="\t", header=False, index=False)
