import sqlite3
import pandas as pd
from pathlib import Path
from datetime import date
import logging

PROC_DIR = Path("data/processed")
DB_DIR = Path("data/database")
DB_PATH = DB_DIR / "cvd_pgx.db"

TODAY = date.today().strftime("%Y%m%d")

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

def build_sqlite_db():
    DB_DIR.mkdir(parents=True, exist_ok=True)
    if DB_PATH.exists():
        DB_PATH.unlink()  # Clean rebuild
    
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute("PRAGMA foreign_keys = ON;")

    logging.info(f"Creating SQLite database at {DB_PATH}...")

    # 1. Create Core Genes Table
    cursor.execute("""
    CREATE TABLE genes (
        gene_symbol TEXT PRIMARY KEY,
        gene_id TEXT,
        gene_name TEXT,
        allele_file TEXT,
        allele_function_source TEXT
    );
    """)

    # 2. Create CPIC Guidelines Table
    cursor.execute("""
    CREATE TABLE cpic_guidelines (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        gene_symbol TEXT NOT NULL,
        drug_name TEXT NOT NULL,
        guideline_name TEXT,
        lookup_method TEXT,
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    );
    """)

    # 3. Create PharmVar Alleles Table
    cursor.execute("""
    CREATE TABLE pharmvar_alleles (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        gene_symbol TEXT NOT NULL,
        star_allele TEXT,
        rsid TEXT,
        ref_seq TEXT,
        pos_grch38 TEXT,
        variant_type TEXT,
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    );
    """)

    # 4. Create ClinVar Variants Table
    cursor.execute("""
    CREATE TABLE clinvar_variants (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        gene_symbol TEXT NOT NULL,
        chrom TEXT NOT NULL,
        pos_grch38 TEXT NOT NULL,
        rsid TEXT,
        ref TEXT,
        alt TEXT,
        clinical_significance TEXT,
        FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
    );
    """)

    conn.commit()

    # Load Processed TSVs into Database
    clinpgx_df = pd.read_csv(PROC_DIR / f"clinpgx_normalized_{TODAY}.tsv", sep="\t").replace({"NA": None})
    cpic_df = pd.read_csv(PROC_DIR / f"cpic_normalized_{TODAY}.tsv", sep="\t").replace({"NA": None})
    pharmvar_df = pd.read_csv(PROC_DIR / f"pharmvar_normalized_{TODAY}.tsv", sep="\t").replace({"NA": None})
    clinvar_df = pd.read_csv(PROC_DIR / f"clinvar_normalized_{TODAY}.tsv", sep="\t").replace({"NA": None})

    # Insert Genes
    genes_data = clinpgx_df[['gene_symbol', 'gene_id', 'gene_name', 'allele_file', 'allele_function_source']].drop_duplicates()
    genes_data.to_sql('genes', conn, if_exists='append', index=False)

    # Insert CPIC
    cpic_data = cpic_df[['gene_symbol', 'drug_name', 'guideline_name', 'lookup_method']]
    cpic_data.to_sql('cpic_guidelines', conn, if_exists='append', index=False)

    # Insert PharmVar
    pv_data = pharmvar_df[['gene_symbol', 'star_allele', 'rsid', 'ref_seq', 'pos_grch38', 'variant_type']]
    pv_data.to_sql('pharmvar_alleles', conn, if_exists='append', index=False)

    # Insert ClinVar
    cv_data = clinvar_df[['gene_symbol', 'chrom', 'pos_grch38', 'rsid', 'ref', 'alt', 'clinical_significance']]
    cv_data.to_sql('clinvar_variants', conn, if_exists='append', index=False)

    # Create Performance Indexes
    logging.info("Creating database indexes...")
    cursor.execute("CREATE INDEX idx_cpic_gene ON cpic_guidelines(gene_symbol);")
    cursor.execute("CREATE INDEX idx_cpic_drug ON cpic_guidelines(drug_name);")
    cursor.execute("CREATE INDEX idx_pv_gene ON pharmvar_alleles(gene_symbol);")
    cursor.execute("CREATE INDEX idx_pv_rsid ON pharmvar_alleles(rsid);")
    cursor.execute("CREATE INDEX idx_cv_gene ON clinvar_variants(gene_symbol);")
    cursor.execute("CREATE INDEX idx_cv_rsid ON clinvar_variants(rsid);")

    conn.commit()
    conn.close()
    logging.info("Database construction complete successfully!")

if __name__ == "__main__":
    build_sqlite_db()
