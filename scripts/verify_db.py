import sqlite3
from pathlib import Path

DB_PATH = Path("data/database/cvd_pgx.db")

def verify_db():
    if not DB_PATH.exists():
        print(f"Error: Database file not found at {DB_PATH}. Run build_database.py first.")
        return

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    print("=== 1. Table Record Counts ===")
    tables = ['genes', 'cpic_guidelines', 'pharmvar_alleles', 'clinvar_variants']
    for t in tables:
        count = cursor.execute(f"SELECT COUNT(*) FROM {t};").fetchone()[0]
        print(f"Table '{t}': {count} rows")

    print("\n=== 2. Relational Benchmark Join (CYP2C19 & Clopidogrel) ===")
    query = """
    SELECT 
        g.gene_symbol,
        c.drug_name,
        c.guideline_name,
        COUNT(DISTINCT p.id) as star_allele_records,
        COUNT(DISTINCT v.id) as clinvar_variant_records
    FROM genes g
    JOIN cpic_guidelines c ON g.gene_symbol = c.gene_symbol
    LEFT JOIN pharmvar_alleles p ON g.gene_symbol = p.gene_symbol
    LEFT JOIN clinvar_variants v ON g.gene_symbol = v.gene_symbol
    WHERE g.gene_symbol = 'CYP2C19' AND c.drug_name = 'clopidogrel'
    GROUP BY g.gene_symbol, c.drug_name;
    """
    row = cursor.execute(query).fetchone()
    if row:
        print(f"Gene: {row[0]} | Drug: {row[1]} | Guideline: {row[2]}")
        print(f"Linked PharmVar Allele Records: {row[3]} | Linked ClinVar Variants: {row[4]}")
    else:
        print("No matching record found for CYP2C19 and clopidogrel.")

    conn.close()

if __name__ == "__main__":
    verify_db()
