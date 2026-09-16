import sqlite3
import json
from pathlib import Path

DB_PATH = Path("data/database/cvd_pgx.db")

class PGxQueryEngine:
    def __init__(self, db_path=DB_PATH):
        self.db_path = db_path

    def get_connection(self):
        return sqlite3.connect(self.db_path)

    def query_gene_profile(self, gene_symbol):
        """Fetch all guideline recommendations, PharmVar star alleles, and ClinVar variants for a gene."""
        conn = self.get_connection()
        cursor = conn.cursor()

        # Guidelines
        guidelines = cursor.execute(
            "SELECT drug_name, guideline_name FROM cpic_guidelines WHERE gene_symbol = ?;", (gene_symbol,)
        ).fetchall()

        # Star Alleles
        alleles = cursor.execute(
            "SELECT star_allele, rsid FROM pharmvar_alleles WHERE gene_symbol = ?;", (gene_symbol,)
        ).fetchall()

        # ClinVar Variants Count
        clinvar_count = cursor.execute(
            "SELECT COUNT(*) FROM clinvar_variants WHERE gene_symbol = ?;", (gene_symbol,)
        ).fetchone()[0]

        conn.close()

        return {
            "gene_symbol": gene_symbol,
            "guideline_count": len(guidelines),
            "drugs": list(set([g[0] for g in guidelines])),
            "star_alleles": [a[0] for a in alleles if a[0]],
            "total_clinvar_variants": clinvar_count
        }

if __name__ == "__main__":
    engine = PGxQueryEngine()
    profile = engine.query_gene_profile("CYP2C19")
    print("=== CYP2C19 Query Profile ===")
    print(json.dumps(profile, indent=2))
