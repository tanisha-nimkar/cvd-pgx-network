import sqlite3
import networkx as nx
import json
from pathlib import Path
from datetime import date
import logging

DB_PATH = Path("data/database/cvd_pgx.db")
PROC_DIR = Path("data/processed")
TODAY = date.today().strftime("%Y%m%d")

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

def sanitize(val, default="NA"):
    """Helper to convert None or 'None' values to string defaults for GraphML compatibility."""
    if val is None or val == "None" or str(val).strip() == "":
        return default
    return str(val)

def build_pgx_network():
    if not DB_PATH.exists():
        logging.error(f"Database not found at {DB_PATH}. Run build_database.py first.")
        return

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    G = nx.MultiDiGraph()
    logging.info("Building Pharmacogenomic Network Graph...")

    # 1. Add Gene Nodes
    genes = cursor.execute("SELECT gene_symbol, gene_name FROM genes;").fetchall()
    for symbol, name in genes:
        G.add_node(symbol, node_type="Gene", name=sanitize(name, symbol))

    # 2. Add CPIC Drug Nodes and Gene-Drug Guideline Edges
    guidelines = cursor.execute("SELECT gene_symbol, drug_name, guideline_name FROM cpic_guidelines;").fetchall()
    for gene, drug, guideline in guidelines:
        if drug and drug != "None":
            clean_drug = sanitize(drug)
            G.add_node(clean_drug, node_type="Drug", name=clean_drug)
            G.add_edge(gene, clean_drug, key=f"cpic_{sanitize(guideline)}", edge_type="has_guideline", source="CPIC", guideline=sanitize(guideline))

    # 3. Add PharmVar Variant Nodes and Edges
    pv_alleles = cursor.execute("SELECT gene_symbol, star_allele, rsid FROM pharmvar_alleles;").fetchall()
    for gene, star, rsid in pv_alleles:
        node_id = star if star and star != "None" else rsid
        if node_id and node_id != "None":
            clean_node = sanitize(node_id)
            G.add_node(clean_node, node_type="Variant", rsid=sanitize(rsid), star_allele=sanitize(star))
            G.add_edge(clean_node, gene, key=f"pv_{clean_node}", edge_type="variant_of", source="PharmVar")

    # 4. Add ClinVar Variant Nodes and Edges
    cv_variants = cursor.execute("SELECT gene_symbol, rsid, clinical_significance FROM clinvar_variants;").fetchall()
    for gene, rsid, clnsig in cv_variants:
        if rsid and rsid != "NA" and rsid != "None":
            clean_rsid = sanitize(rsid)
            G.add_node(clean_rsid, node_type="Variant", rsid=clean_rsid, clinical_significance=sanitize(clnsig))
            G.add_edge(clean_rsid, gene, key=f"cv_{clean_rsid}", edge_type="associated_with", source="ClinVar", significance=sanitize(clnsig))

    conn.close()

    logging.info(f"Graph Construction Complete: {G.number_of_nodes()} Nodes, {G.number_of_edges()} Edges.")

    # Export Network Files
    graphml_path = PROC_DIR / f"cvd_pgx_network_{TODAY}.graphml"
    json_path = PROC_DIR / f"cvd_pgx_network_{TODAY}.json"

    nx.write_graphml(G, graphml_path)
    
    # Export Node-Link JSON format
    data = nx.node_link_data(G)
    with open(json_path, "w") as fp:
        json.dump(data, fp, indent=2)

    logging.info(f"Exported GraphML to {graphml_path}")
    logging.info(f"Exported JSON to {json_path}")

if __name__ == "__main__":
    build_pgx_network()
