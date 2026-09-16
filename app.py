import streamlit as st
import sqlite3
import pandas as pd
import json
import networkx as nx
from pathlib import Path

DB_PATH = Path("data/database/cvd_pgx.db")
PROC_DIR = Path("data/processed")

st.set_page_config(page_title="Cardiovascular PGx Network Explorer", layout="wide")

st.title("🫀 Cardiovascular Pharmacogenomics (PGx) Network Explorer")
st.markdown("Interactive query engine and graph visualizer for heterogeneous PGx datasets.")

# Sidebar Controls
st.sidebar.header("Navigation & Controls")
gene_selection = st.sidebar.selectbox("Select Target Gene", ["CYP2C19", "CYP2C9", "CYP4F2", "SLCO1B1", "ABCG2", "VKORC1"])

conn = sqlite3.connect(DB_PATH)

# Tab 1: Relational Database Query
tab1, tab2 = st.tabs(["📊 Relational Data Explorer", "🕸️ Network Graph Metrics"])

with tab1:
    st.subheader(f"Gene Profile: {gene_selection}")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("### CPIC Guidelines")
        guidelines_df = pd.read_sql_query("SELECT drug_name, guideline_name FROM cpic_guidelines WHERE gene_symbol = ?;", conn, params=(gene_selection,))
        st.dataframe(guidelines_df, use_container_width=True)

    with col2:
        st.write("### PharmVar Star Alleles")
        pv_df = pd.read_sql_query("SELECT star_allele, rsid, variant_type FROM pharmvar_alleles WHERE gene_symbol = ?;", conn, params=(gene_selection,))
        st.dataframe(pv_df, use_container_width=True)

    st.write("### ClinVar Variants Preview")
    cv_df = pd.read_sql_query("SELECT chrom, pos_grch38, rsid, ref, alt, clinical_significance FROM clinvar_variants WHERE gene_symbol = ? LIMIT 20;", conn, params=(gene_selection,))
    st.dataframe(cv_df, use_container_width=True)

with tab2:
    st.subheader("Network Topology Summary")
    json_files = list(PROC_DIR.glob("cvd_pgx_network_*.json"))
    if json_files:
        latest_json = sorted(json_files)[-1]
        with open(latest_json, "r") as fp:
            data = json.load(fp)
        G = nx.node_link_graph(data)
        
        m1, m2, m3 = st.columns(3)
        m1.metric("Total Graph Nodes", G.number_of_nodes())
        m2.metric("Total Graph Edges", G.number_of_edges())
        
        centrality = nx.degree_centrality(G)
        m3.metric(f"{gene_selection} Centrality", f"{centrality.get(gene_selection, 0):.4f}")

conn.close()
