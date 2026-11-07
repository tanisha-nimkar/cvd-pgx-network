#!/usr/bin/env bash
set -e

echo "=== W12 STEP 9: FRESH-CLONE REPLICATION TEST ==="
TEST_DIR="../cvd-pgx-fresh-test"
REF_DIR="../cvd_reference_snapshot"

# 1. Clean previous test clone if present
rm -rf "$TEST_DIR"

# 2. Clone repository from local repo at tag week12-frozen
echo "[*] Cloning fresh repository copy to $TEST_DIR..."
git clone --branch week12-frozen ~/cvd-pgx-network "$TEST_DIR"

# 3. Copy pinned raw data snapshots into fresh clone
echo "[*] Placing archived raw data snapshots into $TEST_DIR/data/raw/..."
mkdir -p "$TEST_DIR/data/raw"
mkdir -p "$TEST_DIR/data/processed"
cp -r data/raw/* "$TEST_DIR/data/raw/" 2>/dev/null || true
cp data/processed/ensembl_pop_cache_single.json "$TEST_DIR/data/processed/" 2>/dev/null || true

# 4. Execute Snakemake pipeline end-to-end
cd "$TEST_DIR"
echo "[*] Executing snakemake --cores 1 in fresh clone..."
snakemake --cores 1

# 5. Run cryptographic and structural byte-by-byte comparison
echo ""
echo "=== VERIFYING REPRODUCIBILITY AGAINST FROZEN REFERENCE ==="
python3 - << 'PYEOF'
import os
import hashlib
import pandas as pd
import networkx as nx

REF_DIR = "../cvd_reference_snapshot"

test_files = [
    "results/tables/variants_annotated.csv",
    "results/tables/population_frequencies.csv",
    "results/tables/centrality.csv",
    "results/tables/disease_module_final.csv",
    "results/tables/network_metrics.csv",
    "results/tables/prioritization_ranked.csv",
    "results/tables/sensitivity_summary.csv",
    "results/tables/validation_results.csv",
    "results/network/network.graphml"
]

def sha256(path):
    if not os.path.exists(path):
        return "MISSING"
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

comparison_results = []
for rel_path in test_files:
    ref_file = os.path.join(REF_DIR, rel_path)
    fresh_file = rel_path
    
    ref_sha = sha256(ref_file)
    fresh_sha = sha256(fresh_file)
    
    if rel_path.endswith(".graphml"):
        # Structural graph comparison (node and edge sets)
        G_ref = nx.read_graphml(ref_file)
        G_fresh = nx.read_graphml(fresh_file)
        same_nodes = set(G_ref.nodes()) == set(G_fresh.nodes())
        same_edges = set(G_ref.edges()) == set(G_fresh.edges())
        status = "EXACT MATCH (Structural)" if same_nodes and same_edges else "MISMATCH"
        detail = f"Nodes: {G_fresh.number_of_nodes()}, Edges: {G_fresh.number_of_edges()}"
    else:
        # Byte identity comparison
        if ref_sha == fresh_sha:
            status = "EXACT MATCH (Byte Identical)"
            detail = f"SHA: {fresh_sha[:12]}..."
        else:
            # Fallback numeric comparison for float tables
            try:
                df_ref = pd.read_csv(ref_file)
                df_fresh = pd.read_csv(fresh_file)
                if df_ref.equals(df_fresh):
                    status = "EXACT MATCH (DataFrame)"
                    detail = f"Shape: {df_fresh.shape}"
                else:
                    status = "FLOAT DISCREPANCY"
                    detail = f"Ref {df_ref.shape} vs Fresh {df_fresh.shape}"
            except Exception as e:
                status = "MISMATCH"
                detail = str(e)

    comparison_results.append({
        "file": rel_path,
        "fresh_sha256": fresh_sha[:16] + "...",
        "status": status,
        "detail": detail
    })

df_comp = pd.DataFrame(comparison_results)
print(df_comp.to_string(index=False))
PYEOF
