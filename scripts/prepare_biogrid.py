#!/usr/bin/env python3
"""
scripts/prepare_biogrid.py
Downloads and filters BioGRID physical interactome strictly per W11.1 & V1.
"""

import os
import sys
import zipfile
import hashlib
import requests
import pandas as pd

RAW_DIR = "data/raw"
PROCESSED_DIR = "data/processed"
os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(PROCESSED_DIR, exist_ok=True)

BIOGRID_ZIP = os.path.join(RAW_DIR, "BIOGRID-ORGANISM-LATEST.tab3.zip")
BIOGRID_PROCESSED = os.path.join(PROCESSED_DIR, "biogrid_edges.tsv")

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

print("=== W11.1: DOWNLOADING & PROCESSING BIOGRID DATA ===")

# 1. Download BioGRID latest release if not already cached
url = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip"

if not os.path.exists(BIOGRID_ZIP):
    print(f"[*] Downloading {url}...")
    headers = {"User-Agent": "Mozilla/5.0"}
    r = requests.get(url, headers=headers, stream=True)
    if r.status_code == 200:
        with open(BIOGRID_ZIP, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024*1024):
                f.write(chunk)
        print(f"[+] Download complete: {BIOGRID_ZIP}")
    else:
        print(f"[!] HTTP error {r.status_code} while downloading BioGRID.")
        sys.exit(1)
else:
    print(f"[*] Found cached BioGRID archive: {BIOGRID_ZIP}")

zip_sha = sha256(BIOGRID_ZIP)
print(f"    Raw Zip SHA-256: {zip_sha}")

# 2. Extract Homo sapiens tab3 file
with zipfile.ZipFile(BIOGRID_ZIP, 'r') as z:
    hs_files = [f for f in z.namelist() if "Homo_sapiens" in f and f.endswith(".tab3.txt")]
    if not hs_files:
        raise FileNotFoundError("Could not find Homo sapiens tab3 file in BioGRID zip archive!")
    hs_file = hs_files[0]
    print(f"[*] Extracting {hs_file} from zip...")
    z.extract(hs_file, RAW_DIR)
    raw_tab3_path = os.path.join(RAW_DIR, hs_file)

# 3. Filter strictly per V1
print(f"[*] Loading and filtering {raw_tab3_path} per V1 specifications...")
df = pd.read_csv(raw_tab3_path, sep="\t", low_memory=False)
print(f"    Raw human interactions: {len(df):,}")

# Filter V1:
# - Organism ID 9606 on both interactors
# - Experimental System Type == 'physical'
# - No self-interactions (A != B)
c_orgA = df["Organism ID Interactor A"] == 9606
c_orgB = df["Organism ID Interactor B"] == 9606
c_phys = df["Experimental System Type"] == "physical"

df_phys = df[c_orgA & c_orgB & c_phys].copy()
print(f"    Human physical interactions: {len(df_phys):,}")

# Remove self-interactions
c_nonself = df_phys["Official Symbol Interactor A"] != df_phys["Official Symbol Interactor B"]
df_phys = df_phys[c_nonself].copy()

# Deduplicate undirected pairs
df_phys["pair"] = df_phys[["Official Symbol Interactor A", "Official Symbol Interactor B"]].apply(
    lambda r: tuple(sorted(r)), axis=1
)
edges = df_phys.drop_duplicates(subset=["pair"]).copy()
print(f"    Deduplicated undirected physical edges: {len(edges):,}")

# 4. Save to data/processed/biogrid_edges.tsv
out_df = pd.DataFrame({
    "symbol1": [p[0] for p in edges["pair"]],
    "symbol2": [p[1] for p in edges["pair"]]
})
out_df.to_csv(BIOGRID_PROCESSED, sep="\t", index=False)
proc_sha = sha256(BIOGRID_PROCESSED)
print(f"[+] Saved {len(out_df):,} physical edges to {BIOGRID_PROCESSED}")
print(f"    Processed SHA-256: {proc_sha}")

