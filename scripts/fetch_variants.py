import requests
import yaml
import json
import time
import logging
from pathlib import Path
from datetime import date

RAW_DIR = Path("data/raw")
RAW_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[
        logging.FileHandler(RAW_DIR / "fetch_log.txt"),
        logging.StreamHandler()
    ]
)

TODAY = date.today().strftime("%Y%m%d")

def load_panel():
    with open("config/config.yaml", "r") as f:
        cfg = yaml.safe_load(f)
    return cfg["gene_panel"]

def _save(obj, path, is_json=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    if is_json:
        with open(path, "w") as f:
            json.dump(obj, f, indent=2)
    else:
        with open(path, "w") as f:
            f.write(obj)
    logging.info(f"Wrote {path} ({path.stat().st_size} bytes)")

def fetch_clinpgx(gene):
    url = "https://api.clinpgx.org/v1/data/gene"
    r = requests.get(url, params={"symbol": gene}, timeout=30)
    time.sleep(0.6)  # Stay under 2 req/sec rate limit
    r.raise_for_status()
    res = r.json()
    data = res.get("data", [])
    if not data:
        logging.warning(f"ZERO RECORDS: clinpgx {gene}")
    _save(data, RAW_DIR / f"clinpgx_{gene}_{TODAY}.json")
    return len(data) if isinstance(data, list) else 1

def fetch_cpic(gene):
    url = "https://api.cpicpgx.org/v1/pair_view"
    r = requests.get(url, params={"genesymbol": f"eq.{gene}"}, timeout=30)
    r.raise_for_status()
    data = r.json()
    if not data:
        logging.warning(f"ZERO RECORDS: cpic {gene}")
    _save(data, RAW_DIR / f"cpic_{gene}_{TODAY}.json")
    return len(data) if isinstance(data, list) else None

def main():
    panel = load_panel()
    counts = {}
    logging.info(f"Starting API retrieval for gene panel: {panel}")
    
    for gene in panel:
        counts[gene] = {}
        
        # ClinPGx Fetch
        try:
            counts[gene]["clinpgx"] = fetch_clinpgx(gene)
        except Exception as e:
            logging.error(f"CLINPGX FAILED {gene}: {e}")
            counts[gene]["clinpgx"] = "FAILED"
            
        # CPIC Fetch
        try:
            counts[gene]["cpic"] = fetch_cpic(gene)
        except Exception as e:
            logging.error(f"CPIC FAILED {gene}: {e}")
            counts[gene]["cpic"] = "FAILED"
            
    summary_path = RAW_DIR / f"fetch_counts_{TODAY}.json"
    _save(counts, summary_path)
    logging.info(f"Fetch completed successfully. Summary: {counts}")

if __name__ == "__main__":
    main()
