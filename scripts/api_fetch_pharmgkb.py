import urllib.request
import json
from pathlib import Path
from datetime import date
import logging

PROC_DIR = Path("data/processed")
TODAY = date.today().strftime("%Y%m%d")

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

TARGET_GENES = ["CYP2C19", "CYP2C9", "CYP4F2", "SLCO1B1", "ABCG2", "VKORC1"]
BASE_URL = "https://api.pharmgkb.org/v1/data/gene"

def fetch_pharmgkb_annotations():
    results = {}
    logging.info("Initiating PharmGKB REST API fetch for target CVD genes...")

    for gene in TARGET_GENES:
        url = f"{BASE_URL}?symbol={gene}"
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as response:
                if response.status == 200:
                    data = json.loads(response.read().decode())
                    results[gene] = {
                        "pharmgkb_id": data['data'][0]['id'] if data.get('data') else "NA",
                        "name": data['data'][0]['name'] if data.get('data') else "NA",
                        "symbol": gene
                    }
                    logging.info(f"Successfully fetched PharmGKB metadata for {gene}")
        except Exception as e:
            logging.warning(f"Failed to fetch {gene} from PharmGKB API: {e}")
            results[gene] = {"symbol": gene, "status": "failed", "error": str(e)}

    out_path = PROC_DIR / f"pharmgkb_api_metadata_{TODAY}.json"
    with open(out_path, "w") as fp:
        json.dump(results, fp, indent=2)

    logging.info(f"PharmGKB API metadata saved to {out_path}")

if __name__ == "__main__":
    fetch_pharmgkb_annotations()
