import requests

def get_clinical_annotations(gene_symbol):
    # Query ClinPGx/PharmGKB clinical annotations using gene symbol search
    url = f"https://api.clinpgx.org/v1/data/clinicalAnnotation?location.genes.symbol={gene_symbol}"
    headers = {"Accept": "application/json"}
    try:
        response = requests.get(url, headers=headers, timeout=30)
        
        # Fallback to PharmGKB public REST endpoint if ClinPGx params differ
        if response.status_code in (400, 404):
            url = f"https://api.pharmgkb.org/v1/data/clinicalAnnotation?location.genes.symbol={gene_symbol}"
            response = requests.get(url, headers=headers, timeout=30)

        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching PharmGKB/ClinPGx data for {gene_symbol}: {e}")
        raise

if __name__ == "__main__":
    test_gene = "CYP2C19"
    result = get_clinical_annotations(test_gene)
    print(f"Successfully fetched annotations for {test_gene}. Record count: {len(result.get('data', result))}")
