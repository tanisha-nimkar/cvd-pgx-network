import requests

def get_cpic_pairs():
    # CPIC gene-drug pair endpoint mirrored via ClinPGx / CPIC API
    url = "https://api.clinpgx.org/v1/cpic/pair"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching CPIC gene-drug pairs: {e}")
        raise

if __name__ == "__main__":
    result = get_cpic_pairs()
    print(f"Successfully fetched CPIC pair data. Total records: {len(result)}")
