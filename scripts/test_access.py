# exploratory connectivity test, not pipeline code
import requests

def test_apis():
    print("--- Testing Connectivity ---")

    # PharmGKB / ClinPGx API test for CYP2C19
    try:
        r = requests.get("https://api.clinpgx.org/v1/clinicalAnnotation?gene=CYP2C19")
        print(f"PharmGKB/ClinPGx API: Status {r.status_code}")
    except Exception as e:
        print(f"PharmGKB/ClinPGx API Error: {e}")

    # DGIdb API test for CYP2C19
    try:
        query = """
        {
          genes(names: ["CYP2C19"]) {
            nodes {
              name
              interactionCategories {
                name
              }
            }
          }
        }
        """
        r = requests.post("https://dgidb.org/api/graphql", json={"query": query})
        print(f"DGIdb API: Status {r.status_code}")
    except Exception as e:
        print(f"DGIdb API Error: {e}")

    # STRING API test for CYP2C19
    try:
        r = requests.get("https://string-db.org/api/json/network?identifiers=CYP2C19&species=9606")
        print(f"STRING API: Status {r.status_code}")
    except Exception as e:
        print(f"STRING API Error: {e}")

    # Open Targets API test
    try:
        r = requests.get("https://api.platform.opentargets.org/api/v4/graphql")
        print(f"Open Targets API: Status {r.status_code}")
    except Exception as e:
        print(f"Open Targets API Error: {e}")

if __name__ == "__main__":
    test_apis()
