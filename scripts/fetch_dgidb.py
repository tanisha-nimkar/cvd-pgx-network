import requests

def get_drug_gene_interactions(gene_symbol):
    url = "https://dgidb.org/api/graphql"
    query = """
    query GetInteractions($gene: String!) {
      genes(names: [$gene]) {
        nodes {
          name
          interactions {
            drug {
              name
            }
            interactionScore
          }
        }
      }
    }
    """
    try:
        response = requests.post(url, json={"query": query, "variables": {"gene": gene_symbol}}, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching DGIdb interactions for {gene_symbol}: {e}")
        raise

if __name__ == "__main__":
    test_gene = "CYP2C19"
    result = get_drug_gene_interactions(test_gene)
    print(f"Successfully fetched DGIdb interactions for {test_gene}.")
