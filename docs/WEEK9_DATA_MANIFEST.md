# Week 9 Data Manifest & Input/Output Checksum Freeze

**Date:** 2026-11-04  
**Project:** cvd-pgx-network  
**Author:** Tanisha  
**Platform:** Ubuntu WSL2 / Python 3.11  

---

## 1. Frozen Pre-Week 9 Artifacts (Immutable)

| File Path | Size (Bytes) | SHA-256 Checksum | Purpose / Provenance |
| :--- | :--- | :--- | :--- |
| `data/processed/population_frequencies_raw.csv` | 328,863 | `102e2ebab6431779a55820cc6b80fa0a410d9a419dbd39d40ead3a14ca5b3892` | Raw population AF extract (Week 5) |
| `results/tables/population_frequencies.csv` | 852,159 | `f0d5dd70837a6b0e1f1c1578da3c63ed9320e81fb235df450518b2c9458d1bc8` | Statistical test results (Week 6 Freeze) |
| `results/network/network.graphml` | 131,295 | `33bcd28abc5583d264ebda841b09cdba2364ceeccd5d4d63b5a6501b068c8b23` | 135-node Week 7 Ego-Network (Retained for local topology) |
| `results/network/string_edges_audit.tsv` | 522,538 | `a480719dd6853295d34b0ea986691d2cd20bd4f708bd4abc83f70a1338a67d05` | Complete edge classification audit (22 TM-only edges excluded) |
| `results/network/string_id_mapping.tsv` | 362 | `cf184d899d9e1a5993ddb6da318c46df0f1414f0bf0e6d70618e6db63fdee57b` | Panel gene exact STRING v12.5 ENSP identifier mappings |
| `results/tables/centrality.csv` | 26,334 | `7ae727a022ac63eb5ff16228c364ee45ecb1fd7f985aedd1a9ff1b78b8e30c2a` | Unweighted degree & betweenness centrality on Week 7 network |
| `data/raw/9606.protein.links.full.v12.5.txt.gz` | 116,855,136 | `c645be3e636a896f82b0acf63d181aa58a1320b8df9b4106ae9348e68e086ad8` | Raw STRING v12.5 human protein links (All channels) |
| `data/raw/9606.protein.info.v12.5.txt.gz` | 1,970,090 | `144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c` | Raw STRING v12.5 human protein info and gene symbols |

---

## 2. Week 9 Reconstructed & Generated Deliverables

| File Path | Description | Generation Script |
| :--- | :--- | :--- |
| `results/tables/pharmacogene_evidence_table.csv` | Field-by-field biological evidence table for 6 pharmacogenes | `scripts/lock_query_sets.py` |
| `results/tables/query_sets_locked.json` | JSON dictionary of Individual, Drug-Specific, and Pooled query sets | `scripts/lock_query_sets.py` |
| `results/tables/disease_module_direct_mondo0005010.csv` | Primary CAD disease module (enableIndirect: false) | `scripts/reconstruct_disease_modules.py` |
| `results/tables/disease_module_indirect_mondo0005010.csv` | Sensitivity CAD disease module (enableIndirect: true) | `scripts/reconstruct_disease_modules.py` |
| `results/tables/disease_modules_threshold_sensitivity.csv` | Summary of candidate modules across score cutoffs (all, >=0.05, 0.10, 0.20, 0.30) | `scripts/reconstruct_disease_modules.py` |
| `data/processed/string_functional_universe_edges.tsv` | Primary STRING v12.5 Functional LCC (all channels >= 0.700) | `scripts/build_interactome_universes.py` |
| `data/processed/string_functional_tm_excluded_edges.tsv` | Sensitivity STRING v12.5 Functional LCC (text-mining-only excluded) | `scripts/build_interactome_universes.py` |
| `data/processed/string_physical_universe_edges.tsv` | Secondary STRING v12.5 Physical LCC (physical subnetwork >= 0.700) | `scripts/build_interactome_universes.py` |
| `results/tables/proximity_results_primary.csv` | Closest-distance, z-scores, p-values (10k permutations) on primary LCC | `scripts/run_network_proximity.py` |
| `results/tables/proximity_results_sensitivity.csv` | Proximity sensitivity matrix across 6 mandatory dimensions | `scripts/run_network_proximity.py` |
