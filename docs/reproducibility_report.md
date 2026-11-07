# Pipeline Reproducibility Report & Execution Inventory

**Audit Date:** 2026-11-07  
**Reference Git Tag:** `pre-week12-reference`  
**Pipeline Framework:** Snakemake (Python 3.11, Conda)  
**Controlling Protocol:** Week 12 Workflow Specifications  

---

## 1. Pipeline Script Inventory & Execution Classification

| Stage | Script | Input Files / Sources | Output Files | Execution Class | Status |
|:---|:---|:---|:---|:---|:---|
| **W2** | `scripts/fetch_variants.py` | ClinPGx & CPIC REST APIs | `data/raw/clinpgx_*.json`, `data/raw/cpic_*.json` | Live API / Pinned Snapshot | Verified (Pinned) |
| **W3** | `scripts/merge_annotate_variants.py` | `data/raw/{cpic,clinvar,pharmvar}_*` | `results/tables/variants_annotated_week3_backup.csv`, `results/tables/annotated_*.tsv` | Deterministic CLI | Verified |
| **W4a**| `scripts/prepare_vep_input.py` | `results/tables/variants_annotated_week3_backup.csv` | `data/processed/vep_input.vcf` | Deterministic CLI | Verified |
| **W4b**| `scripts/annotate_variants_vep.py` | `data/processed/vep_input.vcf`, `data/raw/vep_output_*.txt` | `results/tables/variants_annotated.csv` | Deterministic CLI | Verified (Locked) |
| **W5** | `scripts/prepare_population_frequencies_raw.py` | `data/processed/ensembl_pop_cache_single.json` | `data/processed/population_frequencies_raw.csv` | Deterministic CLI | Verified |
| **W6** | `scripts/frequency_stats.py` | `data/processed/population_frequencies_raw.csv`, `config/config.yaml` | `results/tables/population_frequencies.csv`, `figures/fig2_freq_comparison.png` | Deterministic CLI | Verified (Freeze) |
| **W7** | `scripts/build_network.py` | STRING v12.5 API / `config/config.yaml` | `results/network/network.graphml`, `results/network/string_edges_audit.tsv` | Seeded Graph Builder | Verified (LCC) |
| **W8a**| `scripts/network_centrality.py` | `results/network/network.graphml` | `results/tables/centrality.csv` | Deterministic Graph Metric | Verified |
| **W8b**| `scripts/open_targets_module.py` | Open Targets v4 GraphQL API (`MONDO_0005010`) | `results/tables/disease_module_final.csv`, `results/tables/disease_module_candidates.csv` | Pinned GraphQL Snapshot | Verified |
| **W8c**| `scripts/build_interactome_universes.py` | `data/raw/9606.protein.links.full.v12.5.txt.gz` | `data/processed/string_functional_tm_excluded_edges.tsv` | Deterministic Bulk Filter | Verified |
| **W9** | `scripts/network_proximity.py` | GraphML, Module, STRING Universe | `results/tables/network_metrics.csv`, `results/tables/proximity_null_*.csv` | Stochastic Seeded (B=1000, seed=42) | Verified |
| **W9-P**| `scripts/plot_w9_figures.py` | `network_metrics.csv`, GraphML | `figures/fig3_network.png`, `figures/fig5_proximity.png` | Deterministic Plotter | Verified |
| **W10**| `scripts/prioritization_score.py`| Variants, Frequencies, Network Metrics | `results/tables/prioritization_ranked.csv`, `results/tables/sensitivity_*.csv` | Deterministic Composite | Verified (Tag: week10-ranking) |
| **W10-P**| `scripts/plot_w10_figures.py` | `prioritization_ranked.csv`, `sensitivity_analysis.csv` | `figures/fig6_ranking.png`, `figures/fig7_sensitivity.png` | Deterministic Plotter | Verified |
| **W11a**| `scripts/prepare_biogrid.py` | `data/raw/BIOGRID-ORGANISM-Homo_sapiens-5.0.263.tab3.txt` | `data/processed/biogrid_edges.tsv` | Deterministic Physical Filter | Verified |
| **W11b**| `scripts/biogrid_validation.py` | `biogrid_edges.tsv`, Module, Ranked Variants | `results/tables/biogrid_network_metrics.csv`, `results/tables/validation_detail_biogrid.csv` | Stochastic Seeded (B=1000, seed=42) | Verified |
| **W11c**| `scripts/gnomad_validation.py` | gnomAD v4 GraphQL API | `results/tables/validation_detail_gnomad.csv` | Live API Cross-Check | Verified |
| **W11d**| `scripts/pubmed_validation.py` | NCBI E-Utilities API | `results/tables/validation_detail_pubmed.csv` | Live API Screening | Verified |
| **W11e**| `scripts/pathway_enrichment.py` | `data/raw/ReactomePathways.gmt`, GraphML | `results/tables/pathway_enrichment.csv` | Offline Hypergeometric Test | Verified |
| **W11f**| `scripts/independent_validation.py`| BioGRID, gnomAD, PubMed detail tables | `results/tables/validation_results.csv` | Long-Format Matrix Aggregator | Verified |

---

## 2. P1 Output Collision Resolution
The collision between Week 3 and Week 4 is formally resolved by splitting targets:
- `scripts/merge_annotate_variants.py` targets `results/tables/variants_annotated_week3_backup.csv` (`093fd009...`).
- `scripts/annotate_variants_vep.py` takes the backup plus pinned VEP results to write `results/tables/variants_annotated.csv` (`480d1826...`).

---

## 3. Baseline Checksum Ledger (Frozen Pre-Week 12 Reference)

| Category | File | SHA-256 Checksum |
|:---|:---|:---|
| **GraphML** | `results/network/network.graphml` | `33bcd28abc5583d264ebda841b09cdba2364ceeccd5d4d63b5a6501b068c8b23` |
| **CAD Module**| `results/tables/disease_module_final.csv` | `43d62cbcbc0a05b27f531b21b50b08b22233001fdf63c95e3cab542571f46de4` |
| **STRING Metrics** | `results/tables/network_metrics.csv` | `153b14b3aae83c0486521fdcfe43f98c80c86e4954f49452f6885dabfcf17810` |
| **Ranked Prioritization** | `results/tables/prioritization_ranked.csv` | `c1d0df71491a12195d136d859c999aa79741646e9c8ef48dd507429a3c0b6ba5` |
| **Sensitivity Summary** | `results/tables/sensitivity_summary.csv` | `0b224d2225e67019e59336d575bd5c582ae5eb3cc52f1099ea4763fbc6e7e905` |
| **BioGRID Validation** | `results/tables/validation_detail_biogrid.csv` | `7943612c34dd415ab3915b6ea5a5a3f22203d1964906012eb8c8188c235e2a0f` |
| **gnomAD Validation** | `results/tables/validation_detail_gnomad.csv` | `6c836172de93010cabbf49f905a87f842b56a8d82651758dc7bd6d97226cfcb4` |
| **PubMed Validation** | `results/tables/validation_detail_pubmed.csv` | `068a4840e299074a9d137ab326c668ddb51c9d3f92dded9339034afb2d47ee44` |
| **Reactome Pathways** | `results/tables/pathway_enrichment.csv` | `07a4bcef7cd8009a292798aac6a9e9c49092d9d1d7aff34b8f5763c777dd2ab8` |
| **Integrated Validation** | `results/tables/validation_results.csv` | `1977a490800bb0a46b1f9710367ecdfae9da084068056de8abee200c082f283f` |

## 4. Frozen Results & Deliverables Inventory (SHA-256 Ledger)
**Freeze Date:** 2026-11-07  
**Git Tag:** `week12-frozen`

| Deliverable Path | Size (Bytes) | SHA-256 Checksum | Freeze Status |
|:---|:---|:---|:---|
| `results/network/network.graphml` | 131,295 | `33bcd28abc5583d264ebda841b09cdba2364ceeccd5d4d63b5a6501b068c8b23` | frozen_verified |
| `results/network/string_edges_audit.tsv` | 522,538 | `a480719dd6853295d34b0ea986691d2cd20bd4f708bd4abc83f70a1338a67d05` | frozen_verified |
| `results/tables/centrality.csv` | 26,334 | `7ae727a022ac63eb5ff16228c364ee45ecb1fd7f985aedd1a9ff1b78b8e30c2a` | frozen_verified |
| `data/processed/string_functional_tm_excluded_edges.tsv` | 15,966,251 | `54f6d627731643b284b5eab990cda61ee75619e33f595a4ba5542256fc1160e3` | frozen_verified |
| `results/tables/disease_module_final.csv` | 144,027 | `43d62cbcbc0a05b27f531b21b50b08b22233001fdf63c95e3cab542571f46de4` | frozen_verified |
| `results/tables/network_metrics.csv` | 47,695 | `153b14b3aae83c0486521fdcfe43f98c80c86e4954f49452f6885dabfcf17810` | frozen_verified |
| `results/tables/proximity_null_CYP2C19.csv` | 2,014 | `866925aa20b8864a7f14402cd84715f21484d5f3342b91e0cf3c57ce2affce40` | frozen_verified |
| `results/tables/proximity_null_ADGRG6.csv` | 2,014 | `f379360682f5743eaa5a723674f2abc9414fde1b31c7363818d81c6a7c277d87` | frozen_verified |
| `results/tables/prioritization_ranked.csv` | 165,376 | `c1d0df71491a12195d136d859c999aa79741646e9c8ef48dd507429a3c0b6ba5` | frozen_verified |
| `results/tables/sensitivity_analysis.csv` | 6,413,482 | `a8a66327cb765b6f4acd64c35bb97e97fd9320788b94a88f0b0242cf3896086b` | frozen_verified |
| `results/tables/sensitivity_summary.csv` | 5,484 | `0b224d2225e67019e59336d575bd5c582ae5eb3cc52f1099ea4763fbc6e7e905` | frozen_verified |
| `data/processed/biogrid_edges.tsv` | 12,332,159 | `25cb04a924bed7bd928502e6d222261ce286c09d5ba93ac044cd38585dfe78bc` | frozen_verified |
| `results/tables/biogrid_network_metrics.csv` | 22,545 | `508578a936d3b8c920a1d7bc4ee6c642bb02d0598cd9cd29bfda8a70b560abca` | frozen_verified |
| `results/tables/validation_detail_biogrid.csv` | 1,820 | `7943612c34dd415ab3915b6ea5a5a3f22203d1964906012eb8c8188c235e2a0f` | frozen_verified |
| `results/tables/validation_detail_gnomad.csv` | 5,947 | `6c836172de93010cabbf49f905a87f842b56a8d82651758dc7bd6d97226cfcb4` | frozen_verified |
| `results/tables/validation_detail_pubmed.csv` | 15,693 | `068a4840e299074a9d137ab326c668ddb51c9d3f92dded9339034afb2d47ee44` | frozen_verified |
| `results/tables/pathway_enrichment.csv` | 85,312 | `07a4bcef7cd8009a292798aac6a9e9c49092d9d1d7aff34b8f5763c777dd2ab8` | frozen_verified |
| `results/tables/validation_results.csv` | 22,241 | `1977a490800bb0a46b1f9710367ecdfae9da084068056de8abee200c082f283f` | frozen_verified |
| `figures/fig2_freq_comparison.png` | 680,040 | `8c40f95b1e5a878c461ad1c27b899d592c271bc43733fa24748b98cfe3845861` | frozen_verified |
| `figures/fig3_network.png` | 684,069 | `c0bcff235dc01a794cdb45d09343e32d33a6d38d20ea59a88f6116601420acda` | frozen_verified |
| `figures/fig5_proximity.png` | 428,792 | `6143192494027f13ecdcccff57ab63e591581ab18382d67e34f9f66880800998` | frozen_verified |
| `figures/fig6_ranking.png` | 526,050 | `4dd84c4126d6787ea5f9b8400596c0050a0bea0651f1198e0511e08c9dfdc376` | frozen_verified |
| `figures/fig7_sensitivity.png` | 590,195 | `777ab052d385afb68370deaa0a45645c2a5a3e91f394837a650e4ab47783d73f` | frozen_verified |
