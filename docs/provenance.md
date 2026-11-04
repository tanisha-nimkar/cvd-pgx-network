# Data Provenance & Lineage Record

**Acquisition Date:** 2026-08-30  
**Genome Build:** GRCh38 / hg38  
**Panel Genes:** CYP2C19, CYP2C9, VKORC1, CYP4F2, SLCO1B1, ABCG2  

---

## 1. Primary Data Sources Summary

| Resource | Version / Date | File Format | Total Files | Scope / Description |
| :--- | :--- | :--- | :--- | :--- |
| **ClinPGx** | 2026-08-30 | JSON | 6 | Standardized clinical pharmacogenomic annotations via PharmGKB API |
| **CPIC** | 2026-08-30 | JSON | 6 | Clinical Pharmacogenetics Implementation Consortium recommendations |
| **ClinVar** | GRCh38 (2026-08-30) | VCF | 6 | NCBI ClinVar whole-genome sliced per target gene coordinates |
| **PharmVar** | v6.2.28 (2026-08-30) | TSV | 3 | Core haplotype/star-allele definitions (CYP2C19, CYP2C9, SLCO1B1) |

---

## 2. File-Level Inventory & Variant Counts

### 2.1 ClinPGx & CPIC Data (`data/raw/`)
* **ClinPGx JSONs:**
  * `clinpgx_ABCG2_20260830.json` (1 entry)
  * `clinpgx_CYP2C19_20260830.json` (1 entry)
  * `clinpgx_CYP2C9_20260830.json` (1 entry)
  * `clinpgx_CYP4F2_20260830.json` (1 entry)
  * `clinpgx_SLCO1B1_20260830.json` (1 entry)
  * `clinpgx_VKORC1_20260830.json` (1 entry)
* **CPIC Guidelines JSONs:**
  * `cpic_ABCG2_20260830.json` (1 recommendation record)
  * `cpic_CYP2C19_20260830.json` (27 recommendation records)
  * `cpic_CYP2C9_20260830.json` (25 recommendation records)
  * `cpic_CYP4F2_20260830.json` (3 recommendation records)
  * `cpic_SLCO1B1_20260830.json` (9 recommendation records)
  * `cpic_VKORC1_20260830.json` (1 recommendation record)

### 2.2 ClinVar Gene Slices (`data/raw/`)
* `clinvar_ABCG2_20260830.vcf`: 79 variants (chr4:88090128-88231628)
* `clinvar_CYP2C19_20260830.vcf`: 86 variants (chr10:94762681-94853539)
* `clinvar_CYP2C9_20260830.vcf`: 60 variants (chr10:94938658-94990091)
* `clinvar_CYP4F2_20260830.vcf`: 97 variants (chr19:15878938-15908000)
* `clinvar_SLCO1B1_20260830.vcf`: 73 variants (chr12:21283626-21392658)
* `clinvar_VKORC1_20260830.vcf`: 25 variants (chr16:31093557-31100351)

### 2.3 PharmVar Allele Definitions (`data/raw/`)
* `pharmvar_CYP2C19_20260830.tsv`: ~33.6 KB (Release 6.2.28)
* `pharmvar_CYP2C9_20260830.tsv`: ~27.0 KB (Release 6.2.28)
* `pharmvar_SLCO1B1_20260830.tsv`: ~26.8 KB (Release 6.2.28)

---

## 3. Data Processing Logs
* Operational logs stored in `data/raw/fetch_log.txt`
* Summary metric JSON stored in `data/raw/fetch_counts_20260830.json`

| annotated_cpic_pairs.tsv | v1.0 | 2026-09-17 | results/tables/ | Filtered CPIC levels A/B gene-drug pairs |
| annotated_clinvar_variants.tsv | v1.0 | 2026-09-17 | results/tables/ | Extracted pathogenic/classified variants from ClinVar VCFs |

| master_variant_annotations.tsv | v1.0 | 2026-09-17 | results/tables/ | Unified table linking PharmVar star-alleles, ClinVar clinical annotations, and panel loci |

## Week 6 — Frequency Statistics
| Field | Value |
|---|---|
| Input | `data/processed/population_frequencies_raw.csv`, sha256 `102e2ebab6431779a55820cc6b80fa0a410d9a419dbd39d40ead3a14ca5b3892`, build GRCh38 |
| Population sources/releases | 1000Genomes Phase 3 (GIH, PJL, BEB, STU, ITU), gnomAD v2/v3 (gnomAD_SAS, NFE) |
| Target populations | GIH, PJL, BEB, STU, ITU, gnomAD_SAS |
| Valid comparisons (m) / untested | 688 / 2,102 (1,976 missing target counts, 126 missing reference counts) |
| Test selection | Two-proportion z-test primary if min expected count >= 5 (N=262), else Fisher's exact primary (N=426); Fisher secondary where z-test run |
| Low-power definition | Individuals < 20 (sample_n basis or AN/2); sparse_counts_flag = not z_eligible |
| Odds Ratio | Conditional MLE, target-vs-NFE, exact conditional 95% CI; no pseudocount; sample OR recorded descriptively |
| Multiple testing | Benjamini-Hochberg (BH) + Bonferroni over single family (m=688); alpha = 0.05 |
| Software | scipy 1.17.1, matplotlib 3.11.2, numpy, pandas |
| Output table | `results/tables/population_frequencies.csv` |
| Draft figure | `figures/fig2_freq_comparison.png` (DRAFT) — 6-panel allele frequency comparison across panel genes comparing South Asian cohorts against NFE reference |
| QC results | 688 tested comparisons; 248 BH-significant; 168 Bonferroni-significant; Direction consistency: 100% |
| Known limitations | Allele-count tests assume independent alleles; small 1kGP subpopulation sizes; frequencies integrated across differing sequencing/calling platforms |

## Week 7 — STRING Network (PROVISIONAL until sign-off)
| Field | Value |
|---|---|
| Source / host | STRING API, `https://string-db.org/api` |
| Release/version | 12.5 (stable address: `https://version-12-5.string-db.org`) |
| Access date | 2026-11-04 |
| Organism / network_type | 9606 (*Homo sapiens*) / functional |
| Gene panel + source | CYP2C19, CYP2C9, VKORC1, CYP4F2, SLCO1B1, ABCG2 (from config.yaml) |
| ID mapping | `get_string_ids`, limit=5, exact preferredName match; table: `results/network/string_id_mapping.tsv` (6/6 exact matches) |
| Retrieval | `interaction_partners`, seeds only, no recursion, limit=1000, required_score=0 |
| Score scale | API score 0–1; required_score 0–1000; threshold 700 = 0.700 (verified empirically) |
| Evidence channels used | nscore, fscore, pscore, ascore, escore, dscore, tscore |
| Text-mining-only rule | tscore > 0 and all other channels == 0 (22 edges excluded) |
| Dedup / isolated nodes | Undirected, sorted-pair key / 0 isolated seeds (all 6 connected) |
| Output + checksum | `results/network/network.graphml`, sha256 `33bcd28abc5583d264ebda841b09cdba2364ceeccd5d4d63b5a6501b068c8b23`; audit: `results/network/string_edges_audit.tsv` |
| Manual validation | CYP2C9-PTGS2 (score 0.933, dscore 0.900, PASS); CYP2C9-ALOX15B (score 0.912, dscore 0.900, escore 0.050, PASS) |
| QC results | 135 nodes, 225 edges; min combined score 0.7010; zero text-mining-only edges in graph |
| Approval status | pending (awaiting supervisor sign-off) |

## Week 8 — Centrality + Disease Module (PROVISIONAL)
| Field | Value |
|---|---|
| Input graph / checksum / approval status | `results/network/network.graphml` / `33bcd28abc5583d264ebda841b09cdba2364ceeccd5d4d63b5a6501b068c8b23` / provisional_pending_supervisor |
| NetworkX / graph type | 3.6.1 / undirected simple graph |
| Metrics, normalization | degree ($k/(n-1)$), betweenness (normalized, unweighted); nodes=135; isolated=0 |
| CYP2C19 check | Degree: 53, Degree Centrality: 0.3955, Betweenness Centrality: 0.1886, Connected to 53 interactors (PASS) |
| STRING universe (provisional) | Release 12.5; `9606.protein.links.full.v12.5.txt.gz` (sha256 `c645be3e636a896f82b0acf63d181aa58a1320b8df9b4106ae9348e68e086ad8`); LCC nodes: 17,102 (98.90%), edges: 273,770; all 6 panel genes present |
| Open Targets | Endpoint `https://api.platform.opentargets.org/api/v4/graphql`; Release 26.9.0; retrieval date 2026-11-04; indirect setting: `enableIndirect: true`; page size: 500 |
| Candidates | CAD (`MONDO_0005010`: 6,771 targets, 85 graph overlap, 6,179 in universe LCC); MI (`MONDO_0005068`: 5,613 targets, 63 graph overlap, 5,209 in universe LCC); IHD (`MONDO_0024644`: 3,102 targets, 51 graph overlap, 2,919 in universe LCC) |
| Validation | `results/tables/week8_module_validation.csv`: 15/15 PASS |
| Decision | Pending supervisor meeting (Items D1–D6) |
| Final module file | pending (`results/tables/disease_module_candidates.csv` holds candidates) |

## Week 9 — Proximity Analysis (FINAL LOCKED EXECUTION)
| Field | Value |
|---|---|
| Execution Date | 2026-11-04 |
| Primary Interactome | STRING v12.5 Functional LCC (17,336 nodes, 282,733 edges; sha256 in manifest) |
| Sensitivity Interactomes | TM-Excluded LCC (17,102 nodes, 273,770 edges); Physical LCC (13,996 nodes, 145,498 edges) |
| Primary Disease Module | Coronary artery disorder (`MONDO:0005010`), Direct evidence, cutoff >= 0.10 (1,007 LCC nodes) |
| Permutations | B = 10,000 (Primary, panel-excluded null); B = 1,000 (Guney baseline & sensitivity cells) |
| Degree Binning | 74 bins (min bin size >= 100); seed = 42 |
| Primary Output | `results/tables/proximity_results_primary.csv` |
| Sensitivity Output | `results/tables/proximity_results_sensitivity.csv` (101 cells) |
| Figures Generated | `figures/fig3a_proximity_zscores.png`, `fig3b_distance_distributions.png`, `fig3c_sensitivity_matrix.png` |
| Documentation | `WEEK9_METHODS.md`, `WEEK9_RESULTS.md`, `WEEK9_DECISIONS.md`, `WEEK9_LIMITATIONS.md`, `WEEK9_DATA_MANIFEST.md` |
| Status | Complete, fully reproducible, audited against Part V requirements |
