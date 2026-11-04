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

### W9.5 Toy-Network Validation Results

### W9.5 Toy-Network Validation Results

### W9.5 Toy-Network Validation Results
| test_case | expected | observed | pass_fail | notes |
| --- | --- | --- | --- | --- |
| d_c({A}, S) | 3.0 | 3.0 | PASS | Hand-derived path A-B-C-F |
| d_c({A,D}, S) | 2.0 | 2.0 | PASS | Hand-derived (3+1)/2 = 2.0 |
| d_c({E}, S) | 0.0 | 0.0 | PASS | E in S -> distance 0 |
| T={G} unreachable | inf | inf | PASS | Unreachable pair returns inf |
| Bins min_bin=2 | 2 bins: {A,E,F}, {B,D,C} | 2 bins | PASS | Merged leftover bin |
| draw_matched distinct | 2 distinct matched nodes | ['A', 'C'] | PASS | Without-replacement match |
| z on null=[2,2,3,3] | mu=2.5, sd=0.5, z=-3.0, p=0.2 | mu=2.5, sd=0.5, z=-3.0, p=0.2 | PASS | Hand-derived z and p |
| z on zero variance | NaN (zero_variance_null) | nan | PASS | Zero SD returns NaN |
| toy.sif closest (A,C)->(B,D,E) | 1.0 | 1.0 | PASS | Matches Guney README |
| toy.sif 20k draws | mu ~ 0.4167, sd ~ 0.3436 | mu=0.4169, sd=0.3436 | PASS | Within +-0.02 of hand enumeration |


## Week 9 — Controlling Proximity Run (LOCKED PROTOCOL VERIFIED)

| Field | Value |
| :--- | :--- |
| Inputs + Checksums | `results/network/network.graphml` (`33bcd28abc5583d264ebda841b09cdba2364ceeccd5d4d63b5a6501b068c8b23`); `results/tables/disease_module_final.csv` (`43d62cbcbc0a05b27f531b21b50b08b22233001fdf63c95e3cab542571f46de4`); `data/processed/string_functional_tm_excluded_edges.tsv` (`54f6d627731643b284b5eab990cda61ee75619e33f595a4ba5542256fc1160e3`); STRING v12.5 (Taxon 9606) |
| Sets | Query = CYP2C19, CYP2C9, VKORC1, CYP4F2, SLCO1B1, ABCG2 + Negative Control ADGRG6; Module total = 1,006, in LCC = 1,006; Universe nodes = 17,102, edges = 273,770, LCC = 100% |
| Method | Guney et al. (2016) closest distance measure d_c(S, T); direction: query T -> module S; unweighted shortest paths |
| Reference Code Answers | 1. Replacement: Without replacement inside each set (Guney toolbox uses with-replacement); 2. Random stream: Single seeded NumPy generator stream (`seed=42`); 3. Reappearance: Panel excluded from random query pool; 4. Leftover bins: Merged into previous bin |
| Deviations | Functional network (enzymes/transporters require non-binding functional links); per-gene query (|T|=1); panel exclusion from query null pool; discrete empirical permutation testing |
| Null Specification | Degree bins >= 100 nodes (72 bins, smallest=101); N=1,000 permutations (1,000 completed, 0 failed); `seed=42`; ddof=0 |
| z and p Formulas | z = (d_obs - mu) / sigma (negative = closer); p_emp = (1 + count(d_null <= d_obs)) / (N + 1) (lower-tail) |
| Controls | Positive Control: CYP2C19 (d_c = 0.0000, z = -2.8520, p_emp = 0.099900, in module: True); Negative Control: ADGRG6 (d_c = 1.0000, null mu = 1.7330, z = -1.1895, p_emp = 0.264735) |
| Outputs + Checksums | `results/tables/network_metrics.csv` (`153b14b3aae83c0486521fdcfe43f98c80c86e4954f49452f6885dabfcf17810`); `figures/fig3_network.png` (`c0bcff235dc01a794cdb45d09343e32d33a6d38d20ea59a88f6116601420acda`); `figures/fig5_proximity.png` (`6143192494027f13ecdcccff57ab63e591581ab18382d67e34f9f66880800998`) |
| QC / Approval | W9.1–W9.12 Fully Satisfied / Approved for Supervisor Sign-off |

### W9.8 QC Table
| QC Item | Required Specification | Observed Value / Status |
| :--- | :--- | :--- |
| Week 7 graph and module approvals | Verified in decisions log | Approved in `docs/decisions_log.md` (Gate 1 & Gate 3) |
| Query, module and universe sets | Documented with counts | 6 panel + 1 control; 1,006 module genes; 17,102 universe nodes |
| Identifier matching | Unmapped lists saved | `results/tables/disease_module_unmapped.csv` (550 unmapped targets saved) |
| Network coverage | Query & module genes in LCC | 100% of panel (6/6) and 100% of final module (1,006/1,006) in LCC |
| Disconnected components | LCC applied; no out-of-LCC queries | LCC enforced (17,102 nodes); 0 queries outside LCC |
| Method vs. original | Deviation table complete | Documented in `docs/WEEK9_DECISIONS.md` and provenance |
| Reference code | Inspected; 4 questions answered | Answered in provenance (wrappers.py & selection_utilities.py) |
| Toy validation | All rows pass | 10/10 PASS on Graph 1 and toy.sif 20k draws |
| Degree matching | Bin sizes >= min; capacity check | 72 bins constructed; min size = 101; capacity assertion passed |
| Iterations | Requested=1,000; recorded | 1,000 requested, 1,000 completed, 0 failed |
| Seed; failed iterations | Recorded | `seed=42`, `control_seed=123`, failed = 0 |
| Null distributions | Saved to disk | `results/tables/proximity_null_<gene>.csv` saved for all 7 queries |
| z formula and sign | Checked by hand on one gene | CYP2C19: (0 - 0.912) / 0.3198 = -2.8518 ~ -2.8520 (PASS) |
| NaN or infinite z | 0 unexplained | 0 unexplained NaNs (all 7 finite with status 'ok') |
| Both controls | Investigated and documented | CYP2C19 (positive) and ADGRG6 (negative) documented |
| Required tables & figures | Created | `network_metrics.csv`, `fig3_network.png`, `fig5_proximity.png` generated |
| Provenance | Complete | Synced and locked |

## Week 10 — Composite Variant Prioritization & Sensitivity Analysis (FROZEN)

| Field | Value |
| :--- | :--- |
| Inputs + Checksums | `results/tables/variants_annotated.csv` (480d1826bcef7b31bfc753ec7df725a976364aab21c4b1ad9819224e2d71b617); `results/tables/population_frequencies.csv` (f0d5dd70837a6b0e1f1c1578da3c63ed9320e81fb235df450518b2c9458d1bc8); `results/tables/network_metrics.csv` (153b14b3aae83c0486521fdcfe43f98c80c86e4954f49452f6885dabfcf17810) |
| Ranking Unit | Variant candidates (465 rows total; 116 complete-case ranked, 349 unranked missing GIH counts) |
| Score Weights | Equal baseline: 0.25 Population Enrichment, 0.25 Functional Consequence, 0.25 Network Centrality, 0.25 Disease Proximity (provisional model) |
| Scaler Fit Bounds | P_min_max: (-0.5399, 0.5022); F_min_max: (0.0, 3.0); C_min_max: (0.1040, 0.6639); D_min_max: (-1.6594, 2.8520) |
| Guideline Exclusion Test | Passed (mutation tests confirmed guideline_covered is strictly descriptive; 0 leakage into score or rank) |
| Sensitivity Grid | 593 unique normalized scenarios (+-50% relative weight grid; 5 levels per weight) |
| Baseline Correlations | Frequency-only: Spearman rho = 0.3469 (p = 1.36e-4), Kendall tau = 0.2523 (p = 6.2e-5); Evidence-only: Spearman rho = 0.3735 (p = 3.6e-5), Kendall tau = 0.3064 (p = 6.2e-5) |
| Controls Checked | CYP2C19*3 ranked #1 (Score: 0.7406); Negative control ADGRG6 retained at gene level only (no variants to rank) |
| Output Checksums | `prioritization_ranked.csv` (c1d0df71491a12195d136d859c999aa79741646e9c8ef48dd507429a3c0b6ba5); `sensitivity_analysis.csv` (a8a66327cb765b6f4acd64c35bb97e97fd9320788b94a88f0b0242cf3896086b); `sensitivity_summary.csv` (0b224d2225e67019e59336d575bd5c582ae5eb3cc52f1099ea4763fbc6e7e905); `fig6_ranking.png` (4dd84c4126d6787ea5f9b8400596c0050a0bea0651f1198e0511e08c9dfdc376); `fig7_sensitivity.png` (777ab052d385afb68370deaa0a45645c2a5a3e91f394837a650e4ab47783d73f) |
| Status | Complete, frozen, tagged `week10-ranking` |

### W10.11 QC Table
| QC Item | Required Specification | Observed Result | Status |
| :--- | :--- | :--- | :--- |
| All four inputs verified | Checksums match freeze entries | Verified against provenance | PASS |
| Ranking unit & directions | S1-S10 documented and applied | Variant level; C & D attached; D flipped (-z) | PASS |
| Missing-value policy | Complete-case ranking (S6) | 116 ranked, 349 unranked (status: unranked_missing_components) | PASS |
| Constant columns & duplicates | 0 constant columns; 0 duplicate keys | Checked (P, F, C, D nunique >= 4; 0 dup variant_key) | PASS |
| Weights | Sum to 1.0; 0.25 baseline | Confirmed by code assertion (abs(sum - 1) < 1e-12) | PASS |
| Guideline coverage isolation | Absent from score, scaler, and rank | Passed 3 mutation tests (inverted, all-True, all-False) | PASS |
| Sensitivity grid | +-50% relative; vectors sum to 1.0 | 593 unique normalized scenarios generated | PASS |
| Baseline correlations | Common complete-case set | Spearman (0.347 / 0.374), Kendall (0.252 / 0.306) on 116 rows | PASS |
| Outputs & figures | All files created | Ranked table, sensitivity tables, Figures 6 & 7 on disk | PASS |
