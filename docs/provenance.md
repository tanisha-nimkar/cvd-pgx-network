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
