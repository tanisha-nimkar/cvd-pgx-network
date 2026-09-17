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
