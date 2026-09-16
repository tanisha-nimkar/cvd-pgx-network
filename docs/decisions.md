
## Phase 2: Data Normalization, Standardization & QC
* **Date:** 2026-09-16
* **Status:** Complete
* **Actions Taken:**
  * Created `scripts/normalize_data.py` to ingest 21 raw heterogeneous files and produce 4 standardized TSVs in `data/processed/`.
  * Standardized gene symbols (uppercase HGNC), canonical drug names (lowercase), and GRCh38 genomic coordinates.
  * Generated programmatically verified QC reports (`qc_summary_20260916.json` and `qc_audit_report_20260916.json`).
  * Confirmed cross-source integrity and `CYP2C19` benchmark coverage.

## Phase 2: Data Normalization, Standardization & QC
* **Date:** 2026-09-16
* **Status:** Complete
* **Actions Taken:**
  * Created `scripts/normalize_data.py` to ingest 21 raw heterogeneous files and produce 4 standardized TSVs in `data/processed/`.
  * Standardized gene symbols (uppercase HGNC), canonical drug names (lowercase), and GRCh38 genomic coordinates.
  * Generated programmatically verified QC reports (`qc_summary_20260916.json` and `qc_audit_report_20260916.json`).
  * Confirmed cross-source integrity and `CYP2C19` benchmark coverage.

## Phase 2: Data Normalization, Standardization & QC
* **Date:** 2026-09-16
* **Status:** Complete
* **Actions Taken:**
  * Created `scripts/normalize_data.py` to ingest 21 raw heterogeneous files and produce 4 standardized TSVs in `data/processed/`.
  * Standardized gene symbols (uppercase HGNC), canonical drug names (lowercase), and GRCh38 genomic coordinates.
  * Created `scripts/qc_audit.py` and generated programmatically verified QC reports (`qc_summary_20260916.json` and `qc_audit_report_20260916.json`).
  * Confirmed cross-source integrity and `CYP2C19` benchmark coverage.

## Phase 2: Data Normalization, Standardization & QC
* **Date:** 2026-09-16
* **Status:** Complete
* **Actions Taken:**
  * Created `scripts/normalize_data.py` to ingest 21 raw heterogeneous files and produce 4 standardized TSVs in `data/processed/`.
  * Standardized gene symbols (uppercase HGNC), canonical drug names (lowercase), and GRCh38 genomic coordinates.
  * Created `scripts/qc_audit.py` and generated programmatically verified QC reports (`qc_summary_20260916.json` and `qc_audit_report_20260916.json`).
  * Confirmed cross-source integrity and `CYP2C19` benchmark coverage.

## Phase 3: SQLite Relational Schema & Database Ingestion

- **Database Engine**: SQLite (`data/database/cvd_pgx.db`).
- **Schema Normalization**: Designed a 3NF-compliant schema centered around the `genes` master entity table.
- **Null-Handling Strategy**: Removed strict `NOT NULL` constraints on variable bio-entity fields like `star_allele` to account for incomplete records across upstream PharmVar and ClinVar TSV exports.
- **Performance Optimization**: Created explicit Single-Column B-Tree Indexes on query target fields (`gene_symbol`, `drug_name`, `rsid`) across tables.

## Phase 4: Network Construction & API Integration
* **Date:** 2026-09-16
* **Status:** Complete
* **Actions Taken:**
  * Created `scripts/build_network.py` to construct a 486-node NetworkX graph connecting Genes, Drugs, and Variants.
  * Exported normalized network artifacts (`cvd_pgx_network_20260916.graphml` and `cvd_pgx_network_20260916.json`).
  * Verified graph topology in `scripts/verify_network.py`: Identified `CYP2C19` as the primary hub (degree centrality = 0.2330).
  * Implemented `scripts/api_fetch_pharmgkb.py` for dynamic web service integration with the PharmGKB REST API.
