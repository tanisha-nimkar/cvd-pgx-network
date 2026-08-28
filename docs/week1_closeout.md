# Week 1 Close-Out Summary

## Status of Deliverables
- [x] GitHub repository skeleton and directory hierarchy built matching Part XII.
- [x] Environment configured (`cvd-pgx` conda environment with python=3.11, pandas, requests, numpy, pyyaml).
- [x] Core documentation skeleton created (`scope.md`, `decisions_log.md`, `provenance.md`, `supervisor_questions.md`, `data_access_status.md`).
- [x] Database API access confirmed (PharmGKB/ClinPGx, CPIC, DGIdb, Open Targets, STRING).
- [x] Fetch-script skeletons created (`fetch_pharmgkb.py`, `fetch_cpic_pairs.py`, `fetch_dgidb.py`).
- [x] Configuration locked in `config/config.yaml` with confirmed 6-gene panel (clopidogrel + warfarin + statins).
- [x] Dry-run validation passed across all 6 genes without errors.

## Locked Scope Summary
- **Disease Focus:** Ischemic Heart Disease / Cardiovascular Disease (IHD/CVD).
- **Gene Panel:** CYP2C19, CYP2C9, VKORC1, CYP4F2, SLCO1B1, ABCG2.
- **Population Strata:** India-primary + unpooled 1kGP South Asian strata (GIH, PJL, BEB, STU, ITU).
- **Methodology:** Population frequency + variant annotation + STRING PPI (>=700 confidence) + Open Targets disease module + Reactome pathway enrichment + equal-weight non-ML composite score.
