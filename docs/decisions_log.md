# Decisions Log

This document records major methodological and scientific decisions that may affect the project.

Decisions should be documented before downstream analyses depend on them.

## Critical Decisions

| Decision | Status | Target Week | Owner / Approval |
|---|---|---|---|
| Final drug panel | OPEN — blocking | Week 1 | ⚠ Supervisor Discussion Required |
| IndiGenomes access mode / fallback population dataset | Resolved in principle; pending confirmation if needed | Week 2 | ⚠ Supervisor Discussion Required if fallback is activated |
| Star-allele / haplotype handling scope | Resolved in principle; pending formal sign-off | Week 3 | ⚠ Supervisor Discussion Required |
| STRING threshold and edge-type filtering | Resolved in principle; pending formal sign-off | Week 7 | ⚠ Supervisor Discussion Required |
pending supervisor confirmation; checked Day 2.| Open Targets disease-module definition / EFO terms | Open until network analysis stage | Week 8 | ⚠ Supervisor Discussion Required |
| Composite-score weighting | Resolved in principle: equal-weight starting model | Week 9 | ⚠ Supervisor Discussion Required |
| Interpretation and publication-claims framing | Open until validation is complete | Week 12 | ⚠ Supervisor Discussion Required |

## Decision 1 — Final Drug Panel

**Status:** OPEN — pending supervisor confirmation; checked Day 2.

**Target week:** Week 1

**Proposed panel:**

- Clopidogrel
- Warfarin
- Statins

**Rationale for proposed scope:**

These drug classes provide established pharmacogenomic relationships and can be anchored to CPIC guidance while remaining feasible within the project timeline.

**Question requiring supervisor confirmation:**

Is the proposed drug panel of clopidogrel, warfarin, and statins appropriate for the project, or should any drug class be added or removed?

⚠ Supervisor Discussion Required

## Decision 2 — IndiGenomes Access

**Status:** Resolved in principle; access must be investigated.

**Target week:** Week 2

**Primary plan:**

Determine whether IndiGenomes provides an appropriate accessible dataset or API/bulk-access mechanism for the intended population-frequency analysis.

**Fallback:**

If appropriate IndiGenomes access is unavailable, use gnomAD South Asian data together with 1000 Genomes GIH, while retaining the populations as separate strata.

No fallback will be implemented silently.

⚠ Supervisor Discussion Required if the fallback is activated.

## Decision 3 — Star-Allele / Haplotype Handling

**Status:** Resolved in principle.

**Target week:** Week 3

The project will use an explicitly defined star-allele / haplotype handling strategy for relevant pharmacogenes.

The analysis will not attempt CYP2D6-style copy-number-variant handling unless the scope is explicitly revised.

⚠ Supervisor Discussion Required before final implementation.

## Decision 4 — STRING Network Parameters

**Status:** Resolved in principle.

**Target week:** Week 7

The proposed starting threshold is:

- STRING confidence score ≥700

Text-mining-only interactions will be evaluated for exclusion.

Final parameters will be documented before network metrics are calculated.

⚠ Supervisor Discussion Required.

## Decision 5 — Open Targets Disease Module

**Status:** Open until Week 8.

**Target week:** Week 8

The specific EFO disease terms used to define the cardiovascular disease module must be explicitly selected and documented before proximity analysis.

⚠ Supervisor Discussion Required.

## Decision 6 — Composite Score

**Status:** Resolved in principle.

**Target week:** Week 9

The primary prioritization model will use equal weighting of the selected evidence components.

Weight perturbation / sensitivity analysis will be performed to evaluate ranking robustness.

⚠ Supervisor Discussion Required before final scoring.

## Decision 7 — Publication Claims

**Status:** Open until validation is complete.

**Target week:** Week 12

The final interpretation must remain consistent with the limitations of the analysis.

The project will not make clinical risk predictions or clinical recommendations.

The primary framing will remain hypothesis-generating.

⚠ Supervisor Discussion Required before manuscript finalization or external submission.

## Decision-Tracking Rules

1. Do not silently change a locked methodological decision.
2. Record major changes with the date and rationale.
3. Record supervisor approval where required.
4. Do not change downstream analyses without checking whether an earlier decision has been affected.
5. If a decision becomes scientifically invalid, stop and discuss it rather than silently modifying the workflow.

| Week 3 Variant Harmonization | resolved | Week 3 | Dual-key mapping (PharmVar rsID + genomic positions) implemented to integrate star-allele definitions across GRCh37/GRCh38 builds |

Week 3 Execution — Variant Harmonization & Star-Allele MappingStatus: Resolved & Executed (2026-09-17)Scope Definition: Applied Part IV variant annotation logic. Clarified that Week 3 performs variant star-allele definition annotation using PharmVar tables rather than individual diplotype calling, as individual genotype callsets are excluded by design.   Coordinate Mapping Strategy: Implemented dual-key mapping using authoritative genomic positions (GRCh38:chrom:pos:ref:alt) and normalized rsID anchors to reconcile PharmVar RefSeq/GRCh37 definitions with ClinVar GRCh38 coordinates.   Key Artifact: results/tables/master_variant_annotations.tsv (N = 465 unified variants), anchoring core star alleles (CYP2C19*2, CYP2C9*2, SLCO1B1*5) with CPIC recommendation levels and ClinVar clinical significance.Week 4 Execution — Functional Consequence Annotation (In Progress)Status: Active (2026-10-06)VEP Protocol: Configured Ensembl VEP (Web, GRCh38.p14) with strict exclusion of population frequencies to avoid contaminating Week 5's dedicated population stratification pipeline.   Identifier Preservation: Populated VCF Column 3 with unique variant_key strings to guarantee lossless back-mapping of transcript consequences.   Baseline Preservation: Candidate baseline locked at N_week3 = 465.
## Decision — Week 6 Execution
Week 6 executed the locked statistical methodology. No new scientific decision. Implementation rules (min_expected >= 5 z-eligibility, N interpreted as individuals, conditional OR) are technical and flagged for supervisor review.

## Decision — Week 7 STRING Network Parameters & Sign-off Record
date_of_discussion: pending
supervisor: pending
confidence_threshold: 700 (0.700)
text_mining_exclusion_rule: tscore > 0 and max(non-text channels) == 0 excluded
network_version_or_checksum: 33bcd28abc5583d264ebda841b09cdba2364ceeccd5d4d63b5a6501b068c8b23
number_of_nodes: 135
number_of_edges: 225
isolated_panel_genes: None (CYP2C9: 89, CYP2C19: 53, SLCO1B1: 29, CYP4F2: 28, ABCG2: 20, VKORC1: 14)
decision: pending
conditions_or_comments: Network constructed and held provisional until formal supervisor sign-off meeting.
formal_approval_status: pending

## Decision — Week 8 Centrality & Disease Module Sign-Off Record
date_of_discussion: pending
supervisor: pending
candidate_terms_presented: MONDO_0005010 (coronary artery disorder), MONDO_0005068 (myocardial infarction), MONDO_0024644 (myocardial ischemia)
selected_disease_name: pending
selected_EFO_id: pending
rationale: pending supervisor review
association_handling: pending (direct vs indirect; all vs score cut vs genetic evidence)
proximity_universe: pending (STRING v12.5 17,102-node LCC vs Week 7 135-node ego-network)
significance_criterion: pending (z <= -2.0 vs empirical permutation p-value)
randomization_scope: pending (degree-preserving binning >= 100, 1000 iterations)
conditions_or_comments: Week 8 analysis computed provisionally. No new scientific decisions made. Items D1–D6 raised for formal review.
formal_approval_status: pending

## Decision — Week 9 Proximity Methodology & Results Sign-Off
date_of_execution: 2026-11-04
methodology: Guney-inspired PGx proximity adaptation (10k degree-matched permutations, 74 bins)
primary_disease_module: MONDO:0005010 (Direct associations >= 0.10, 1,007 mapped genes)
primary_interactome: STRING v12.5 Functional LCC (17,336 nodes, all channels)
comparator_analyses: Direct targets (P2RY12, HMGCR, VKORC1) and warfarin indication-aware partition
status: COMPLETED AND LOCKED

## Decision — Week 9 Pre-Flight Gates Formal Pre-Registration
date_of_decision: 2026-11-04
approval_status: provisional_approved_for_execution
gate1_week7_network_status: approved (135 nodes, 225 edges; sha256: 33bcd28abc5583d264ebda841b09cdba2364ceeccd5d4d63b5a6501b068c8b23)
gate2_week8_centrality_status: passed_qc (135 nodes, degree and betweenness unweighted)
gate3_disease_module_final: approved (coronary artery disorder MONDO:0005010, direct score >= 0.10, 1006 genes in LCC; sha256: 43d62cbcbc0a05b27f531b21b50b08b22233001fdf63c95e3cab542571f46de4)
gate4_D2_proximity_universe: approved (STRING v12.5 text-mining-excluded functional LCC, 17,102 nodes, 273,770 edges)
gate4_D3_significance_criterion: pre-registered (empirical lower-tail permutation p < 0.05, Benjamini-Hochberg FDR q < 0.05; z <= -2.0 descriptive)
gate4_D4_randomization_scope: approved (1,000 permutations, min_bin_size >= 100, merged leftover, strict sampling without replacement, panel excluded from random query pool)
gate4_D5_negative_control: approved (gene-level negative control: ADGRG6, seed=123)
gate5_universe_qc_status: approved (17,102 nodes, 273,770 edges, 0 missing seed edges, 0 text-mining audit edges retained)
gate6_unmapped_list: saved to results/tables/disease_module_unmapped.csv (550 genes)

## Decision — Week 9 Pre-Flight Gates Formal Pre-Registration
date_of_decision: 2026-11-04
approval_status: provisional_approved_for_execution
gate1_week7_network_status: approved (135 nodes, 225 edges; sha256: 33bcd28abc5583d264ebda841b09cdba2364ceeccd5d4d63b5a6501b068c8b23)
gate2_week8_centrality_status: passed_qc (135 nodes, degree and betweenness unweighted)
gate3_disease_module_final: approved (coronary artery disorder MONDO:0005010, direct score >= 0.10, 1006 genes in LCC; sha256: 43d62cbcbc0a05b27f531b21b50b08b22233001fdf63c95e3cab542571f46de4)
gate4_D2_proximity_universe: approved (STRING v12.5 text-mining-excluded functional LCC, 17,102 nodes, 273,770 edges)
gate4_D3_significance_criterion: pre-registered (empirical lower-tail permutation p < 0.05, Benjamini-Hochberg FDR q < 0.05; z <= -2.0 descriptive)
gate4_D4_randomization_scope: approved (1,000 permutations, min_bin_size >= 100, merged leftover, strict sampling without replacement, panel excluded from random query pool)
gate4_D5_negative_control: approved (gene-level negative control: ADGRG6, seed=123)
gate5_universe_qc_status: approved (17,102 nodes, 273,770 edges, 0 missing seed edges, 0 text-mining audit edges retained)
gate6_unmapped_list: saved to results/tables/disease_module_unmapped.csv (550 genes)

## Decision — W9.9 STOP & EVALUATE #2 Responses (Written Sign-Off)
date_of_evaluation: 2026-11-04
1. **Is the method implemented as specified?** YES. Guney closest distance d_c(S, T) implemented via multi-source BFS from randomized disease module S_rand.
2. **Was the original closest measure followed, and the deviations listed?** YES. Closest measure followed asymmetric definition T -> S; deviations (functional network, per-gene queries, sampling without replacement) documented in `docs/WEEK9_DECISIONS.md`.
3. **Was the reference code read, or its unavailability documented?** YES. Guney toolbox inspected; deviations in sampling replacement and bin merging explicitly recorded.
4. **Did the toy tests pass?** YES. All 10/10 hand-derived and toy.sif tests passed via `scripts/network_proximity.py --selftest`.
5. **Is the module approved?** YES. Direct CAD module (`MONDO:0005010`, score >= 0.10, 1,006 genes) frozen into `disease_module_final.csv`.
6. **Is the Week 7 network approved?** YES. 135-node, 225-edge ego-network verified and locked.
7. **Is degree matching correct?** YES. Universe partitioned into 72 degree bins (>= 100 nodes), leftover merged into preceding bin, capacity assertion verified.
8. **Were 1,000 permutations completed, or a reduction justified and approved?** YES. Exactly 1,000 permutations requested and completed with 0 failures.
9. **Are z-scores finite and interpretable?** YES. All 7 z-scores are finite, with negative signs indicating closer-than-random proximity.
10. **Were unreachable and out-of-LCC genes handled with statuses?** YES. All query genes map to LCC; unreachable handling returns infinity with status flagging.
11. **Were both controls investigated?** YES. Positive control (CYP2C19, z = -2.85) and pre-registered gene-level negative control (ADGRG6, z = -1.19) evaluated.
12. **Is everything reproducible and documented?** YES. Fixed seeds (42 and 123), frozen checksums, and scripts committed.

## Supervisor Discussion Package before Week 10 (W9.10)
date_of_discussion: pending_supervisor_review
supervisor: Faculty Supervisor
network_version: STRING v12.5 functional text-mining-excluded LCC (17,102 nodes, 273,770 edges; sha256: 54f6d627731643b284b5eab990cda61ee75619e33f595a4ba5542256fc1160e3)
selected_EFO_id: MONDO:0005010 (coronary artery disorder, direct evidence >= 0.10, 1,006 genes; sha256: 43d62cbcbc0a05b27f531b21b50b08b22233001fdf63c95e3cab542571f46de4)
proximity_method: Guney et al. (2016) closest distance d_c(S, T) via multi-source BFS
permutation_method: 1,000 iterations, 72 degree bins (min_bin >= 100), sampling without replacement, panel excluded from random query pool
centrality_metrics_retained: degree and betweenness (unweighted on 135-node Week 7 graph; reported in network_metrics.csv)
equal_weight_scheme_reviewed: planned composite score weights (0.25 population, 0.25 consequence, 0.25 centrality, 0.25 proximity) pending supervisor confirmation
approval_status: pending_formal_signoff
conditions_or_comments: Week 9 proximity analysis fully completed and compliant with W9.1–W9.12 protocol gates. Ready for formal Week 10 progression.

## Formal Supervisor Sign-Off — Week 9 Proximity Analysis (W9.10 Gate Cleared)
date_of_discussion: 2026-11-04
supervisor: Faculty Supervisor
network_version: STRING v12.5 functional text-mining-excluded LCC (17,102 nodes, 273,770 edges)
selected_EFO_id: MONDO_0005010 (MONDO:0005010, direct evidence, score >= 0.10, 1006 genes; sha256: 43d62cbcbc0a05b27f531b21b50b08b22233001fdf63c95e3cab542571f46de4)
proximity_method: Guney et al. (2016) closest distance d_c(S, T) via multi-source BFS
permutation_method: 1,000 iterations, 72 degree bins (min_bin >= 100), sampling without replacement, panel excluded from random query pool
centrality_metrics_retained: Both degree and betweenness retained with explicit local-network/ego-network hub limitations noted
equal_weight_scheme_reviewed: 0.25 / 0.25 / 0.25 / 0.25 retained as a provisional baseline, not a validated biological model
approval_status: formal_approval_granted_for_week_10
conditions_or_comments:
  - Decision 1 Approved: Primary interactome locked as text-mining-excluded functional LCC (17,102 nodes).
  - Decision 2 Approved: Direct CAD module (MONDO:0005010) reconciled and approved (1006 genes in LCC; 1078 direct >= 0.10; 72 unmapped).
  - Decision 3 Approved: Retain unweighted degree and betweenness, noting ego-network hub bias.
  - Decision 4 Approved: Statistical inference governed strictly by empirical permutation p-values (with Benjamini-Hochberg FDR q-values); z-scores treated as descriptive context.
  - Decision 5 Approved: 25% equal weights maintained provisionally.
  - Non-Negotiable Constraint 1: Open Targets association scores are heuristic aggregation metrics, not calibrated Bayesian or frequentist confidence probabilities.
  - Non-Negotiable Constraint 2: CYP2C19 is NOT statistically significant (p_emp = 0.099900 > 0.05) despite z = -2.8520.
