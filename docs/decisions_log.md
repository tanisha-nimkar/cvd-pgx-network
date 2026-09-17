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
