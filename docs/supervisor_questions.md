# Supervisor Questions

This document contains questions requiring supervisor input or approval.

Questions should be resolved before the relevant downstream analysis begins.

## Week 1 — Drug Panel

### Question 1 — Final Drug Panel

⚠ Supervisor Discussion Required

The current proposed panel is:

- Clopidogrel
- Warfarin
- Statins

Is this drug panel appropriate for the project's IHD/CVD pharmacogenomic scope?

Should any drug class be added or removed?

**Decision deadline:** End of Week 1

## Week 2 — Population Data

### Question 2 — IndiGenomes Fallback

⚠ Supervisor Discussion Required if required

If appropriate IndiGenomes access is unavailable, is the proposed fallback of gnomAD South Asian data + 1000 Genomes GIH acceptable?

The populations would remain unpooled.

**Decision deadline:** Week 2

## Week 3 — Pharmacogene Representation

### Question 3 — Star-Allele / Haplotype Scope

⚠ Supervisor Discussion Required

Is the planned star-allele / haplotype handling sufficient for the selected pharmacogenes?

**Decision deadline:** Week 3

## Week 7 — PPI Network

### Question 4 — STRING Parameters

⚠ Supervisor Discussion Required

Confirm:

- STRING confidence threshold ≥700
- Treatment of text-mining-only interactions
- Final network construction parameters

**Decision deadline:** Week 7

## Week 8 — Disease Module

### Question 5 — CVD/IHD Disease Definition

⚠ Supervisor Discussion Required

Which Open Targets EFO terms should define the disease module?

Possible scope needs to be explicitly selected rather than changing during analysis.

**Decision deadline:** Week 8

## Week 9 — Composite Score

### Question 6 — Prioritization Weighting

⚠ Supervisor Discussion Required

Confirm equal weighting as the primary prioritization model.

Confirm that sensitivity analysis is sufficient to assess robustness.

**Decision deadline:** Week 9

## Week 11–12 — Interpretation

### Question 7 — Validation and Publication Claims

⚠ Supervisor Discussion Required

Review whether the validation results support the planned interpretation.

Confirm that the manuscript remains hypothesis-generating and does not make clinical risk or treatment claims.

**Decision deadline:** Week 12

## General Rule

If a result requires a scientific interpretation that could materially change the project's methodology or claims, pause and discuss it with the supervisor rather than making an unrecorded decision.

## Week 6 Statistical Methodology Review
1. **Z-test vs. Fisher's Exact Priority**: The technical implementation applied two-proportion z-test when min expected count >= 5 and Fisher's exact when min expected count < 5. Confirm whether Fisher's exact is preferred as the primary test for all comparisons.
2. **IndiGenomes vs. gnomAD Fallback**: Confirm formal acceptance of gnomAD South Asian data alongside 1000 Genomes strata as the population frequency baseline in lieu of direct IndiGenomes API access.

## Week 7 STRING Network Sign-Off Checkpoints
1. **Confidence Threshold & Text-Mining Filter**: Formal sign-off on combined score >= 0.700 with exclusion of text-mining-only edges (22 edges excluded).
2. **Network Type**: Confirm functional (default) vs. physical interactome mode.
3. **Edge Scope**: Confirm whether to evaluate seed-to-partner edges only (current 225 edges) or query edges among interactors within the 135-node set.
4. **Disease Proximity Network**: For Week 9 shortest path analysis, evaluate whether to map disease modules onto this 135-node ego-network or an expanded background network.
