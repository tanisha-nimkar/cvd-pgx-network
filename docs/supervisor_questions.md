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

## Week 8 Supervisor Sign-Off Agenda (D1–D6)
- **D1. Disease Module Selection & Inclusion Rule**: 
  - Candidate 1: *Coronary Artery Disorder* (`MONDO_0005010`, 6,771 targets). Comprehensive phenotypic match for clopidogrel, statin, and warfarin CPIC guidelines.
  - Candidate 2: *Myocardial Infarction* (`MONDO_0005068`, 5,613 targets). Acute ischemic endpoint.
  - Candidate 3: *Myocardial Ischemia* (`MONDO_0024644`, 3,102 targets). Broader parent ontology grouping.
  - Inclusion Handling: Propose filtering to `score >= 0.1` (1,210 genes for CAD) or `score >= 0.3` (476 genes) or genetic-evidence datatypes to prevent dilution by literature text-mining.
- **D2. Proximity Universe**: Confirm using the 17,102-node STRING v12.5 high-confidence LCC for shortest-path distances and degree-binned permutations (ego-network paths are artificially bounded at <= 2).
- **D3. Significance Criterion**: Establish pre-registered cutoff for significant proximity: $z \le -2.0$ or empirical $p < 0.05$.
- **D4. Randomization Scope**: Preserve degree distribution using bins of $\ge 100$ nodes with 1,000 iterations per Guney et al. (2016). Query set size = 1 (per panel gene). Confirm whether panel genes are excluded from the random background pool.
- **D5. Negative-Control Selection**: Pre-register negative control disease phenotype (e.g., non-vascular, non-metabolic condition such as an unlinked sensory/neurodevelopmental trait) to benchmark false-positive proximity.
- **D6. Centrality Metric Representation**: Confirm reporting unweighted degree and betweenness centrality on the 135-node first-degree interactome, acknowledging that panel gene hub status reflects ego-network topology.
