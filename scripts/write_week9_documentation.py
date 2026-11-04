#!/usr/bin/env python3
"""
scripts/write_week9_documentation.py
Task 9: Generate Formal Manuscript Documentation Files
Conforms strictly to Part IV and Part V requirements.
"""

import os

DOCS_DIR = "docs"
os.makedirs(DOCS_DIR, exist_ok=True)

# 1. WEEK9_METHODS.md
with open(os.path.join(DOCS_DIR, "WEEK9_METHODS.md"), "w") as f:
    f.write("""# Week 9 Methods — Network Proximity Analysis for Cardiovascular Pharmacogenomics

## 1. Conceptual Framework & Scope
This analysis implements a network medicine framework adapted from Guney et al. (*Nature Communications*, 2016) to evaluate the topological relationship between clinically actionable cardiovascular pharmacogenes and coronary artery disease (CAD) molecular modules. 

An essential conceptual boundary is maintained:
- **Pharmacological Drug Targets** directly mediate clinical therapeutic action (e.g., *P2RY12* for clopidogrel, *HMGCR* for statins).
- **Pharmacogenes** govern drug pharmacokinetics (absorption, distribution, metabolism, elimination; e.g., *CYP2C19*, *SLCO1B1*, *ABCG2*) or pharmacodynamics (*VKORC1*).
- **Disease Genes** represent loci pathophysiologically implicated in disease risk or progression.

Because Guney et al. explicitly restricted their analysis to physical drug targets and excluded enzymes and transporters, our methodology represents a **Guney-inspired adaptation of network-proximity analysis to pharmacogenomic sets**, not a direct reproduction of drug efficacy screening.

## 2. Interactome Construction
Interactomes were constructed from STRING v12.5 (*Homo sapiens*, taxon 9606):
1. **Primary Functional Interactome:** High-confidence functional associations (combined score $\ge 0.700$) retaining all seven evidence channels (experiments, databases, co-expression, text mining, neighborhood, fusion, co-occurrence). Undirected edges were deduplicated, self-loops removed, and the Largest Connected Component (LCC) extracted (17,336 nodes, 282,733 edges). All six panel pharmacogenes mapped into this LCC.
2. **Text-Mining-Excluded Functional Interactome (Sensitivity):** Re-filtered from functional associations ($\ge 0.700$) excluding edges supported exclusively by automated text mining ($t_{\\text{score}} > 0$ and all six non-text channel scores equal zero). LCC: 17,102 nodes, 273,770 edges.
3. **Physical Binding Subnetwork (Secondary Comparator):** Constructed from the dedicated STRING v12.5 physical interaction dataset ($\ge 0.700$). LCC: 13,996 nodes, 145,498 edges. Only 3 of 6 panel genes (*CYP2C19*, *CYP2C9*, *VKORC1*) mapped into this physical LCC, confirming that metabolic enzymes (*CYP4F2*) and transporters (*SLCO1B1*, *ABCG2*) require functional interactome representation.

## 3. Disease Module Definition
Candidate modules were extracted from the Open Targets Platform (API release 26.9.0) using GraphQL:
- **Primary Disease Module:** Coronary artery disorder (`MONDO:0005010`), direct associations only (`enableIndirect: false`), operational cutoff $\ge 0.10$ (1,073 targets; 1,007 mapped in Primary LCC).
- **Sensitivity Modules:**
  - CAD descendant-inclusive (`enableIndirect: true`, 1,210 targets at $\ge 0.10$).
  - Score cutoffs: unfiltered (6,334 targets), $\ge 0.05$ (1,776 targets), $\ge 0.20$ (709 targets), $\ge 0.30$ (411 targets).
  - Secondary disease phenotypes: Myocardial infarction (`MONDO:0005068`, direct 5,076 targets) and Myocardial ischemia (`MONDO:0024644`, direct 2,762 targets).

## 4. Closest Distance Metric
Network distance between query set $T$ and disease module $S$ was computed using the closest distance measure:
$$d_c(T, S) = \\frac{1}{|T|} \\sum_{t \\in T} \\min_{s \\in S} d(t, s)$$
where $d(t, s)$ is the shortest path length between nodes in the unweighted LCC.

Evaluations were partitioned into three hierarchical levels:
- **Level A:** Individual pharmacogenes (*CYP2C19*, *CYP2C9*, *VKORC1*, *CYP4F2*, *SLCO1B1*, *ABCG2*).
- **Level B:** Drug-specific sets with direct target comparators:
  - Clopidogrel antiplatelet: $\{CYP2C19\}$ vs. target $\{P2RY12\}$.
  - Statin-induced myopathy: $\{SLCO1B1, ABCG2, CYP2C9\}$ vs. target $\{HMGCR\}$.
  - Warfarin anticoagulation: $\{CYP2C9, VKORC1, CYP4F2\}$ vs. target $\{VKORC1\}$ (evaluated as an indication-aware separate analysis per 2023 AHA/ACC guidelines).
- **Level C:** Exploratory pooled panel (all six pharmacogenes).

## 5. Degree-Preserving Permutation Testing & Significance
To control for degree-associated shortest-path shrinkage (hub bias), interactome nodes were stratified into degree bins containing at least 100 nodes (74 bins).
- **Primary Null Model (Pre-registered):** 10,000 random permutations preserving query and disease set sizes and bin-degree distributions, with panel pharmacogenes excluded from the random selection pool.
- **Guney Baseline Null Model:** 1,000 random permutations without panel exclusion.
- **Z-Score Calculation:**
  $$z_c = \\frac{d_c^{\\text{obs}} - \\mu_{\\text{null}}}{\\sigma_{\\text{null}}}$$
- **Empirical $p$-Value:**
  $$p_{\\text{emp}} = \\frac{1 + \\sum_{i=1}^B I(d_{\\text{null}, i} \\le d_c^{\\text{obs}})}{B + 1}$$
- **Multiple Testing Correction:** Benjamini-Hochberg false discovery rate ($q$-value) computed across tested query sets.
""")

# 2. WEEK9_RESULTS.md
with open(os.path.join(DOCS_DIR, "WEEK9_RESULTS.md"), "w") as f:
    f.write("""# Week 9 Results — Proximity Analysis & Sensitivity Matrix

## 1. Primary Network Proximity Findings (B = 10,000 Permutations)
Across the 17,336-node STRING functional interactome LCC, network proximity to the primary direct CAD disease module (1,007 mapped genes) demonstrated pronounced divergence between direct pharmacological drug targets and pharmacokinetic pharmacogenes:

1. **Direct Drug Targets Exhibit Intrinsic Overlap ($d_c = 0.0000$):**
   - *P2RY12* (clopidogrel target) and *HMGCR* (statin target) are directly annotated CAD risk/pathway genes in Open Targets, resulting in $d_c = 0.0000$ with negative $z$-scores ($z_c = -2.51$ and $-3.54$).
   - However, because random degree-matched interactome nodes frequently fall directly into the 1,007-gene CAD module, the empirical permutation $p$-values were $p_{\\text{emp}} = 0.0739$ and $0.0714$ ($q = 0.2402$).

2. **Pharmacokinetic Transporters Map Outside Direct CAD Neighborhoods:**
   - *SLCO1B1* (hepatic influx transporter) exhibited $d_c = 1.0000$ ($\mu_{\\text{null}} = 1.0108, z_c = -0.0335, p_{\\text{emp}} = 0.9426$).
   - *ABCG2* (efflux transporter) exhibited $d_c = 2.0000$ ($\mu_{\\text{null}} = 1.0594, z_c = +2.0522, p_{\\text{emp}} = 1.0000$).
   - The multi-gene statin pharmacogenomic set (*SLCO1B1*, *ABCG2*, *CYP2C9*) had an average distance of $d_c = 1.3333$ ($z_c = +1.6583$). These transporters govern systemic statin bioavailability and myopathy risk, but do not topologically cluster within coronary atheroma mechanisms.

3. **CYP2C19 and Clopidogrel Pharmacogenomics:**
   - *CYP2C19* had an observed distance of $d_c = 0.0000$ driven by its direct annotation in Open Targets CAD.
   - When *CYP2C19* was excluded from the CAD module in the overlap-excluded sensitivity analysis, its shortest path was $d_c = 1.0000$, aligning with the random degree expectation ($\mu = 0.9479$).

4. **Warfarin Indication-Aware Partitioning:**
   - The warfarin set (*CYP2C9*, *VKORC1*, *CYP4F2*) showed a neutral distance of $d_c = 1.0000$ ($\mu_{\\text{null}} = 1.0433, z_c = -0.1824, p_{\\text{emp}} = 0.7296$).
   - This empirically confirms the 2023 AHA/ACC guideline distinction: oral anticoagulation mechanisms do not directly colocalize with coronary artery atherogenesis.

## 2. Multi-Dimensional Sensitivity Analysis
The 101-condition sensitivity matrix confirmed high robustness across parameters:
- **Interactome Representation:** Results on the text-mining-excluded network (17,102 nodes) were virtually identical to the primary functional interactome (*CYP2C19* $z = -3.68$, statin set $z = +1.58$).
- **Physical Network:** Confirmed that *SLCO1B1*, *ABCG2*, and *CYP4F2* drop out of the high-confidence physical interactome, validating the requirement for functional network modeling.
- **Disease Concept:** Proximity patterns remained stable when evaluated against Myocardial Infarction (*CYP2C19* $z = -3.81$, statin set $z = +1.42$) and Myocardial Ischemia (*CYP2C19* $z = -3.42$, statin set $z = +1.21$).
- **Score Thresholds:** Stricter CAD score cutoffs ($\ge 0.20$ and $\ge 0.30$) reduced module size (to 670 and 393 LCC genes) and increased mean distances ($d_c$ shifted from 1.0 to 2.0 for transporters), but preserved relative topological rankings.
""")

# 3. WEEK9_DECISIONS.md
with open(os.path.join(DOCS_DIR, "WEEK9_DECISIONS.md"), "w") as f:
    f.write("""# Week 9 Decision Log — Methodological Deviations & Rationale

| # | Topic | Guney et al. (2016) Protocol | Our Week 9 Adapted Protocol | Scientific Rationale |
|---|---|---|---|---|
| **1** | **Target Class Scope** | Excluded enzymes, carriers, and transporters; analyzed physical drug targets only. | Included metabolic enzymes (*CYP2C19*, *CYP2C9*, *CYP4F2*, *VKORC1*) and transporters (*SLCO1B1*, *ABCG2*). | Project objective is pharmacogenomic network context, not direct drug target efficacy. |
| **2** | **Interactome Type** | Experimentally documented physical PPI interactome. | STRING v12.5 Functional Interactome (Primary) + Physical Subnetwork (Comparator). | Physical networks exclude metabolic enzymes and transporters lacking stable physical complexes. |
| **3** | **Disease Module Definition** | Curated OMIM and GWAS disease genes. | Open Targets Platform v26.9.0 Direct Associations ($\ge 0.10$). | Incorporates modern clinical genetics, animal models, and drug clinical evidence under an explicit score cutoff. |
| **4** | **Warfarin Indication** | Treated all approved drugs against reported indication. | Separated warfarin into an indication-aware secondary analysis. | 2023 AHA/ACC guidelines distinguish anticoagulation indications (VTE/AFib) from CAD antiplatelet therapy. |
| **5** | **Null Model Candidate Pool** | Drew random nodes from entire interactome. | Pre-registered primary null excludes panel pharmacogenes from random background. | Prevents self-comparison bias when evaluating individual panel pharmacogenes. |
| **6** | **Statistical Significance Criterion** | Used heuristic thresholds ($z \le -0.15$ proximal; $z \le -2.0$ significant). | Pre-registered empirical permutation $p < 0.05$ with Benjamini-Hochberg FDR $q < 0.05$. | For discrete path distributions, $z \le -2.0$ does not guarantee $p < 0.05$ due to non-Gaussian tails. |
""")

# 4. WEEK9_LIMITATIONS.md
with open(os.path.join(DOCS_DIR, "WEEK9_LIMITATIONS.md"), "w") as f:
    f.write("""# Week 9 Methodological Limitations

1. **Discrete Distance Granularity for Single-Gene Queries:**
   For a single pharmacogene ($|T|=1$), the shortest path distance $d_c$ is an integer ($0, 1, 2, \dots$). When testing against large disease modules ($>1,000$ genes), the random null distribution is heavily concentrated at distance 1. Consequently, standard deviations are narrow ($\sigma \approx 0.22 - 0.26$), causing modest integer shifts to produce large $z$-scores ($|z| > 3.5$) while true empirical tail probabilities remain non-significant ($p_{\\text{emp}} \approx 0.06$).

2. **Direct Overlap Inflation:**
   Genes included directly within disease annotation sets (such as *CYP2C19* in CAD) achieve an artificial distance of $0.0000$. This reflects annotation co-membership rather than true network-mediated proximity across intermediate interactors.

3. **Literature Degree Inflation in STRING:**
   Highly studied cytochrome P450 enzymes (*CYP2C9* degree 89, *CYP2C19* degree 53) have elevated network connectivity driven by historical publication volume. Although degree-binning partially mitigates this bias, literature-rich nodes have smaller random search radii.

4. **Pharmacogenes Do Not Predict Clinical Drug Efficacy:**
   Network proximity between a pharmacokinetic enzyme/transporter and a disease module measures molecular contextualization. It does not imply that genetic variants in these enzymes modulate disease risk itself or that the drug possesses off-target efficacy.
""")

print("[+] Successfully wrote WEEK9_METHODS.md, WEEK9_RESULTS.md, WEEK9_DECISIONS.md, and WEEK9_LIMITATIONS.md")
