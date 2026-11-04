# Week 9 Methods — Network Proximity Analysis for Cardiovascular Pharmacogenomics

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
2. **Text-Mining-Excluded Functional Interactome (Sensitivity):** Re-filtered from functional associations ($\ge 0.700$) excluding edges supported exclusively by automated text mining ($t_{\text{score}} > 0$ and all six non-text channel scores equal zero). LCC: 17,102 nodes, 273,770 edges.
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
$$d_c(T, S) = \frac{1}{|T|} \sum_{t \in T} \min_{s \in S} d(t, s)$$
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
  $$z_c = \frac{d_c^{\text{obs}} - \mu_{\text{null}}}{\sigma_{\text{null}}}$$
- **Empirical $p$-Value:**
  $$p_{\text{emp}} = \frac{1 + \sum_{i=1}^B I(d_{\text{null}, i} \le d_c^{\text{obs}})}{B + 1}$$
- **Multiple Testing Correction:** Benjamini-Hochberg false discovery rate ($q$-value) computed across tested query sets.
