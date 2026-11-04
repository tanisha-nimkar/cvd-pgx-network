# Week 9 Results — Proximity Analysis & Sensitivity Matrix

## 1. Primary Network Proximity Findings (B = 10,000 Permutations)
Across the 17,336-node STRING functional interactome LCC, network proximity to the primary direct CAD disease module (1,007 mapped genes) demonstrated pronounced divergence between direct pharmacological drug targets and pharmacokinetic pharmacogenes:

1. **Direct Drug Targets Exhibit Intrinsic Overlap ($d_c = 0.0000$):**
   - *P2RY12* (clopidogrel target) and *HMGCR* (statin target) are directly annotated CAD risk/pathway genes in Open Targets, resulting in $d_c = 0.0000$ with negative $z$-scores ($z_c = -2.51$ and $-3.54$).
   - However, because random degree-matched interactome nodes frequently fall directly into the 1,007-gene CAD module, the empirical permutation $p$-values were $p_{\text{emp}} = 0.0739$ and $0.0714$ ($q = 0.2402$).

2. **Pharmacokinetic Transporters Map Outside Direct CAD Neighborhoods:**
   - *SLCO1B1* (hepatic influx transporter) exhibited $d_c = 1.0000$ ($\mu_{\text{null}} = 1.0108, z_c = -0.0335, p_{\text{emp}} = 0.9426$).
   - *ABCG2* (efflux transporter) exhibited $d_c = 2.0000$ ($\mu_{\text{null}} = 1.0594, z_c = +2.0522, p_{\text{emp}} = 1.0000$).
   - The multi-gene statin pharmacogenomic set (*SLCO1B1*, *ABCG2*, *CYP2C9*) had an average distance of $d_c = 1.3333$ ($z_c = +1.6583$). These transporters govern systemic statin bioavailability and myopathy risk, but do not topologically cluster within coronary atheroma mechanisms.

3. **CYP2C19 and Clopidogrel Pharmacogenomics:**
   - *CYP2C19* had an observed distance of $d_c = 0.0000$ driven by its direct annotation in Open Targets CAD.
   - When *CYP2C19* was excluded from the CAD module in the overlap-excluded sensitivity analysis, its shortest path was $d_c = 1.0000$, aligning with the random degree expectation ($\mu = 0.9479$).

4. **Warfarin Indication-Aware Partitioning:**
   - The warfarin set (*CYP2C9*, *VKORC1*, *CYP4F2*) showed a neutral distance of $d_c = 1.0000$ ($\mu_{\text{null}} = 1.0433, z_c = -0.1824, p_{\text{emp}} = 0.7296$).
   - This empirically confirms the 2023 AHA/ACC guideline distinction: oral anticoagulation mechanisms do not directly colocalize with coronary artery atherogenesis.

## 2. Multi-Dimensional Sensitivity Analysis
The 101-condition sensitivity matrix confirmed high robustness across parameters:
- **Interactome Representation:** Results on the text-mining-excluded network (17,102 nodes) were virtually identical to the primary functional interactome (*CYP2C19* $z = -3.68$, statin set $z = +1.58$).
- **Physical Network:** Confirmed that *SLCO1B1*, *ABCG2*, and *CYP4F2* drop out of the high-confidence physical interactome, validating the requirement for functional network modeling.
- **Disease Concept:** Proximity patterns remained stable when evaluated against Myocardial Infarction (*CYP2C19* $z = -3.81$, statin set $z = +1.42$) and Myocardial Ischemia (*CYP2C19* $z = -3.42$, statin set $z = +1.21$).
- **Score Thresholds:** Stricter CAD score cutoffs ($\ge 0.20$ and $\ge 0.30$) reduced module size (to 670 and 393 LCC genes) and increased mean distances ($d_c$ shifted from 1.0 to 2.0 for transporters), but preserved relative topological rankings.
