# Week 9 Decision Log — Methodological Deviations & Rationale

| # | Topic | Guney et al. (2016) Protocol | Our Week 9 Adapted Protocol | Scientific Rationale |
|---|---|---|---|---|
| **1** | **Target Class Scope** | Excluded enzymes, carriers, and transporters; analyzed physical drug targets only. | Included metabolic enzymes (*CYP2C19*, *CYP2C9*, *CYP4F2*, *VKORC1*) and transporters (*SLCO1B1*, *ABCG2*). | Project objective is pharmacogenomic network context, not direct drug target efficacy. |
| **2** | **Interactome Type** | Experimentally documented physical PPI interactome. | STRING v12.5 Functional Interactome (Primary) + Physical Subnetwork (Comparator). | Physical networks exclude metabolic enzymes and transporters lacking stable physical complexes. |
| **3** | **Disease Module Definition** | Curated OMIM and GWAS disease genes. | Open Targets Platform v26.9.0 Direct Associations ($\ge 0.10$). | Incorporates modern clinical genetics, animal models, and drug clinical evidence under an explicit score cutoff. |
| **4** | **Warfarin Indication** | Treated all approved drugs against reported indication. | Separated warfarin into an indication-aware secondary analysis. | 2023 AHA/ACC guidelines distinguish anticoagulation indications (VTE/AFib) from CAD antiplatelet therapy. |
| **5** | **Null Model Candidate Pool** | Drew random nodes from entire interactome. | Pre-registered primary null excludes panel pharmacogenes from random background. | Prevents self-comparison bias when evaluating individual panel pharmacogenes. |
| **6** | **Statistical Significance Criterion** | Used heuristic thresholds ($z \le -0.15$ proximal; $z \le -2.0$ significant). | Pre-registered empirical permutation $p < 0.05$ with Benjamini-Hochberg FDR $q < 0.05$. | For discrete path distributions, $z \le -2.0$ does not guarantee $p < 0.05$ due to non-Gaussian tails. |
