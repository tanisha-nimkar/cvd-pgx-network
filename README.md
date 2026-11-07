# Cardiovascular Pharmacogenomics Prioritization Network (cvd-pgx-network)

An end-to-end computational biology workflow for population-stratified variant prioritization across a 6-gene cardiovascular pharmacogenomics panel (CYP2C19, CYP2C9, VKORC1, CYP4F2, SLCO1B1, ABCG2).

## Pipeline Execution
```bash
# Activate environment
conda activate cvd-pgx

# Reproduce entire pipeline
snakemake --cores 1
```

## Key Publication Figures
- Fig 2: South Asian vs. NFE allele frequency comparison (`figures/fig2_freq_comparison.png`)
- Fig 3 & 4: STRING v12.5 first-degree interactome network (`figures/fig4_cytoscape_network.png`)
- Fig 5: Network proximity null distributions (`figures/fig5_proximity.png`)
- Fig 6: Top 30 prioritized variants stacked scores (`figures/fig6_ranking.png`)
- Fig 7: Sensitivity analysis across 593 scenarios (`figures/fig7_sensitivity.png`)

## Status & Governance
- Workflow: Snakemake 9.22.0 (20 rules, fully deterministic)
- Verification: Clean-clone replication passed (100% top-10 candidate variant recovery)
- Git Tag: week12-frozen
- License: MIT Open Academic License
