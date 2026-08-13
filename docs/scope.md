# Project Scope

## Project

Network-Based Pharmacogenomic Risk Mapping for Ischemic Heart Disease / Cardiovascular Disease in South Asian Populations

## Current Status

Week 1 — scope lock in progress.

## Disease Scope

The project focuses on ischemic heart disease (IHD) within the broader cardiovascular disease (CVD) context.

The disease definition and overall disease scope are locked and should not be reopened without explicit supervisor approval.

## Population Scope

The primary population framing is India-focused.

South Asian populations from the 1000 Genomes Project will be retained as separate, unpooled population strata rather than being combined into a single South Asian population.

Relevant population strata include:

- GIH
- PJL
- BEB
- STU
- ITU

No silent pooling of these populations will be performed.

## Methodological Scope

The project will integrate:

1. Pharmacogenomic variant identification and annotation
2. Population-level allele-frequency analysis
3. Protein-protein interaction network analysis using STRING
4. Disease-module analysis using Open Targets
5. Pathway analysis using Reactome
6. Network centrality and disease-module proximity analysis
7. An equal-weight composite prioritization score
8. Independent validation and sensitivity analysis

The project is computational and hypothesis-generating.

No machine-learning model will be used as part of the core methodology.

## Drug and Gene Panel

The final drug/gene panel is pending supervisor confirmation.

Decision #1:

- Proposed scope: clopidogrel, warfarin, and statins
- Final inclusion requires supervisor confirmation
- Target decision week: Week 1

No additional drug classes will be added without explicit scientific justification and supervisor approval.

## Primary Research Questions

### Population Pharmacogenomics

Do pharmacogenomic variants relevant to the selected cardiovascular drug classes show meaningful allele-frequency differences across Indian and unpooled South Asian population strata?

### Network-Based Prioritization

Can integration of pharmacogenomic evidence, protein-protein interaction topology, disease-module proximity, and pathway information prioritize genes/variants that warrant further investigation in cardiovascular pharmacogenomics?

## Core Methodological Principles

- Population groups will remain unpooled.
- Publicly accessible aggregate data will be used.
- Individual-level clinical or genotype data will not be used.
- The analysis will be hypothesis-generating rather than clinically predictive.
- Statistical significance will not be treated as equivalent to clinical significance.
- Conflicting or uncertain evidence will be reported rather than silently resolved.
- All datasets, versions, access dates, and relevant licensing information will be documented.
- Analysis decisions and methodological changes will be recorded in the project decision log.

## Planned Computational Framework

The project will ultimately use a reproducible computational workflow incorporating:

- Python
- pandas
- NumPy
- requests
- NetworkX
- Snakemake
- STRING
- Open Targets
- Reactome
- Cytoscape where appropriate

Additional software will be introduced only when required by the relevant project stage.

## Scope Boundaries

The project does not currently include:

- Wet-lab genotyping
- Individual-level genotype analysis
- Individual patient clinical data
- Electronic health record analysis
- Clinical decision-support software development
- A machine-learning prediction model
- Raw sequencing variant calling
- Clinical validation of predicted pharmacogenomic risk

