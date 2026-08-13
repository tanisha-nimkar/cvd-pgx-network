# Data Access Status

This document tracks access to external databases and resources required by the project.

No external data has been downloaded for analysis as of Day 1.

## Database Status

| Database / Resource | Purpose | Access Status | Action Required |
|---|---|---|---|
| PharmGKB | Pharmacogenomic evidence and variants | Not yet checked | Check in Week 2 |
| CPIC | Clinical pharmacogenomic guidelines | Not yet checked | Check in Week 2 |
| ClinVar | Variant clinical significance | Not yet checked | Check in Week 3 |
| IndiGenomes | India-specific population genomic data | Needs investigation | Check access mode in Week 1–2 |
| 1000 Genomes Project | South Asian population frequencies | Not yet checked | Check in Week 2–5 |
| gnomAD | Population allele frequencies / validation | Not yet checked | Check in Week 2–5 |
| STRING | Protein-protein interaction network | Not yet checked | Week 7 |
| Open Targets | Disease-gene associations / disease module | Not yet checked | Week 8 |
| Reactome | Pathway analysis | Not yet checked | Week 11 |
| BioGRID | Independent PPI validation | Not yet checked | Week 11 |

## IndiGenomes

**Status:** Needs investigation.

The immediate goal is to determine:

1. Whether the required data are publicly accessible.
2. Whether bulk download or API access is available.
3. What population data are available.
4. What genome build is used.
5. What licensing or usage restrictions apply.
6. Whether the data can be used for this project.
7. Whether derived results can be publicly released.

Do not assume access or redistribution rights without checking.

## Access Rules

For each resource, record:

- Access date
- Version/release
- URL
- API version where relevant
- License/usage restrictions
- Download method
- Any authentication requirement

No API keys, passwords, or access tokens should ever be committed to GitHub.
