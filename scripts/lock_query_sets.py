#!/usr/bin/env python3
"""
scripts/lock_query_sets.py
Task 2: Lock Biological Query Sets & Generate Pharmacogene Classification Table
Conforms strictly to W9-1 SOP guidelines.
"""

import os
import json
import pandas as pd

OUT_DIR = "results/tables"
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Complete Evidence Table conforming to W9-1 schema
evidence_data = [
    {
        "HGNC symbol": "CYP2C19",
        "Ensembl ID": "ENSG00000165841",
        "STRING Protein ID": "9606.ENSP00000360372",
        "Drug association": "Clopidogrel",
        "Molecular role": "Cytochrome P450 monooxygenase (Enzyme)",
        "PGx role": "Pharmacokinetic (Bioactivation of prodrug to active thiol metabolite)",
        "CPIC evidence": "CPIC Level A (Lee et al., 2022 Guideline update)",
        "PharmGKB/ClinPGx evidence": "Level 1A (ClinPGx PA124; PA166251443)",
        "Clinical indication": "Antiplatelet therapy for acute coronary syndromes (ACS) and percutaneous coronary intervention (PCI)",
        "Target vs PGx distinction": "Pharmacogene (Metabolic enzyme). Direct pharmacological target is P2RY12.",
        "Source": "Lee et al. Clin Pharmacol Ther 2022;112(5):959-967; CPIC; PharmGKB"
    },
    {
        "HGNC symbol": "CYP2C9",
        "Ensembl ID": "ENSG00000138109",
        "STRING Protein ID": "9606.ENSP00000260682",
        "Drug association": "Warfarin; Fluvastatin",
        "Molecular role": "Cytochrome P450 monooxygenase (Enzyme)",
        "PGx role": "Pharmacokinetic (Primary clearance of active S-warfarin; fluvastatin metabolism)",
        "CPIC evidence": "CPIC Level A (Johnson et al., 2017 Warfarin; Cooper-DeHoff et al., 2022 Statins)",
        "PharmGKB/ClinPGx evidence": "Level 1A (ClinPGx PA126; PA166251465; PA166264281)",
        "Clinical indication": "Oral anticoagulation (Warfarin for VTE/AFib/mechanical valves); Dyslipidemia (Fluvastatin)",
        "Target vs PGx distinction": "Pharmacogene (Metabolic clearance enzyme). Direct target of warfarin is VKORC1.",
        "Source": "Johnson et al. Clin Pharmacol Ther 2017;102(3):397-404; Cooper-DeHoff et al. CPT 2022"
    },
    {
        "HGNC symbol": "VKORC1",
        "Ensembl ID": "ENSG00000167397",
        "STRING Protein ID": "9606.ENSP00000378426",
        "Drug association": "Warfarin",
        "Molecular role": "Vitamin K epoxide reductase complex subunit 1 (Reductase)",
        "PGx role": "Pharmacodynamic (Sensitivity to vitamin K cycle inhibition)",
        "CPIC evidence": "CPIC Level A (Johnson et al., 2017 Warfarin guideline)",
        "PharmGKB/ClinPGx evidence": "Level 1A (ClinPGx PA133787052; PA166251465)",
        "Clinical indication": "Oral anticoagulation (VTE prophylaxis/treatment, non-valvular AFib)",
        "Target vs PGx distinction": "Dual role: Pharmacological drug target of warfarin AND pharmacodynamic pharmacogene.",
        "Source": "Johnson et al. Clin Pharmacol Ther 2017;102(3):397-404; PharmGKB"
    },
    {
        "HGNC symbol": "CYP4F2",
        "Ensembl ID": "ENSG00000186529",
        "STRING Protein ID": "9606.ENSP00000221700",
        "Drug association": "Warfarin",
        "Molecular role": "Cytochrome P450 monooxygenase (Enzyme)",
        "PGx role": "Pharmacokinetic/Pharmacodynamic (Vitamin K1 oxidase; clearance of vitamin K)",
        "CPIC evidence": "CPIC Level A (Johnson et al., 2017 Warfarin guideline)",
        "PharmGKB/ClinPGx evidence": "Level 1A (ClinPGx PA27121; PA166251465)",
        "Clinical indication": "Oral anticoagulation (Modulates dose requirement via vitamin K catabolism)",
        "Target vs PGx distinction": "Pharmacogene (Catabolic enzyme in vitamin K pathway).",
        "Source": "Johnson et al. Clin Pharmacol Ther 2017;102(3):397-404; McDonald et al. Blood 2009"
    },
    {
        "HGNC symbol": "SLCO1B1",
        "Ensembl ID": "ENSG00000134538",
        "STRING Protein ID": "9606.ENSP00000256958",
        "Drug association": "Statins (Simvastatin, Atorvastatin, Rosuvastatin, Pravastatin, Pitavastatin)",
        "Molecular role": "Solute carrier organic anion transporter 1B1 / OATP1B1 (Transporter)",
        "PGx role": "Pharmacokinetic (Hepatic uptake transporter; determines systemic statin exposure)",
        "CPIC evidence": "CPIC Level A (Cooper-DeHoff et al., 2022 Statin guideline update)",
        "PharmGKB/ClinPGx evidence": "Level 1A (ClinPGx PA134865839; PA166264281)",
        "Clinical indication": "Atherosclerotic cardiovascular disease (ASCVD) primary/secondary prevention; lipid-lowering",
        "Target vs PGx distinction": "Pharmacogene (Influx transporter). Direct pharmacological target of statins is HMGCR.",
        "Source": "Cooper-DeHoff et al. Clin Pharmacol Ther 2022;111(5):1007-1021; CPIC"
    },
    {
        "HGNC symbol": "ABCG2",
        "Ensembl ID": "ENSG00000118777",
        "STRING Protein ID": "9606.ENSP00000498246",
        "Drug association": "Rosuvastatin, Atorvastatin",
        "Molecular role": "ATP-binding cassette transporter / BCRP (Efflux Transporter)",
        "PGx role": "Pharmacokinetic (Intestinal/hepatic efflux transporter; modulates systemic statin bioavailability)",
        "CPIC evidence": "CPIC Level A (Cooper-DeHoff et al., 2022 Statin guideline update)",
        "PharmGKB/ClinPGx evidence": "Level 1A (ClinPGx PA390; PA166264281)",
        "Clinical indication": "Atherosclerotic cardiovascular disease lipid management",
        "Target vs PGx distinction": "Pharmacogene (Efflux transporter). Direct pharmacological target of statins is HMGCR.",
        "Source": "Cooper-DeHoff et al. Clin Pharmacol Ther 2022;111(5):1007-1021; CPIC"
    }
]

df_evidence = pd.DataFrame(evidence_data)
df_evidence.to_csv(os.path.join(OUT_DIR, "pharmacogene_evidence_table.csv"), index=False)
print(f"[+] Saved pharmacogene evidence table to {OUT_DIR}/pharmacogene_evidence_table.csv")

# 2. Structured Query Sets conforming to W9-1
query_sets = {
    "level_A_individual": {
        "CYP2C19": ["CYP2C19"],
        "CYP2C9": ["CYP2C9"],
        "VKORC1": ["VKORC1"],
        "CYP4F2": ["CYP4F2"],
        "SLCO1B1": ["SLCO1B1"],
        "ABCG2": ["ABCG2"]
    },
    "level_B_drug_specific": {
        "clopidogrel_antiplatelet": {
            "genes": ["CYP2C19"],
            "indication": "Coronary artery disease / ACS / Post-PCI antiplatelet therapy",
            "direct_targets_comparator": ["P2RY12"]
        },
        "statins_lipid_lowering": {
            "genes": ["SLCO1B1", "ABCG2", "CYP2C9"],
            "indication": "Coronary artery disease / Atherosclerosis prevention",
            "direct_targets_comparator": ["HMGCR"]
        },
        "warfarin_anticoagulation_separate": {
            "genes": ["CYP2C9", "VKORC1", "CYP4F2"],
            "indication": "Venous thromboembolism / Atrial fibrillation anticoagulation (Separate from primary CAD antiplatelet therapy)",
            "direct_targets_comparator": ["VKORC1"]
        }
    },
    "level_C_exploratory_whole_panel": {
        "all_6_pharmacogenes": ["CYP2C19", "CYP2C9", "VKORC1", "CYP4F2", "SLCO1B1", "ABCG2"]
    }
}

with open(os.path.join(OUT_DIR, "query_sets_locked.json"), "w") as f:
    json.dump(query_sets, f, indent=2)

print(f"[+] Saved structured query sets to {OUT_DIR}/query_sets_locked.json")
print("\n=== SUMMARY OF LOCKED QUERY SETS ===")
for level, content in query_sets.items():
    print(f"\n{level}:")
    if level == "level_B_drug_specific":
        for k, v in content.items():
            print(f"  - {k}: Genes={v['genes']} | Direct Target Comparator={v['direct_targets_comparator']}")
    else:
        for k, v in content.items():
            print(f"  - {k}: {v}")
