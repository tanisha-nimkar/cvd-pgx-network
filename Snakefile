# =============================================================================
# Snakefile — Cardiovascular Pharmacogenomics Prioritization Pipeline
# Framework: Snakemake >= 8.0 (Tested on 9.22.0)
# Controlling Protocol: Week 12 Reproducibility Specifications
# =============================================================================

configfile: "config/config.yaml"

rule all:
    input:
        # Week 3-4 Variant Annotation Pipeline
        "results/tables/variants_annotated_week3_backup.csv",
        "data/processed/vep_input.vcf",
        "results/tables/variants_annotated.csv",
        # Week 5-6 Population Stratification
        "data/processed/population_frequencies_raw.csv",
        "results/tables/population_frequencies.csv",
        "figures/fig2_freq_comparison.png",
        # Week 7-8 Interactome & Disease Module Extraction
        "results/network/network.graphml",
        "results/tables/centrality.csv",
        "data/processed/string_functional_tm_excluded_edges.tsv",
        "results/tables/disease_module_final.csv",
        # Week 9 Network Proximity Analysis & Figures
        "results/tables/network_metrics.csv",
        "figures/fig3_network.png",
        "figures/fig5_proximity.png",
        # Week 10 Composite Prioritization & Sensitivity Grid
        "results/tables/prioritization_ranked.csv",
        "results/tables/sensitivity_analysis.csv",
        "results/tables/sensitivity_summary.csv",
        "figures/fig6_ranking.png",
        "figures/fig7_sensitivity.png",
        # Week 11 Multi-Stream Independent Validation
        "data/processed/biogrid_edges.tsv",
        "results/tables/biogrid_network_metrics.csv",
        "results/tables/validation_detail_biogrid.csv",
        "results/tables/validation_detail_gnomad.csv",
        "results/tables/validation_detail_pubmed.csv",
        "results/tables/pathway_enrichment.csv",
        "results/tables/validation_results.csv"

# -----------------------------------------------------------------------------
# Stage 1: Variant Harmonization (Week 3)
# -----------------------------------------------------------------------------
rule merge_variants:
    input:
        manifest = "data/raw_manifest.tsv",
        clinpgx = expand("data/raw/clinpgx_{gene}_20260830.json", gene=config["gene_panel"]),
        cpic = expand("data/raw/cpic_{gene}_20260830.json", gene=config["gene_panel"]),
        clinvar = expand("data/raw/clinvar_{gene}_20260830.vcf", gene=config["gene_panel"])
    output:
        w3_backup = "results/tables/variants_annotated_week3_backup.csv",
        cpic_tsv = "results/tables/annotated_cpic_pairs.tsv",
        clinvar_tsv = "results/tables/annotated_clinvar_variants.tsv",
        master_tsv = "results/tables/master_variant_annotations.tsv"
    log:
        "logs/merge_variants.log"
    shell:
        "python scripts/merge_annotate_variants.py > {log} 2>&1"

# -----------------------------------------------------------------------------
# Stage 2: Functional Annotation Preparation & VEP Integration (Week 4)
# -----------------------------------------------------------------------------
rule prepare_vep_input:
    input:
        w3_backup = "results/tables/variants_annotated_week3_backup.csv"
    output:
        vep_vcf = "data/processed/vep_input.vcf"
    log:
        "logs/prepare_vep_input.log"
    shell:
        "python scripts/prepare_vep_input.py > {log} 2>&1"

rule annotate_vep:
    input:
        vep_vcf = "data/processed/vep_input.vcf",
        vep_output = "data/raw/vep_output_20261031.txt",
        w3_backup = "results/tables/variants_annotated_week3_backup.csv"
    output:
        annotated = "results/tables/variants_annotated.csv"
    log:
        "logs/annotate_vep.log"
    shell:
        "python scripts/annotate_variants_vep.py > {log} 2>&1"

# -----------------------------------------------------------------------------
# Stage 3: Population Allele Frequencies & Statistics (Weeks 5-6)
# -----------------------------------------------------------------------------
rule prepare_pop_raw:
    input:
        cache = "data/processed/ensembl_pop_cache_single.json",
        variants = "results/tables/variants_annotated.csv"
    output:
        pop_raw = "data/processed/population_frequencies_raw.csv"
    log:
        "logs/prepare_pop_raw.log"
    shell:
        "python scripts/prepare_population_frequencies_raw.py > {log} 2>&1"

rule frequency_stats:
    input:
        pop_raw = "data/processed/population_frequencies_raw.csv",
        cfg = "config/config.yaml"
    output:
        pop_freq = "results/tables/population_frequencies.csv",
        fig2 = "figures/fig2_freq_comparison.png"
    log:
        "logs/frequency_stats.log"
    shell:
        "python scripts/frequency_stats.py > {log} 2>&1"

# -----------------------------------------------------------------------------
# Stage 4: Network Construction & Interactome Universes (Weeks 7-8)
# -----------------------------------------------------------------------------
rule build_network:
    input:
        cfg = "config/config.yaml"
    output:
        graphml = "results/network/network.graphml",
        audit = "results/network/string_edges_audit.tsv",
        mapping = "results/network/string_id_mapping.tsv"
    log:
        "logs/build_network.log"
    shell:
        "python scripts/build_network.py > {log} 2>&1"

rule network_centrality:
    input:
        graphml = "results/network/network.graphml"
    output:
        centrality = "results/tables/centrality.csv"
    log:
        "logs/network_centrality.log"
    shell:
        "python scripts/network_centrality.py > {log} 2>&1"

rule build_string_universe:
    input:
        string_links = "data/raw/9606.protein.links.full.v12.5.txt.gz",
        string_info = "data/raw/9606.protein.info.v12.5.txt.gz"
    output:
        universe_edges = "data/processed/string_functional_tm_excluded_edges.tsv"
    log:
        "logs/build_string_universe.log"
    shell:
        "python scripts/build_interactome_universes.py > {log} 2>&1"

rule open_targets_module:
    input:
        cfg = "config/config.yaml",
        universe_edges = "data/processed/string_functional_tm_excluded_edges.tsv"
    output:
        candidates = "results/tables/disease_module_candidates.csv",
        comparison = "results/tables/disease_module_comparison.csv"
    log:
        "logs/open_targets_module.log"
    shell:
        "python scripts/open_targets_module.py > {log} 2>&1"

# -----------------------------------------------------------------------------
# Stage 5: Proximity Analysis & Visualization (Week 9)
# -----------------------------------------------------------------------------
rule network_proximity:
    input:
        graphml = "results/network/network.graphml",
        centrality = "results/tables/centrality.csv",
        module = "results/tables/disease_module_final.csv",
        universe = "data/processed/string_functional_tm_excluded_edges.tsv",
        cfg = "config/config.yaml"
    output:
        metrics = "results/tables/network_metrics.csv",
        null_cyp2c19 = "results/tables/proximity_null_CYP2C19.csv",
        null_cyp2c9 = "results/tables/proximity_null_CYP2C9.csv",
        null_vkorc1 = "results/tables/proximity_null_VKORC1.csv",
        null_cyp4f2 = "results/tables/proximity_null_CYP4F2.csv",
        null_slco1b1 = "results/tables/proximity_null_SLCO1B1.csv",
        null_abcg2 = "results/tables/proximity_null_ABCG2.csv",
        null_adgrg6 = "results/tables/proximity_null_ADGRG6.csv"
    log:
        "logs/network_proximity.log"
    shell:
        "python scripts/network_proximity.py > {log} 2>&1"

rule plot_w9_figures:
    input:
        graphml = "results/network/network.graphml",
        metrics = "results/tables/network_metrics.csv"
    output:
        fig3 = "figures/fig3_network.png",
        fig5 = "figures/fig5_proximity.png"
    log:
        "logs/plot_w9_figures.log"
    shell:
        "python scripts/plot_w9_figures.py > {log} 2>&1"

# -----------------------------------------------------------------------------
# Stage 6: Composite Prioritization & Sensitivity Analysis (Week 10)
# -----------------------------------------------------------------------------
rule prioritization_score:
    input:
        variants = "results/tables/variants_annotated.csv",
        pop_freq = "results/tables/population_frequencies.csv",
        metrics = "results/tables/network_metrics.csv",
        cpic = "results/tables/annotated_cpic_pairs.tsv"
    output:
        ranked = "results/tables/prioritization_ranked.csv",
        sens_long = "results/tables/sensitivity_analysis.csv",
        sens_summ = "results/tables/sensitivity_summary.csv"
    log:
        "logs/prioritization_score.log"
    shell:
        "python scripts/prioritization_score.py > {log} 2>&1"

rule plot_w10_figures:
    input:
        ranked = "results/tables/prioritization_ranked.csv",
        sens_long = "results/tables/sensitivity_analysis.csv"
    output:
        fig6 = "figures/fig6_ranking.png",
        fig7 = "figures/fig7_sensitivity.png"
    log:
        "logs/plot_w10_figures.log"
    shell:
        "python scripts/plot_w10_figures.py > {log} 2>&1"

# -----------------------------------------------------------------------------
# Stage 7: Independent Validation & Plausibility Enrichment (Week 11)
# -----------------------------------------------------------------------------
rule prepare_biogrid:
    input:
        raw_biogrid = "data/raw/BIOGRID-ORGANISM-Homo_sapiens-5.0.263.tab3.txt"
    output:
        bg_edges = "data/processed/biogrid_edges.tsv"
    log:
        "logs/prepare_biogrid.log"
    shell:
        "python scripts/prepare_biogrid.py > {log} 2>&1"

rule biogrid_validation:
    input:
        bg_edges = "data/processed/biogrid_edges.tsv",
        module = "results/tables/disease_module_final.csv",
        string_metrics = "results/tables/network_metrics.csv",
        ranked = "results/tables/prioritization_ranked.csv"
    output:
        bg_metrics = "results/tables/biogrid_network_metrics.csv",
        bg_detail = "results/tables/validation_detail_biogrid.csv"
    log:
        "logs/biogrid_validation.log"
    shell:
        "python scripts/biogrid_validation.py > {log} 2>&1"

rule gnomad_validation:
    input:
        ranked = "results/tables/prioritization_ranked.csv",
        variants = "results/tables/variants_annotated.csv"
    output:
        gnomad_detail = "results/tables/validation_detail_gnomad.csv"
    log:
        "logs/gnomad_validation.log"
    shell:
        "python scripts/gnomad_validation.py > {log} 2>&1"

rule pubmed_validation:
    output:
        pubmed_detail = "results/tables/validation_detail_pubmed.csv"
    log:
        "logs/pubmed_validation.log"
    shell:
        "python scripts/pubmed_validation.py > {log} 2>&1"

rule pathway_enrichment:
    input:
        graphml = "results/network/network.graphml",
        universe = "data/processed/string_functional_tm_excluded_edges.tsv",
        gmt = "data/raw/ReactomePathways.gmt"
    output:
        pathways = "results/tables/pathway_enrichment.csv"
    log:
        "logs/pathway_enrichment.log"
    shell:
        "python scripts/pathway_enrichment.py > {log} 2>&1"

rule independent_validation:
    input:
        bg = "results/tables/validation_detail_biogrid.csv",
        gn = "results/tables/validation_detail_gnomad.csv",
        pm = "results/tables/validation_detail_pubmed.csv"
    output:
        val_results = "results/tables/validation_results.csv"
    log:
        "logs/independent_validation.log"
    shell:
        "python scripts/independent_validation.py > {log} 2>&1"
