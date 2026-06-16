Analysis workflow for taxonomic community-level analysis:
```mermaid
flowchart TD;
    A[create_phyloseq_obj]-->B[collect_omnicrobe_metadata];
    A-->C[raw_data_plots];
    B-->D[assess_contamination];
    D-->E[normalisations];
    D-->F[sample_summary_plots];
    D-->G[core_microbiomes];
    D-->G[alpha_div];
    E-->H[multivariate_stats];
    E-->I[microbial_phenotypes];
    E-->J[models];
    E-->K[machine_learning];
    E-->L[diff_abund];
```

Analysis workflow for functional community-level analysis:
```mermaid
flowchart TD;
    A[prepare_data]-->B[pathway_completeness];
    B-->C[diff_abund_path];
    B-->D[functional_core];
    A-->F[distillR];
    B,F-->G[func_tax_procrustes];
    B,F-->H[multivariate_stats_func];
    B,F-->I[func_models];
    A-->J[nitrate_reduction];
```

Analysis workflow for metagenomes-assembled genomes analysis:

```mermaid # Not sure if accurate, double check
flowchart TD;
    A[get_mag_tree_and_metadata]-->B[collect_omnicrobe_metadata]
    A-->C[mag_mapping_stats]
    A-->H[gsea_da_taxa];
    B-->D[plot_mag_tree]
    B-->E[codiversification_test]
    E-->F[mag_annotations]
    E-->G[clades_for_popgen]
```
