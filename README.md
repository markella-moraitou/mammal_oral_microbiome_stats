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
    E-->J[phylosymbiosis_test];
    E-->K[microbial_phenotypes];
    E-->L[models];
    E-->L[machine_learning];
```

Analysis workflow for functional community-level analysis:
```mermaid
flowchart TD;
    A[Prepare data]-->B[Plot ordinations];
    A-->C[Exploratory plots];
    A-->D[Pathway analysis];
    A-->E[DistillR];
    D-->F[Functional core microbiome];
```

Analysis workflow for metagenomes-assembled genomes analysis:

```mermaid # Not sure if accurate, double check
flowchart TD;
    A[get_mag_tree_and_metadata]-->B[collect_omnicrobe_metadata]
    A-->C[mag_mapping_stats]
    B-->D[plot_mag_tree]
    B-->E[codiversification_test]
    A-->F[mag_annotations]
```
