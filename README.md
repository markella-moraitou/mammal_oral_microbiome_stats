Analysis workflow for taxonomic community-level analysis:
```mermaid
flowchart TD;
    A[Create phyloseq object]-->B[Access omnicrobe database];
    A-->C[Raw dataset plots];
    B-->D[Assess and remove contaminants];
    D-->E[Normalisation];
    D-->F[Sample summary plots];
    D-->G[Core microbiome];
    D-->G[Alpha diversity];
    E-->H[Filtered dataset plots];
    E-->I[RDA analysis];
    E-->J[Phylosymbiosis tests];
    E-->K[Microbial phenotypes - host diet hypotheses];
    E-->L[Statistical tests];
    E-->L[Machine learning];
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
    A-->D[codiversification_test]
```
