##### DIFFERENTIAL ABUNDANCE OF FUNCTIONS AND FUNCTIONS OF DIFFERENTIALLY ABUNDANT TAXA #####

#### Differential abundance analysis of KEGG pathways
#### and functional enrichment of the differentially abundant taxa identified with MCMCglmm and PGLMM/Pagel's lambda

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(phyloseq)
library(phyr)
library(phytools)
library(ggplot2)
library(clusterProfiler)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "function")) # subdirectory for the output of this script
taxdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Directory with taxonomy analysis
pathdir <- normalizePath(file.path("..", "..", "output", "function", "pathway_completeness")) # Directory with pathway analysis output
datadir <- normalizePath(file.path(outdir, "data")) # Directory with data files
subdir <- normalizePath(file.path(outdir, "diff_abund_taxa")) # subdirectory for the output of this script

dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

# Get ordination functions
source(file.path("..","ordination_functions.R"))

#######################
#####  LOAD INPUT #####
#######################

# Phyloseq objects
phy_gene_f <- readRDS(file.path(datadir, "phy_gene_f.RDS"))
phy_gene_f_clr <- readRDS(file.path(datadir, "phy_gene_f_clr.RDS"))

phy_sp_f <- readRDS(file.path(taxdir, "phyloseq_objects", "phy_sp_f.RDS"))

phy_pathway <- readRDS(file.path(pathdir, "phy_pathway.RDS"))
phy_pathway_clr <- readRDS(file.path(pathdir, "phy_pathway_clr.RDS"))

# Differentially abundant taxa
pglmm_taxa_res <- read.csv(file.path(taxdir, "differential_abundance", "combined_phylo_ecology_results.csv"))

# Stratified sample data
gene_str <- read.table(file.path(datadir, "gene_abundance_stratified_modified.tsv"),
                      quote = "", comment.char = "", header = TRUE, sep = "\t")

# Host phylogeny
host_consensus <- read.tree(file.path(taxdir, "host_consensus.tre"))

#########################
#### PREP TEST INPUT ####
#########################

# Combine Sirenia/Proboscidea in one clades
phy_pathway_clr@sam_data$Order <- ifelse(phy_pathway_clr@sam_data$Order %in% c("Sirenia", "Proboscidea"), "Sirenia/Proboscidea", as.character(phy_pathway_clr@sam_data$Order))

# Turn habitat to factor
phy_pathway_clr@sam_data$habitat.general <- factor(phy_pathway_clr@sam_data$habitat.general, levels = c("Terrestrial", "Marine"))

# Test for ruminant differences
phy_pathway_clr@sam_data$ruminant <- factor(ifelse(phy_pathway_clr@sam_data$digestion == "Ruminant", "Ruminant", "Other"), levels = c("Other", "Ruminant"))

# Test for hypsodont differences
phy_pathway_clr@sam_data$hypsodont <- factor(ifelse(grepl("hyps", phy_pathway_clr@sam_data$molar_category), "Hypsodont", "Other"), levels = c("Other", "Hypsodont"))

# Melt
data <- psmelt(phy_pathway_clr) %>%
        select(OTU, Abundance, Sample, Species, Order, diet.general, habitat.general, ruminant, Fruit, Animal, Fruit, Seed)

data <- data %>%
    # Turn S. scrofa domesticus to S. scrofa to match tree
    mutate(Species = case_when(Species == "Sus domesticus" ~ "Sus scrofa",
                               TRUE ~ Species))

# Phylogeny
host_consensus$node.label <- paste0("node", c(1:length(host_consensus$node.label)))
host_consensus$tip.label <- gsub("_", " ", host_consensus$tip.label)
Ainv <- inverseA(host_consensus)$Ainv

#############################
#### PATHWAY DIFF ABUND #####
#############################

# Approach inspired by Youngblut et al. 2019

# list of all unique genera
paths <- taxa_names(phy_pathway_clr)

#### PGLMM: Test for ecology after accounting for phylogeny
pglmm_res <- data.frame(term = character(), 
                        coef = numeric(),
                        pval = numeric(),
                        pathway = character())

for (p in paths) {
    cat("Running PGLMM for OTU", which(paths == p), "of", length(paths), "--", p, "\n")
    # Subset data for the current OTU
    data_filt <- data %>% filter(OTU == p)
    # Run PGLMM
    model <- pglmm(Abundance ~ Animal + Fruit + habitat.general + ruminant +
                   (1 | Species__), 
                   data = data_filt, 
                   cov_ranef = list(Species = host_consensus),
                   family = "gaussian")
    # Extract results
    res <- cbind(model$B, model$B.pvalue) %>% as.data.frame %>%
            rownames_to_column %>% filter(rowname != "(Intercept)") %>%
            mutate(rowname = str_remove(str_remove(rowname, "habitat.general"), "ruminant"))
    colnames(res) <- c("term", "coef", "pval")
    # Add random phylogenetic effect
    #res <- rbind(res, data.frame(term = "phylogeny", coef = unname(model$s2r[2]), pval = NA))
    res$pathway <- p
    # Combine results
    pglmm_res <- rbind(pglmm_res, res)
}

# Make wide
pglmm_wide <- pglmm_res %>%
    pivot_wider(names_from = term, values_from = c(coef, pval))

write.csv(pglmm_wide, file = file.path(subdir, "pglmm_pathway_results.csv"), quote = FALSE, row.names = FALSE)

#### Phylogenetic signal after regressing out ecology

phy_res <- data.frame(pathway = character(),
                       lambda = numeric(),
                       pval = numeric())

for (p in paths) {
    cat("Running phylogenetic signal for OTU", which(paths == p), "of", length(paths), "--", p, "\n")
    # Subset data for the current OTU
    data_filt <- data %>% filter(OTU == p)
    # Run PGLMM
    model <- lm(Abundance ~ Animal + Fruit + habitat.general + ruminant,
                   data = data_filt)
    # Extract residuals
    resids_df <- data.frame(Sample = data_filt$Sample,
                            Species = data_filt$Species,
                            residuals = residuals(model))
    resids_by_species <- resids_df %>%
                        group_by(Species) %>%
                        summarise(residuals = mean(residuals)) %>% as.data.frame %>%
                        column_to_rownames("Species")
    # Make vector
    resids_vec <- resids_by_species$residuals
    names(resids_vec) <- rownames(resids_by_species)
    # Calculate Pagel's lamda
    pagel <- phylosig(host_consensus, resids_vec, method="lambda", test = TRUE)
    res <- data.frame(pathway = p,
                      lambda = pagel$lambda,
                      pval = pagel$P)
    phy_res <- rbind(phy_res, res)
}

write.csv(phy_res, file = file.path(subdir, "phylogenetic_results.csv"), quote = FALSE, row.names = FALSE)

#### Combine results ####
combined_res <- phy_res %>% rename(coef_lambda = lambda,
                                   pval_lambda = pval) %>%
    full_join(pglmm_wide, by = "pathway") %>%
    # Add metadata
    left_join(unique(rownames_to_column(data.frame(phy_pathway@tax_table), "pathway")))

# Adjust p-values
combined_res <- combined_res %>%
    mutate(padj_Animal = p.adjust(pval_Animal, method = "holm"),
           padj_Fruit = p.adjust(pval_Fruit, method = "holm"),
           padj_Marine = p.adjust(pval_Marine, method = "holm"),
           padj_Ruminant = p.adjust(pval_Ruminant, method = "holm"),
           padj_lambda = p.adjust(pval_lambda, method = "holm"))

write.csv(combined_res, file = file.path(subdir, "combined_phylo_ecology_pathway_results.csv"), quote = TRUE, row.names = FALSE)

combined_coef <- select(pivot_longer(combined_res, cols = starts_with("coef_"), names_to = "term", values_to = "coefficient"), "pathway", "term", "coefficient") %>%
    mutate(term = str_remove(term, "coef_"))

combined_pval <- select(pivot_longer(combined_res, cols = starts_with("pval_"), names_to = "term", values_to = "pval"), "pathway", "term", "pval") %>%
    mutate(term = str_remove(term, "pval_"))

combined_padj <- select(pivot_longer(combined_res, cols = starts_with("padj_"), names_to = "term", values_to = "padj"), "pathway", "term", "padj", "path_name", "path_class") %>%
    mutate(term = str_remove(term, "padj_"))

combined_long <- full_join(combined_coef, combined_pval, by = c("pathway", "term")) %>%
    full_join(combined_padj, by = c("pathway", "term")) %>%
    mutate(term = factor(term, levels = c("lambda", "Animal", "Fruit", "Ruminant", "Marine")),
           path_class = str_remove(path_class, ".*; "),
           significant = ifelse(pval < 0.05, ifelse(padj < 0.05, "yes", "pre-adjustment"), "no"))

# Identify classes with the most significant pathways
most_sign <- combined_long %>%
    group_by(path_class) %>% summarise(n_signif = sum(significant == "yes")) %>%
    arrange(desc(n_signif)) %>% slice_head(n = 6) %>% pull(path_class)

combined_long$path_class <- factor(ifelse(combined_long$path_class %in% most_sign,
                                          combined_long$path_class,
                                          "Other"),
                                    levels = rev(c(most_sign, "Other")))

# Plot
p <- ggplot(combined_long, aes(x = coefficient, y = path_class, colour = path_class, shape = significant)) +
    geom_point(alpha = 0.7, size = 1.5, position = position_jitter(height = 0.2)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    facet_grid(cols = vars(term), scales = "free_x") +
    scale_shape_manual(values = c("yes" = 19, "pre-adjustment" = 6, "no" = 4), name = "Significant (padj < 0.05)") +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 9), legend.title = element_text(size = 10),
          axis.title.y = element_blank(), axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 6)) +
    guides(colour = "none",
           shape = guide_legend(nrow = 2, byrow = TRUE))

ggsave(p, filename = file.path(subdir, "combined_phylo_ecology_pathways.png"), width = 7, height = 5)

####################
#### FOCAL TAXA ####
####################

# Identify differentially abundant taxa and their association

# From both outputs, extract significant taxa, and what their association is (positive or negative)
da_taxa <- pglmm_taxa_res %>% filter(padj_Animal < 0.05 | padj_Fruit < 0.05 | padj_Marine < 0.05 | padj_Ruminant < 0.05) %>%
            # label association (positive and negative)
            mutate(label = case_when(padj_Animal < 0.05 & coef_Animal > 0 ~ "Animal +",
                                     padj_Animal < 0.05 & coef_Animal < 0 ~ "Animal -",
                                     padj_Fruit < 0.05 & coef_Fruit > 0 ~ "Fruit +",
                                     padj_Fruit < 0.05 & coef_Fruit < 0 ~ "Fruit -",
                                     padj_Marine < 0.05 & padj_Marine > 0 ~ "Marine +",
                                     padj_Marine < 0.05 & padj_Marine < 0 ~ "Marine -",
                                     padj_Ruminant < 0.05 & padj_Ruminant > 0 ~ "Ruminant +",
                                     padj_Ruminant < 0.05 & padj_Ruminant < 0 ~ "Ruminant -")) %>%
            separate(label, into = c("term", "association"), sep = " ") %>% select(OTU, term, association)

#mcmc_sig <- mcmc_res %>% filter(pMCMC < 0.05) %>%
#    rename(taxon = OTU) %>%
#    mutate(association = ifelse(post.mean > 0, "+", "-")) %>%
#    select(taxon, term, association) %>%
#    mutate(term = str_remove(term, "ruminant|habitat.general")) %>%
#    mutate(method = "MCMCglmm")

#da_taxa <- full_join(ancom_sig, mcmc_sig) %>% arrange(taxon)

tax_info <- data.frame(phy_sp_f@tax_table) %>% select(phylum, family, genus) %>% unique %>%
    filter(genus %in% da_taxa$OTU)

da_taxa$phylum <- tax_info$phylum[match(da_taxa$OTU, tax_info$genus)]

# Keep only positive associations at least for now
da_pos <- da_taxa %>% filter(association == "+") %>% select(OTU, term, association) %>% unique

# Add taxa from RDA results manually
da_rda <- rbind(c("Desulfovibrio desulfuricans_D", "Ruminant", "+"),
                c("Actinomyces dentalis", "Ruminant", "+"),
                c("Actinomyces glycerinitolerans", "Ruminant", "+"),
                c("Actinomyces lilanjuaniae", "Ruminant", "+"),
                c("Propionibacterium australiense", "Ruminant", "+"),
                c("Actinomyces qiguomingii", "Ruminant", "+"),
                c("Actinomyces glycerinitolerans", "Ruminant", "+"),
                c("CAJPSE01 sp003860125", "Animal", "+"),
                c("Tannerella forsythia", "Animal", "+"),
                c("Tannerella forsythia_A", "Animal", "+"),
                c("Johnsonella sp03052905", "Animal", "+"),
                c("F0058 sp000768855", "Animal", "+")) %>%
            as.data.frame() %>%
            setNames(c("OTU", "term", "association"))

da_pos_all <- bind_rows(da_pos, da_rda) %>% unique

##############################
#### ENRICHMENT ANALYSIS #####
##############################

# As background, use all genes in the dataset

K_all <- gene_str %>% filter(database == "KEGG" & genus %in% da_pos_all$OTU) %>% pull(gene_id) %>% unique

taxa <- da_pos_all$OTU

enrichment_res <- data.frame()

for (t in taxa) {
    gene_filtered <- gene_str %>% filter(genus == t | species == t)
    if (nrow(gene_filtered) == 0) {
        cat("No genes found for taxon", t, "\n")
        next
    }
    # Genes in differentially abundant taxa
    K_da <- gene_filtered %>% filter(database == "KEGG") %>% pull(gene_id) %>% unique
    ek <- enrichKEGG(gene     = K_da,
                     universe = K_all,
                     organism = "ko",        # KEGG orthology space
                     pvalueCutoff = 0.1,
                     pAdjustMethod = "holm")
    
    # Save results
    ek_df <- as.data.frame(ek)
    if (nrow(ek_df) == 0) {
        cat("No enriched terms found for taxon", t, "\n")
        next
    }
    write.csv(ek_df, file = file.path(subdir, paste0("kegg_enrichment_", t, ".csv")), row.names = FALSE)
    
    # Add to main results
    enrichment_res <- rbind(enrichment_res,
                            ek_df %>% select(Description, category, subcategory, FoldEnrichment, p.adjust) %>% mutate(taxon = t))
    
    # Keep 20 most significant
    ek_df_f <- ek_df %>% select(-geneID) %>%
        arrange(p.adjust) %>%
        slice_head(n = 20) %>%
        arrange(desc(FoldEnrichment)) %>%
        mutate(Description = factor(Description, levels = rev(Description)))

    # Plot
    p <- ggplot(ek_df_f, aes(x = FoldEnrichment, y = Description, size = RichFactor, colour = FoldEnrichment)) +
        geom_point() +
        ggtitle(paste0("KEGG enrichment dotplot for ", t)) +
        scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1)
    
    ggsave(p, filename = file.path(subdir, paste0("enrichment_dotplot_", t, ".png")), width = 8, height = 6)
    
    p <- cnetplot(ek, showCategory = 5) + ggtitle(paste0("KEGG enrichment cnetplot for ", t)) +
        theme(plot.background = element_rect(fill = "white"))
    
    ggsave(p, filename = file.path(subdir, paste0("enrichment_cnetplot_", t, ".png")), width = 8, height = 6)
}

# Combine with taxa information
enrichment_res <- enrichment_res %>% left_join(da_pos_all, by = c("taxon" = "OTU")) %>%
        # Filter to positive diff. abund. associations and significantly enriched terms
        filter(association == "+") %>% filter(p.adjust < 0.05)

# Get term ordering based on fold enrichment
term_order <- enrichment_res %>% group_by(Description) %>%
    summarise(n_clusters = n_distinct(taxon),
              mean_FE = mean(FoldEnrichment)) %>%
    arrange(n_clusters, mean_FE)

enrichment_res$Description <- factor(enrichment_res$Description, levels = term_order$Description)
enrichment_res$term <- factor(enrichment_res$term, levels = c("Ruminant", "Fruit", "Animal", "Marine"))
# Get facet labels for plotting
enrichment_res$label <- ifelse(enrichment_res$category == "Metabolism", enrichment_res$subcategory, enrichment_res$category) %>%
                            str_replace_all(., " ", "\n")

p <- ggplot(enrichment_res, aes(x = taxon, y = Description, size = FoldEnrichment)) +
        geom_point(colour = "darkgreen") +
        facet_grid(cols = vars(paste0(term, association)),
                   rows = vars(label), scales = "free", space = "free") +
        theme(legend.position = "top", axis.text.x = element_text(hjust = 1, vjust = 0.5),
              axis.title = element_blank(), strip.text.y = element_text(angle = 0, size = 8))

ggsave(p, filename = file.path(subdir, "comparison_dotplot.png"), width = 10, height = 10)
