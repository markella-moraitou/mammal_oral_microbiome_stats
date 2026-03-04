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
subdir <- normalizePath(file.path(outdir, "gsea_diff_abund_taxa")) # subdirectory for the output of this script

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
all_da <- read.csv(file.path(taxdir, "differential_abundance", "all_diffabund_results.csv"))

# Stratified sample data
gene_str <- read.table(file.path(datadir, "gene_abundance_stratified_modified.tsv"),
                      quote = "", comment.char = "", header = TRUE, sep = "\t")

####################
#### FOCAL TAXA ####
####################

# Identify differentially abundant taxa and their association

# From both outputs, extract significant taxa, and what their association is (positive or negative)
# From both outputs, extract significant taxa, and what their association is (positive or negative)
da_taxa <- all_da %>%
  filter(padj_t < 0.05 & pval_b < 0.05) %>%
  group_by(OTU, term, superkingdom, family) %>%
  summarise(association = ifelse(all(across(starts_with("coeff"), ~ . > 0)), "+", "-"))

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

K_all <- gene_str %>% filter(database == "KEGG") %>% pull(gene_id) %>% unique

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

# Keep only metabolism-associated terms for plotting
enrichment_res_filt <- enrichment_res %>%
    filter(!is.na(term) & category %in% c("Metabolism", "Cellular Processes"))

enrichment_res_filt$subcategory <- factor(enrichment_res_filt$subcategory, levels = unique(c("Global and overview maps", enrichment_res_filt$subcategory)))

p <- ggplot(enrichment_res_filt, aes(x = taxon, y = Description, size = FoldEnrichment)) +
        geom_point(colour = "darkgreen") +
        facet_grid(cols = vars(paste0(term, association)),
                   rows = vars(subcategory), scales = "free", space = "free") +
        theme(legend.position = "top", axis.text.x = element_text(hjust = 0, vjust = 0.5, angle = -45),
              axis.title = element_blank(), strip.text.y = element_text(angle = 0, size = 8), strip.text.x = element_text(angle = 90))

ggsave(p, filename = file.path(subdir, "comparison_dotplot.png"), width = 12, height = 10)
