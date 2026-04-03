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
library(readr)
library(stringr)
library(phyloseq)
library(phyr)
library(phytools)
library(ggplot2)
library(clusterProfiler)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "mags")) # subdirectory for the output of this script
taxdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Directory with taxonomy analysis
funcdir <- normalizePath(file.path("..", "..", "output", "function")) # Directory with taxonomy analysis
pathdir <- normalizePath(file.path(outdir, "pathway_completeness")) # Directory with pathway analysis output
subdir <- normalizePath(file.path(outdir, "gsea_diff_abund_mags")) # subdirectory for the output of this script

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
phy_sp_f <- readRDS(file.path(taxdir, "phyloseq_objects", "phy_sp_f.RDS"))

# Differentially abundant taxa
all_da <- read.csv(file.path(taxdir, "differential_abundance", "all_diffabund_results.csv"))

# MAG annotations
genes <- read_tsv(file.path(indir, "MAG_annotations.tsv.gz"), quote = "", comment = "")
colnames(genes)[1] <- "contig"

# All genes in data
gene_str <- read.table(file.path(funcdir, "data", "gene_abundance_stratified_modified.tsv"),
                      quote = "", comment.char = "", header = TRUE, sep = "\t")

# MAG metadata
bac_meta <- read.table(file.path(outdir, "bac_meta.tsv"), sep="\t", header=TRUE)
ar_meta <- read.table(file.path(outdir, "ar_meta.tsv"), sep="\t", header=TRUE)

####################
#### FOCAL TAXA ####
####################

# Identify differentially abundant taxa and their association

# Extract significant taxa, and what their association is (positive or negative)
# For ruminants, consider taxa only identified using PGLMM as the two methods didn't converge on many taxa
da_taxa <- all_da %>%
  filter((padj_t < 0.05 & pval_b < 0.05) | (padj_t < 0.05 & term == "Ruminant")) %>%
  mutate(OTU = gsub("\\.", "*", OTU)) %>% # Fix names
  group_by(OTU, term) %>%
  summarise(association = ifelse(all(across(starts_with("coeff"), ~ . > 0)), "+", "-"))

tax_info <- data.frame(phy_sp_f@tax_table) %>% select(phylum, order, family, genus) %>% unique %>%
    filter(genus %in% da_taxa$OTU)

da_taxa <- da_taxa %>% left_join(tax_info, by = c("OTU" = "genus"))

# Keep only positive associations at least for now
da_pos <- da_taxa %>% filter(association == "+") %>% unique

##############################
#### ENRICHMENT ANALYSIS #####
##############################

# Identify HQ MAGs corresponding to the differentially abundant taxa
hq_meta <- rbind(bac_meta, ar_meta) %>%
    filter((Completeness >= 90 & Contamination <= 5) | (phylum == "Patescibacteria" & Contamination <= 5)) %>%
    select(label, bin, Sample, host_species, species, genus, family, order, phylum) %>%
    mutate(matching_tax = str_remove(label, "[sgfo]__") %>% str_remove(., "_.*"))

# Match diff abund OTUs to MAGs
da_pos_mags <- da_pos %>%
        mutate(matching_tax = case_when(OTU %in% hq_meta$matching_tax ~ OTU,
                                        family %in% hq_meta$matching_tax ~ family,
                                        order %in% hq_meta$matching_tax ~ order,
                                        TRUE ~ NA)) %>%
        filter(!is.na(matching_tax))

da_pos_mags <- hq_meta  %>% select(-c(family, order, phylum)) %>%
    right_join(da_pos_mags, by = c("matching_tax"))

write.csv(da_pos_mags, file = file.path(subdir, "mags_matched_diff_abund_taxa.csv"), row.names = FALSE, quote = FALSE)

##############################
#### ENRICHMENT ANALYSIS #####
##############################

# As background, use all KEGG genes in the dataset

K_all <- gene_str %>% filter(database == "KEGG") %>% pull(gene_id)

genomes <- da_pos_mags$bin

enrichment_res <- data.frame()

for (g in genomes) {
    label <- hq_meta$label[hq_meta$bin == g]
    cat("Processing genome", label, "\n")
    gene_filtered <- genes %>% filter(fasta == g)
    if (nrow(gene_filtered) == 0) {
        cat("No genes found for genome", label, "\n")
        next
    }
    # Genes in differentially abundant taxa
    K_da <- gene_filtered %>% filter(!is.na(ko_id)) %>% pull(ko_id) %>% unique
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
    write.csv(ek_df, file = file.path(subdir, paste0("kegg_enrichment_", label, ".csv")), row.names = FALSE)
    
    # Add to main results
    enrichment_res <- rbind(enrichment_res,
                            ek_df %>% select(Description, category, subcategory, FoldEnrichment, p.adjust) %>% mutate(bin = g, label = label))
    
    # Keep 20 most significant
    ek_df_f <- ek_df %>% select(-geneID) %>%
        arrange(p.adjust) %>%
        slice_head(n = 20) %>%
        arrange(desc(FoldEnrichment)) %>%
        mutate(Description = factor(Description, levels = rev(Description)))

    # Plot
    p <- ggplot(ek_df_f, aes(x = FoldEnrichment, y = Description, size = RichFactor, colour = FoldEnrichment)) +
        geom_point() +
        ggtitle(paste0("KEGG enrichment dotplot for ", label)) +
        scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1) +
        theme(axis.text.y = element_text(size = 8))
    
    ggsave(p, filename = file.path(subdir, paste0("enrichment_dotplot_", label, ".png")), width = 8, height = 6)
    
    p <- cnetplot(ek, showCategory = 5) + ggtitle(paste0("KEGG enrichment cnetplot for ", label)) +
        theme(plot.background = element_rect(fill = "white"))
    
    ggsave(p, filename = file.path(subdir, paste0("enrichment_cnetplot_", label, ".png")), width = 8, height = 6)
}

# Combine with taxa information
enrichment_res <- enrichment_res %>% left_join(da_pos_mags, by = c("bin" = "bin", "label" = "label")) %>%
        # Filter to positive diff. abund. associations and significantly enriched terms
        filter(association == "+") %>% filter(p.adjust < 0.05)

# Get term ordering based on number of occurences and fold enrichment
term_order <- enrichment_res %>% group_by(Description) %>%
    summarise(n_clusters = n_distinct(label),
              mean_FE = mean(FoldEnrichment)) %>%
    arrange(n_clusters, mean_FE)

enrichment_res$Description <- factor(enrichment_res$Description, levels = term_order$Description)
enrichment_res$term <- factor(enrichment_res$term, levels = c("Ruminant", "Fruit", "Animal", "Marine"))

# Keep only metabolism-associated terms for plotting
enrichment_res_filt <- enrichment_res %>%
    filter(!is.na(term) & category %in% c("Metabolism", "Cellular Processes"))

enrichment_res_filt$subcategory <- factor(enrichment_res_filt$subcategory, levels = unique(c("Global and overview maps", enrichment_res_filt$subcategory)))

p <- ggplot(enrichment_res_filt, aes(x = label, y = Description, size = FoldEnrichment)) +
        geom_point(colour = "darkgreen") +
        facet_grid(cols = vars(paste0(term, association)),
                   rows = vars(subcategory), scales = "free", space = "free") +
        theme(legend.position = "top",
              axis.text.x = element_text(size = 10, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 9, hjust = 1, vjust = 0.5),
              axis.title = element_blank(),
              strip.text.y = element_text(angle = 0, size = 8), strip.text.x = element_text(angle = 90))

ggsave(p, filename = file.path(subdir, "comparison_dotplot.png"), width = 12, height = 10)
