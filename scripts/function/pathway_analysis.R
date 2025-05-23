##### PATHWAYS ANALYSIS #####

#### Visualisation of KEGG pathways

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(phyloseq)
library(tibble)
library(KEGGREST)
library(pathview)
library(stringr)
library(viridis)
library(rlist)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "function")) # subdirectory for the output of this script
subdir <- normalizePath(file.path(outdir, "pathway_plots"))

dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

# Get ordination functions
source(file.path("..","ordination_functions.R"))

setwd(subdir)

#######################
#####  LOAD INPUT #####
#######################

# Phyloseq object
phy_gene_clr <- readRDS(file.path(outdir, "phy_gene_clr.RDS"))
phy_function <- readRDS(file.path(outdir, "phy_function.RDS"))

# ANCOMBC results
mcmc_results <- read.csv(file.path(outdir, "statistical_tests", "mcmcglmm_results.csv"))

low_content <- read.csv(file.path(outdir, "low_content_samples.txt"), header = TRUE)

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

# Remove low content samples
phy_gene_clr <- phy_gene_clr %>% prune_samples(!sample_names(phy_gene_clr) %in% low_content$x, .)

# Contig annotations
contig_annot <- read.table("/cfs/klemming/projects/snic/sllstore2017021/MAMMALIAN_DC/output/T1_contig_taxonomy/combined_contig_annotations.tsv",
                            header = TRUE, sep = "\t", quote = "", comment = "")

####################
#### PREP INPUT ####
####################

# Get only significant results that pass sensitivity test
mcmc_res <- mcmc_results %>%
            filter(!term %in% c("Species", "Total_abundance")) %>%
            # Consider as differentially abundant the taxa that have qval below 0.05 and passed the sensitivity test
            mutate(diff_abund = (qval < 0.05)) %>%
            # If not diff abund, consider lfc to be 0 (no change)
            mutate(post.mean = case_when(diff_abund ~ post.mean,
                                   TRUE ~ 0)) %>%
            # label association (positive and negative)
            mutate(term = str_remove(term, "Order|habitat_general|ruminant")) %>%
            mutate(assoc = case_when(post.mean < 0 ~ paste0(term, "- "),
                                     post.mean > 0 ~ paste0(term, "+ "),
                                     post.mean == 0 ~ "")) %>%
            # Then summarise all associations per taxon
            group_by(gene) %>%
            mutate(label = paste(assoc, collapse=""), assoc = NULL, qval = NULL, sensitivity = NULL, pval = NULL)

# Add more info on genes
mcmc_res <- data.frame(phy_gene_clr@tax_table) %>% rownames_to_column("gene") %>% right_join(mcmc_res)

# Keep only significant results
mcmc_signif <- mcmc_res %>% filter(diff_abund)

#### Get unique pathways that contain differentially abundant genes

# Get pathways for each KO number
if (file.exists(file.path(subdir, "pathways_and_genes.yml"))) {
  cat("Loading pathways_and_genes.yml")
  pathways <- list.load(file.path(subdir, "pathways_and_genes.yml"))
} else {
  pathways <- list()
  for (ko in unique(mcmc_signif$gene)) {
    # Identify pathway ids where this KO is found
    pathids <- keggLink("pathway", ko) %>% str_remove("path:map") %>% str_remove("path:ko") %>% unique
    # Get pathway names names
    pathnames <- lapply(unique_pathids, function(x) {keggFind(database = "pathway", query = x)}) %>% unlist
    # Path classes
    pathclasses <- lapply(unique_pathids, function(x) {keggGet(paste0("ko", x))[[1]]$CLASS}) %>% unlist
    pathways[[ko]] <- paths
  }
  list.save(pathways, file.path(subdir, "pathways_and_genes.yml"))
}

unique_pathids <- unlist(pathways) %>% unique

unique_pathways <- lapply(unique_pathids, function(x) {keggFind(database = "pathway", query = x)}) %>% unlist
names(unique_pathways) <- unique_pathids

write.table(unique_pathways, file.path(subdir, "pathway_ids_and_names.txt"), quote = FALSE, col.names = FALSE)

#######################
#### PLOT PATHWAYS ####
#######################

for (i in 1:length(unique_pathways)) {
  pathname <- unique_pathways[i]
  pathid <- names(pathname)
  
  cat("GETTING INFO ON ", pathname, pathid, "\n")
  # Get pathways genes
  genes <- keggGet(paste0("ko", pathid))[[1]]$ORTHOLOGY %>% names
  if (is.null(genes)) {
    cat("No genes slot found for", pathname, ". Skipping...\n")
    next
  }
  # Keep only lfc for terms that have at least on differentially abundant gene in this pathway
  data_to_plot <- mcmc_res %>% filter(gene %in% genes) %>%
      # Identify if a term has differentially abundant genes associated with it
      group_by(term) %>% mutate(has_diff_abund = (sum(post.mean != 0) > 0)) %>%
      filter(has_diff_abund) %>%
      # Make wider and use gene as rownames
      select(gene, term, post.mean) %>%
      pivot_wider(names_from = term, values_from = post.mean) %>%
      column_to_rownames("gene")
  
  # Plot a separate plot for every host group that shows an association
  for (j in 1:ncol(data_to_plot)) {
    term = colnames(data_to_plot)[j]
    # Plot pathway
    path <- pathview(
      gene.data = select(data_to_plot, j),
      pathway.id = pathid,
      species = "ko",
      out.suffix = paste(make.names(pathname), make.names(term), sep = "_"),
      # We are plotting number of hits, so specify TRUE for this
      # If plotting, say, gene/transcript abundance, set this to FALSE
      descrete = list(
        gene = TRUE,
        cpd = TRUE),
      # Tally colours
      low = "red",
      mid = "yellow",
      high = "green",
    )
  }
}

##############################
#### PLOT GENE ABUNDANCES ####
##############################

pathid="00130"
genes <- keggGet(paste0("ko", pathid))[[1]]$ORTHOLOGY %>% names

path_abundances <- phy_gene_clr %>% subset_taxa(taxa_names(phy_gene_clr) %in% genes) %>%
          psmelt %>% select(OTU, Abundance, Sample, Common.name, Species, Order, diet.general, gene_description) %>%
          rename("gene" = "OTU")

p <- ggplot(path_abundances, aes(x = Common.name, y = Abundance, colour = Order, fill = diet.general)) +
    geom_boxplot() +
    scale_fill_manual(values = diet_palette) +
    scale_colour_manual(values = order_palette) +
    facet_wrap(~ gene_description, ncol = 3)

ggsave(p, filename = file.path(subdir, paste0("ko", pathid, "_gene_abundances.png")), width = 15, height = 15)

# Find taxonomy of these genes
filter(contig_annot, grepl("K00355", paste(.)))