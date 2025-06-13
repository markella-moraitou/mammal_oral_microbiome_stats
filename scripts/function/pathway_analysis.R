##### PATHWAYS ANALYSIS #####

#### Visualisation of KEGG pathways

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(phyloseq)
library(microbiome)
library(tibble)
library(KEGGREST)
library(pathview)
library(stringr)
library(rlist)
library(purrr)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths
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
phy_gene <- readRDS(file.path(outdir, "phy_gene.RDS"))
phy_function <- readRDS(file.path(outdir, "phy_function.RDS"))

# MCMCglmm results
mcmc_results <- read.csv(file.path(outdir, "statistical_tests", "mcmcglmm_results.csv"))

# Big gene table
all_genes <- read.csv(file.path(outdir, "big_gene_table.csv"))

low_content <- read.csv(file.path(outdir, "low_content_samples.txt"), header = TRUE)

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

# DRAM and CAT output (contig annotations)
annot <- read.table(file.path(indir, "sample_annotations_counts.tsv"),
                        sep = "\t", header = TRUE, check.names=FALSE, quote = "", comment = "")

####################
#### PREP INPUT ####
####################

# Remove low content samples
phy_gene <- phy_gene %>% prune_samples(!sample_names(phy_gene) %in% low_content$x, .)

# Get only significant results that pass sensitivity test
mcmc_res <- mcmc_results %>%
            filter(!term %in% c("Species", "Total_abundance")) %>%
            # Consider as differentially abundant the taxa that have qval below 0.05 and passed the sensitivity test
            mutate(diff_abund = (pMCMC < 0.05)) %>%
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
mcmc_res <- data.frame(phy_gene@tax_table) %>% rownames_to_column("gene") %>% right_join(mcmc_res)

# Keep only significant results
mcmc_signif <- mcmc_res %>% filter(diff_abund)

# The big gene table contains genes that were not used for MCMCglmm and statistical tests
# But I should still show them (as non differentially abundant) in pathways, if they are there

# First remove low content samples
all_genes_filt <- all_genes[,!colnames(all_genes) %in% low_content$x]

# From the remaining, keep those in at least 10% of samples
all_genes_filt <- all_genes_filt[prevalence(all_genes_filt) > 0.1,]

#### Get unique pathways that contain differentially abundant genes

# Get pathways for each KO number
if (file.exists(file.path(subdir, "pathways_and_genes.yml"))) {
  cat("Loading pathways_and_genes.yml\n")
  pathways <- list.load(file.path(subdir, "pathways_and_genes.yml"))
} else {
  pathways <- list()
  for (ko in unique(mcmc_signif$gene)) {
    # Identify pathway ids where this KO is found
    pathids <- keggLink("pathway", ko) %>% str_remove("path:map") %>% str_remove("path:ko") %>% unique %>% as.character
    # Get pathway names names
    pathnames <- lapply(pathids, function(x) {keggFind(database = "pathway", query = x)}) %>% unlist
    # Path classes
    pathclasses <- lapply(pathids, 
                         function(x) { result <- keggGet(paste0("ko", x))[[1]]$CLASS
                                       if (is.null(result)) {
                                        return("no class")
                                       } else {
                                        return(result)
                                       }}) %>% unlist()
    df <- data.frame(pathid = pathids, pathname = pathnames, pathclass = pathclasses)
    rownames(df) <- NULL
    df <- df %>% column_to_rownames("pathid") %>% t %>% as.data.frame
    pathways[[ko]] <- df
  }
  list.save(pathways, file.path(subdir, "pathways_and_genes.yml"))
}

unique_pathways <- data.frame(pathid = character(), pathname = character(), pathclass = character())
# Get unique pathways, names, classes
for (i in 1:length(pathways)) {
  ko <- names(pathways)[i]
  df <- pathways[[ko]] %>% t %>% as.data.frame %>% rownames_to_column("pathid")
  unique_pathways <- unique_pathways %>% bind_rows(df)
}

unique_pathways <- unique_pathways %>%
  group_by(pathid, pathname, pathclass) %>%
  summarise(n_genes = n())

write.table(unique_pathways, file.path(subdir, "pathway_ids_and_names.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

#######################
#### PLOT PATHWAYS ####
#######################

# Remove top level pathways that cannot be plotted
# as well as human disease
unique_pathways <- unique_pathways %>% filter(pathclass != "no class" & !grepl("Human Diseases", pathclass))

for (i in 1:nrow(unique_pathways)) {
  pathname <- unique_pathways$pathname[i]
  pathid <- unique_pathways$pathid[i]
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

selected_pathids <- c("00130", # Ubiquinone biosynthesis
                      "00860", # Porphyrin metabolism
                      "02040", # Flagella assembly
                      #"00910", # Nitrogen metabolism
                      #"00920" # Sulfur metabolism
                      )

for (pathid in selected_pathids) {
  if(file.exists(file.path(subdir, paste0(pathid, "_gene_info.tsv")))) {
    gene_info <- read.table(file.path(subdir, paste0(pathid, "_gene_info.tsv")), sep = "\t", header = TRUE)
  } else {
    # Get gene id and names
    cat("Getting gene ids for", pathid, "...\n")
    gene_ids <- keggGet(paste0("ko", pathid))[[1]]$ORTHOLOGY %>% names
    gene_symbols <- lapply(gene_ids, function(x) {keggGet(x)[[1]]$SYMBOL}) %>% unlist
    gene_info <- data.frame(gene_id = gene_ids, gene_symbol = gene_symbols)
    print(gene_info)
    write.table(gene_info, file.path(subdir, paste0(pathid, "_gene_info.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  cat("Gene prevalences and sample metadata...\n")
  path_prevalence <- phy_gene %>% subset_taxa(taxa_names(phy_gene) %in% gene_info$gene_id) %>%
            psmelt %>% select(OTU, Abundance, Sample, Common.name, Species, Order, diet.general, gene_description) %>%
            rename("gene_id" = "OTU") %>%
            # Calculate prevalence per species
            group_by(gene_id, gene_description, Common.name, Species, Order, diet.general) %>%
            summarise(prevalence = sum(Abundance > 0)/n_distinct(Sample)) %>% ungroup %>%
            # Add gene symbols
            left_join(gene_info)
  # Reorder host species
  species_levels <- phy_gene@sam_data %>% data.frame %>% arrange(as.character(Order), as.character(digestion), Common.name) %>% select(Order, Common.name) %>% unique
  path_prevalence$Common.name <- factor(path_prevalence$Common.name, levels = species_levels$Common.name)
  cat("Plotting CLR prevalences per gene...\n")
  p <- ggplot(path_prevalence, aes(x = Common.name, y = prevalence, colour = Order, fill = diet.general)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = diet_palette) +
      scale_colour_manual(values = order_palette) +
      facet_wrap(~ gene_symbol, ncol = 3) +
      theme(axis.text.x = element_text(hjust = 1))
    
  ggsave(p, filename = file.path(subdir, paste0("ko", pathid, "_gene_prevalences.png")), width = 15, height = floor(length(unique(gene_ids))/3))
  
  # Find taxonomy of these genes
  tax_annot <- annot %>% right_join(gene_info) %>%
              # Clean up taxonomy
              mutate(across(superkingdom:species, ~ str_remove(str_remove(., ": .*$"), "^.__"))) %>%
              filter(!is.na(fasta) & phylum != "no support") %>%
              # Group rare phyla together
              mutate(phylum_grouped = factor(case_when(str_remove(phylum, "_.") %in% names(phylum_palette) ~ str_remove(phylum, "_."),
                                                        superkingdom == "Bacteria" ~ "Other Bacteria",
                                                        superkingdom == "Archaea" ~ "Other Archaea"), levels = names(phylum_palette)))
  
  # Add host species
  tax_annot$Common.name <- phy_gene@sam_data$Common.name[match(str_remove(tax_annot$fasta, "_final_contigs"),  phy_gene@sam_data$Ext.ID)]
  tax_annot$Common.name <- factor(tax_annot$Common.name, levels = species_levels$Common.name)
  tax_annot <- tax_annot[!is.na(tax_annot$Common.name),]
  # Get average count per species
  tax_annot <- tax_annot %>%
              # First sum genes per phylum per sample
              group_by(fasta, Common.name, phylum_grouped, gene_symbol) %>%
              summarise(count = sum(count)) %>% ungroup %>%
              group_by(Common.name, phylum_grouped, gene_symbol) %>%
              summarise(mean.count = mean(count))
  
  cat("Plotting microbial origin per gene...\n")
  p <- ggplot(tax_annot, aes(x = Common.name, y = mean.count, fill = phylum_grouped)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = phylum_palette, name = "Genus") +
      facet_wrap(~ gene_symbol, ncol = 3, scales = "free_y") +
      theme(axis.text.x = element_text(hjust = 1))

  ggsave(p, filename = file.path(subdir, paste0("ko", pathid, "_gene_taxa.png")), width = 15, height = floor(length(unique(gene_ids))/3))
}
