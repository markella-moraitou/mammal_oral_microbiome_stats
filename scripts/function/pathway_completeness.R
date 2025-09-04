##### PATHWAYS COMPLETENESS #####

#### Collect information on KEGG pathways and calculate their completeness

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
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "function")) # subdirectory for the output of this script
datadir <- normalizePath(file.path(outdir, "data")) # Directory with data files
subdir <- normalizePath(file.path(outdir, "pathway_completeness"))

here <- getwd() # Get current working directory

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

# Phyloseq object
phy_gene_f <- readRDS(file.path(outdir, "data", "phy_gene_f.RDS"))

# Gene abundance stratified
gene_str <- read.table(file.path(datadir, "gene_abundance_stratified_modified.tsv"),
                      quote = "", comment.char = "", header = TRUE, sep = "\t")

##################################
#### GET PATHWAY COMPLETENESS ####
##################################

ko_list <- phy_gene_f %>% subset_taxa(database == "KEGG") %>% taxa_names

#### Map KOs to pathways ####

# Check if file exists
if (file.exists(file.path(subdir, "ko_to_pathways.csv"))) {
  cat("Loading ko_to_pathways.csv\n")
  pathway_hits_df <- read.csv(file.path(subdir, "ko_to_pathways.csv"), stringsAsFactors = FALSE)
  colnames(pathway_hits_df) <- c("path", "ko")
} else {
  pathway_hits_df <- NULL
  write.table(data.frame(path = character(), ko = character()), sep = ",", file.path(subdir, "ko_to_pathways.csv"), row.names = FALSE, quote = FALSE, col.names = TRUE)
}

# Go through each KO and get pathways, skipping those already in the file
# and appending any new ones

lapply(ko_list, function(ko) {
  ko_in_file = FALSE
  # Check if the KO exists in the pathways list file and get info
  if (exists("pathway_hits_df")) {
    if (ko %in% str_remove(pathway_hits_df$ko, "ko:")) {
      cat("KO", which(ko_list == ko), "of", length(ko_list), "already in file...\n")
      ko_in_file = TRUE
    }
  }
  if (ko_in_file == FALSE) {
    # if not, query KEGG
    cat("Mapping KO", which(ko_list == ko), "of", length(ko_list), "to pathways...\n")
    path_hits <- tryCatch({
      Sys.sleep(1)
      keggLink("pathway", ko)
    }, error = function(e) {
    cat("Error for KO:", ko, "\n")})
    # Flatten and clean
    path_df <- data.frame(path = unlist(path_hits),
                          ko = names(path_hits),
                          stringsAsFactors = FALSE) %>%
              filter(grepl("ko", path))
    write.table(path_df, sep = ",", file.path(subdir, "ko_to_pathways.csv"), row.names = FALSE, quote = FALSE,
            append = TRUE, col.names = FALSE)
  }
})

pathway_hits_df <- read.csv(file.path(subdir, "ko_to_pathways.csv"), stringsAsFactors = FALSE) %>% unique

#### Get info on pathways ####
unique_pathways <- unique(pathway_hits_df$path)

path_info <- sapply(unique_pathways, function(pw) {
  cat("Getting info for pathway", pw, "...\n")
  tryCatch({
    Sys.sleep(0.4)
    keggGet(pw)
    }, error= function(e) pw)})

# Extract pathways class and name
path_info_df <- data.frame(path = unique_pathways)

for (i in 1:length(path_info)) {
    # Path name
    name <- paste(path_info[[i]]$NAME, collapse = "; ")
    path_info_df$path_name[i] <- name
    # Path class
    class <- paste(path_info[[i]]$CLASS, collapse = "; ")
    path_info_df$path_class[i] <- ifelse(is.null(class), NA, class)
}

# Add to main table
pathway_hits_df <- pathway_hits_df %>%
                    left_join(path_info_df, by = "path")

#### Get components of all identified pathways ####
pathway_kos <- lapply(unique_pathways, function(pw) {
  cat("Getting KOs for pathway", pw, which(unique_pathways == pw), "of", length(unique_pathways), "...\n")
  tryCatch({
    keggLink("ko", pw)
  }, error = function(e) NULL)
})

pathway_kos_df <- data.frame(pathway = names(unlist(pathway_kos)),
                             kos = unlist(pathway_kos))

write.csv(pathway_kos_df, file.path(subdir, "pathways_to_kos.csv"), row.names = FALSE, quote = FALSE)

#### Calculate pathway completeness in the entire dataset ####
completeness <- sapply(pathway_kos, function(kos) {
  if (is.null(kos)) return(NA)
  total_kos <- length(kos)
  matched_kos <- sum(str_remove(kos, "ko:") %in% ko_list)
  matched_kos / total_kos
})

compl_df <- data.frame(path = unique_pathways,
                      completeness_in_dataset = completeness)

# Add to pathway info
path_info_df <- path_info_df %>%
                  left_join(compl_df, by = "path")

write.csv(path_info_df, file.path(subdir, "pathway_info.csv"), row.names = FALSE, quote = TRUE)

# Keep only relevant pathways
path_info_filt <- path_info_df %>% filter(!grepl("Human Diseases;", path_class)) %>% filter(!grepl("Organismal Systems;", path_class))

# Same in pathway_kos_df
pathway_kos_filt <- pathway_kos_df %>% filter(pathway %in% path_info_filt$path)

#### Calculate pathways completeness per sample ####

ko_list_samples <- list()

# Create a list with all KOs present in each sample
for (sample in sample_names(phy_gene_f)) {
  cat("Extracting KOs present in", sample, "...\n")
  phy_sub <- prune_samples(sample, phy_gene_f) %>% prune_taxa(taxa_sums(.) > 0, .)
  # Get KOs present in this sample
  ko_present <- taxa_names(phy_sub)[taxa_names(phy_sub) %in% ko_list]
  ko_list_samples[[sample]] <- ko_present
}

unique_pathways <- unique(pathway_kos_filt$pathway)

pathway_compl_sample <- matrix(nrow = length(unique_pathways), ncol = length(sample_names(phy_gene_f)))
rownames(pathway_compl_sample) <- unique_pathways
colnames(pathway_compl_sample) <- sample_names(phy_gene_f)

for (i in 1:length(unique_pathways)) {
  pw <- unique_pathways[i]
  kos <- pathway_kos[[i]] %>% str_remove("ko:")
  if (is.null(kos)) {
    cat("No KOs found for pathway", pw, "\n")
    pathway_compl_sample[pw, ] <- NA
    next
  }
  for (j in 1:length(sample_names(phy_gene_f))) {
    sample <- sample_names(phy_gene_f)[j]
    kos_in_sample <- ko_list_samples[[sample]]
    matched_kos <- sum(kos %in% kos_in_sample)
    pathway_compl_sample[pw, sample] <- matched_kos / length(kos)
  }
}

write.csv(pathway_compl_sample, file.path(subdir, "pathway_completeness_per_sample.csv"), row.names = TRUE, quote = FALSE)

## Plot pathways completeness
pathway_compl_l <- as.data.frame(pathway_compl_sample) %>%
                    rownames_to_column("path") %>%
                    pivot_longer(-path, names_to = "Sample", values_to = "completeness") %>%
                    left_join(path_info_filt %>% select(path, path_name, path_class), by = "path") %>%
                    # Keep only average completeness above 0.2
                    group_by(path) %>%
                    filter(mean(completeness) > 0.1)

pathway_compl_l$Common.name <- phy_gene_f@sam_data$Common.name[match(pathway_compl_l$Sample, sample_names(phy_gene_f))]
pathway_compl_l$class <- str_remove(pathway_compl_l$path_class, ".*; ")

p <- ggplot(aes(x = Sample, y = path, fill = completeness), data = pathway_compl_l) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma", na.value = "grey90") +
      labs(x = "Sample", y = "KEGG Pathway", fill = "Completeness") +
      facet_grid(cols = vars(Common.name), rows = vars(class), scales = "free", space = "free") +
      theme(axis.text.x = element_blank(),
            strip.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0))

ggsave(p, filename = file.path(subdir, "pathway_completeness_per_sample_heatmap.png"), width = 20, height = 20)

#### Calculate pathways completeness per microbial genus per sample ####
ko_per_taxon_df <- gene_str %>% filter(database == "KEGG") %>%
                    # At least an average coverage of 1
                    filter(totalAvgDepth >= 1) %>% filter(genus != "no support") %>%
                    select(sample, genus, gene_id) %>%
                    # Keep genera with at least 1000 genes
                    group_by(genus, sample) %>%
                    filter(n_distinct(gene_id) >= 1000)

# Get a list of KOs per genus per sample
ko_list_s_g <- list()

for (s_g in unique(paste(ko_per_taxon_df$sample, ko_per_taxon_df$genus, sep = "&"))) {
  cat("Extracting KOs present in", s_g, "...\n")
  ko_present <- ko_per_taxon_df %>% filter(sample == str_split(s_g, "&")[[1]][1] & genus == str_split(s_g, "&")[[1]][2]) %>%
                pull(gene_id) %>% unique
  ko_list_s_g[[s_g]] <- ko_present
}

pathway_compl_s_g <- matrix(nrow = length(unique_pathways), ncol = length(names(ko_list_s_g)))
rownames(pathway_compl_s_g) <- unique_pathways
colnames(pathway_compl_s_g) <- names(ko_list_s_g)

for (i in 1:length(unique_pathways)) {
  pw <- unique_pathways[i]
  kos <- pathway_kos[[i]] %>% str_remove("ko:")
  if (is.null(kos)) {
    cat("No KOs found for pathway", pw, "\n")
    pathway_compl_s_g[pw, ] <- NA
    next
  }
  for (j in 1:length(names(ko_list_s_g))) {
    s_g <- names(ko_list_s_g)[j]
    kos_in_s_g <- ko_list_s_g[[s_g]]
    matched_kos <- sum(kos %in% kos_in_s_g)
    pathway_compl_s_g[pw, s_g] <- matched_kos / length(kos)
  }
}

pathway_compl_s_g <- pathway_compl_s_g %>% t %>% data.frame %>% rownames_to_column("sample_genus") %>%
                      separate(sample_genus, into = c("Sample", "genus"), sep = "&", remove = TRUE)

write.csv(pathway_compl_s_g, file.path(subdir, "pathway_completeness_per_sample_per_micrgenus.csv"), row.names = FALSE, quote = FALSE)

# How different are the same taxa between species and across species

# Calculate pairwise euclidean distances
#dist_matrix <- dist(pathway_compl_s_g %>% select(-Sample, -genus), method = "manhattan") %>% as.matrix
cor_matrix <- cor(pathway_compl_s_g %>% select(-Sample, -genus) %>% t, method = "pearson") %>% as.matrix
s_g <- paste(pathway_compl_s_g$Sample, pathway_compl_s_g$genus, sep = "&")

rownames(cor_matrix) <- s_g
colnames(cor_matrix) <- s_g

cor_df <- as.data.frame(cor_matrix) %>% rownames_to_column("term1") %>% pivot_longer(-term1, names_to = "term2", values_to = "correlation") %>%
      separate(term1, into = c("Sample1", "genus1"), sep = "&", remove = TRUE) %>%
      separate(term2, into = c("Sample2", "genus2"), sep = "&", remove = TRUE) %>%
      left_join(data.frame(phy_gene_f@sam_data) %>% select(Ext.ID, Species) %>% rename(Sample1 = Ext.ID, Species1 = Species), by = "Sample1") %>%
      left_join(data.frame(phy_gene_f@sam_data) %>% select(Ext.ID, Species) %>% rename(Sample2 = Ext.ID, Species2 = Species), by = "Sample2")

cor_df <- cor_df %>% filter(!is.na(Species1) & !is.na(Species2)) %>%
      filter(Sample1 != Sample2) %>%
      mutate(relationship = case_when(
        Species1 == Species2 & genus1 == genus2 ~ "same host species & genus",
        Species1 == Species2 & genus1 != genus2 ~ "same host species, different genus",
        Species1 != Species2 & genus1 == genus2 ~ "different host species, same genus",
        Species1 != Species2 & genus1 != genus2 ~ "different host species & genus",
        TRUE ~ "other"
      )) %>%
      mutate(relationship = factor(relationship, levels = c("same host species & genus",
                                                            "different host species, same genus",
                                                            "same host species, different genus",
                                                            "different host species & genus")))

# Plot
p1 <- ggplot(aes(y = relationship, x = correlation), data = cor_df) +
    geom_boxplot() + 
    geom_vline(xintercept = 0.90, linetype = "dashed", colour = "red") +
    xlim(c(min(cor_df$correlation), max(cor_df$correlation)))
    
# Plot
p2 <- ggplot(aes(y = genus1, x = correlation), data = filter(cor_df, genus1 == genus2)) +
    geom_boxplot() + 
    geom_vline(xintercept = 0.90, linetype = "dashed", colour = "red") +
    xlim(c(min(cor_df$correlation), max(cor_df$correlation)))

p <- plot_grid(p1, p2, ncol = 2, align = "h")

ggsave(p, filename = file.path(subdir, "pathway_completeness_dist.png"), width = 10, height = 10)

#### Calculate pathways completeness per microbial genus across dataset ####
## Genera seem to have more or less stable pathway profiles across host species,
## so we can consider them across the dataset
ko_per_taxon_df <- ko_per_taxon_df %>% ungroup %>% select(-sample) %>% unique

# Get a list of KOs per genus per sample
ko_list_genus <- list()

for (g in unique(ko_per_taxon_df$genus)) {
  cat("Extracting KOs present in", g, "...\n")
  ko_present <- ko_per_taxon_df %>% filter(genus == g) %>%
                pull(gene_id) %>% unique
  ko_list_genus[[g]] <- ko_present
}

pathway_compl_genus <- matrix(nrow = length(unique_pathways), ncol = length(names(ko_list_genus)))
rownames(pathway_compl_genus) <- unique_pathways
colnames(pathway_compl_genus) <- names(ko_list_genus)

for (i in 1:length(unique_pathways)) {
  pw <- unique_pathways[i]
  kos <- pathway_kos[[i]] %>% str_remove("ko:")
  if (is.null(kos)) {
    cat("No KOs found for pathway", pw, "\n")
    pathway_compl_genus[pw, ] <- NA
    next
  }
  for (j in 1:length(names(ko_list_genus))) {
    g <- names(ko_list_genus)[j]
    kos_in_g <- ko_list_genus[[g]]
    matched_kos <- sum(kos %in% kos_in_g)
    pathway_compl_genus[pw, g] <- matched_kos / length(kos)
  }
}

write.csv(pathway_compl_genus, file.path(subdir, "pathway_completeness_per_micrgenus.csv"), row.names = TRUE, quote = FALSE)

## Plot pathways completeness
pathway_compl_l <- as.data.frame(pathway_compl_genus) %>%
                    rownames_to_column("path") %>%
                    pivot_longer(-path, names_to = "genus", values_to = "completeness") %>%
                    left_join(path_info_filt %>% select(path, path_name, path_class), by = "path") %>%
                    # Keep only average completeness above 0.1
                    group_by(path) %>%
                    filter(mean(completeness) > 0.1)

tax_info <- gene_str %>% select(genus, family, order, class, phylum) %>% unique

pathway_compl_l$phylum <- tax_info$phylum[match(pathway_compl_l$genus, tax_info$genus)]
pathway_compl_l$order <- tax_info$order[match(pathway_compl_l$genus, tax_info$genus)]

pathway_compl_l$path_class <- str_remove(pathway_compl_l$path_class, ".*; ")

p <- ggplot(aes(x = genus, y = path, fill = completeness), data = pathway_compl_l) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma", na.value = "grey90") +
      labs(x = "Microbial genus", y = "KEGG Pathway", fill = "Completeness") +
      facet_grid(cols = vars(order), rows = vars(path_class), scales = "free", space = "free") +
      theme(axis.text.x = element_blank(),
            strip.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0))

ggsave(p, filename = file.path(subdir, "pathway_completeness_per_taxon_heatmap.png"), width = 20, height = 20)

#######################
#### PLOT PATHWAYS ####
#######################

#### Plot all relevant pathways ####
setwd(subdir)
dir.create("pathway_plots", showWarnings = FALSE)
setwd("pathway_plots")

path_info_filt <- path_info_filt %>% filter(completeness_in_dataset >= 0.2 & path_class != "")

for (i in 1:nrow(path_info_filt)) {
  pathname <- path_info_filt$path_name[i]
  pathid <- path_info_filt$path[i]
  cat("GETTING INFO ON ", pathname, pathid, "\n")

  # Plot pathway
  tryCatch({
    path <- pathview(
      gene.data = ko_list,
      pathway.id = str_remove(pathid, "path:"),
      species = "ko",
      out.suffix = make.names(pathname))
    }, error = function(e) {
    cat("!! Couldn't plot", pathname, "\n")})
}
