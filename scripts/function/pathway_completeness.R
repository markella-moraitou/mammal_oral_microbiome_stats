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
phy_gene_f_clr <- readRDS(file.path(outdir, "data", "phy_gene_f_clr.RDS"))

# Gene abundance stratified
gene_str <- read.table(file.path(datadir, "gene_abundance_stratified_modified.tsv"),
                      quote = "", comment.char = "", header = TRUE, sep = "\t")

##################################
#### GET PATHWAY COMPLETENESS ####
##################################

ko_list <- phy_gene_f %>% subset_taxa(database == "KEGG") %>% taxa_names

#### Map KOs to pathways ####

add_paths <- FALSE # Change this to append more pathways
# Check if file exists
if (file.exists(file.path(subdir, "ko_to_pathways.csv"))) {
  cat("Loading ko_to_pathways.csv\n")
  pathway_hits_df <- read.csv(file.path(subdir, "ko_to_pathways.csv"), stringsAsFactors = FALSE)
  colnames(pathway_hits_df) <- c("path", "ko")
} 
if (add_paths == TRUE | !file.exists(file.path(subdir, "ko_to_pathways.csv"))) {
  if(!file.exists(file.path(subdir, "ko_to_pathways.csv"))) {
      pathway_hits_df <- data.frame(path = character(),
                                    ko = character(),
                                    stringsAsFactors = FALSE)
    } else {
      pathway_hits_df <- read.csv(file.path(subdir, "ko_to_pathways.csv"), stringsAsFactors = FALSE)
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
}

#### Get components of all identified pathways ####
unique_pathways <- unique(pathway_hits_df$path)

if (file.exists(file.path(subdir, "pathways_to_kos.csv"))) {
  cat("Loading pathways_to_kos.csv\n")
  pathway_kos_df <- read.csv(file.path(subdir, "pathways_to_kos.csv"), stringsAsFactors = FALSE)
  pathway_kos <- split(pathway_kos_df$kos, pathway_kos_df$pathway)
} else {
  pathway_kos <- lapply(unique_pathways, function(pw) {
    cat("Getting KOs for pathway", pw, which(unique_pathways == pw), "of", length(unique_pathways), "...\n")
    tryCatch({
      keggLink("ko", pw)
    }, error = function(e) NULL)
  })
  pathway_kos_df <- data.frame(pathway = names(unlist(pathway_kos)),
                               kos = unlist(pathway_kos))
  write.csv(pathway_kos_df, file.path(subdir, "pathways_to_kos.csv"), row.names = FALSE, quote = FALSE)
}

#### Get info on pathways ####

if (file.exists(file.path(subdir, "pathway_info.csv"))) {
  cat("Loading pathway_info.csv\n")
  path_info_df <- read.csv(file.path(subdir, "pathway_info.csv"), stringsAsFactors = FALSE)
} else {
  path_info <- sapply(unique_pathways, function(pw) {
    cat("Getting info for pathway", pw, which(unique_pathways == pw), "of", length(unique_pathways), "...\n")
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
  path_info_df <- path_info_df %>%
                    mutate(path_class = ifelse(path_class == "", "Overview Pathway", path_class))
  # Add to main table
  pathway_hits_df <- pathway_hits_df %>%
                      left_join(path_info_df, by = "path")
  #### Calculate pathway completeness and abundance in the entire dataset ####
  completeness <- sapply(pathway_kos, function(kos) {
    if (is.null(kos)) return(NA)
    total_kos <- length(kos)
    matched_kos <- sum(str_remove(kos, "ko:") %in% ko_list)
    matched_kos / total_kos
  })
  abundance <- sapply(pathway_kos, function(kos) {
    if (is.null(kos)) return(NA)
    kos_clean <- str_remove(kos, "ko:")
    ko_abundances <- data.frame(phy_gene_f@otu_table) %>% filter(rownames(.) %in% kos_clean) %>%
                        rowSums
    mean(ko_abundances)
  })
  compl_df <- data.frame(path = unique_pathways,
                         completeness_in_dataset = completeness,
                          mean_reads = mean_reads)
  # Add to pathway info
  path_info_df <- path_info_df %>%
                    left_join(compl_df, by = "path")
  write.csv(path_info_df, file.path(subdir, "pathway_info.csv"), row.names = FALSE, quote = TRUE)
}

# Keep only relevant pathways
path_info_filt <- path_info_df %>% filter(!grepl("Human Diseases;", path_class)) %>%
    filter(!grepl("Organismal Systems;", path_class)) %>%
    filter(!grepl("viruses", path_class)) %>%
    filter(completeness_in_dataset > 0)

# Same in pathway_kos_df
pathway_kos_filt <- pathway_kos_df %>% filter(pathway %in% path_info_filt$path)

## Plot
p <- ggplot(aes(y = completeness_in_dataset, x = abundance), data = path_info_filt) +
      geom_point() +
      labs(y = "Pathway completeness in dataset", x = "Mean pathway abundance\n(mean mapped reads)")

ggsave(p, filename = file.path(subdir, "pathway_completeness_vs_abundance.png"), width = 5, height = 5)

#### Calculate pathways completeness and abundance per sample ####

ko_list_samples <- list()

# Create a list with all KOs present in each sample
for (sample in sample_names(phy_gene_f)) {
  cat("Extracting KOs present in", sample, "...\n")
  phy_sub <- prune_samples(sample, phy_gene_f) %>% prune_taxa(taxa_sums(.) > 0, .)
  # Get KOs present in this sample
  ko_present <- taxa_names(phy_sub)[taxa_names(phy_sub) %in% ko_list]
  ko_list_samples[[sample]] <- ko_present
}

# Collect info only on the filtered pathways
filt_pathways <- path_info_filt$path

if (file.exists(file.path(subdir, "pathway_completeness_per_sample.csv")) &
    file.exists(file.path(subdir, "pathway_abundance_per_sample.csv")) &
    file.exists(file.path(subdir, "pathway_abundance_per_sample_clr.csv"))) {
  cat("Loading pathway completeness and abundance per sample...\n")
  pathway_compl_sample <- read.csv(file.path(subdir, "pathway_completeness_per_sample.csv"), row.names = 1, check.names = FALSE)
  pathway_abund_sample <- read.csv(file.path(subdir, "pathway_abundance_per_sample.csv"), row.names = 1, check.names = FALSE)
  pathway_abund_sample_clr <- read.csv(file.path(subdir, "pathway_abundance_per_sample_clr.csv"), row.names = 1, check.names = FALSE)
} else {
  pathway_compl_sample <- matrix(nrow = length(filt_pathways), ncol = length(sample_names(phy_gene_f)))
  rownames(pathway_compl_sample) <- filt_pathways
  colnames(pathway_compl_sample) <- sample_names(phy_gene_f)

  pathway_abund_sample <- matrix(nrow = length(filt_pathways), ncol = length(sample_names(phy_gene_f)))
  rownames(pathway_abund_sample) <- filt_pathways
  colnames(pathway_abund_sample) <- sample_names(phy_gene_f)

  for (i in 1:length(filt_pathways)) {
    cat("Calculating completeness and abundance for pathway", i, "of", length(filt_pathways), "...\n")
    pw <- filt_pathways[i]
    kos <- pathway_kos[[pw]] %>% str_remove("ko:")
    if (is.null(kos)) {
      cat("No KOs found for pathway", pw, "\n")
      pathway_compl_sample[pw, ] <- NA
      next
    }
    for (j in 1:length(sample_names(phy_gene_f))) {
      sample <- sample_names(phy_gene_f)[j]
      kos_in_sample <- ko_list_samples[[sample]]
      matched_kos <- intersect(kos, kos_in_sample)
      pathway_compl_sample[pw, sample] <- length(matched_kos) / length(kos)
      if (length(matched_kos) == 0) {
        abund <- 0
      } else { 
        abund <- prune_samples(sample, phy_gene_f) %>% prune_taxa(matched_kos,.) %>% otu_table %>% sum
      }
    pathway_abund_sample[pw, sample] <- abund
    }
  }

  # CLR-normalise
  pathway_abund_sample_clr <- microbiome::transform(phyloseq(otu_table(pathway_abund_sample, taxa_are_rows = TRUE)), "clr") %>%
                          otu_table %>% data.frame

  write.csv(pathway_compl_sample, file.path(subdir, "pathway_completeness_per_sample.csv"), row.names = TRUE, quote = FALSE)
  write.csv(pathway_abund_sample, file.path(subdir, "pathway_abundance_per_sample.csv"), row.names = TRUE, quote = FALSE)
  write.csv(pathway_abund_sample_clr, file.path(subdir, "pathway_abundance_per_sample_clr.csv"), row.names = TRUE, quote = FALSE)
}

## Plot pathways completeness and abundance
pathway_compl_l <- as.data.frame(pathway_compl_sample) %>%
                    rownames_to_column("path") %>%
                    pivot_longer(-path, names_to = "Sample", values_to = "completeness")

pathway_abund_l <- as.data.frame(pathway_abund_sample_clr) %>%
                    rownames_to_column("path") %>%
                    pivot_longer(-path, names_to = "Sample", values_to = "abundance")

pathway_l <- full_join(pathway_compl_l, pathway_abund_l,
                        by = c("path", "Sample")) %>%
                        left_join(path_info_filt %>% select(path, path_name, path_class), by = "path") %>%
                        filter(!is.na(path_name))

pathway_l$Common.name <- phy_gene_f@sam_data$Common.name[match(pathway_l$Sample, sample_names(phy_gene_f))]
pathway_l$class <- str_remove(pathway_l$path_class, ".*; ")
pathway_l$class <- factor(pathway_l$class, levels = c("Overview Pathway", unique(pathway_l$class[-which(pathway_l$class == "Overview Pathway")])))
pathway_l$path_name_short <- str_trunc(pathway_l$path_name, 30, "right")

p <- ggplot(aes(x = Sample, y = path_name_short, fill = abundance), data = pathway_l) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma", na.value = "grey90") +
      labs(x = "Sample", y = "KEGG Pathway", fill = "CLR abundance") +
      facet_grid(cols = vars(Common.name), rows = vars(class), scales = "free", space = "free") +
      theme(axis.text.x = element_blank(),
            strip.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0),
            panel.background = element_rect(fill = "black", color = "black"))

ggsave(p, filename = file.path(subdir, "pathway_abundance_per_sample_heatmap.png"), width = 20, height = 30)

p <- ggplot(aes(x = Sample, y = path_name_short, fill = completeness), data = pathway_l) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma", na.value = "grey90") +
      labs(x = "Sample", y = "KEGG Pathway", fill = "Completeness") +
      facet_grid(cols = vars(Common.name), rows = vars(class), scales = "free", space = "free") +
      theme(axis.text.x = element_blank(),
            strip.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0),
            panel.background = element_rect(fill = "black", color = "black"))

ggsave(p, filename = file.path(subdir, "pathway_completeness_per_sample_heatmap.png"), width = 20, height = 30)

# Add completeness and abundance to pathway info
path_info_filt <- path_info_filt %>%
                  left_join(pathway_l[c("path", "completeness")] %>%
                            group_by(path) %>%
                            summarise(mean_completeness = mean(completeness, na.rm = TRUE)),
                            by = "path")

#### Plot gene abundance per species per pathway ####

# Write function
gene_abundances <- function(phyloseq, pathway, kos) {
  # Subset to KOs in pathway
  phy_sub <- prune_taxa(taxa_names(phyloseq) %in% str_remove(kos, "ko:"), phyloseq)
  cat(paste("Number of KOs in pathway:", length(kos), "\n"))
  cat(paste("Number of KOs present in dataset:", length(taxa_names(phy_sub)), "\n"))
  # Get abundance per sample
  abundances <- psmelt(phy_sub) %>% select(OTU, Sample, Abundance, Common.name, Species, diet.general, Order_grouped)
  # Plot
  p <- ggplot(aes(y = Sample, x = OTU, fill = Abundance), data = abundances) +
        geom_tile(stat = "identity") +
        scale_fill_viridis_c(option = "magma", na.value = "grey90") +
        facet_grid(rows = vars(Species), scales = "free", space = "free") +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(),
              strip.text.y = element_text(angle = 0),
              panel.background = element_rect(fill = "black", color = "black")) +
        labs(title = paste(pathway, path_info_df$path_name[path_info_df$path == pathway]))
  return(p)
}

# Only filtered pathways except for overview pathways
pathway_to_plot <- path_info_filt %>%
        filter(path_class != "Overview Pathway" & mean_completeness > 0.2 & path_class != "") %>%
        arrange(desc(mean_completeness)) %>% pull(path)

pdf(file.path(subdir, "gene_abundance_per_pathway.pdf"), width = 10, height = 10)
i <- 1
for (pathway in pathway_to_plot) {
  kos <- pathway_kos[[pathway]]
  cat(paste("Plotting pathway", i ,"-", pathway, "...", "\n"))
  if (is.null(kos)) {
    cat("No KOs found for pathway", pathway, "\n")
    next
  }
  p <- gene_abundances(phy_gene_f_clr, pathway, kos)
  plot(p)
  i <- i + 1
}
dev.off()

#### Calculate pathways completeness per microbial genus per sample ####
ko_per_taxon_df <- gene_str %>% filter(database == "KEGG") %>%
                    # At least an average coverage of 1
                    filter(mapped_reads >= 1) %>% filter(genus != "no support") %>%
                    select(Sample, genus, gene_id) %>%
                    # Keep genera with at least 1000 genes
                    group_by(genus, Sample) %>%
                    filter(n_distinct(gene_id) >= 1000)

# Get a list of KOs per genus per sample
if (file.exists(file.path(subdir, "pathway_completeness_per_sample_per_micrgenus.csv"))) {
  cat("Loading pathway_completeness_per_sample_per_micrgenus.csv\n")
  pathway_compl_s_g <- read.csv(file.path(subdir, "pathway_completeness_per_sample_per_micrgenus.csv"), stringsAsFactors = FALSE)
} else {
    ko_list_s_g <- list()

  for (s_g in unique(paste(ko_per_taxon_df$Sample, ko_per_taxon_df$genus, sep = "&"))) {
    cat("Extracting KOs present in", s_g, "...\n")
    ko_present <- ko_per_taxon_df %>% filter(Sample == str_split(s_g, "&")[[1]][1] & genus == str_split(s_g, "&")[[1]][2]) %>%
                  pull(gene_id) %>% unique
    ko_list_s_g[[s_g]] <- ko_present
  }

  pathway_compl_s_g <- matrix(nrow = length(filt_pathways), ncol = length(names(ko_list_s_g)))
  rownames(pathway_compl_s_g) <- filt_pathways
  colnames(pathway_compl_s_g) <- names(ko_list_s_g)

  for (i in 1:length(filt_pathways)) {
    pw <- filt_pathways[i]
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
}

# How different are the same taxa between species and across species

# Calculate correlation between pathways completeness
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

ggsave(p, filename = file.path(subdir, "pathway_completeness_cor.png"), width = 10, height = 10)

#### Calculate pathways completeness per microbial genus across dataset ####
## Genera seem to have more or less stable pathway profiles across host species,
## so we can consider them across the dataset
ko_per_taxon_df <- ko_per_taxon_df %>% ungroup %>% select(-Sample) %>% unique

# Get a list of KOs per genus per sample
ko_list_genus <- list()

for (g in unique(ko_per_taxon_df$genus)) {
  cat("Extracting KOs present in", g, "...\n")
  ko_present <- ko_per_taxon_df %>% filter(genus == g) %>%
                pull(gene_id) %>% unique
  ko_list_genus[[g]] <- ko_present
}

pathway_compl_genus <- matrix(nrow = length(filt_pathways), ncol = length(names(ko_list_genus)))
rownames(pathway_compl_genus) <- filt_pathways
colnames(pathway_compl_genus) <- names(ko_list_genus)

for (i in 1:length(filt_pathways)) {
  pw <- filt_pathways[i]
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
                    left_join(path_info_filt %>% select(path, path_name, path_class), by = "path")

tax_info <- gene_str %>% select(genus, family, order, class, phylum) %>% unique

pathway_compl_l$phylum <- tax_info$phylum[match(pathway_compl_l$genus, tax_info$genus)]
pathway_compl_l$order <- tax_info$order[match(pathway_compl_l$genus, tax_info$genus)]

pathway_compl_l$path_class <- str_remove(pathway_compl_l$path_class, ".*; ")

p <- ggplot(aes(x = genus, y = path, fill = completeness), data = pathway_compl_l) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma", na.value = "grey90") +
      labs(x = "Microbial genus", y = "KEGG Pathway", fill = "Completeness") +
      facet_grid(cols = vars(order), rows = vars(path_class), scales = "free", space = "free") +
      theme(strip.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0))

ggsave(p, filename = file.path(subdir, "pathway_completeness_per_taxon_heatmap.png"), width = 20, height = 20)

#######################
####  GET PHYLOSEQ ####
#######################

# Create a pathway phyloseq

otu <- otu_table(pathway_abund_sample, taxa_are_rows = TRUE)
tax <- tax_table(as.matrix(column_to_rownames(select(path_info_filt, c(path, path_name, path_class)), "path")))
sam <- sample_data(data.frame(phy_gene_f@sam_data))

phy_pathway <- phyloseq(otu, tax, sam)

# Keep pathways at least 20% complete in the dataset
compl_pathways <- path_info_filt$path[which(path_info_filt$mean_completeness > 0.2)]

phy_pathway <- prune_taxa(taxa_names(phy_pathway) %in% compl_pathways, phy_pathway)
phy_pathway_clr <- microbiome::transform(phy_pathway, "clr")

# Save
saveRDS(phy_pathway, file.path(subdir, "phy_pathway.RDS"))
saveRDS(phy_pathway_clr, file.path(subdir, "phy_pathway_clr.RDS"))

# Keep only metabolism pathways
phy_metabolism <- subset_taxa(phy_pathway, grepl("Metabolism;", path_class))
phy_metabolism_clr <- microbiome::transform(phy_metabolism, "clr")

# Save
saveRDS(phy_metabolism, file.path(subdir, "phy_metabolism.RDS"))
saveRDS(phy_metabolism_clr, file.path(subdir, "phy_metabolism_clr.RDS"))

#######################
#### PLOT PATHWAYS ####
#######################

#### Plot all relevant pathways ####
setwd(subdir)
dir.create("pathway_plots", showWarnings = FALSE)
setwd("pathway_plots")

path_info_filt <- path_info_filt %>% filter(mean_completeness > 0.2 & path_class != "Overview Pathway")

# Keep only KOs present in the data and average abundance in animalivores and herbivores
kos <- phy_gene_f_clr %>% subset_taxa(taxa_names(phy_gene_f) %in% str_remove(pathway_kos_filt$kos, "ko:")) %>%
    subset_samples(diet.general %in% c("Animalivore", "Herbivore")) %>% psmelt() %>%
    group_by(OTU, diet.general) %>% summarise(median_abundance = median(Abundance)) %>%
    pivot_wider(names_from = diet.general, values_from = median_abundance, values_fill = 0)

# Add pathway info
kos <- pathway_hits_df[c("path", "ko")] %>% mutate(ko = str_remove(ko, "ko:")) %>%
      right_join(kos, by = c("ko" = "OTU")) %>%
      filter(path %in% path_info_filt$path)

for (i in 1:nrow(path_info_filt)) {
  pathname <- path_info_filt$path_name[i]
  pathid <- path_info_filt$path[i]
  cat("GETTING INFO ON ", pathname, pathid, "\n")
  kos_path <- kos %>% filter(path == pathid) %>% column_to_rownames("ko") %>% select(-path) %>% as.matrix
  # Plot pathway
  tryCatch({
    path <- pathview(
      gene.data = kos_path,
      pathway.id = str_remove(pathid, "path:"),
      species = "ko",
      out.suffix = make.names(pathname),
      limit = list(gene = c(min(kos_path[,1], na.rm = TRUE), max(kos_path[,1], na.rm = TRUE)),
                   gene2 = c(min(kos_path[,2], na.rm = TRUE), max(kos_path[,2], na.rm = TRUE))),
      low = "white", mid = "pink", high = "red")
    }, error = function(e) {
    cat("!! Couldn't plot", pathname, "\n")})
}
