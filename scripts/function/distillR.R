##### INPUT FOR FUNCTIONAL ANALYSIS #####

#### Prepares tables and phyloseq objects for analysis of functional annotations

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(distillR)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "function")) # subdirectory for the output of this script
taxdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Directory with taxonomy analysis
datadir <- normalizePath(file.path(outdir, "data")) # Directory with data files
subdir <- normalizePath(file.path(outdir, "distillR")) # subdirectory for the output of this script

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

phy_gene_f <- readRDS(file.path(datadir, "phy_gene_f.RDS"))

# Stratified sample data
gene_str <- read.table(file.path(datadir, "gene_abundance_stratified_modified.tsv"),
                      quote = "", comment.char = "", header = TRUE, sep = "\t")

# Differentially abundant taxa
ancom_res <- read.csv(file.path(taxdir, "statistical_tests", "ancomb_long.csv"))

################################
#### DISTILL COMMUNITY WIDE ####
################################

# Identify functions in entire communities

if(file.exists(file.path(subdir, "GIFTs_community.csv"))) {
  cat("Loading GIFTs from file\n")
  GIFTs <- read.csv(file.path(subdir, "GIFTs_community.csv"), row.names = 1)
} else {
  # Prep input for distillR
  data <- psmelt(phy_gene_f) %>% select(OTU, Sample, gene_name, Abundance, Total_abundance) %>%
      # Remove zero abundances and samples with low read depth
      filter(Abundance > 0)
  cat("Running distillR to get GIFTs\n")
  GIFTs <- distill(data, GIFT_db, genomecol=2, annotcol=c(1, 3))
  write.csv(GIFTs, file.path(subdir, "GIFTs_community.csv"), row.names = TRUE)
}

#Aggregate bundle-level GIFTs into the compound level
GIFTs_elements <- to.elements(GIFTs, GIFT_db)
GIFTs_elements <- GIFTs_elements[, colSums(GIFTs_elements) > 0] # Remove empty columns

GIFTs_elements_long <- 
  GIFTs_elements %>% data.frame %>% rownames_to_column("Sample") %>%
  pivot_longer(cols = -c(Sample), names_to = "Code_element", values_to = "Completeness") %>%
  left_join(select(rownames_to_column(data.frame(phy_gene_f@sam_data), "Sample"), c(Sample, Common.name, Order, diet.general))) %>%
  left_join(unique(select(GIFT_db, c("Code_element", "Element", "Function", "Domain")))) %>%
  # Turn function into factor
  arrange(Domain, Function) %>% mutate(Function = factor(Function, levels = unique(Function)))

# Remove communities with low total abundance
low_depth <- data.frame(phy_gene_f@sam_data) %>% filter(Total_abundance < 10^5) %>%
                  rownames

GIFTs_elements_long <- GIFTs_elements_long %>% filter(!Sample %in% low_depth) %>% unique

# Match this on wide table
GIFTs_elements <- GIFTs_elements[!rownames(GIFTs_elements) %in% low_depth, ]

p <- ggplot(GIFTs_elements_long, aes(y=Sample, x=Element)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Compound", fill="Completeness") +
  facet_grid(rows = vars(Common.name), cols = vars(Function), scales = "free", space = "free")
  
ggsave(p, filename = file.path(subdir, "GIFTs_elements_community.png"), width = 25, height = 20)

#Aggregate element-level GIFTs into the function level
GIFTs_functions <- to.functions(GIFTs_elements, GIFT_db)
GIFTs_functions <- GIFTs_functions[, colSums(GIFTs_functions) > 0] # Remove empty columns

GIFTs_functions_long <- 
  GIFTs_functions %>% data.frame %>% rownames_to_column("Sample") %>%
  pivot_longer(cols = -c(Sample), names_to = "Code_function", values_to = "Completeness") %>%
  left_join(select(rownames_to_column(data.frame(phy_gene_f@sam_data), "Sample"), c(Sample, Common.name, Order, diet.general))) %>%
  left_join(unique(select(GIFT_db, c("Code_function", "Function", "Domain")))) %>%
  # Remove samples in incomplete communities
  filter(!Sample %in% low_depth)

p <- ggplot(GIFTs_functions_long, aes(y=Sample, x=Function)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Function", fill="Completeness") +
  facet_grid(rows = vars(Common.name), cols = vars(Domain), scales = "free", space = "free")

ggsave(p, filename = file.path(subdir, "GIFTs_functions_community.png"), width = 10, height = 15)

#### PCA ####
ord <- prcomp(GIFTs_elements)

# Get some info for plotting
var_explained <- round(ord$sdev^2 * 100 / sum(ord$sdev^2), 1) # Variance explained

loadings <- data.frame(Code_element = rownames(ord$rotation[,c(1,2)]), ord$rotation[,c("PC1", "PC2")]) %>%
  left_join(unique(select(GIFT_db, c("Code_element", "Element", "Function", "Domain")))) %>%
  mutate(Variables = paste(Element, Domain, sep = " ")) %>%
  # Keep the longest arrows
  arrange(desc(sqrt(PC1^2 + PC2^2))) %>%
  select(PC1, PC2, Variables, Function)

metadata <- data.frame(phy_gene_f@sam_data)[rownames(ord$x),]

# Color by order, shape by diet
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

# Plot
pca <- ggplot(aes(x = PC1, y = PC2, colour = metadata$Order), data = data.frame(ord$x)) +
  geom_point(aes(shape = metadata$diet.general, size = metadata$Total_abundance)) +
  scale_colour_manual(values = order_palette, name = "") +
  scale_shape_manual(values = diet_shape_scale, name = "") +
  scale_size_continuous(trans = "log10") +
  xlab(paste("PC1 -", var_explained[1], "%")) +
  ylab(paste("PC2 -", var_explained[2], "%")) +
  geom_segment(data = loadings[1:8,], aes(x = 0, y = 0, xend = (PC1*10),
                                       yend = (PC2*10)), arrow = arrow(length = unit(0.5, "picas")),
               color = "black") +
  geom_label(data = loadings[1:8,], aes(x = (PC1*10), y = (PC2*10), label = Variables),
            size = 2, hjust = 0.5, vjust = -0.5, color = "black", alpha = 0.7) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow =3), shape = guide_legend(nrow = 3))

ggsave(pca, filename = file.path(subdir, "PCA_GIFTs_community.png"), width = 10, height = 10)

# Extract info on functions explaining the most variance (top 20)
var_gifts <- loadings %>% head(20) %>%
    separate(Variables, into = c("Element", "Domain"), sep = " ") %>%
    # Get genes implicated in these orders
    left_join(GIFT_db) %>%
    select(Element, Domain, Function, Definition)

write.csv(var_gifts, file.path(subdir, "GIFTs_explaining_variance.csv"), row.names = FALSE)

##########################
#### DISTILL BY TAXON ####
##########################

# Identify diff abundant taxa
#diff_taxa <- select(gene_str, c(species, genus)) %>% unique %>%
#      inner_join(ancom_res, by = c("genus" = "taxon")) %>%
#      filter(pval < 0.05 & sensitivity == TRUE & lfc > 0) %>%
#      select(species, genus, term) %>% unique %>%
#      group_by(species, genus) %>%
#      summarise(term = paste(term, collapse = ", "))

if(file.exists(file.path(subdir, "GIFTs_by_taxon.csv"))) {
  cat("Loading GIFTs from file\n")
  GIFTs <- read.csv(file.path(subdir, "GIFTs_by_taxon.csv"), row.names = 1)
} else {
  cat("Running distillR to get GIFTs\n")
  data <- gene_str %>% select(gene_id, database, species, sample, totalAvgDepth) %>%
    # Keep entries with a totalAvgDepth > 1
    filter(totalAvgDepth > 1) %>%
    # Keep only taxa identified at the species level
    filter(species != "no support") %>%
    #filter(species %in% diff_taxa$species) %>%
    # Keep genes that are present in phy_gene_f and are KEGG
    filter(gene_id %in% taxa_names(phy_gene_f)) %>%
    filter(database == "KEGG") %>%
    # Combine species and sample into a "genome" column
    mutate(genome = paste(species, sample, sep = " - ")) %>%
    # Remove genomes with less than 200 genes
    group_by(genome) %>% %>% mutate(n_genes = n_distinct(gene_id)) %>% filter(n_genes > 400) %>%
    # Keep microbial species found in at least 5 samples
    group_by(species) %>% filter(n_distinct(sample) > 5) %>% ungroup %>%
    select(genome, gene_id)
  GIFTs <- distill(data, GIFT_db, genomecol=1, annotcol=2)
  write.csv(GIFTs, file.path(subdir, "GIFTs_by_taxon.csv"), row.names = TRUE)
}

#Aggregate bundle-level GIFTs into the compound level
GIFTs_elements <- to.elements(GIFTs, GIFT_db)
GIFTs_elements <- GIFTs_elements[, colSums(GIFTs_elements) > 0] # Remove empty columns

GIFTs_elements_long <- 
  GIFTs_elements %>% data.frame %>% rownames_to_column("genome") %>%
  separate(genome, into = c("species", "Ext.ID"), sep = " - ", remove = FALSE) %>%
  pivot_longer(cols = -c(species, Ext.ID, genome), names_to = "Code_element", values_to = "Completeness") %>%
  left_join(select(rownames_to_column(data.frame(phy_gene_f@sam_data), "Sample"), c(Ext.ID, Sample, Common.name, Order, diet.general))) %>%
  left_join(unique(select(GIFT_db, c("Code_element", "Element", "Function", "Domain")))) %>%
  mutate(genus = str_remove(species, " .*")) %>%
  # Turn function into factor
  arrange(Domain, Function) %>% mutate(Function = factor(Function, levels = unique(Function)))

p <- ggplot(GIFTs_elements_long, aes(y=genome, x=Element)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Compound", fill="Completeness") +
  facet_grid(rows = vars(genus), cols = vars(Function), scales = "free", space = "free")
  
ggsave(p, filename = file.path(subdir, "GIFTs_elements_by_taxon.png"), width = 25, height = 20)

#Aggregate element-level GIFTs into the function level
GIFTs_functions <- to.functions(GIFTs_elements, GIFT_db)
GIFTs_functions <- GIFTs_functions[, colSums(GIFTs_functions) > 0] # Remove empty columns

GIFTs_functions_long <- 
  GIFTs_functions %>% data.frame %>% rownames_to_column("genome") %>%
  separate(genome, into = c("species", "Ext.ID"), sep = " - ", remove = FALSE) %>%
  pivot_longer(cols = -c(species, Ext.ID, genome), names_to = "Code_function", values_to = "Completeness") %>%
  left_join(select(rownames_to_column(data.frame(phy_gene_f@sam_data), "Sample"), c(Ext.ID, Sample, Common.name, Order_grouped, diet.general))) %>%
  left_join(unique(select(GIFT_db, c("Code_function", "Function", "Domain")))) %>%
  mutate(genus = str_remove(species, " .*"))

p <- ggplot(GIFTs_functions_long, aes(y=genome, x=Function)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Function", fill="Completeness") +
  facet_grid(rows = vars(genus), cols = vars(Domain), scales = "free", space = "free")

ggsave(p, filename = file.path(subdir, "GIFTs_functions_by_taxon.png"), width = 10, height = 15)

#### PCA ####
ord <- prcomp(GIFTs_elements)

# Get some info for plotting
var_explained <- round(ord$sdev^2 * 100 / sum(ord$sdev^2), 1) # Variance explained

loadings <- data.frame(Code_element = rownames(ord$rotation[,c(1,2)]), ord$rotation[,c("PC1", "PC2")]) %>%
  left_join(unique(select(GIFT_db, c("Code_element", "Element", "Function", "Domain")))) %>%
  mutate(Variables = paste(Element, Domain, sep = " ")) %>%
  # Keep the longest arrows
  arrange(desc(sqrt(PC1^2 + PC2^2))) %>%
  slice(1:12) %>% select(PC1, PC2, Variables, Function)

metadata <- GIFTs_functions_long[match(rownames(ord$x), GIFTs_functions_long$genome),] %>%
  select(Sample, Common.name, Order_grouped, diet.general, genome, genus) %>%
  # Group genomes to highlight most abundant genera in plot
  group_by(genus) %>%
  mutate(n_genomes = n_distinct(genome))

top_genera <- metadata %>% arrange(desc(n_genomes)) %>% pull(genus) %>% unique %>% head(8)

metadata <- metadata %>%
  mutate(genus_grouped = case_when(genus %in% top_genera ~ genus,
                                   TRUE ~ "Other")) %>%
  mutate(genus_grouped = factor(genus_grouped, levels = c(top_genera, "Other")))

# Get palette with RColorBrewer
genus_palette <- brewer.pal(n = length(unique(metadata$genus_grouped))-1, name = "Dark2")
names(genus_palette) <- unique(metadata$genus_grouped)[-which(unique(metadata$genus_grouped) == "Other")]
genus_palette["Other"] <- "grey"

# Plot
pca <- ggplot(aes(x = PC1, y = PC2, colour=metadata$genus_grouped), data = data.frame(ord$x)) +
  geom_point(size=2, alpha = 0.8) +
  scale_colour_manual(values = genus_palette, name = "") +
  xlab(paste("PC1 -", var_explained[1], "%")) +
  ylab(paste("PC2 -", var_explained[2], "%")) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*6),
                                       yend = (PC2*6)), arrow = arrow(length = unit(0.5, "picas")),
               color = "black") +
  geom_label(data = loadings, aes(x = (PC1*6), y = (PC2*6), label = Variables),
            size = 2, hjust = 0.5, vjust = -0.5, color = "black", alpha = 0.7) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))

ggsave(pca, filename = file.path(subdir, "PCA_GIFTs_by_taxon.png"), width = 10, height = 10)

# Plot
pca <- ggplot(aes(x = PC1, y = PC2, colour=metadata$Order_grouped, shape = metadata$diet.general), data = data.frame(ord$x)) +
  geom_point(size=1, alpha = 0.8) +
  scale_colour_manual(values = order_palette, name = "") +
  scale_shape_manual(values = diet_shape_scale, name = "") +
  xlab(paste("PC1 -", var_explained[1], "%")) +
  ylab(paste("PC2 -", var_explained[2], "%")) +
  facet_wrap(~metadata$genus_grouped, ncol = 2) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3), shape = guide_legend(nrow = 3))

ggsave(pca, filename = file.path(subdir, "PCA_GIFTs_by_taxon_faceted.png"), width = 10, height = 10)
