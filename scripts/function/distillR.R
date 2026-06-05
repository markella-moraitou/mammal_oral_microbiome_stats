##### INPUT FOR FUNCTIONAL ANALYSIS #####

#### Uses distillR R package to obtain genome-inferred functional traits for communities and taxa

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

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

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
  # Turn function and common name into factors
  arrange(Domain, Function, diet.general) %>%
  mutate(Function = factor(Function, levels = unique(Function)),
         Common.name = factor(Common.name, levels = unique(Common.name)))

p <- ggplot(GIFTs_elements_long, aes(y=Sample, x=Element)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90)) +
  labs(y="Sample", x="Compound", fill="Completeness") +
  facet_grid(rows = vars(Common.name),
             cols = vars(gsub(
                gsub(Function, pattern = " biosynthesis", replacement = "\nbiosynthesis"),
                pattern = " degradation", replacement = "\ndegradation")),
             scales = "free", space = "free")
  
ggsave(p, filename = file.path(subdir, "GIFTs_elements_community.png"), width = 25, height = 20)

write.csv(GIFTs_elements_long, file.path(subdir, "GIFTs_elements_long.csv"), row.names = FALSE)

# Create phyloseq
phy_gifts_el <- phyloseq(otu_table(GIFTs_elements, taxa_are_rows = FALSE),
                        sample_data(phy_gene_f))

saveRDS(phy_gifts_el, file.path(datadir, "phy_gifts_el.RDS"))

#Aggregate element-level GIFTs into the function level
GIFTs_functions <- to.functions(GIFTs_elements, GIFT_db)
GIFTs_functions <- GIFTs_functions[, colSums(GIFTs_functions) > 0] # Remove empty columns

GIFTs_functions_long <- 
  GIFTs_functions %>% data.frame %>% rownames_to_column("Sample") %>%
  pivot_longer(cols = -c(Sample), names_to = "Code_function", values_to = "Completeness") %>%
  left_join(select(rownames_to_column(data.frame(phy_gene_f@sam_data), "Sample"), c(Sample, Common.name, Order, diet.general))) %>%
  left_join(unique(select(GIFT_db, c("Code_function", "Function", "Domain"))))

p <- ggplot(GIFTs_functions_long, aes(y=Sample, x=Function)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Function", fill="Completeness") +
  facet_grid(rows = vars(Common.name), cols = vars(Domain), scales = "free", space = "free")

ggsave(p, filename = file.path(subdir, "GIFTs_functions_community.png"), width = 10, height = 15)

write.csv(GIFTs_functions_long, file.path(subdir, "GIFTs_functions_long.csv"), row.names = FALSE)

##########################
#### DISTILL BY TAXON ####
##########################

if(file.exists(file.path(subdir, "GIFTs_by_taxon.csv"))) {
  cat("Loading GIFTs from file\n")
  GIFTs <- read.csv(file.path(subdir, "GIFTs_by_taxon.csv"), row.names = 1)
} else {
  cat("Running distillR to get GIFTs\n")
  data <- gene_str %>% select(gene_id, database, species, Sample, mapped_reads) %>%
    # Keep entries with a totalAvgDepth > 1
    filter(mapped_reads > 1) %>%
    # Keep only taxa identified at the species level
    filter(species != "no support") %>%
    #filter(species %in% diff_taxa$species) %>%
    # Keep genes that are present in phy_gene_f and are KEGG
    filter(gene_id %in% taxa_names(phy_gene_f)) %>%
    filter(database == "KEGG") %>%
    # Combine species and sample into a "genome" column
    mutate(genome = paste(species, Sample, sep = " - ")) %>%
    # Remove genomes with less than 200 genes
    group_by(genome) %>% mutate(n_genes = n_distinct(gene_id)) %>% filter(n_genes > 400) %>%
    # Keep microbial species found in at least 5 samples
    group_by(species) %>% filter(n_distinct(Sample) > 5) %>% ungroup %>%
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

write.csv(GIFTs_elements_long, file.path(subdir, "GIFTs_elements_by_taxon.csv"), row.names = FALSE)

# Plot genera with at least 10 taxa
GIFTs_elements_long <- GIFTs_elements_long %>% group_by(genus) %>%
  filter(n_distinct(genome) > 10) %>% ungroup()

p <- ggplot(GIFTs_elements_long, aes(x=genome, y=Element)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90)) +
  labs(y="Compound", x="Sample", fill="Completeness") +
  facet_grid(cols = vars(genus), rows = vars(Function), scales = "free", space = "free")

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

write.csv(GIFTs_functions_long, file.path(subdir, "GIFTs_functions_by_taxon.csv"), row.names = FALSE)
