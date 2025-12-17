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
library(rphylopic)
library(ggExtra)
library(RColorBrewer)
library(distillR)
library(microViz)

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
  # Turn function into factor
  arrange(Domain, Function) %>% mutate(Function = factor(Function, levels = unique(Function)))

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

#### RDA ####

phy_gifts_el <- phyloseq(otu_table(GIFTs_elements, taxa_are_rows = FALSE),
                        sample_data(phy_gene_f))

# Recode order and habitat as TRUE and FALSE
phy_gifts_el <- phy_gifts_el %>%
        ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  Rodentia = (Order == "Rodentia"),
                  ruminant = (digestion == "Ruminant"),
                  marine = (habitat.general == "Marine"),
                  animalivore = (diet.general == "Animalivore"))

# Species traits to use as constraints
species_traits <- c("Artiodactyla", "Perissodactyla", "Primates", "Rodentia",
                    "ruminant", "marine", "Fruit", "Animal")

# Ordinate using all data
ord <- ord_calc(phy_gifts_el, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "ord_community_screeplot.png"), p, width=3, height=3)

#### Get shape scales for plotting ####
order_shape_scale <- c("Carnivora" = 4, "Primates" = 19, "Artiodactyla" = 5, "Perissodactyla" = 2, "Rodentia" = 1, "Rest" = 12, "Proboscidea_Sirenia" = 12)

p <- ord_plot(ord, colour="diet.general", shape="Order_grouped", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=order_shape_scale, name = "Order") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  geom_phylopic(data = centroids(ord@ord, phy_gifts_el), aes(colour = diet.general), uuid = centroids(ord@ord, phy_gifts_el)$uid, width = 0.1, fill = "transparent") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(shape = guide_legend(ncol = 2), colour = guide_legend(ncol = 2))

p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(p, filename = file.path(subdir, "ord_community_diet.png"), width=6, height=6)

# Plot arrows
# Get loading arrows coordinaties
arrows <- arrow_coord(ord@ord, axes = c(1, 2))

element_info <- GIFTs_elements_long %>% select(Code_element, Element, Function, Domain) %>%
                  distinct() %>% column_to_rownames("Code_element")

# Get gene category
arrows <- arrows %>% cbind(element_info[rownames(arrows),])

arrows$plot_label <- (rownames(arrows) %in% head(rownames(arrows), 10))

# Save
write.csv(arrows, file.path(subdir, "ord_community_arrows.csv"), quote = TRUE, row.names = FALSE)

# Keep only strongest associations
arrows_sum <- arrows %>% group_by(Function, Domain) %>%
              summarise(RDA1 = mean(RDA1), RDA2 = mean(RDA2)) %>%
              arrange(desc(sqrt(RDA1^2 + RDA2^2))) %>% head(10)

p <- ggplot(data = arrows) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = Domain), linewidth = 0.5, alpha = 0.3, linetype = "dashed") +
  geom_label(data = arrows[arrows$plot_label,], aes(x = RDA1*1.1, y = RDA2*1.1, label = Element, colour = Domain), size = 2, fill = "white", alpha = 0.7) +
  geom_segment(data = arrows_sum, aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = Domain), linewidth = 1, alpha = 0.7) +
  geom_label(data = arrows_sum, aes(x = RDA1*1.1, y = RDA2*1.1, label = str_remove(Function, tolower(Domain)), fill = Domain), size = 2, colour = "white", alpha = 0.7) +
  scale_colour_manual(values = c("Degradation" = "#D55E00", "Biosynthesis" = "#0072B2", "Structure" = "#6A6A6A"), name = "") +
  scale_fill_manual(values = c("Degradation" = "#D55E00", "Biosynthesis" = "#0072B2", "Structure" = "#6A6A6A"), name = "") +
  xlab("RDA1") + ylab("RDA2") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1))

ggsave(p, filename = file.path(subdir, "ord_community_arrows.png"), width=6, height=5)

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
