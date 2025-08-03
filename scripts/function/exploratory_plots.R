##### EXPLORATORY PLOTS #####

#### Plots exploratory plots such as ordinations and heatmaps

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(microViz)
library(microbiome)
library(rphylopic)
library(ggplot2)
library(ggnewscale)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "function")) # subdirectory for the output of this script
subdir <- normalizePath(file.path(outdir, "exploratory_plots"))

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
phy_gene <- readRDS(file.path(outdir, "phy_gene.RDS"))
phy_gene_clr <- readRDS(file.path(outdir, "phy_gene_clr.RDS"))

low_content_samples <- read.csv(file.path(outdir, "low_content_samples.txt"), header = TRUE)

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

###################
#### RDA PLOTS ####
###################

#### GENE ABUNDANCE (CLR) ####

# Filter out low content samples
phy_gene_filt <- phy_gene_clr %>% subset_samples(!sample_names(phy_gene_clr) %in% low_content_samples$x)

# Recode order and habitat as TRUE and FALSE
# Also scale protein, fiber and carbohydrate content
phy_gene_filt <- phy_gene_filt %>%
        ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  Rodentia = (Order == "Rodentia"),
                  ruminant = (digestion == "Ruminant"),
                  marine = (habitat.general == "Marine"),
                  herbivore = (diet.general == "Herbivore"),
                  frugivore = (diet.general == "Frugivore"),
                  omnivore = (diet.general == "Omnivore"),
                  animalivore = (diet.general == "Animalivore"))


# Species traits to use as constraints
species_traits <- c("Artiodactyla", "Perissodactyla", "Primates", "Rodentia",
                    "ruminant", "marine", "herbivore", "frugivore", "omnivore", "animalivore")

# Ordinate using all data
ord <- ord_calc(phy_gene_filt, constraints = species_traits, method = "RDA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_genes.png"), p, width=8, height=6)

# Color by order
#### Get shape scales for plotting ####
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

# Get loading arrows coordinaties
arrows <- arrow_coord(ord@ord, phy_gene_filt)

# Get gene module
arrows$module <- as.character(phy_gene_filt@tax_table[match(rownames(arrows),  rownames(phy_gene_clr@tax_table)), "module"])

arrows$to_plot <- (arrows$r > quantile(arrows$r, 0.50))

# Save
write.csv(arrows, file.path(subdir, "gene_ordination_arrows.txt"), quote = FALSE, row.names = FALSE)

# Keep only strongest associations
arrows_filt <- arrows %>% filter(to_plot) %>%
              select(contains(c("1", "2", "3", "4")), module)

# Group uncommon modules
common_modules <- table(arrows_filt$module) %>% sort(decreasing = TRUE) %>% head(6) %>% names

arrows_filt <- arrows_filt %>%
    mutate(module_grouped = factor(case_when(module %in% common_modules ~ module,
                                            TRUE ~ "Other"), levels = c(common_modules, "Other")))

arrow_colours <- c("Flagellar Assembly" = "#FF5454",
                   "Ribosome, archaea & Ribosome, eukaryotes" = "#FFD454",
                   "RNA polymerase, archaea" = "#FFC003",
                   "Menaquinone biosynthesis" = "#50F250",
                   "Ubiquinone biosynthesis" = "#6565EB",
                   "Siroheme biosynthesis, glutamate => siroheme" = "#F276E7",
                   "Other" = "grey80")

p <- ord_plot(ord, colour="Order_grouped", shape="diet.general", alpha = 0.5, size = "Total_abundance") +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  scale_size(name = "# Features") +
  geom_phylopic(data = centroids(ord@ord, phy_gene_filt), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_gene_filt)$uid, width = 0.3, fill = "transparent") +
  new_scale_colour() +
  geom_segment(data = arrows_filt, aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = module_grouped), linewidth = 0.8, alpha = 0.5) +
  scale_color_manual(values = arrow_colours, name = "Module") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(ncol = 1, size = 1, byrow = TRUE))

ggsave(p, filename = file.path(subdir, "gene_ordination.png"), width=8, height=10)

#############################
#### NON KEGG ANNOTATIOS ####
#############################

# Plot peptidase diversity

phy_merops <- phy_gene %>% subset_taxa(database == "MEROPS") %>%
    subset_samples(!(sample_names(phy_gene) %in% low_content_samples$x))

merops_richness <- 
    estimate_richness(rarefy_even_depth(phy_merops, 5000), measures = "Observed") %>%
    rownames_to_column("Sample") %>%
    left_join(sample_data(phy_merops) %>% as.data.frame() %>% rownames_to_column("Sample"))

p <-
    ggplot(merops_richness, aes(x = Species, y = Observed, fill = diet.general)) +
    geom_boxplot() +
    scale_fill_manual(values = diet_palette) +
    labs(x = "Order", y = "Peptidase richness (Observed)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(subdir, "peptidase_richness.png"), width=8, height=6)

# Plot peptidase abundance

merops_abundance <- 
    data.frame(transform(phy_merops, "compositional")@otu_table) %>% rownames_to_column("OTU") %>%
    pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Abundance") %>%
    left_join(sample_data(phy_merops) %>% data.frame() %>% rownames_to_column("Sample") %>% select(Sample, Species, diet.general)) %>%
    group_by(Sample, Species, diet.general) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

p <-
    ggplot(merops_abundance, aes(x = Species, y = Abundance, fill = diet.general)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = diet_palette) +
    labs(x = "Order", y = "Peptidase abundance (CLR)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(subdir, "peptidase_abundance.png"), width=8, height=6)

# Plot CAZy enzyme diversity

phy_cazy <- phy_gene %>% subset_taxa(database == "CAZY") %>%
    subset_samples(!(sample_names(phy_gene) %in% low_content_samples$x))

cazy_richness <- 
    estimate_richness(rarefy_even_depth(phy_cazy, 500), measures = "Observed") %>%
    rownames_to_column("Sample") %>%
    left_join(sample_data(phy_merops) %>% as.data.frame() %>% rownames_to_column("Sample"))

p <-
    ggplot(cazy_richness, aes(x = Species, y = Observed, fill = diet.general)) +
    geom_boxplot() +
    scale_fill_manual(values = diet_palette) +
    labs(x = "Order", y = "CAZy richness (Observed)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(subdir, "cazy_richness.png"), width=8, height=6)

# Plot CAZy enzyme abundance
cazy_abundance <- 
    data.frame(transform(phy_cazy, "compositional")@otu_table) %>% rownames_to_column("OTU") %>%
    pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Abundance") %>%
    left_join(sample_data(phy_cazy) %>% data.frame() %>% rownames_to_column("Sample") %>% select(Sample, Species, diet.general)) %>%
    group_by(Sample, Species, diet.general) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

p <-
  ggplot(cazy_abundance, aes(x = Species, y = Abundance, fill = diet.general)) +
  geom_boxplot() +
  scale_fill_manual(values = diet_palette) +
  labs(x = "Order", y = "CAZy abundance (CLR)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(subdir, "cazy_abundance.png"), width=8, height=6)

