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
library(RColorBrewer)
library(ggnewscale)
library(scales)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
datadir <- normalizePath(file.path("..", "..", "output", "function", "data"))
subdir <- normalizePath(file.path("..", "..", "output", "function", "exploratory_plots")) # subdirectory for the output of this script

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

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

###################
#### RDA PLOTS ####
###################

#### GENE ABUNDANCE (CLR) ####

# Recode order and habitat as TRUE and FALSE
# Also scale protein, fiber and carbohydrate content
phy_gene_f_clr <- phy_gene_f_clr %>%
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
                    "ruminant", "marine", "cf", "cp")

# Ordinate using all data
ord <- ord_calc(phy_gene_f_clr, constraints = species_traits, method = "RDA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_genes.png"), p, width=8, height=6)

# Color by order

#### Get shape scales for plotting ####
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

p <- ord_plot(ord, colour="Order_grouped", shape="diet.general", alpha = 0.5, size = "Total_abundance") +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  scale_size(name = "# Features") +
  geom_phylopic(data = centroids(ord@ord, phy_gene_f), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_gene_f)$uid, width = 0.3, fill = "transparent") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(ncol = 1, size = 1, byrow = TRUE))

ggsave(p, filename = file.path(subdir, "gene_ordination.png"), width=8, height=10)

# Plot arrows
# Get loading arrows coordinaties
arrows <- arrow_coord(ord@ord, axes = c(1, 2))

# Get gene category
arrows$category <- as.character(phy_gene_f_clr@tax_table[match(rownames(arrows),  rownames(phy_gene_f_clr@tax_table)), "category"])

arrows$to_plot <- (rownames(arrows) %in% head(rownames(arrows), 500))

# Save
write.csv(rownames_to_column(arrows, "gene"), file.path(subdir, "gene_ordination_arrows.txt"), quote = FALSE, row.names = FALSE)

# Keep only strongest associations
arrows_filt <- arrows %>% filter(to_plot) %>%
              select(contains(c("1", "2")), category)

# Group uncommon categories
common_categories <- table(arrows_filt$category) %>% sort(decreasing = TRUE) %>% head(6) %>% names

arrows_filt <- arrows_filt %>%
    mutate(category_grouped = factor(case_when(category %in% common_categories ~ category,
                                            TRUE ~ "Other"), levels = c(common_categories, "Other")))

# Set colours for categories using colour brewer
arrow_colours <- brewer.pal(n = length(unique(arrows_filt$category_grouped))-1, name = "Dark2")
names(arrow_colours) <- unique(arrows_filt$category_grouped)[-length(unique(arrows_filt$category_grouped))] # Remove "Other" from names
arrow_colours["Other"] <- "grey90" # Set "Other" to grey

p <- ggplot(data = arrows_filt) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = category_grouped), linewidth = 0.5, alpha = 0.5) +
  scale_color_manual(values = arrow_colours, name = "Module") +
  xlab("RDA1") + ylab("RDA2")

ggsave(p, filename = file.path(subdir, "gene_ordination_arrows.png"), width=8, height=6)

#############################
#### NON KEGG ANNOTATIOS ####
#############################

# Plot peptidase diversity

phy_merops <- phy_gene_f %>% subset_taxa(database == "MEROPS")

phy_merops_comp <- phy_gene_f %>% transform("compositional") %>% subset_taxa(database == "MEROPS")

merops_richness <- 
    estimate_richness(rarefy_even_depth(phy_merops, 3000), measures = "Observed") %>%
    rownames_to_column("Sample") %>%
    left_join(select(data.frame(sample_data(phy_merops)), Common.name, cp, diet.general, Order_grouped) %>%
                as.data.frame() %>% rownames_to_column("Sample"))

p <-
    ggplot(merops_richness, aes(x = cp/100, y = Observed, colour = diet.general)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, alpha = 0.5, colour = "black") +
    scale_colour_manual(values = diet_palette, name = "Diet") +
    labs(x = "crude protein", y = "Peptidase richness (Observed)") +
    scale_x_continuous(labels = percent) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(subdir, "peptidase_richness.png"), width=8, height=6)

# Plot peptidase abundance

merops_abundance <- 
    data.frame(phy_merops_comp@otu_table) %>% rownames_to_column("OTU") %>%
    pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Abundance") %>%
    left_join(sample_data(phy_merops_comp) %>% data.frame() %>% rownames_to_column("Sample") %>%
                select(Sample, Common.name, cp, diet.general, Order_grouped)) %>%
    group_by(Common.name, cp, diet.general, Order_grouped) %>%
    # Get mean and quartile abundance of peptidases per species
    summarise(Abundance = mean(Abundance, na.rm = TRUE),
              q1 = quantile(Abundance, probs = 0.25),
              q3 = quantile(Abundance, probs = 0.75), .groups = "drop")

p <-
    ggplot(merops_abundance, aes(x = cp/100, y = Abundance, colour = diet.general)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, alpha = 0.5, colour = "black") +
    scale_colour_manual(values = diet_palette, name = "Diet") +
    labs(x = "crude protein", y = "Peptidase relative abundance") +
    scale_x_continuous(labels = percent) +
    scale_y_continuous(labels = percent) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(subdir, "peptidase_abundance.png"), width=8, height=6)

# Plot CAZy enzyme diversity

phy_cazy <- phy_gene_f %>% subset_taxa(database == "CAZY")

phy_cazy_comp <- phy_gene_f %>% transform("compositional") %>% subset_taxa(database == "CAZY")

cazy_richness <- 
    estimate_richness(rarefy_even_depth(phy_cazy, 1000), measures = "Observed") %>%
    rownames_to_column("Sample") %>%
    left_join(select(data.frame(sample_data(phy_cazy)), Common.name, cp, diet.general, Order_grouped) %>%
                as.data.frame() %>% rownames_to_column("Sample"))

p <-
    ggplot(cazy_richness, aes(x = cp/100, y = Observed, colour = diet.general)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, alpha = 0.5, colour = "black") +
    scale_colour_manual(values = diet_palette, name = "Diet") +
    labs(x = "crude protein", y = "CAZy richness (Observed)") +
    scale_x_continuous(labels = percent) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(subdir, "cazy_richness.png"), width=8, height=6)

# Plot CAZy abundance
cazy_abundance <- 
    data.frame(phy_cazy_comp@otu_table) %>% rownames_to_column("OTU") %>%
    pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Abundance") %>%
    left_join(sample_data(phy_cazy_comp) %>% data.frame() %>% rownames_to_column("Sample") %>%
                select(Sample, Common.name, cp, diet.general, Order_grouped)) %>%
    group_by(Common.name, cp, diet.general, Order_grouped) %>%
    # Get mean and quartile abundance of CAZy enzymes per species
    summarise(Abundance = mean(Abundance, na.rm = TRUE),
              q1 = quantile(Abundance, probs = 0.25),
              q3 = quantile(Abundance, probs = 0.75), .groups = "drop")

p <-
    ggplot(cazy_abundance, aes(x = cp/100, y = Abundance, colour = diet.general)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, alpha = 0.5, colour = "black") +
    scale_colour_manual(values = diet_palette, name = "Diet") +
    labs(x = "crude protein", y = "CAZy relative abundance") +
    scale_x_continuous(labels = percent) +
    scale_y_continuous(labels = percent) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(subdir, "cazy_abundance.png"), width=8, height=6)
