##### EXPLORATORY PLOTS #####

#### Plots exploratory plots such as ordinations and heatmaps

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(phyloseq)
library(microViz)
library(microbiome)
library(rphylopic)
library(ggplot2)
library(RColorBrewer)
library(ggnewscale)
library(scales)
library(ggExtra)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
datadir <- normalizePath(file.path("..", "..", "output", "function", "data"))
pathdir <- normalizePath(file.path("..", "..", "output", "function", "pathway_completeness")) # Directory with pathway analysis output
subdir <- normalizePath(file.path("..", "..", "output", "function", "multivariate_func")) # subdirectory for the output of this script

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

phy_pathway <- readRDS(file.path(pathdir, "phy_pathway.RDS"))
phy_pathway_clr <- readRDS(file.path(pathdir, "phy_pathway_clr.RDS"))

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

###################
#### RDA PLOTS ####
###################

#### GENE ABUNDANCE (CLR) ####

# Recode order and habitat as TRUE and FALSE
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
                    "ruminant", "marine", "Fruit", "Animal")

# Ordinate using all data
ord <- ord_calc(phy_gene_f_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_genes.png"), p, width=3, height=3)

# Color by order

#### Get shape scales for plotting ####
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

p <- ord_plot(ord, colour="Order_grouped", shape="diet.general", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  geom_phylopic(data = centroids(ord@ord, phy_gene_f), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_gene_f)$uid, width = 0.3, fill = "transparent") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(shape = guide_legend(ncol = 2), colour = guide_legend(ncol = 2))
  
p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(p, filename = file.path(subdir, "gene_ordination_order.png"), width=6, height=6)

# Color by diet

#### Get shape scales for plotting ####
order_shape_scale <- c("Carnivora" = 4, "Primates" = 19, "Artiodactyla" = 5, "Perissodactyla" = 2, "Rodentia" = 1, "Rest" = 12, "Proboscidea_Sirenia" = 12)

p <- ord_plot(ord, colour="diet.general", shape="Order_grouped", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=order_shape_scale, name = "Order") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  geom_phylopic(data = centroids(ord@ord, phy_gene_f), aes(colour = diet.general), uuid = centroids(ord@ord, phy_gene_f)$uid, width = 0.3, fill = "transparent") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(shape = guide_legend(ncol = 2), colour = guide_legend(ncol = 2))

p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(p, filename = file.path(subdir, "gene_ordination_diet.png"), width=6, height=6)

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
  xlab("RDA1") + ylab("RDA2") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 2))

ggsave(p, filename = file.path(subdir, "gene_ordination_arrows.png"), width=6, height=6)

#### PATH ABUNDANCE (CLR) ####

phy_pathway_clr <- phy_pathway_clr %>% 
        ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  Rodentia = (Order == "Rodentia"),
                  ruminant = (digestion == "Ruminant"),
                  marine = (habitat.general == "Marine"),
                  animalivore = (diet.general == "Animalivore"))

# Ordinate using all data
ord <- ord_calc(phy_pathway_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_pathways.png"), p, width=3, height=3)

# Color by order
p <- ord_plot(ord, colour="Order_grouped", shape="diet.general", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  geom_phylopic(data = centroids(ord@ord, phy_pathway_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_pathway_clr)$uid, width = 0.2, fill = "transparent") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(shape = guide_legend(ncol = 2), colour = guide_legend(ncol = 2))

p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(p, filename = file.path(subdir, "pathway_ordination_order.png"), width=6, height=6)

# Color by diet
p <- ord_plot(ord, colour="diet.general", shape="Order_grouped", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=order_shape_scale, name = "Order") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  geom_phylopic(data = centroids(ord@ord, phy_pathway_clr), aes(colour = diet.general), uuid = centroids(ord@ord, phy_pathway_clr)$uid, width = 0.2, fill = "transparent") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(shape = guide_legend(ncol = 2), colour = guide_legend(ncol = 2))

p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(p, filename = file.path(subdir, "pathway_ordination_diet.png"), width=6, height=6)

# Plot arrows
# Get loading arrows coordinaties
arrows <- arrow_coord(ord@ord, axes = c(1, 2))

# Get gene category
arrows$name <- as.character(phy_pathway_clr@tax_table[match(rownames(arrows), rownames(phy_pathway_clr@tax_table)), "path_name"])
arrows$category <- as.character(phy_pathway_clr@tax_table[match(rownames(arrows), rownames(phy_pathway_clr@tax_table)), "path_class"])
arrows$category <- str_remove(arrows$category, ".*; ")

arrows$to_plot <- (rownames(arrows) %in% head(rownames(arrows), nrow(arrows)))

# Save
write.csv(rownames_to_column(arrows, "gene"), file.path(subdir, "gene_ordination_arrows.txt"), quote = FALSE, row.names = FALSE)

# Keep only strongest associations
arrows_filt <- arrows %>% filter(to_plot) %>%
              select(contains(c("1", "2")), name, category)

# Group uncommon categories
common_categories <- table(arrows_filt$category) %>% sort(decreasing = TRUE) %>% head(8) %>% names

arrows_filt <- arrows_filt %>%
    mutate(category_grouped = factor(case_when(category %in% common_categories ~ category,
                                            TRUE ~ "Other"), levels = c(common_categories, "Other")))

# Set colours for categories using colour brewer
arrow_colours <- brewer.pal(n = length(unique(arrows_filt$category_grouped))-1, name = "Dark2")
names(arrow_colours) <- unique(arrows_filt$category_grouped)[-length(unique(arrows_filt$category_grouped))] # Remove "Other" from names
arrow_colours["Other"] <- "grey90" # Set "Other" to grey

p <- ggplot(data = arrows_filt) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = category_grouped), linewidth = 0.5, alpha = 0.5) +
  #geom_text(aes(x = RDA1, y = RDA2, label = name, hjust = ifelse(RDA1 < 0, 1, 0)), size = 1.5, vjust = 1) +
  scale_color_manual(values = arrow_colours, name = "Pathway class") +
  xlim(c(min(arrows_filt$RDA1)*2, max(arrows_filt$RDA1)*2)) +
  xlab("RDA1") + ylab("RDA2") +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.text = element_text(size = 7), legend.title = element_text(size = 10)) +
  guides(colour = guide_legend(ncol = 2))

ggsave(p, filename = file.path(subdir, "pathway_ordination_arrows.png"), width=6, height=5)

###################
#### PERMANOVA ####
###################

# Explanatory variables
order <- sample_data$Order
diet <- sample_data$diet.general
habitat <- sample_data$habitat.general
ruminant <- (sample_data$digestion == "Ruminant")
hypsodont <- grepl("hyps", sample_data$molar_category)
species <- sample_data$Species

#### GENE DATA ####
set.seed(123)

otu_table <- as.data.frame(phy_gene_f_clr@otu_table)

perm <- adonis2(otu_table ~ order + diet + habitat + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_clr_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(otu_table ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_clr_onlyspecies.csv"), row.names = TRUE, quote = TRUE)


#### PATHWAY DATA ####
set.seed(123)

otu_table <- as.data.frame(phy_gene_f_clr@otu_table)

perm <- adonis2(otu_table ~ order + diet + habitat + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_clr_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(otu_table ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_clr_onlyspecies.csv"), row.names = TRUE, quote = TRUE)
