##### EXPLORATORY PLOTS #####

#### Plots exploratory plots such as ordinations and heatmaps

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
#library(stringr)
#library(tibble)
library(phyloseq)
library(microViz)
#library(cowplot)
library(microbiome)
library(rphylopic)
#library(vegan)
library(ggplot2)
library(ggnewscale)
#library(RColorBrewer)

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

# Phyloseq object
phy_gene_clr <- readRDS(file.path(outdir, "phy_gene_clr.RDS"))
phy_function <- readRDS(file.path(outdir, "phy_function.RDS"))

low_content_samples <- read.csv(file.path(outdir, "low_content_samples.txt"), header = TRUE)

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

###################
#### PCA PLOTS ####
###################

#### GENE ABUNDANCE (CLR) ####

# Filter out low content samples
phy_gene_filt <- phy_gene_clr %>% subset_samples(!sample_names(phy_gene_clr) %in% low_content_samples$x)

ord <- ord_calc(phy_gene_filt, method = "PCA")

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

# Group uncommon headers
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
  geom_phylopic(data = centroids(ord@ord, phy_gene_filt), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_gene_filt)$uid, width = 0.2, fill = "transparent") +
  new_scale_colour() +
  geom_segment(data = arrows_filt, aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = module_grouped), linewidth = 0.8, alpha = 0.5) +
  scale_color_manual(values = arrow_colours, name = "Header") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(ncol = 1, size = 1, byrow = TRUE))

ggsave(p, filename = file.path(subdir, "gene_ordination.png"), width=8, height=10)

#### FUNCTION PRESENCE-ABSENCE ####

# Filter out low content samples
phy_function <- phy_function %>% subset_samples(sample_sums(phy_function) > 0 | Order == "Control/blank")
phy_function <- phy_function %>% subset_samples(!(sample_names(phy_function) %in% low_content_samples$x))

ord <- ord_calc(phy_function, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_functions.png"), p, width=8, height=6)

# Color by order
#### Get shape scales for plotting ####
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

# Get loading arrows coordinaties
arrows <- arrow_coord(ord@ord, phy_function)

# Get gene module
arrows$category <- as.character(phy_function@tax_table[match(rownames(arrows),  rownames(phy_function@tax_table)), "category"])

arrows$to_plot <- (arrows$r > quantile(arrows$r, 0.50))

# Save
write.csv(arrows, file.path(subdir, "function_ordination_arrows.txt"), quote = TRUE)

# Keep only strongest associations
arrows_filt <- arrows %>% filter(to_plot) %>%
              select(contains(c("1", "2", "3", "4")), category)

# Group uncommon function categories
common_categories <- table(arrows_filt$category) %>% sort(decreasing = TRUE) %>% head(4) %>% names

arrows_filt <- arrows_filt %>%
    mutate(category_grouped = factor(case_when(category %in% common_categories ~ category,
                                            TRUE ~ "Other"), levels = c(common_categories, "Other")))

arrow_colours <- c("Methanogenesis and methanotrophy" = "#FF5F17",
                   "Nitrogen metabolism" = "#FFEB17",
                   "Sulfur metabolism" = "#13D38A",
                   "CAZy" = "#7123D4",
                   "Other" = "grey80")

p <- ord_plot(ord, colour="Order_grouped", shape="diet.general", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  geom_phylopic(data = centroids(ord@ord, phy_function), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_function)$uid, width = 0.05, fill = "transparent") +
  new_scale_colour() +
  geom_segment(data = arrows_filt, aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = category_grouped), linewidth = 0.8, alpha = 0.5) +
  scale_color_manual(values = arrow_colours, name = "Header") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(ncol = 1, size = 1, byrow = TRUE))

ggsave(p, filename = file.path(subdir, "function_ordination.png"), width=8, height=10)

###############################
#### PLOT FUNCTION HEATMAP ####
###############################

func_completeness <- psmelt(phy_function) %>% rename("Completeness" = "Abundance", "Function" = "OTU") %>% select(Completeness, Function, category, Sample, Species)

p <- ggplot(func_completeness, aes(x = Sample, y = Function, fill = Completeness)) +
    geom_tile() +
    scale_fill_gradient(high = "#146AB7", low = "#FFE5C0") +
    facet_grid(cols = vars(Species), rows = vars(category), scales = "free", space = "free") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 10),
         strip.text.y = element_text(angle = 0), strip.text.x = element_text(angle = 90))

ggsave(p, filename = file.path(subdir, "function_heatmap.png"), width=25, height=16)
