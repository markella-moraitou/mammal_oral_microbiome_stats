##### ALPHA DIVERSITY #####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(microbiomeutilities)
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "alpha_diversity")) # subdirectory for the output of this script
phydir <- normalizePath(file.path("..", "..", "output", "community_analysis", "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Load all phyloseq objects in phydir
for (phy_file in list.files(phydir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(phydir, phy_file)))
}

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

#########################
#### ALPHA DIVERSITY ####
#########################

# Rarefy to min depth
phy_sp_rarefied <- rarefy_even_depth(subset_samples(phy_sp_f, sample_sums(phy_sp_f) > 1000), sample.size = 100, rngseed = 1)

# Create a named vector for the relabeling
order_labels <- c("Rodentia" = "Rod.", "Carnivora" = "Carniv.")

# Plot
p <- plot_richness(phy_sp_f, x="Species", measures=c("Observed")) +
  geom_boxplot(aes(fill=diet.general)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=diet_palette, name = "Species") +
  scale_x_discrete(labels = setNames(phy_sp_f@sam_data$Common.name, phy_sp_f@sam_data$Species)) +
  facet_grid(Order_grouped ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(Order_grouped = as_labeller(order_labels, default = label_value))) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Observed species richness") +
  coord_flip()

ggsave(file.path(subdir, "alpha_diversity.png"), p, width=8, height=10)

# Plot
p <- plot_richness(phy_sp_rarefied, x="Species", measures=c("Observed")) +
  geom_boxplot(aes(fill=diet.general)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=diet_palette, name = "Species") +
  scale_x_discrete(labels = setNames(phy_sp_f@sam_data$Common.name, phy_sp_f@sam_data$Species)) +
  facet_grid(Order_grouped ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(Order_grouped = as_labeller(order_labels, default = label_value))) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Observed species richness (after rarefaction)") +
  coord_flip()

ggsave(file.path(subdir, "alpha_diversity_rarefied.png"), p, width=8, height=10)

#### RAREFACTION CURVES ####
# get it separate per tax. order
subsamples <- seq(0, 1000, by=100)[-1]
plot_list <- list()
for (ord in unique(phy_sp_f@sam_data$Order_grouped)) { 
  phy_subset <- phy_sp_f %>% subset_samples(Order_grouped == ord)
  p <- plot_alpha_rcurve(phy_subset, index="observed", subsamples=subsamples,
                         type = "SD",
                         group="Species", label.color = "brown3",
                         label.size = 3, label.min = TRUE) +
    scale_fill_manual(values=species_palette, name = "Species") +
    scale_color_manual(values=species_palette, name = "Species") +
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3)))
  plot_list[[ord]] <- p
}

p <- plot_grid(plotlist = plot_list, ncol = 1, align = "v")

ggsave(file.path(subdir, "rarefaction_curves.png"), p, width=8, height=16)
