#### RDA analysis ####
## Use redundancy analysis to investigate how taxonomic order, diet and habitat affect microbiome composition

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(microbiome)
library(reshape2)
library(ggplot2)
library(vegan)
library(microViz)
library(rphylopic)
library(ggnewscale)
library(RColorBrewer)
library(colorspace)
library(ggExtra)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Directory with output from the community analysis
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "RDA_analysis")) # subdirectory for the output of this script
phydir <- normalizePath(file.path("..", "..", "output", "community_analysis", "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

source(file.path("..","ordination_functions.R"))

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

#######################
#####  LOAD INPUT #####
#######################

# Load all phyloseq objects in phydir
for (phy_file in list.files(phydir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(phydir, phy_file)))
}

########################
#### RDA ORDINATION ####
########################

#### PREP SCALE FOR PLOTTING ####
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

uniq_species <- unique(subset_samples(phy_sp_f_clr, Order %in% c("Primates", "Carnivora", "Artiodactyla") | habitat.general == "Aquatic")@sam_data$Common.name)
species_shape_scale <- c(1:25, 35:35+26-length(uniq_species))
names(species_shape_scale) <- uniq_species

order_shape_scale <- c("Carnivora" = 4, "Primates" = 19, "Artiodactyla" = 5, "Perissodactyla" = 2, "Rodentia" = 1, "Rest" = 12)

# Recode order and habitat as TRUE and FALSE
# Also scale protein, fiber and carbohydrate content
phy_sp_f_clr <- phy_sp_f_clr %>%
        ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  Rodentia = (Order == "Rodentia"),
                  ruminant = (digestion == "Ruminant"),
                  marine = (habitat.general == "Marine"))

# Species traits to use as constraints
species_traits <- c("Artiodactyla", "Perissodactyla", "Primates", "Rodentia",
                    "ruminant", "marine", "cf", "cp")

# Ordinate using all data
ord <- ord_calc(phy_sp_f_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:10))

ggsave(file.path(subdir, "rda_screeplot_all.png"), p, width=8, height=6)

#### SAMPLE PLOTS ####
# Color by diet
# Axes 1,2 -- colour by diet
p <- ord_plot(ord, colour="diet.general", shape="Order_grouped", alpha = 0.5, size = 1) +
  custom_theme() +
  scale_shape_manual(values=order_shape_scale, name = "Order") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_f_clr), aes(colour = diet.general), uuid = centroids(ord@ord, phy_sp_f_clr)$uid, width = 0.3, alpha = 0.6) +
  theme(legend.position = "bottom", legend.direction = "vertical")

p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(file.path(subdir, "RDA_all_1_2_diet.png"), p, width=5, height=7)

# Axes 2,3 -- colour by diet
p <- ord_plot(ord, colour="diet.general", shape="Order_grouped", alpha = 0.5, size = 1, axes = c(2, 3)) +
  custom_theme() +
  scale_shape_manual(values=order_shape_scale[names(order_shape_scale) %in% phy_sp_f@sam_data$Order_grouped], name = "Order") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_f_clr), aes(colour = diet.general), uuid = centroids(ord@ord, phy_sp_f_clr)$uid, width = 0.3, alpha = 0.6) +
  theme(legend.position = "bottom", legend.direction = "vertical")

p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(file.path(subdir, "RDA_all_2_3_diet.png"), p, width=5, height=7)

# Axes 1,2 -- colour by order
p <- ord_plot(ord, colour="Order_grouped", shape="diet.general", alpha = 0.5, size = 1) +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette[names(order_palette) %in% phy_sp_f@sam_data$Order_grouped], name = "Order") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_f_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_sp_f_clr)$uid, width = 0.3, alpha = 0.6) +
  theme(legend.position = "bottom", legend.direction = "vertical")

p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(file.path(subdir, "RDA_all_1_2_order.png"), p, width=5, height=7)

# Axes 2,3 -- colour by order
p <- ord_plot(ord, colour="Order_grouped", shape="diet.general", alpha = 0.5, axes = c(2, 3)) +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette[names(order_palette) %in% phy_sp_f@sam_data$Order_grouped], name = "Order") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_f_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_sp_f_clr)$uid, width = 0.3, alpha = 0.6) +
  theme(legend.position = "bottom", legend.direction = "vertical")

p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)

ggsave(file.path(subdir, "RDA_all_2_3_order.png"), p, width=5, height=7)

#### TAXA PLOTS ####

# Get taxa scores and add taxonomic info
taxa_rda <- data.frame(scores(ord@ord, display="species", choices=1:3)) %>%
            cbind(tax_table(phy_sp_f_clr)[, c("superkingdom", "phylum", "genus", "species")]) %>%
            rownames_to_column("OTU") %>%
            # Get grouped phylum
            mutate(phylum_grouped = factor(case_when(phylum %in% names(phylum_palette) ~ phylum,
                                      superkingdom == "Bacteria" ~ "Other Bacteria",
                                      superkingdom == "Archaea" ~ "Other Archaea"), levels = names(phylum_palette)))

# Add mean relative abundance log10-transformed with pseudocount
taxa_rda$mean_abund <- rowMeans(phy_sp_f_clr@otu_table)[taxa_rda$OTU]

# Identify 15 genera with most abundance, group remaining genera to "Other"
top_genera <- taxa_rda %>% group_by(genus) %>% summarise(total = sum(mean_abund, na.rm = TRUE)) %>% top_n(15, total) %>% pull(genus)

taxa_rda <- taxa_rda %>% mutate(genus = ifelse(genus %in% top_genera, genus, "Other")) %>%
            mutate(genus = factor(genus, levels = top_genera)) %>%
            filter(genus != "Other") 

genus_centroids <- taxa_rda %>% select(OTU, RDA1, RDA2, RDA3, genus, phylum_grouped) %>% ungroup %>%
                group_by(genus, phylum_grouped) %>%
                summarise(RDA1 = mean(RDA1),
                          RDA2 = mean(RDA2),
                          RDA3 = mean(RDA3)) %>%
                filter(!grepl("Other", genus))

# Get colours
genus_palette <- expand_palette(select(genus_centroids, c(phylum_grouped, genus)), phylum_palette)

# Axes 1,2
set.seed(245)
p <- ggplot(taxa_rda, aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = genus)) +
  geom_segment(linewidth = 0.2, alpha = 0.8) +
  scale_color_manual(values = genus_palette, name = "Microbial genus") +
  geom_label(data = genus_centroids, aes(label = genus, x = RDA1, y = RDA2),
            size = 2, alpha = 0.8, position = position_jitter(height = 0.02)) +
  scale_size_continuous(range = c(0.01, 2), name = "Mean CLR-abundance") +
  custom_theme() + xlab("RDA1 scores") + ylab("RDA2 scores") +
  theme(legend.position = "none")

ggsave(file.path(subdir, "RDA_all_1_2_taxa.png"), p, width=5, height=5)

# Axes 2, 3
p <- ggplot(taxa_rda, aes(x = 0, y = 0, xend = RDA2, yend = RDA3, colour = genus)) +
  geom_segment(linewidth = 0.2, alpha = 0.8) +
  scale_color_manual(values = genus_palette, name = "Microbial genus") +
  geom_label(data = genus_centroids, aes(label = genus, x = RDA2, y = RDA3),
            size = 2, alpha = 0.8, position = position_jitter(height = 0.02)) +
  scale_size_continuous(range = c(0.01, 2), name = "Mean CLR-abundance") +
  custom_theme() + xlab("RDA1 scores") + ylab("RDA2 scores") +
  theme(legend.position = "none")

ggsave(file.path(subdir, "RDA_all_2_3_taxa.png"), p, width=5, height=5)

#### Test
anova(ord@ord, by = "margin", perm = 500)
anova(ord@ord, by = "axis", perm = 500)

#### VIOLIN PLOTS ####

rda_res <- data.frame(phy_deep_clr@sam_data) %>% select(new_name, Species, Order, Order_grouped, diet.general) %>%
           left_join(rownames_to_column(data.frame(scores(ord@ord, display="sites", choices=1:3))), by=c("new_name"="rowname")) %>%
           pivot_longer(cols = c(RDA1, RDA2, RDA3), names_to = "Axis", values_to = "Value")

p <- ggplot(rda_res, aes(y = Order, x = Value, fill = diet.general, colour = diet.general)) +
      geom_boxplot() +
      scale_fill_manual(values = diet_palette, name = "Diet category") +
      scale_color_manual(values = darken(diet_palette, amount = 0.3), name = "Diet category") +
      facet_grid(cols = vars(Axis), scales = "free") +
      theme(legend.position = "none", axis.title = element_blank())

ggsave(file.path(subdir, "RDA_axis_comparison.png"), p, width=8, height=8)
