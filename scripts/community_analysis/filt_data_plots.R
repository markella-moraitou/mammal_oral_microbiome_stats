##### EXPLORE FILTERED & DECONTAMINATED DATA #####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(microbiome)
library(microViz)
library(rphylopic)
library(RColorBrewer)
library(wesanderson)
library(ggplot2)
library(vegan)
library(microbiomeutilities)
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "filtered_data")) # subdirectory for the output of this script
phydir <- normalizePath(file.path("..", "..", "output", "community_analysis", "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

# Get ordination functions
source(file.path("ordination_functions.R"))

#######################
#####  LOAD INPUT #####
#######################

# Load all phyloseq objects in phydir
for (phy_file in list.files(phydir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(phydir, phy_file)))
}

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

#########################
#####  COMPOSITION  #####
#########################

#### PHYLUM COMPOSITION ####
phy_phylum <- tax_glom(phy_sp_f, taxrank = "phylum")
taxa_names(phy_phylum) <- as.vector(phy_phylum@tax_table[,"phylum"])

# Get only 10 most common phyla and turn rest to other
phylum_grouped = data.frame(abundance = taxa_sums(phy_phylum), superkingdom = phy_phylum@tax_table[,"superkingdom"]) %>% rownames_to_column("phylum") %>% arrange(-abundance) %>%
  mutate(phylum_grouped = ifelse(row_number() > 5, paste("Other", superkingdom, sep = " "), phylum))

# Change phylum names to grouped names and reaggregate
phy_phylum@tax_table[,"phylum"] <- phylum_grouped$phylum_grouped[match(phy_phylum@tax_table[,"phylum"], phylum_grouped$phylum)]

# Aggregate again
phy_phylum <- tax_glom(phy_phylum, taxrank = "phylum")
taxa_names(phy_phylum) <- phy_phylum@tax_table[,"phylum"] 

# Melt and turn phyla into a factor and reorder
phy_phylum_melt <- psmelt(transform(phy_phylum, "compositional"))
phy_phylum_melt$OTU <- factor(phy_phylum_melt$OTU , levels=names(phylum_palette))

# Order by fusobacteriota
sample_levels <- select(phy_phylum_melt, c(Sample, Species, Order_grouped, OTU, Abundance)) %>% filter(OTU == "Pseudomonadota") %>%
  arrange(Order_grouped, Species, desc(Abundance)) %>% pull(Sample)

phy_phylum_melt$Sample <- factor(phy_phylum_melt$Sample, levels=sample_levels)

p = ggplot(data = phy_phylum_melt, aes(x = Abundance, y = Sample, fill = OTU)) +
  geom_bar(stat = "identity") +
  facet_grid(Order_grouped~., space = "free_y", scales = "free_y", switch = "y") +
  scale_fill_manual(values=phylum_palette, name = "Phylum") +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title.position = "top", legend.key.spacing.x = unit(0.5, "cm"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  xlab("")

## Get sample_metadata
species_bar <- phy_phylum_melt %>% select(Sample, Species, Common.name, Order_grouped) %>% unique %>%
               arrange(Sample) %>%
               # Choose one label per species (in the middle of the species group)
               group_by(Species) %>% mutate(order = row_number()) %>%
               mutate(label = ifelse(order == floor(mean(order)), Common.name, "")) %>%
               mutate(label = ifelse(is.na(label), as.character(Species), label))

p_bar <-
  ggplot(data = species_bar, aes(y = Sample, x=1, fill = Species, group = Common.name)) +
  geom_tile() + scale_fill_manual(values = species_palette, name = "") +
  facet_grid(rows = vars(Order_grouped), scales = "free", space = "free", switch = "y") +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_blank()) + scale_y_discrete(position = "right", label = setNames(species_bar$label, species_bar$Sample)) +
    xlab("") + ylab("")

ggsave(filename = file.path(subdir, "phy_sp_f_composition.png"), device="png", width=8, height=12,
       plot_grid(p, p_bar, ncol = 2, align = "h", rel_widths = c(4, 2)))

####################################
#### PCA ORDINATION ON CLR DATA ####
####################################

#### Get shape scales for plotting ####
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

uniq_species <- unique(subset_samples(phy_sp_f_clr, Order %in% c("Primates", "Carnivora", "Artiodactyla") | habitat.general == "Aquatic")@sam_data$Common.name)
species_shape_scale <- c(1:25, 35:35+26-length(uniq_species))
names(species_shape_scale) <- uniq_species

order_shape_scale <- c("Carnivora" = 4, "Primates" = 19, "Artiodactyla" = 5, "Perissodactyla" = 2, "Rodentia" = 1, "Rest" = 12)

#### All data ####
ord <- ord_calc(phy_sp_f_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_all.png"), p, width=8, height=6)

# Color by order
p <- ord_plot(ord, colour="Order_grouped", shape="diet.general", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_f_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_sp_f_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(1, 2), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_order_1_2.png"), p, width=8, height=10)

# Axes 3 & 4
p <- ord_plot(ord, colour="Order_grouped", shape="diet.general", axes = c(3, 4)) +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_f_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_sp_f_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(3, 4), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_order_3_4.png"), p, width=8, height=10)

# Color by diet
p <- ord_plot(ord, colour="diet.general", shape="Order_grouped") +
  custom_theme() +
  scale_shape_manual(values=order_shape_scale, name = "Order") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_f_clr), aes(colour = diet.general), uuid = centroids(ord@ord, phy_sp_f_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(1, 2), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_diet_1_2.png"), p, width=8, height=10)

p <- ord_plot(ord, colour="diet.general", shape="Order_grouped", axes = 3:4) +
  custom_theme() +
  scale_shape_manual(values=order_shape_scale, name = "Order") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_f_clr), aes(colour = diet.general), uuid = centroids(ord@ord, phy_sp_f_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(3, 4), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_diet_3_4.png"), p, width=8, height=10)

#### Deep dataset ####
ord <- ord_calc(phy_deep_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_deep.png"), p, width=8, height=6)

p <- ord_plot(ord, colour="Order", shape="diet.general") +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_deep_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_deep_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(1, 2), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_subset_deep_1_2.png"), p, width=8, height=10)

# Axes 3 & 4
p <- ord_plot(ord, colour="Order", shape="diet.general", axes = c(3, 4)) +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_deep_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_deep_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(3, 4), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_subset_deep_3_4.png"), p, width=8, height=10)

#### Artiodactyla ####
ord <- ord_calc(phy_artio_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_artio.png"), p, width=8, height=6)

p <- ord_plot(ord, colour="diet.general", shape="Common.name") +
  custom_theme() +
  scale_shape_manual(values = species_shape_scale, name = "Species") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_artio_clr), aes(colour = diet.general), uuid = centroids(ord@ord, phy_artio_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(1, 2), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_subset_artio_1_2.png"), p, width=8, height=10)

#### Carnivora ####
ord <- ord_calc(phy_carni_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_carni.png"), p, width=8, height=6)

p <- ord_plot(ord, colour="habitat.general", shape="Common.name") +
  custom_theme() +
  scale_shape_manual(values = species_shape_scale) +
  scale_color_manual(values=habitat_palette, name = "Habitat") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_carni_clr), aes(colour = habitat.general), uuid = centroids(ord@ord, phy_carni_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(1, 2), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_subset_carni_1_2.png"), p, width=8, height=10)

#### Primates ####
ord <- ord_calc(phy_prim_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_prim.png"), p, width=8, height=6)

p <- ord_plot(ord, colour="diet.general", shape="Common.name") +
  custom_theme() +
  scale_shape_manual(values = species_shape_scale, name = "Species") +
  scale_color_manual(values = diet_palette, name = "Estimated diet") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_prim_clr), aes(colour = diet.general), uuid = centroids(ord@ord, phy_prim_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(1, 2), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_subset_prim_1_2.png"), p, width=8, height=10)

#### Habitat ####
ord <- ord_calc(phy_habitat_clr, method = "PCA")

p <- ord_plot(ord, colour="habitat.general", shape = "Common.name") +
  custom_theme() +
  scale_shape_manual(values = species_shape_scale, name = "Species") +
  scale_color_manual(values = habitat_palette, name = "Habitat") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_habitat_clr), aes(colour = habitat.general), uuid = centroids(ord@ord, phy_habitat_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(1, 2), top_taxa = 30)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_subset_habitat_1_2.png"), p, width=8, height=10)

#### Print phylopics legend ####
phylopics$Common.name <- phy_sp_f@sam_data$Common.name[match(phylopics$Species, as.vector(phy_sp_f@sam_data$Species))]

p <- ggplot(phylopics[!is.na(phylopics$Common.name),], aes(y = Common.name, x = 1)) +
    geom_phylopic(aes(uuid = uid), height = 0.8) +
    theme_void() + theme(axis.title = element_blank(), axis.text.x = element_blank(),
                         plot.background = element_rect(fill = "white", color = "white"),
                         axis.text.y = element_text(size = 5, hjust = 1))

ggsave(file.path(subdir, "phylopics_legend.png"), p, width=1.5, height=4)

######################################
#### PCA ORDINATION ON PHILR DATA ####
######################################

#### All data ####
ord <- ord_calc(phy_sp_philr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "screeplot_all_philr.png"), p, width=8, height=6)

# Color by order
p <- ord_plot(ord, colour="Order_grouped", shape="diet.general", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_philr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_sp_philr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(1, 2), top_taxa = 20)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_order_philr_1_2.png"), p, width=8, height=10)

# Axes 3 & 4
p <- ord_plot(ord, colour="Order_grouped", shape="diet.general", axes = c(3, 4)) +
  custom_theme() +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_philr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_sp_philr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(3, 4), top_taxa = 20)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_order_philr_3_4.png"), p, width=8, height=10)

# Color by diet
p <- ord_plot(ord, colour="diet.general", shape="Order_grouped") +
  custom_theme() +
  scale_shape_manual(values=order_shape_scale, name = "Order") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_f_clr), aes(colour = diet.general), uuid = centroids(ord@ord, phy_sp_f_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(1, 2), top_taxa = 20)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_diet_philr_1_2.png"), p, width=8, height=10)

p <- ord_plot(ord, colour="diet.general", shape="Order_grouped", axes = 3:4) +
  custom_theme() +
  scale_shape_manual(values=order_shape_scale, name = "Order") +
  scale_color_manual(values=diet_palette, name = "Estimated diet") +
  theme(legend.position = "left") +
  geom_phylopic(data = centroids(ord@ord, phy_sp_f_clr), aes(colour = diet.general), uuid = centroids(ord@ord, phy_sp_f_clr)$uid, width = 0.2, alpha = 0.8)

# Plot loadings alongside
p_l <- loadings_plot(ord@ord, axes = c(3, 4), top_taxa = 20)
p <- plot_grid(p + theme(legend.position = "right"),
               p_l + theme(legend.position = "right"),
               nrow = 2, rel_heights = c(4, 3))

ggsave(file.path(subdir, "PCA_diet_philr_3_4.png"), p, width=8, height=10)

#################
#### HEATMAP ####
#################

grad_palette <- colorRampPalette(c("#2D627B","#FFF7A4", "#E7C46E","#C24141"))
grad_palette <- grad_palette(10)

png(filename = file.path(subdir, "heatmap_all_taxa.png"), width=16, height=20, units="in", res=300)
plot_taxa_heatmap(phy_sp_f, subset.top=ntaxa(phy_sp_f), transformation="clr",
                  VariableA=c("Order", "diet.general", "unmapped_count"),
                  annotation_colors = list("Order" = order_palette, "diet.general" = diet_palette),
                  show_rownames = FALSE,
                  show_colnames = FALSE, heatcolors = grad_palette)$plot
dev.off()
