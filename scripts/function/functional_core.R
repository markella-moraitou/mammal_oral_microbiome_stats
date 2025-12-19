##### IDENTIFY FUNCTIONAL CORE MICROBIOME #####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(ggVennDiagram)
library(microViz)
library(vegan)
library(rphylopic)
library(wesanderson)
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "function", "functional_core")) # subdirectory for the output of this script
datadir <- normalizePath(file.path("..", "..", "output", "function", "data")) # Directory with data files
pathdir <- normalizePath(file.path("..", "..", "output", "function", "pathway_completeness")) # Directory with pathway completeness results
dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

# Get ordination functions
source(file.path("..", "ordination_functions.R"))

set.seed(123)

#######################
#####  LOAD INPUT #####
#######################

phy_pathway <- readRDS(file.path(pathdir, "phy_pathway.RDS"))
phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

#################
#### PROCESS ####
#################

# Keep host species with at least 4 samples
sample_counts <- data.frame(phy_pathway@sam_data) %>%
        group_by(Species) %>%
        summarise(n_samples = n())

phy_pathway <- phy_pathway %>% subset_samples(Species %in% sample_counts$Species[sample_counts$n_samples >= 4])
taxa_names(phy_pathway) <- phy_pathway@tax_table[, "path_name"]

## Calculate prevalence of each pathway in each species

species <- phy_pathway@sam_data$Species %>% levels
prevalence <- data.frame()

for (spe in species) {
  subset <- phy_pathway %>% subset_samples(Species == spe)
  temp <- prevalence(subset, sort = TRUE) %>% data.frame() %>% t
  rownames(temp) <- spe
  prevalence <- rbind(prevalence, temp)
}

prevalence <- t(prevalence) %>% data.frame %>% arrange(desc(rowSums(.))) %>%
                rownames_to_column("pathway")

write.csv(prevalence, file = file.path(subdir, "pathway_prevalence.csv"), quote = TRUE, row.names = FALSE)

#### Get core genes per species ####
## 75% prevalence

core_pathways <- data.frame(core_genes = character(), host_species = character())
core_plots <- list()

prev = 0.75
det = 10^-3

for (sp in species) {
  phy_sub <- phy_pathway %>% subset_samples(Species == sp) %>% transform("compositional")
  phy_sub <- phy_sub %>% subset_taxa(taxa_sums(phy_sub) > 0)
  # Plot core genes for this species
  p <- plot_core(phy_sub, prevalences = c(0.5, 0.6, 0.7, 0.8, 0.9, 1), detections = 10^(c(-8, -7, -6, -5, -4, -3)), 
          plot.type = "heatmap") +
          scale_fill_viridis_c(limits = c(0.5, 1)) +
          theme(legend.position = "none", axis.text.y = element_blank()) +
          xlab("") + ylab("") + labs(title = sp)
  if (sp == species[length(species)]) {p <- p + theme(legend.position = "right")} # Add legend to last plot
  core_plots[[sp]] <- p
  # Get core microbiota for this species
  core <- phy_sub %>% core_members(prevalence = prev, detection = det)
  # If there are no core genes, print warning and skip
  # otherwise, add to df
  if (length(core) == 0) {
    warning(paste("No core pathways found for species:", sp))
    next
  }
  df <- data.frame(core_pathways = core,
                   host_species = sp)
  core_pathways <- rbind(core_pathways, df)
}

# Save core pathways plots
p <- plot_grid(plotlist = core_plots, ncol = 4, align = "hv", axis = "tb")
ggsave(p, file = file.path(subdir, "core_thresholds.png"), width = 20, height = 20)

# Add extra information
# on the pathways
core_pathways$path_class <- as.vector(phy_pathway@tax_table[match(core_pathways$core_pathways, taxa_names(phy_pathway)), "path_class"])

# and on the samples
core_pathways$host_order <- phy_pathway@sam_data$Order[match(core_pathways$host_species, phy_pathway@sam_data$Species)]

core_pathways <- core_pathways %>% arrange(host_order, host_species, core_pathways) %>% 
        select(host_order, host_species, path_class, core_pathways)

# Calculate number of core pathways per species
core_pathways_summary <- 
        core_pathways %>% group_by(host_species, host_order) %>% 
        summarise(n_pathways = n(),
                  n_classes = n_distinct(path_class))

# Add number of samples
nsamples <- data.frame(phy_pathway@sam_data) %>% group_by(Species) %>% summarise(n_samples = n())

core_pathways_summary <- core_pathways_summary %>%
        left_join(nsamples, by = c("host_species" = "Species")) %>%
        arrange(host_order, host_species)

# Save tables
write.csv(core_pathways, file = file.path(subdir, "core_mb_per_host_species.csv"), quote = TRUE, row.names = FALSE)
write.csv(core_pathways_summary, file = file.path(subdir, "core_mb_per_host_species_summary.csv"), quote = TRUE, row.names = FALSE)

cor.test(core_pathways_summary$n_pathways, core_pathways_summary$n_samples, method = "spearman")

#### Core pathways per host order ####
# For orders with more than two species, get core genes per order
# defined as genes that are considered core in at least half of species in that order

core_pathways_per_order <- core_pathways %>% group_by(host_order) %>%
        # Calculate number of species in order
        mutate(n_species = n_distinct(host_species)) %>%
        filter(n_distinct(host_species) > 2) %>%
        # Calculate prevalence in order
        group_by(core_pathways, path_class, host_order) %>%
        summarise(prevalence_in_order = n_distinct(host_species)/first(n_species)) %>%
        filter(prevalence_in_order >= 2/3) %>%
        arrange(host_order, core_pathways) %>%
        # How many orders is this gene core in?
        group_by(core_pathways, path_class) %>%
        mutate(n_classes = n_distinct(path_class))

# Calculate number of core pathways per order
core_pathways_per_order_summary <- core_pathways_per_order %>% group_by(host_order) %>%
        summarise(n_pathways = n(),
                n_classes = n_distinct(path_class))

# Save tables
write.csv(core_pathways_per_order, file = file.path(subdir, "core_mb_per_host_order.csv"), quote = TRUE, row.names = FALSE)
write.csv(core_pathways_per_order_summary, file = file.path(subdir, "core_mb_per_host_order_summary.csv"), quote = FALSE, row.names = FALSE)

#### Plot Venn Diagrams ####

# Host order level
# Collect core genes in a list
core <- list()

for (ord in unique(core_pathways_per_order$host_order)) {
    core[[ord]] <-
        core_pathways_per_order %>% filter(host_order == ord) %>%
        pull(core_pathways) %>% unique
}

core_pathways_venn <- 
  ggVennDiagram(core, set_color = order_palette[names(core)], label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "grey40", name = "N. genera") +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "bottom")
ggsave(core_pathways_venn, file=file.path(subdir, "core_pathways_venn.png"),
       device = "png", width = 5, height = 5)

#### Core genera per host species heatmap ####
core_df <- core_pathways %>% mutate(is.core = 1) %>%
        ggplot(aes(x = host_species, y = str_trunc(core_pathways, 50), fill = is.core)) +
        facet_grid(cols = vars(host_order), rows = vars(str_remove(path_class, ".*; ")), scales = "free", space = "free") +
        geom_tile() +
        theme(strip.text.y = element_text(angle = 0,),
              strip.text.x = element_text(angle = 90),
              axis.text.x = element_text(hjust = 1, vjust = 0.5),
              legend.position = "none")

ggsave(core_df, file = file.path(subdir, "core_pathways_heatmap.png"),
       device = "png", width = 15, height = 20)

#### PCA with only core order pathways ####

phy_core <- phy_pathway %>% subset_taxa(taxa_names(phy_pathway) %in% unique(unlist(core))) %>%
    subset_samples(Order %in% core_pathways_per_order$host_order) 

phy_core_clr <- transform(phy_core, "clr")

# PCA ordination
ord <- ord_calc(phy_core_clr, method = "PCA") 

# Get shape scales for plotting
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

# Axes 1 & 2
p <- 
  ord_plot(ord, colour="Order", shape = "diet.general", alpha = 0.5) +
  custom_theme() +
  scale_shape_manual(values = diet_shape_scale, name = "Diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  geom_phylopic(data = centroids(ord@ord, phy_core_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_core_clr)$uid, width = 0.07, alpha = 0.8)

ggsave(p, file = file.path(subdir, "pca_core_order_1_2.png"), width = 5, height = 6)

# Axes 3 & 4
p <- 
  ord_plot(ord, colour="Order", shape = "diet.general", alpha = 0.5, axes = c(3, 4)) +
  custom_theme() +
  scale_color_manual(values=order_palette, name = "Order") +
  scale_shape_manual(values = diet_shape_scale, name = "Diet") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  geom_phylopic(data = centroids(ord@ord, phy_core_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_core_clr)$uid, width = 0.07, alpha = 0.8)

ggsave(p, file = file.path(subdir, "pca_core_order_3_4.png"), width = 5, height = 6)

#### Abundance heatmap ####

# Get more info on gene category
heat_data <- transform(phy_core, "compositional") %>% psmelt %>%
        select(OTU, Sample, Species, Order, Abundance) %>%
        # Turn abundances below detection limit to NA
        mutate(Abundance = ifelse(Abundance <= det, NA, Abundance)) %>%
        filter(OTU %in% taxa_names(phy_core)) %>%
        filter(Sample %in% sample_names(phy_core)) %>%
        mutate(pathway_type = case_when(OTU %in% Reduce(intersect, core) ~ "Mammalian core species",
                                      OTU %in% setdiff(core[["Primates"]], unlist(core[!names(core) %in% "Primates"])) ~ "Primates core only",
                                      OTU %in% setdiff(core[["Perissodactyla"]], unlist(core[!names(core) %in% "Perissodactyla"])) ~ "Perissodactyla core only",
                                      OTU %in% setdiff(core[["Carnivora"]], unlist(core[!names(core) %in% "Carnivora"])) ~ "Carnivora core only",
                                      TRUE ~ "Other")) %>%
        mutate(pathway_type = factor(pathway_type, levels = c("Mammalian core species", "Primates core only", "Perissodactyla core only", "Carnivora core only", "Other"))) %>%
        # Shorten order names for plotting
        mutate(Order = recode(Order, "Perissodactyla" = "Peris.", "Carnivora" = "Carn.")) %>%
        # Arrange sample by species name
        arrange(Order, Species, Sample) %>%
        mutate(Sample = factor(Sample, levels = unique(Sample)))

grad_palette <- colorRampPalette(c("#2D627B","#FFF7A4", "#E7C46E","#C24141"))
grad_palette <- grad_palette(10)

p <- ggplot(heat_data, aes(x = Sample, y = str_trunc(OTU, 50), fill = Abundance)) +
        geom_tile() +
        facet_grid(cols = vars(Order), scales = "free", space = "free",
                   rows = vars(pathway_type), switch = "y") +
        scale_fill_gradientn(colors = grad_palette, name = "Relative abundance", transform = "log10", breaks = c(0.01, 0.1)) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 8, angle = 0, hjust = 1, vjust = 0),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              legend.position = "top",
              legend.text = element_text(size = 8, angle=45),
              panel.spacing.y = unit(1.5, "lines"),
              strip.placement = "outside",
              strip.clip = "off",
              strip.text.y.left = element_text(size = 12, angle=0, vjust=1, hjust = 1, face = "bold"),
              strip.text.y = element_text(margin = margin(t=-15, r = -60)),
              strip.background.y = element_blank())

# Add species bar
species_bar <- data.frame(Sample = heat_data$Sample,
                          Species = phy_core@sam_data$Species[match(heat_data$Sample, rownames(phy_core@sam_data))],
                          Order = phy_core@sam_data$Order[match(heat_data$Sample, rownames(phy_core@sam_data))]) %>%
                left_join(phylopics, by = "Species") %>%
                # Keep uid for middle sample
                group_by(Species) %>%
                mutate(plot_phylopic = ifelse(row_number() == ceiling(n()/2), TRUE, FALSE)) %>%
                mutate(uid = ifelse(plot_phylopic, uid, NA)) %>%
                mutate(label = ifelse(plot_phylopic, Species, ""))

# Create a named vector for mapping Sample to label
sample_to_label <- setNames(species_bar$label, species_bar$Sample)

p_bar <- ggplot(species_bar, aes(x = Sample, y = 1)) +
        geom_tile(aes(fill = Species)) +
        scale_fill_manual(values = species_palette, name = "Species") +
        geom_phylopic(aes(x = Sample, y = 1, uuid = uid), height = 0.1, alpha = 1, fill = "white") +
        facet_grid(cols = vars(Order), scales = "free", space = "free") +
        scale_x_discrete(labels = sample_to_label) +
        theme(axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(hjust = 0, vjust = 0, size = 8),
              axis.title.x = element_blank(),
              legend.position = "none",
              panel.background = element_rect(fill = "white", color = "white"),
              strip.background.x = element_blank(),
              strip.text.x = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

# Combine heatmap and species bar
p <- plot_grid(p, p_bar, ncol = 1, rel_heights = c(1, 0.2), align = "v", axis = "lr")

ggsave(file.path(subdir, "heatmap_core_genes.png"), p, width=8, height=14)

#### Percentage of abundance consisting of core pathways ####
# This reflects how representative the above PCA is of the whole community

# Calculate the proportion of abundance that is made up of core pathways
phy_melt <- transform(phy_pathway, "compositional") %>% psmelt %>% filter(Abundance > 0)

core_pathways_per_order$is.core <- "Core in this order"

core_pathways_abund <-
        phy_melt %>%
        select(OTU, path_class, Sample, Species, Common.name, Order, Abundance) %>%
        # Keep only the same orders we use for the PCA
        filter(Order %in% core_pathways_per_order$host_order) %>%
        # Get core pathways per order
        left_join(core_pathways_per_order, by = c("OTU" = "core_pathways", "Order" = "host_order", "path_class" = "path_class"),
                  relationship = "many-to-one") %>%
        # Indicate if a gene is part of the mammalian core mb.
        # Alternatively, if it isnt core in all mammals or that specific order, check if it is core in another order
        mutate(is.core = case_when(OTU %in% Reduce(intersect, core_pathways) ~ "Mammalian core",
                                   !is.na(is.core) ~ is.core,
                                   OTU %in% core_pathways_per_order$core_pathways ~ "Core in other order")) %>%
        mutate(is.core = factor(is.core, levels = c("Core in other order", "Core in this order", "Mammalian core"))) %>%
        filter(!is.na(is.core)) %>%
        # Get relative abundance per sample
        group_by(Sample, Species, Common.name, Order, is.core) %>%
        summarise(Abundance = sum(Abundance)) %>%
        # Summarize by species
        group_by(Species, Common.name, Order, is.core) %>%
        summarise(mean_abundance = mean(Abundance),
                  sd = sd(Abundance, na.rm = TRUE)) %>%
        # Shorten "Perissodactyla" for plotting
        mutate(Order = recode(Order, "Perissodactyla" = "Periss.", "Carnivora" = "Carn."))

# Plot
p1 <- ggplot(core_pathways_abund, aes(y = Common.name, x = mean_abundance, fill = is.core)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("#FFDC7C", "#DA8A3D", "#DA4C3D"), name = "") +
        facet_grid(rows = vars(Order), scales = "free_y", space = "free_y") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75)) +
        labs(x = "Rel. abundance\nof core genes") +
        theme(legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank(),
              axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
              strip.background.y = element_blank(), strip.text.y = element_blank()) +
        guides(fill = guide_legend(reverse = TRUE))

#### Abundance of core genera ####
mamm_core_abund <- phy_melt %>%
        select(OTU, path_class, Sample, Species, Common.name, Order, Abundance) %>%
        # Keep only the same orders we use for the PCA
        filter(Order %in% core_pathways$host_order) %>%
        # Get pathways that are the mammalian core pathways
        filter(OTU %in% unique(unlist(core_pathways))) %>%
        # Sum by genus
        group_by(Sample, Species, Common.name, Order, path_class) %>%
        summarise(Abundance = sum(Abundance)) %>%
        # Calculate mean abundance per species
        group_by(Species, Common.name, Order, path_class) %>%
        summarise(mean_abundance = mean(Abundance)) %>%
        # Shorten "Perissodactyla" for plotting
        mutate(Order = recode(Order, "Perissodactyla" = "Periss.", "Carnivora" = "Carn."))

# Identify the 5 most abundant classes, group the rest as "Other"
top_classes <- mamm_core_abund %>% group_by(path_class) %>% summarise(mean_abundance = mean(mean_abundance)) %>%
        arrange(desc(mean_abundance)) %>% top_n(5, mean_abundance) %>% pull(path_class)

mamm_core_abund <- mamm_core_abund %>% mutate(path_class = ifelse(path_class %in% top_classes, path_class, "Other")) %>%
                  mutate(path_class = str_remove(path_class, ".*; ")) %>%
                  mutate(path_class = factor(path_class, levels = rev(c(str_remove(top_classes, ".*; "), "Other"))))

# Plot
colours <- rev(c("lightgrey", c(wes_palette("Darjeeling1", 5, type = "discrete"))))
names(colours) <- c(str_remove(top_classes, ".*; "), "Other")

p2 <- ggplot(mamm_core_abund, aes(y = Common.name, x = mean_abundance, fill = path_class)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = colours) +
        facet_grid(rows = vars(Order), scales = "free_y", space = "free_y") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75)) +
        labs(x = "Rel. abundance of\n core path classes") +
        theme(legend.position = "bottom", legend.direction = "vertical",
              legend.title = element_blank(),  
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        guides(fill = guide_legend(reverse = TRUE, ncol = 1))

p <- plot_grid(p1, p2, ncol = 2, align = "h", axis = "tb", rel_widths = c(1.5, 1))

ggsave(p, file = file.path(subdir, "core_pathway_abundance.png"), width = 8, height = 8)
