##### IDENTIFY FUNCTIONAL CORE MICROBIOME #####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
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

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

#################
#### PROCESS ####
#################

# Since we are dealing with presence absence, I will rarefy the same number of reads
phy_gene_rarefied <- phy_gene_f %>% rarefy_even_depth(rngseed=123, sample.size = 50000, replace = FALSE)

## Calculate prevalence of each genus in each species

# Prevalence per species
species <- phy_gene_rarefied@sam_data$Species %>% levels
prevalence <- data.frame()

for (spe in species) {
  subset <- phy_gene_rarefied %>% subset_samples(Species == spe)
  temp <- prevalence(subset, sort = TRUE) %>% data.frame() %>% t
  rownames(temp) <- spe
  prevalence <- rbind(prevalence, temp)
}

prevalence <- t(prevalence) %>% data.frame %>% arrange(desc(rowSums(.))) %>%
                rownames_to_column("gene")

write.csv(prevalence, file = file.path(subdir, "gene_prevalence.csv"), quote = FALSE, row.names = FALSE)

#### Get core genes per species ####
## 75% prevalence

core_genes <- data.frame(core_genes = character(), host_species = character())
core_plots <- list()

for (sp in species) {
  phy_sub <- phy_gene_rarefied %>% subset_samples(Species == sp) %>% transform("compositional")
  # Plot core genes for this species
  p <- plot_core(phy_sub, prevalences = c(0.5, 0.6, 0.7, 0.8, 0.9, 1), detections = 10^(c(-8, -7, -6, -5, -4, -3)), 
          plot.type = "heatmap") +
          scale_fill_viridis_c(limits = c(0.5, 1)) +
          theme(legend.position = "none", axis.text.y = element_blank()) +
          xlab("") + ylab("") + labs(title = sp)
  if (sp == species[length(species)]) {p <- p + theme(legend.position = "right")} # Add legend to last plot
  core_plots[[sp]] <- p
  # Get core microbiota for this species
  core <- phy_sub %>% core_members(prevalence = 0.75, detection = 0)
  # If there are no core genes, print warning and skip
  # otherwise, add to df
  if (length(core) == 0) {
    warning(paste("No core genes found for species:", sp))
    next
  }
  df <- data.frame(core_genes = core,
                   host_species = sp)
  core_genes <- rbind(core_genes, df)
}

# Save core genes plots
p <- plot_grid(plotlist = core_plots, ncol = 4, align = "hv", axis = "tb")
ggsave(p, file = file.path(subdir, "core_thresholds.png"), width = 20, height = 20)

# Add extra information
# on the genes
core_genes$gene_category <- as.vector(phy_gene_rarefied@tax_table[match(core_genes$core_genes, taxa_names(phy_gene_rarefied)), "category"])

# and on the samples
core_genes$host_order <- phy_gene_rarefied@sam_data$Order[match(core_genes$host_species, phy_gene_rarefied@sam_data$Species)]

core_genes <- core_genes %>% arrange(host_order, host_species, core_genes) %>% 
        select(host_order, host_species, gene_category, core_genes)

# Calculate number of core genes per species
core_genes_summary <- 
        core_genes %>% group_by(host_species, host_order) %>% 
        summarise(n_core_genes = n(),
                  n_categories = n_distinct(gene_category))

# Add number of samples
nsamples <- data.frame(phy_gene_rarefied@sam_data) %>% group_by(Species) %>% summarise(n_samples = n())

core_genes_summary <- core_genes_summary %>%
        left_join(nsamples, by = c("host_species" = "Species")) %>%
        arrange(host_order, host_species)

# Save tables
write.csv(core_genes, file = file.path(subdir, "core_mb_per_host_species.csv"), quote = FALSE, row.names = FALSE)
write.csv(core_genes_summary, file = file.path(subdir, "core_mb_per_host_species_summary.csv"), quote = FALSE, row.names = FALSE)

# For orders with more than two species, get core genes per order
# defined as genes that are considered core in at least half of species in that order

core_genes_per_order <- core_genes %>% group_by(host_order) %>%
        filter(n_distinct(host_species) > 2) %>%
        # Calculate number of species in order
        mutate(n_species = n_distinct(host_species)) %>%
        # Calculate prevalence in order
        group_by(core_genes, gene_category, host_order) %>%
        summarise(prevalence_in_order = n_distinct(host_species)/first(n_species)) %>%
        filter(prevalence_in_order >= 0.5) %>%
        arrange(host_order, core_genes) %>%
        # How many orders is this gene core in?
        group_by(core_genes, gene_category) %>%
        mutate(n_categories = n_distinct(gene_category))

# Calculate number of core genes per order
core_genes_per_order_summary <- core_genes_per_order %>% group_by(host_order) %>%
        summarise(n_core_genes = n(),
                n_categories = n_distinct(gene_category))

# Save tables
write.csv(core_genes_per_order, file = file.path(subdir, "core_mb_per_host_order.csv"), quote = FALSE, row.names = FALSE)
write.csv(core_genes_per_order_summary, file = file.path(subdir, "core_mb_per_host_order_summary.csv"), quote = FALSE, row.names = FALSE)

#### Plot Venn Diagrams ####

# Host order level
# Collect core genes in a list
core <- list()

for (ord in unique(core_genes_per_order$host_order)) {
    core[[ord]] <-
        core_genes_per_order %>% filter(host_order == ord) %>%
        pull(core_genes) %>% unique
}

core_genes_venn <- 
  ggVennDiagram(core, set_color = order_palette[names(core)], label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "grey40", name = "N. genera") +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "bottom")

ggsave(core_genes_venn, file=file.path(subdir, "core_genes_venn.png"),
       device = "png", width = 5, height = 5)

#### Core genera per host species heatmap ####
core_df <- core_genes %>% mutate(is.core = 1) %>%
        ggplot(aes(x = host_species, y = core_genes, fill = is.core)) +
        facet_grid(cols = vars(host_order), rows = vars(gene_category), scales = "free", space = "free") +
        geom_tile() +
        theme(strip.text.y = element_text(angle = 0,),
              strip.text.x = element_text(angle = 90),
              axis.text.y = element_blank(),
              axis.text.x = element_text(hjust = 1, vjust = 0.5),
              legend.position = "none")

ggsave(core_df, file = file.path(subdir, "core_genes_heatmap.png"),
       device = "png", width = 10, height = 16)

#### PCA with only core order genes ####

phy_core <- phy_gene_f %>% subset_taxa(taxa_names(phy_gene_f) %in% unique(unlist(core))) %>%
    subset_samples(Order %in% core_genes_per_order$host_order) 

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
        select(OTU, Sample, Order, Abundance) %>%
        # Log transform relative abundance
        mutate(Abundance = log10(Abundance + 0.0001)) %>%
        mutate(gene_type = case_when(OTU %in% Reduce(intersect, core) ~ "Mammalian core species",
                                      OTU %in% setdiff(core[["Primates"]], unlist(core[!names(core) %in% "Primates"])) ~ "Primates core only",
                                      OTU %in% setdiff(core[["Perissodactyla"]], unlist(core[!names(core) %in% "Perissodactyla"])) ~ "Perissodactyla core only",
                                      OTU %in% setdiff(core[["Carnivora"]], unlist(core[!names(core) %in% "Carnivora"])) ~ "Carnivora core only",
                                      TRUE ~ "Other")) %>%
        mutate(gene_type = factor(gene_type, levels = c("Mammalian core species", "Primates core only", "Perissodactyla core only", "Carnivora core only", "Other"))) %>%
        # Shorten order names for plotting
        mutate(Order = recode(Order, "Perissodactyla" = "Peris.", "Carnivora" = "Carn."))

grad_palette <- colorRampPalette(c("#2D627B","#FFF7A4", "#E7C46E","#C24141"))
grad_palette <- grad_palette(10)

p <- ggplot(heat_data, aes(x = Sample, y = OTU, fill = Abundance)) +
        geom_tile() +
        facet_grid(cols = vars(Order), scales = "free", space = "free",
                   rows = vars(gene_type), switch = "y") +
        scale_fill_gradientn(colors = grad_palette, name = "Rel.abund.\n(log10)") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "bottom",
              panel.spacing.y = unit(1.5, "lines"),
              strip.placement = "outside",
              strip.clip = "off",
              strip.text.y.left = element_text(angle=0, vjust=1),
              strip.text.y = element_text(margin = margin(t=-15, r = -160),
                                          size = 15, hjust = 1),
              strip.background.y = element_blank())

ggsave(file.path(subdir, "heatmap_core_genes.png"), p, width=8, height=11)

#### Percentage of abundance consisting of core genes ####
# This reflects how representative the above PCA is of the whole community

# Calculate the proportion of abundance that is made up of core genes
phy_melt <- transform(phy_gene_rarefied, "compositional") %>% psmelt %>% filter(Abundance > 0)

core_genes_per_order$is.core <- "Core in this order"

core_genes_abund <-
        phy_melt %>%
        select(OTU, category, Sample, Species, Common.name, Order, Abundance) %>%
        # Keep only the same orders we use for the PCA
        filter(Order %in% core_genes_per_order$host_order) %>%
        # Get core genes per order
        left_join(core_genes_per_order, by = c("OTU" = "core_genes", "Order" = "host_order", "category" = "gene_category"),
                  relationship = "many-to-one") %>%
        # Indicate if a gene is part of the mammalian core mb.
        # Alternatively, if it isnt core in all mammals or that specific order, check if it is core in another order
        mutate(is.core = case_when(OTU %in% Reduce(intersect, core_genes) ~ "Mammalian core",
                                   !is.na(is.core) ~ is.core,
                                   OTU %in% core_genes_per_order$core_genes ~ "Core in other order")) %>%
        mutate(is.core = factor(is.core, levels = c("Core in other order", "Core in this order", "Mammalian core"))) %>%
        filter(!is.na(is.core)) %>%
        # Get relative abunddance of host genes per sample
        group_by(Sample, Species, Common.name, Order, is.core) %>%
        summarise(Abundance = sum(Abundance)) %>%
        # Summarize by species
        group_by(Species, Common.name, Order, is.core) %>%
        summarise(mean_abundance = mean(Abundance),
                  sd = sd(Abundance, na.rm = TRUE)) %>%
        # Shorten "Perissodactyla" for plotting
        mutate(Order = recode(Order, "Perissodactyla" = "Periss.", "Carnivora" = "Carn."))

# Plot
p1 <- ggplot(core_genes_abund, aes(y = Common.name, x = mean_abundance, fill = is.core)) +
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
        select(OTU, category, Sample, Species, Common.name, Order, Abundance) %>%
        # Keep only the same orders we use for the PCA
        filter(Order %in% core_genes_per_order$host_order) %>%
        # Get genes that are the mammalian core genes
        filter(OTU %in% unique(unlist(core_genes))) %>%
        # Sum by genus
        group_by(Sample, Species, Common.name, Order, category) %>%
        summarise(Abundance = sum(Abundance)) %>%
        # Calculate mean abundance per species
        group_by(Species, Common.name, Order, category) %>%
        summarise(mean_abundance = mean(Abundance)) %>%
        # Shorten "Perissodactyla" for plotting
        mutate(Order = recode(Order, "Perissodactyla" = "Periss.", "Carnivora" = "Carn."))

# Identify the 5 most abundant categories, group the rest as "Other"
top_categories <- mamm_core_abund %>% group_by(category) %>% summarise(mean_abundance = mean(mean_abundance)) %>%
        arrange(desc(mean_abundance)) %>% top_n(5, mean_abundance) %>% pull(category)

mamm_core_abund <- mamm_core_abund %>% mutate(category = ifelse(category %in% top_categories, category, "Other")) %>%
                  mutate(category = factor(category, levels = rev(c(top_categories, "Other"))))

# Plot
colours <- rev(c("transparent", c(wes_palette("Darjeeling1", 5, type = "discrete"))))
names(colours) <- c(top_categories, "Other")

p2 <- ggplot(mamm_core_abund, aes(y = Common.name, x = mean_abundance, fill = category)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = colours) +
        facet_grid(rows = vars(Order), scales = "free_y", space = "free_y") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75)) +
        labs(x = "Rel. abundance of\n core categories") +
        theme(legend.position = "bottom", legend.direction = "vertical",
              legend.title = element_blank(),  
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        guides(fill = guide_legend(reverse = TRUE, ncol = 2))

p <- plot_grid(p1, p2, ncol = 2, align = "h", axis = "tb", rel_widths = c(1.5, 1))

ggsave(p, file = file.path(subdir, "core_genes_abundance.png"), width = 8, height = 8)
