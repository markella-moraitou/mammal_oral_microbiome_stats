##### IDENTIFY CORE MICROBIOMES #####

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
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "core_microbiomes")) # subdirectory for the output of this script
phydir <- normalizePath(file.path("..", "..", "output", "community_analysis", "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

source(file.path("..", "ordination_functions.R"))

set.seed(123)

#######################
#####  LOAD INPUT #####
#######################

# Load all phyloseq objects in phydir
for (phy_file in list.files(phydir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(phydir, phy_file)))
}

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

#################
#### PROCESS ####
#################

# Aggregate to genus level
phy_gen <- phy_sp_f %>%
    tax_glom(taxrank = "genus")

taxa_names(phy_gen) <- make.unique(phy_gen@tax_table[,"genus"])

# Keep host species with at least 4 samples
sample_counts <- data.frame(phy_gen@sam_data) %>%
        group_by(Species) %>%
        summarise(n_samples = n())

phy_gen <- phy_gen %>% subset_samples(Species %in% sample_counts$Species[sample_counts$n_samples >= 4])

## Calculate prevalence of each genus in each species

species <- phy_gen@sam_data$Species %>% levels
prevalence <- data.frame()

prev = 0.75
det = 10^-3

for (spe in species) {
  subset <- phy_gen %>% subset_samples(Species == spe) %>% transform("compositional")
  temp <- prevalence(subset, detection = det, sort = TRUE) %>% data.frame() %>% t
  rownames(temp) <- spe
  prevalence <- rbind(prevalence, temp)
}

prevalence <- t(prevalence) %>% data.frame %>% arrange(desc(rowSums(.))) %>%
                rownames_to_column("taxon")

write.csv(prevalence, file = file.path(subdir, "taxa_prevalence.csv"), quote = FALSE, row.names = FALSE)

#### Get core taxa per species ####
## 75% prevalence

## Try difference prevalence thresholds

core_genera <- data.frame(core_taxa = character(), host_species = character())
core_plots <- list()

for (sp in species) {
  phy_sub <- phy_gen %>% subset_samples(Species == sp) %>% transform("compositional")
  phy_sub <- phy_sub %>% subset_taxa(taxa_sums(phy_sub) > 0)
  # Plot core taxa for this species
  p <- plot_core(phy_sub, prevalences = c(0.5, 0.6, 0.7, 0.8, 0.9, 1), detections = 10^(c(-6, -5, -4, -3)), 
          plot.type = "heatmap") +
          scale_fill_viridis_c(limits = c(0.5, 1)) +
          theme(legend.position = "none", axis.text.y = element_blank()) +
          xlab("") + ylab("") + labs(title = sp)
  if (sp == species[length(species)]) {p <- p + theme(legend.position = "right")} # Add legend to last plot
  core_plots[[sp]] <- p
  # Get core microbiota for this species
  core_taxa <- phy_sub %>% core_members(prevalence = prev, detection = det)
  # If there are no core taxa, print warning and skip
  # otherwise, add to df
  if (length(core_taxa) == 0) {
    warning(paste("No core taxa found for species:", sp))
    next
  }
  df <- data.frame(core_taxa = core_taxa,
                   host_species = sp)
  core_genera <- rbind(core_genera, df)
}

# Save core taxa plots
p <- plot_grid(plotlist = core_plots, ncol = 4, align = "hv", axis = "tb")
ggsave(p, file = file.path(subdir, "core_thresholds.png"), width = 20, height = 20)

# Add extra information
# on the taxa
core_genera$core_phylum <- as.vector(phy_gen@tax_table[match(core_genera$core_taxa, taxa_names(phy_gen)), "phylum"])

# and on the samples
core_genera$host_order <- phy_gen@sam_data$Order[match(core_genera$host_species, phy_gen@sam_data$Species)]

core_genera <- core_genera %>% arrange(host_order, host_species, core_taxa) %>% 
        select(host_order, host_species, core_phylum, core_taxa)

# Calculate number of core taxa per species
core_genera_summary <- 
        core_genera %>% group_by(host_species, host_order) %>% 
        summarise(n_core_taxa = n(),
                  n_core_phyla = n_distinct(core_phylum))

# Add number of samples
nsamples <- data.frame(phy_gen@sam_data) %>% group_by(Species) %>% summarise(n_samples = n())

core_genera_summary <- core_genera_summary %>%
        left_join(nsamples, by = c("host_species" = "Species")) %>%
        arrange(host_order, host_species)

# Save tables
write.csv(core_genera, file = file.path(subdir, "core_mb_per_host_species.csv"), quote = FALSE, row.names = FALSE)
write.csv(core_genera_summary, file = file.path(subdir, "core_mb_per_host_species_summary.csv"), quote = FALSE, row.names = FALSE)

cor.test(core_genera_summary$n_core_taxa, core_genera_summary$n_samples, method = "spearman")

#### Core taxa per host order ####
# For orders with more than two species, get core taxa per order
# defined as taxa that are considered core in at least half of species in that order

core_genera_per_order <- core_genera %>% group_by(host_order) %>%
        # Calculate number of species in order
        mutate(n_species = n_distinct(host_species)) %>%
        filter(n_species > 2) %>%
        # Calculate prevalence in order
        group_by(core_taxa, core_phylum, host_order) %>%
        summarise(prevalence_in_order = n_distinct(host_species)/first(n_species)) %>%
        filter(prevalence_in_order >= 2/3) %>%
        arrange(host_order, core_taxa) %>%
        # How many orders is this taxon core in?
        group_by(core_taxa, core_phylum) %>%
        mutate(n_orders = n_distinct(host_order))

# Calculate number of core taxa per order
core_genera_per_order_summary <- core_genera_per_order %>% group_by(host_order) %>%
        summarise(n_core_taxa = n(),
                n_core_phyla = n_distinct(core_phylum))

# Save tables
write.csv(core_genera_per_order, file = file.path(subdir, "core_mb_per_host_order.csv"), quote = FALSE, row.names = FALSE)
write.csv(core_genera_per_order_summary, file = file.path(subdir, "core_mb_per_host_order_summary.csv"), quote = FALSE, row.names = FALSE)

#### Plot Venn Diagrams ####

# Host order level
# Collect core taxa in a list
core_taxa <- list()

for (ord in unique(core_genera_per_order$host_order)) {
    core_taxa[[ord]] <-
        core_genera_per_order %>% filter(host_order == ord) %>%
        pull(core_taxa) %>% unique
}

core_genera_venn <- 
  ggVennDiagram(core_taxa, set_color = order_palette[names(core_taxa)], set_size = 0, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "grey40", name = "N. genera") +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "bottom", plot.margin = margin(t = 0, r = 20, b = 0, l = 20))

ggsave(core_genera_venn, file=file.path(subdir, "core_genera_venn.png"),
       device = "png", width = 5, height = 5)

#### Does number of species per order affect core size ####

# Permute which species belong to which order to understand if a small number of species
# leads to a larger core microbiome size

set.seed(123)

n_perm <- 100

species_order_map <- core_genera %>% select(host_species, host_order) %>% unique

perm_results <- data.frame()

for (i in 1:n_perm) {
        # Permute host orders across species
        species_order_map_perm <- species_order_map
        species_order_map_perm$host_order <- sample(species_order_map_perm$host_order)
        
        # Add permuted host order
        perm_core_genera <- core_genera
        perm_core_genera$host_order <- species_order_map_perm$host_order[match(perm_core_genera$host_species,
                                                                 species_order_map_perm$host_species)]
        
        # Calculate number of core taxa per order
        perm_core_genera_per_order <- 
                perm_core_genera %>% group_by(host_order) %>%
                # Calculate number of species in order
                mutate(n_species = n_distinct(host_species)) %>%
                filter(n_species > 2) %>%
                # Calculate prevalence in order
                group_by(core_taxa, core_phylum, host_order) %>%
                summarise(prevalence_in_order = n_distinct(host_species)/first(n_species)) %>%
                filter(prevalence_in_order >= 2/3) %>%
                arrange(host_order, core_taxa) %>%
                # How many orders is this taxon core in?
                group_by(core_taxa, core_phylum) %>%
                mutate(n_orders = n_distinct(host_order))
        
        # Calculate number of core taxa per order
        perm_core_genera_per_order_summary <- perm_core_genera_per_order %>%
                        group_by(host_order) %>%
                        summarise(n_core_taxa = n(),
                        n_core_phyla = n_distinct(core_phylum))
        
        # Add results to df
        perm_results <- rbind(perm_results, perm_core_genera_per_order_summary)
}

#### Plot ####

p <- ggplot(data = perm_results, aes(y = n_core_taxa, x = host_order)) +
        geom_boxplot(aes(fill = host_order)) +
        scale_fill_manual(values = order_palette, guide = "none") +
        geom_point(data = core_genera_per_order_summary, aes(y = n_core_taxa, x = host_order),
                color = "red", size = 3, shape = 18) +
        xlab("") + ylab("Number of core genera") +
        theme(axis.text.x = element_text(hjust = 1))

ggsave(p, file = file.path(subdir, "core_size_permutations.png"), width = 3, height = 5)

#### Core genera per host species heatmap ####
core_df <- core_genera %>% mutate(is.core = 1) %>%
        ggplot(aes(x = host_species, y = core_taxa, fill = is.core)) +
        facet_grid(cols = vars(host_order), rows = vars(core_phylum), scales = "free", space = "free") +
        geom_tile() +
        theme(strip.text.y = element_text(angle = 0,),
              strip.text.x = element_text(angle = 90),
              axis.text.x = element_text(hjust = 1, vjust = 0.5),
              legend.position = "none")

ggsave(core_df, file = file.path(subdir, "core_genera_heatmap.png"),
       device = "png", width = 10, height = 20)

#### PCA with only core order taxa ####

phy_core <- phy_gen %>% subset_taxa(genus %in% unique(unlist(core_taxa))) %>%
    subset_samples(Order %in% core_genera_per_order$host_order) 

phy_core_clr <- transform(phy_core, "clr")

# PCA ordination
ord <- ord_calc(phy_core_clr, method = "PCA") 

# Get shape scales for plotting
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

# Axes 1 & 2
p <- 
  ord_plot(ord, colour="Order", shape = "diet.general", alpha = 0.5, plot_taxa = 1:10) +
  custom_theme() +
  scale_shape_manual(values = diet_shape_scale, name = "Diet") +
  scale_color_manual(values=order_palette, name = "Order") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  geom_phylopic(data = centroids(ord@ord, phy_core_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_core_clr)$uid, width = 0.07, alpha = 0.8)

ggsave(p, file = file.path(subdir, "pca_core_order_1_2.png"), width = 5, height = 6)

# Axes 3 & 4
p <- 
  ord_plot(ord, colour="Order", shape = "diet.general", alpha = 0.5, plot_taxa = 1:10, axes = c(3, 4)) +
  custom_theme() +
  scale_color_manual(values=order_palette, name = "Order") +
  scale_shape_manual(values = diet_shape_scale, name = "Diet") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  geom_phylopic(data = centroids(ord@ord, phy_core_clr), aes(colour = Order_grouped), uuid = centroids(ord@ord, phy_core_clr)$uid, width = 0.07, alpha = 0.8)

ggsave(p, file = file.path(subdir, "pca_core_order_3_4.png"), width = 5, height = 6)

#### Abundance heatmap ####

# Get more info on taxon category
heat_data <- phy_gen %>% transform("compositional") %>% psmelt %>%
        select(OTU, genus, Sample, Species, Order, Abundance) %>%
        # Turn abundances below detection limit to NA
        mutate(Abundance = ifelse(Abundance <= det, NA, Abundance)) %>%
        filter(OTU %in% taxa_names(phy_core)) %>%
        filter(Sample %in% sample_names(phy_core)) %>%
        mutate(taxon_type = case_when(OTU %in% Reduce(intersect, core_taxa) ~ "Mammalian core",
                                      OTU %in% setdiff(core_taxa[["Primates"]], unlist(core_taxa[!names(core_taxa) %in% "Primates"])) ~ "Primates core",
                                      OTU %in% setdiff(core_taxa[["Perissodactyla"]], unlist(core_taxa[!names(core_taxa) %in% "Perissodactyla"])) ~ "Perissodactyla core",
                                      OTU %in% setdiff(core_taxa[["Carnivora"]], unlist(core_taxa[!names(core_taxa) %in% "Carnivora"])) ~ "Carnivora core",
                                      TRUE ~ "Other")) %>%
        mutate(taxon_type = factor(taxon_type, levels = c("Mammalian core", "Primates core", "Perissodactyla core", "Carnivora core", "Other"))) %>%
        # Shorten order names for plotting
        mutate(Order = recode(Order, "Perissodactyla" = "Peris.", "Carnivora" = "Carn.")) %>%
        # Arrange sample by species name
        arrange(Order, Species, Sample) %>%
        mutate(Sample = factor(Sample, levels = unique(Sample)))

grad_palette <- colorRampPalette(c("#2D627B","#FFF7A4", "#E7C46E","#C24141"))
grad_palette <- grad_palette(10)

p <- ggplot(heat_data, aes(x = Sample, y = OTU, fill = Abundance)) +
        geom_tile() +
        facet_grid(cols = vars(Order), scales = "free", space = "free",
                   rows = vars(taxon_type), switch = "y") +
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
              strip.text.y.left = element_text(angle=0, vjust=1, size = 12, hjust = 1, face = "bold"),
              strip.text.y = element_text(margin = margin(t=-15, r = -60)),
              strip.background.y = element_blank())

# Add species bar
species_bar <- data.frame(Sample = unique(heat_data$Sample),
                          Species = phy_core@sam_data$Species[match(unique(heat_data$Sample), rownames(phy_core@sam_data))],
                          Order = phy_core@sam_data$Order[match(unique(heat_data$Sample), rownames(phy_core@sam_data))]) %>%
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
        geom_phylopic(aes(x = Sample, y = 1, uuid = uid), width = 15, alpha = 1, fill = "white",
                      position = position_jitter(height = 0.1), hjust = 0.5) +
        facet_grid(cols = vars(Order), scales = "free", space = "free") +
        scale_x_discrete(labels = sample_to_label) +
        theme(axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(hjust = 1, vjust = 0.5),
              axis.title.x = element_blank(),
              legend.position = "none",
              panel.background = element_rect(fill = "white", color = "white"),
              strip.background.x = element_blank(),
              strip.text.x = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 20, l = 0))

# Combine heatmap and species bar
p_grid <- plot_grid(p, p_bar, ncol = 1, rel_heights = c(1, 0.5), align = "v", axis = "lr")

ggsave(file.path(subdir, "core_abundance_heatmap.png"), p_grid, width=8, height=8)

#### Percentage of abundance consisting of core taxa ####
# This reflects how representative the above PCA is of the whole community

# Calculate the proportion of abundance that is made up of core taxa
phy_melt <- transform(phy_gen, "compositional") %>% psmelt %>% filter(Abundance > 0)

core_genera_per_order$is.core <- "Core in this order"

core_taxa_abund <-
        phy_melt %>%
        select(OTU, genus, Sample, Species, Common.name, Order, Abundance) %>%
        # Keep only the same orders we use for the PCA
        filter(Order %in% core_genera_per_order$host_order) %>%
        # Get core taxa per order
        left_join(core_genera_per_order, by = c("OTU" = "core_taxa", "Order" = "host_order", "genus" = "core_phylum"),
                  relationship = "many-to-one") %>%
        # Indicate if a taxon is part of the mammalian core mb.
        # Alternatively, if it isnt core in all mammals or that specific order, check if it is core in another order
        mutate(is.core = case_when(OTU %in% Reduce(intersect, core_taxa) ~ "Mammalian core",
                                   !is.na(is.core) ~ is.core,
                                   OTU %in% core_genera_per_order$core_taxa ~ "Core in other order")) %>%
        mutate(is.core = factor(is.core, levels = c("Core in other order", "Core in this order", "Mammalian core"))) %>%
        filter(!is.na(is.core)) %>%
        # Get relative abunddance per sample
        group_by(Sample, Species, Common.name, Order, is.core) %>%
        summarise(Abundance = sum(Abundance)) %>%
        # Summarize by species
        group_by(Species, Common.name, Order, is.core) %>%
        summarise(mean_abundance = mean(Abundance),
                  sd = sd(Abundance, na.rm = TRUE)) %>%
        # Shorten "Perissodactyla" for plotting
        mutate(Order = recode(Order, "Perissodactyla" = "Periss.", "Carnivora" = "Carn."))

# Plot
p1 <- ggplot(core_taxa_abund, aes(y = Common.name, x = mean_abundance, fill = is.core)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("#FFDC7C", "#DA8A3D", "#DA4C3D"), name = "") +
        facet_grid(rows = vars(Order), scales = "free_y", space = "free_y") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75)) +
        labs(x = "Rel. abundance\nof core taxa") +
        theme(legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank(),
              axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
              strip.background.y = element_blank(), strip.text.y = element_blank()) +
        guides(fill = guide_legend(reverse = TRUE))

#### Abundance of core genera ####
mamm_core_abund_all <- phy_melt %>%
        select(OTU, genus, Sample, Species, Common.name, Order, Abundance) %>%
        # Sum by genus
        group_by(Sample, Species, Common.name, Order, genus) %>%
        summarise(Abundance = sum(Abundance)) %>%
        # Calculate mean abundance per species
        group_by(Species, Common.name, Order, genus) %>%
        summarise(mean_abundance = mean(Abundance))

mamm_core_abund <- mamm_core_abund_all %>%
        # Keep only the same orders we use for the PCA
        filter(Order %in% core_genera_per_order$host_order) %>%
        # Get taxa that are the mammalian core taxa
        filter(OTU %in% unique(unlist(core_taxa))) %>%
        # Shorten "Perissodactyla" for plotting
        mutate(Order = recode(Order, "Perissodactyla" = "Periss.", "Carnivora" = "Carn."))

# Identify the 5 most abundant genera, group the rest as "Other"
top_genera <- mamm_core_abund %>% group_by(genus) %>% summarise(mean_abundance = mean(mean_abundance)) %>%
        arrange(desc(mean_abundance)) %>% top_n(5, mean_abundance) %>% pull(genus)

mamm_core_abund <- mamm_core_abund %>% mutate(genus = ifelse(genus %in% top_genera, genus, "Other")) %>%
                  mutate(genus = factor(genus, levels = rev(c(top_genera, "Other"))))

# Plot
colours <- rev(c("transparent", c(wes_palette("Darjeeling1", 5, type = "discrete"))))
names(colours) <- c(top_genera, "Other")

p2 <- ggplot(mamm_core_abund, aes(y = Common.name, x = mean_abundance, fill = genus)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = colours) +
        facet_grid(rows = vars(Order), scales = "free_y", space = "free_y") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75)) +
        labs(x = "Rel. abundance of\n core genera") +
        theme(legend.position = "bottom", legend.direction = "vertical",
              legend.title = element_blank(),  
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        guides(fill = guide_legend(reverse = TRUE, ncol = 2))

p <- plot_grid(p1, p2, ncol = 2, align = "h", axis = "tb", rel_widths = c(1.5, 1))

ggsave(p, file = file.path(subdir, "core_taxa_abundance.png"), width = 8, height = 8)

######################
#### ORCA CORE MB ####
######################

#### Compare orca core mb to Artiodactyla and Carnivora ####

carni_core <- core_genera_per_order %>% filter(host_order == "Carnivora" & n_orders == 1) %>% pull(core_taxa)

orca_core_comparison <- core_genera %>%
        filter(host_species %in% c("Capreolus capreolus", "Orcinus orca", "Sus scrofa", "Hippopotamus amphibius", "Meles meles", "Otaria flavescens")) %>%
        # Get abundances 
        left_join(mamm_core_abund_all, by = c("core_taxa" = "genus", "host_species" = "Species", "host_order" = "Order")) %>%
        group_by(host_species, Common.name, host_order) %>%
        summarise(Carnivora_taxa = sum(core_taxa %in% carni_core)) %>%
        ungroup() %>%
        mutate(Common.name = ifelse(host_species == "Otaria flavescens", "Sea lion",
                                 ifelse(host_species == "Meles meles", "Badger", Common.name))) %>%
        mutate(Common.name = factor(Common.name, levels = c("Roe deer", "Hippo", "Wild boar", "Orca", "Sea lion", "Badger")))

write.csv(orca_core_comparison, file = file.path(subdir, "orca_core_comparison.csv"), quote = FALSE, row.names = FALSE)

p <- ggplot(data = orca_core_comparison, aes(x = Common.name, y = Carnivora_taxa, fill = host_order)) +
        geom_bar(stat = "identity", position = position_dodge2()) +
        scale_fill_manual(values = order_palette, name = "Host order") +
        theme(axis.text.x = element_text(hjust = 1), axis.title.x = element_blank(), legend.title = element_blank(),
                legend.position = "bottom", legend.direction = "vertical") +
        labs(y = "Genera shared with\nCarnivora core exclusive")

ggsave(file.path(subdir, "orca_core_comparison.png"), p, width = 3, height = 5)
