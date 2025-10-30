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
library(ape)
library(reshape2)
library(microViz)
library(rphylopic)
library(RColorBrewer)
library(wesanderson)
library(ggplot2)
library(ggnewscale)
library(ggtree)
library(ggtreeExtra)
library(vegan)
library(microbiomeutilities)
library(cowplot)
library(ggExtra)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # subdirectory for the output of this script
subdir <- normalizePath(file.path(outdir, "multivariate_stats")) # subdirectory for the output of this script
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

# Get functions
source(file.path("..", "ordination_functions.R"))
source(file.path("..", "phylo_functions.R"))

#######################
#####  LOAD INPUT #####
#######################

# Load all phyloseq objects in phydir
for (phy_file in list.files(phydir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(phydir, phy_file)))
}

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

host_consensus <- read.tree(file.path(outdir, "host_consensus.tre"))

#########################
#####  BETA DISPER  #####
#########################

# Calculate beta disper for species, order and diet, and run tukey's test

tukey_results <- data.frame()
plot_list <- list()

# Compare species
for (var in c("Species", "Order", "diet.general")) {
  disp <- betadisper(vegdist(t(otu_table(phy_sp_f_clr)), method = "euclidean"), group = phy_sp_f_clr@sam_data[[var]])
  
  disp_tukey <- TukeyHSD(disp, which = "group", ordered = FALSE)$group %>% data.frame %>% rownames_to_column("Comparison") %>%
    separate(Comparison, into = c("Group1", "Group2"), sep = "-")
  
  disp_tukey$Variable <- var
    
  tukey_results <- rbind(tukey_results, disp_tukey)
  
  disp_df <- data.frame(Sample = sample_names(phy_sp_f_clr),
                        Species = phy_sp_f_clr@sam_data$Species,
                        diet.general = phy_sp_f_clr@sam_data$diet.general,
                        Order_grouped = phy_sp_f_clr@sam_data$Order_grouped,
                        Order = phy_sp_f_clr@sam_data$Order,
                        Distance = disp$distances)
  if (var %in% c("Species", "diet.general")) {
    p <- ggplot(data = disp_df, aes(x = !!sym(var), y = Distance)) +
      geom_boxplot(outlier.shape = NA, aes(fill = diet.general)) +
      geom_jitter(alpha = 0.5, width = 0.2) +
      scale_fill_manual(values = diet_palette, name = "Estimated diet") +
      theme(legend.position = "none", axis.text.x = element_text(hjust = 1))
  } else {
    p <- ggplot(data = disp_df, aes(x = !!sym(var), y = Distance)) +
     geom_boxplot(outlier.shape = NA, aes(fill = Order)) +
     geom_jitter(alpha = 0.5, width = 0.2) +
     scale_fill_manual(values = order_palette, name = "Order") +
     theme(legend.position = "none", axis.text.x = element_text(hjust = 1))
  }
  plot_list[[var]] <- p
}

write.csv(tukey_results, file = file.path(subdir, "betadisper_tukey_results.csv"), row.names = FALSE, quote = FALSE)

p <- plot_grid(plotlist = plot_list, ncol = 1, align = "v")

ggsave(file.path(subdir, "betadisper.png"), p, width = 8, height = 12)

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

# Order by Pseudomonadota
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

##### Composition on phylogenetic tree ####

# Create a copy of the tree
sample_tree <- host_consensus
sample_tree$tip.label <- gsub("_", " ", sample_tree$tip.label)

# Create a list of samples per species
sample_map <- split(rownames(phy_sp_f@sam_data), phy_sp_f@sam_data$Species)

# For each species, add samples as zero-length branches
for (species in names(sample_map)) {
  if (species %in% sample_tree$tip.label) {
    cat("Creating sample tips for", species, "\n")
    samples <- sample_map[[species]]
    # Create a mini tree (polytomy) for the samples
    polytomy <- stree(length(samples), type = "star")
    polytomy$tip.label <- samples
    polytomy$node.label <- species
    # Bind the polytomy at the species tip
    polytomy$edge.length <- rep(0.5, nrow(polytomy$edge))
    cat("Polytomy:\n")
    print(polytomy)
    sample_tree <- bind.tree(sample_tree, polytomy, where = which(sample_tree$tip.label == species))
    # Drop the original species tip
    sample_tree <- drop.tip(sample_tree, species)
  }
}

# Get traits for tips
host_traits <- data.frame(tip = sample_tree$tip.label,
                          node=nodeid(sample_tree, sample_tree$tip.label),
                          Order = phy_sp_f@sam_data[sample_tree$tip.label, "Order"],
                          Species = phy_sp_f@sam_data[sample_tree$tip.label, "Species"],
                          Common.name = phy_sp_f@sam_data[sample_tree$tip.label, "Common.name"],
                          diet.general = phy_sp_f@sam_data[sample_tree$tip.label, "diet.general"],
                          PlantO = phy_sp_f@sam_data[sample_tree$tip.label, "PlantO"],
                          Animal = phy_sp_f@sam_data[sample_tree$tip.label, "Animal"],
                          habitat.general = phy_sp_f@sam_data[sample_tree$tip.label, "habitat.general"])

# Add order trait to nodes
node_traits <- tips_to_nodes(sample_tree, host_traits, "Order") %>% left_join(host_traits) %>% arrange(node)

# Plot tree
p <- 
  ggtree(sample_tree, aes(colour = Order), layout = "circular", open.angle = 180, size=1) %<+% node_traits +
  # Colour branches by order
  scale_color_manual(values = order_palette) +
  new_scale_colour() +
  # Colour tips by species
  geom_tippoint(aes(colour = Species)) +
  scale_color_manual(values = species_palette) +
  theme(legend.position = "none")

# Plot
p_fruit <- p +
  geom_fruit(
    data = rename(phy_phylum_melt, "label" = "Sample"),
    geom = geom_bar(),
    mapping = aes(y = label, x = Abundance, fill = OTU),
    stat = "identity",
    axis.params = list(axis = "x", text.size = 1, hjust = 1, vjust = 0., nbreak = 3),
    offset = 0.05,
    pwidth = 0.3
  ) +
  scale_fill_manual(values = phylum_palette)

#### Add phylopic at the middle sample of each species
tree_phylopics <- host_traits %>% left_join(phylopics) %>% select(tip, Species, Common.name, uid) %>% group_by(Species, Common.name) %>%
            mutate(is.middle = (row_number() == floor(n_distinct(tip)/2))) %>% ungroup %>%
            mutate(uid = case_when(is.middle ~ uid),
                  species_lab = case_when(is.middle ~ Common.name)) %>%
            select(tip, uid, species_lab) %>% rename(label = tip)

p_fruit <- p_fruit %<+% tree_phylopics +
           new_scale_colour() +
           geom_phylopic(aes(uuid=uid, colour = diet.general), width = 8, position=position_nudge(x=70)) +
           geom_text(aes(label = species_lab), size = 3.5, position = position_nudge(x=110)) +
           scale_colour_manual(values = diet_palette)

ggsave(file.path(subdir, "phy_sp_f_composition_tree.png"), p_fruit, width = 10, height = 10)

##############################################
#### ORDINATIONS & PERMANOVA FULL DATASET ####
##############################################

#### Get shape scales for plotting ####
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)
uniq_species <- unique(subset_samples(phy_sp_f_clr, Order %in% c("Primates", "Carnivora", "Artiodactyla") | habitat.general == "Marine")@sam_data$Common.name)
species_shape_scale <- c(1:25, 35:35+26-length(uniq_species))
names(species_shape_scale) <- uniq_species
order_shape_scale <- c("Carnivora" = 4, "Primates" = 19, "Artiodactyla" = 5, "Perissodactyla" = 2, "Rodentia" = 1, "Rest" = 12, "Proboscidea_Sirenia" = 12)

#### PCA ####

ord <- ord_calc(phy_sp_f_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "PCA_all_screeplot.png"), p, width=8, height=6)

# Color by order
p <- custom_ord_plot(phy_sp_f_clr, ord, colour="Order_grouped", shape="diet.general", arrows_scaling = 2, type = "PCA")

ggsave(file.path(subdir, "PCA_all_clr_order_1_2.png"), p, width=6, height=6)

# Color by diet
p <- custom_ord_plot(phy_sp_f_clr, ord, colour="diet.general", shape="Order_grouped", arrows_scaling = 2, type = "PCA")

ggsave(file.path(subdir, "PCA_all_clr_diet_1_2.png"), p, width=6, height=6)
 
#### RDA ####

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
                    "ruminant", "marine", "Fruit", "Animal")

# Ordinate using all data
ord <- ord_calc(phy_sp_f_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:10))

ggsave(file.path(subdir, "RDA_all_screeplot.png"), p, width=8, height=6)

## SAMPLE PLOTS

# Color by diet
p <- custom_ord_plot(phy_sp_f_clr, ord, colour="diet.general", shape="Order_grouped", type = "RDA")

ggsave(file.path(subdir, "RDA_all_clr_diet_1_2.png"), p, width=6, height=6)

# Colour by order
p <- custom_ord_plot(phy_sp_f_clr, ord, colour="Order_grouped", shape="diet.general", type = "RDA")

ggsave(file.path(subdir, "RDA_all_clr_order_1_2.png"), p, width=6, height=6)

## TAXA PLOT
p <- taxa_plot(ord, phy_sp_f_clr)
ggsave(file.path(subdir, "RDA_all_clr_taxa_1_2.png"), p, width=5, height=5)

#### PERMANOVA ####

# Explanatory variables
sample_data <- as.data.frame(phy_sp_f_clr@sam_data)

order <- sample_data$Order
diet <- sample_data$diet.general
habitat <- sample_data$habitat.general
ruminant <- (sample_data$digestion == "Ruminant")
hypsodont <- grepl("hyps", sample_data$molar_category)
species <- sample_data$Species

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(t(otu_table(phy_sp_f_clr)) ~ order + diet + habitat + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_all_clr_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(t(otu_table(phy_sp_f_clr)) ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_all_clr_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

##############################################
#### ORDINATIONS & PERMANOVA DEEP DATASET ####
##############################################

#### PCA ####

ord <- ord_calc(phy_deep_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "PCA_deep_screeplot.png"), p, width=8, height=6)

# Color by order
p <- custom_ord_plot(phy_deep_clr, ord, colour="Order_grouped", shape="diet.general", arrows_scaling = 2, type = "PCA")

ggsave(file.path(subdir, "PCA_deep_clr_order_1_2.png"), p, width=6, height=6)

# Color by diet
p <- custom_ord_plot(phy_deep_clr, ord, colour="diet.general", shape="Order_grouped", arrows_scaling = 2, type = "PCA")

ggsave(file.path(subdir, "PCA_deep_clr_diet_1_2.png"), p, width=6, height=6)
 
#### RDA ####

# Recode order and habitat as TRUE and FALSE
# Also scale protein, fiber and carbohydrate content
phy_deep_clr <- phy_deep_clr %>%
        ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  ruminant = (digestion == "Ruminant"),
                  marine = (habitat.general == "Marine"))

# Species traits to use as constraints
species_traits <- c("Artiodactyla", "Perissodactyla", "Primates",
                    "ruminant", "marine", "Fruit", "Animal")

# Ordinate using all data
ord <- ord_calc(phy_deep_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:10))

ggsave(file.path(subdir, "RDA_deep_screeplot.png"), p, width=8, height=6)

## SAMPLE PLOTS

# Color by diet
p <- custom_ord_plot(phy_deep_clr, ord, colour="diet.general", shape="Order_grouped", type = "RDA")

ggsave(file.path(subdir, "RDA_deep_clr_diet_1_2.png"), p, width=6, height=6)

# Colour by order
p <- custom_ord_plot(phy_deep_clr, ord, colour="Order_grouped", shape="diet.general", type = "RDA")

ggsave(file.path(subdir, "RDA_deep_clr_order_1_2.png"), p, width=6, height=6)

## TAXA PLOT
p <- taxa_plot(ord, phy_deep_clr)
ggsave(file.path(subdir, "RDA_deep_clr_taxa_1_2.png"), p, width=5, height=5)

## RDA axis violin plots

rda_res <- data.frame(phy_deep_clr@sam_data) %>% select(new_name, Species, Order, Order_grouped, diet.general) %>%
           left_join(rownames_to_column(data.frame(vegan::scores(ord@ord, display="sites", choices=1:3))), by=c("new_name"="rowname")) %>%
           pivot_longer(cols = c(RDA1, RDA2, RDA3), names_to = "Axis", values_to = "Value")

p <- ggplot(rda_res, aes(y = Order, x = Value, fill = diet.general, colour = diet.general)) +
      geom_violin() +
      scale_fill_manual(values = diet_palette, name = "Diet category") +
      scale_color_manual(values = darken(diet_palette, amount = 0.3), name = "Diet category") +
      facet_grid(cols = vars(Axis), scales = "free") +
      theme(legend.position = "none", axis.title = element_blank())

ggsave(file.path(subdir, "RDA_deep_clr_axis_comparison.png"), p, width=8, height=8)

## Test

#anova(ord@ord, by = "margin", perm = 500)
#anova(ord@ord, by = "axis", perm = 500)

#### PERMANOVA ####

# Explanatory variables
sample_data <- as.data.frame(phy_deep_clr@sam_data)

order <- sample_data$Order
diet <- sample_data$diet.general
habitat <- sample_data$habitat.general
ruminant <- (sample_data$digestion == "Ruminant")
hypsodont <- grepl("hyps", sample_data$molar_category)
species <- sample_data$Species

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(t(otu_table(phy_deep_clr)) ~ order + diet + habitat + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_deep_clr_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(t(otu_table(phy_deep_clr)) ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_deep_clr_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

##############################################
#### ORDINATIONS & PERMANOVA ARTIODACTYLA ####
##############################################

#### PCA ####

ord <- ord_calc(phy_artio_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "PCA_artio_screeplot.png"), p, width=8, height=6)

# Color by diet
p <- custom_ord_plot(phy_artio_clr, ord, colour="diet.general", shape="digestion", arrows_scaling = 2, type = "PCA")

ggsave(file.path(subdir, "PCA_artio_clr_1_2.png"), p, width=6, height=6)

#### RDA ####

# Recode order and habitat as TRUE and FALSE
# Also scale protein, fiber and carbohydrate content
phy_artio_clr <- phy_artio_clr %>%
        ps_mutate(ruminant = (digestion == "Ruminant"),
                  marine = (habitat.general == "Marine"))

# Species traits to use as constraints
species_traits <- c("ruminant", "marine", "Fruit")

# Ordinate using all data
ord <- ord_calc(phy_artio_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:10))

ggsave(file.path(subdir, "RDA_artio_screeplot.png"), p, width=8, height=6)

## SAMPLE PLOTS

# Color by diet
p <- custom_ord_plot(phy_artio_clr, ord, colour="diet.general", shape="digestion", type = "RDA")

ggsave(file.path(subdir, "RDA_artio_clr_1_2.png"), p, width=6, height=6)

## TAXA PLOT 
p <- taxa_plot(ord, phy_artio_clr)
ggsave(file.path(subdir, "RDA_artio_clr_taxa_1_2.png"), p, width=5, height=5)

#### PERMANOVA ####

# Explanatory variables
sample_data <- as.data.frame(phy_artio_clr@sam_data)

diet <- sample_data$diet.general
ruminant <- (sample_data$digestion == "Ruminant")
hypsodont <- grepl("hyps", sample_data$molar_category)
species <- sample_data$Species

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(t(otu_table(phy_artio_clr)) ~ diet + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_artio_clr_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(t(otu_table(phy_artio_clr)) ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_artio_clr_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

###########################################
#### ORDINATIONS & PERMANOVA CARNIVORA ####
###########################################

#### PCA ####

ord <- ord_calc(phy_carni_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "PCA_carni_screeplot.png"), p, width=8, height=6)

# Color by diet
p <- custom_ord_plot(phy_carni_clr, ord, colour="habitat.general", shape="Common.name", arrows_scaling = 2, type = "PCA")

ggsave(file.path(subdir, "PCA_carni_clr_1_2.png"), p, width=6, height=6)

#### RDA ####

# Recode order and habitat as TRUE and FALSE
# Also scale protein, fiber and carbohydrate content
phy_carni_clr <- phy_carni_clr %>%
        ps_mutate(marine = (habitat.general == "Marine"))

# Species traits to use as constraints
species_traits <- c("marine", "Animal")

# Ordinate using all data
ord <- ord_calc(phy_carni_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:10))

ggsave(file.path(subdir, "RDA_carni_screeplot.png"), p, width=8, height=6)

## SAMPLE PLOTS

# Color by diet
p <- custom_ord_plot(phy_carni_clr, ord, colour="habitat.general", shape="Common.name", type = "RDA")

ggsave(file.path(subdir, "RDA_carni_clr_1_2.png"), p, width=6, height=6)

## TAXA PLOT 
p <- taxa_plot(ord, phy_carni_clr)
ggsave(file.path(subdir, "RDA_carni_clr_taxa_1_2.png"), p, width=5, height=5)

#### PERMANOVA ####

# Explanatory variables
sample_data <- as.data.frame(phy_carni_clr@sam_data)

diet <- sample_data$diet.general
habitat <- sample_data$habitat.general
species <- sample_data$Species

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(t(otu_table(phy_carni_clr)) ~ diet + habitat,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_carni_clr_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(t(otu_table(phy_carni_clr)) ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_carni_clr_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

##########################################
#### ORDINATIONS & PERMANOVA PRIMATES ####
##########################################

#### PCA ####

ord <- ord_calc(phy_prim_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "PCA_prim_screeplot.png"), p, width=8, height=6)

# Color by diet
p <- custom_ord_plot(phy_prim_clr, ord, colour="diet.general", shape="Common.name", arrows_scaling = 2, type = "PCA")

ggsave(file.path(subdir, "PCA_prim_clr_1_2.png"), p, width=6, height=6)

#### RDA ####

# Species traits to use as constraints
species_traits <- c("PlantO", "Animal")

# Ordinate using all data
ord <- ord_calc(phy_prim_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:10))

ggsave(file.path(subdir, "RDA_prim_screeplot.png"), p, width=8, height=6)

## SAMPLE PLOTS

# Color by diet
p <- custom_ord_plot(phy_prim_clr, ord, colour="diet.general", shape="Common.name", type = "RDA")

ggsave(file.path(subdir, "RDA_prim_clr_1_2.png"), p, width=6, height=6)

## TAXA PLOT 
p <- taxa_plot(ord, phy_prim_clr)
ggsave(file.path(subdir, "RDA_prim_clr_taxa_1_2.png"), p, width=5, height=5)

#### PERMANOVA ####

# Explanatory variables
sample_data <- as.data.frame(phy_prim_clr@sam_data)

diet <- sample_data$diet.general
species <- sample_data$Species

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(t(otu_table(phy_prim_clr)) ~ diet,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_prim_clr_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(t(otu_table(phy_prim_clr)) ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_prim_clr_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

#########################################
#### ORDINATIONS & PERMANOVA HABITAT ####
#########################################

#### PCA ####

ord <- ord_calc(phy_habitat_clr, method = "PCA")

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

ggsave(file.path(subdir, "PCA_prim_screeplot.png"), p, width=8, height=6)

# Color by diet
p <- custom_ord_plot(phy_habitat_clr, ord, colour="habitat.general", shape="Order_grouped", arrows_scaling = 2, type = "PCA")

ggsave(file.path(subdir, "PCA_habitat_clr_1_2.png"), p, width=6, height=6)

#### RDA ####
# Recode order and habitat as TRUE and FALSE
phy_habitat_clr <- phy_habitat_clr %>%
        ps_mutate(marine = (habitat.general == "Marine"),
                  Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Order_grouped = case_when(Order %in% c("Proboscidea", "Sirenia") ~ "Proboscidea_Sirenia",
                                        TRUE ~ Order_grouped))

# Species traits to use as constraints
species_traits <- c("marine", "Artiodactyla", "Carnivora")

# Ordinate using all data
ord <- ord_calc(phy_habitat_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:10))

ggsave(file.path(subdir, "RDA_habitat_screeplot.png"), p, width=8, height=6)

## SAMPLE PLOTS

# Color by habitat
p <- custom_ord_plot(phy_habitat_clr, ord, colour="habitat.general", shape="Order_grouped", type = "RDA")

ggsave(file.path(subdir, "RDA_habitat_clr_1_2.png"), p, width=6, height=6)

## TAXA PLOT 
p <- taxa_plot(ord, phy_habitat_clr)
ggsave(file.path(subdir, "RDA_habitat_clr_taxa_1_2.png"), p, width=5, height=5)

#### PERMANOVA ####

# Explanatory variables
sample_data <- as.data.frame(phy_habitat_clr@sam_data)

habitat <- sample_data$habitat.general
order <- sample_data$Order_grouped
species <- sample_data$Species

set.seed(123)

# Run PERMANOVA with all factors and only species
perm <- adonis2(t(otu_table(phy_habitat_clr)) ~ order + habitat,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_habitat_clr_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(t(otu_table(phy_habitat_clr)) ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_habitat_clr_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

##########################
#### PHYLOPICS LEGEND ####
##########################

phylopics$Common.name <- phy_sp_f@sam_data$Common.name[match(phylopics$Species, as.vector(phy_sp_f@sam_data$Species))]

p <- ggplot(phylopics[!is.na(phylopics$Common.name),], aes(y = Common.name, x = 1, colour = Species)) +
    geom_phylopic(aes(uuid = uid), height = 0.8) +
    scale_color_manual(values = species_palette, name = "Species") +
    theme_void() + theme(axis.title = element_blank(), axis.text.x = element_blank(),
                         axis.text.y = element_text(size = 5, hjust = 1),
                         legend.position = "none",
                         plot.background = element_rect(fill = "transparent", color = "transparent"),
                         panel.background = element_rect(fill = "transparent", color = "transparent"))

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

set.seed(123)
png(filename = file.path(subdir, "heatmap_all_taxa.png"), width=16, height=20, units="in", res=300)
plot_taxa_heatmap(phy_sp_f, subset.top=ntaxa(phy_sp_f), transformation="clr",
                  VariableA=c("Order", "diet.general", "unmapped_count"),
                  annotation_colors = list("Order" = order_palette, "diet.general" = diet_palette),
                  show_rownames = FALSE, cluster_cols = TRUE,
                  show_colnames = FALSE, heatcolors = grad_palette)$plot
dev.off()
