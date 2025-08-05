##### PHYLOSYMBIOSIS TEST #####
## This script compares host phylogenetic distances, with microbiome distances and dietary distances
## to evaluate if microbiome evolution follows a phylosymbiotic pattern

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(reshape2)
library(vegan)
library(ape)
library(phytools)
library(phyloseq)
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Directory with output from the community analysis
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "phylosymbiosis")) # subdirectory for the output of this script
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

# Host phylogeny
host_trees <- read.nexus(file.path(indir, "mammal_vertlife.nex"))

########################
#### PHYLOSYMBIOSIS ####
########################

#### Distances ####
## Do host distances correlate with microbiome distances?

# Get host distances
# Get consensus tree and fix tip labels
host_consensus <- consensus.edges(host_trees, method="least.squares")

# Change tip labels to match metadata
host_consensus$tip.label <- host_consensus$tip.label %>%
                            gsub(pattern="Equus_quagga", replacement="Equus_burchellii") %>%
                            gsub(pattern="Procolobus_badius", replacement="Piliocolobus_foai") %>%
                            gsub(pattern="Otaria_bryonia", replacement="Otaria_byronia") %>%
                            gsub(pattern="_", replacement=" ", fixed=TRUE)

# Species distances
host_dist_melt <- cophenetic(host_consensus) %>% as.data.frame() %>% rownames_to_column("Item1") %>%
                melt(idvar = "Item1") %>% rename(Item2 = variable, host_distance = value)

# Samples to species
samples_to_species <- phy_sp_f@sam_data %>% data.frame %>% select(Species) %>% 
                      mutate(Species = recode(Species, "Sus scrofa domesticus" ="Sus scrofa")) %>%
                      rownames_to_column("Sample")

# Add to dist
samples_dist_melt <- host_dist_melt %>%
                    right_join(samples_to_species, by=c("Item1"="Species"), relationship = "many-to-many") %>%
                    right_join(samples_to_species, by=c("Item2"="Species"), relationship = "many-to-many", suffix = c("1","2")) %>%
                    select(Sample1, Sample2, host_distance)

samples_dist <- dcast(samples_dist_melt, Sample1 ~ Sample2, value = "host_distance") %>% column_to_rownames("Sample1") %>% as.dist

# Microbiome distances
mb_dist_clr <- vegdist(t(phy_sp_f_clr@otu_table), method="euclidean")
mb_dist_philr <- vegdist(phy_sp_philr@otu_table, method="euclidean")

# Diet distances after scaling
diet_scaled <- data.frame(phy_sp_f@sam_data) %>% select("cp", "cf", "nfe", "ee", "ash") %>% scale(center = TRUE, scale = TRUE) %>% as.data.frame
diet_dist <- vegdist(diet_scaled, method="euclidean")

## Mantel test
# Microbiome vs host phylogeny accounting for diet 
# CLR
mantel_mvp_clr <- mantel.partial(xdis = mb_dist_clr, ydis = samples_dist, zdis = diet_dist, method="spearman", permutations=500)
# PHILR
mantel_mvp_philr <- mantel.partial(xdis = mb_dist_philr, ydis = samples_dist, zdis = diet_dist, method="spearman", permutations=500)

#### Plot ####
# Melt microbiome distances
mb_dist_melt_clr <- mb_dist_clr %>% as.matrix %>% as.data.frame %>% rownames_to_column("Sample1") %>%
                melt(idvar = "Sample1") %>% rename(Sample2 = variable, microbiome_distance_clr = value)

mb_dist_melt_philr <- mb_dist_philr %>% as.matrix %>% as.data.frame %>% rownames_to_column("Sample1") %>%
                melt(idvar = "Sample1") %>% rename(Sample2 = variable, microbiome_distance_philr = value)

mb_dist_melt <- mb_dist_melt_clr %>% full_join(mb_dist_melt_philr, by=c("Sample1"="Sample1", "Sample2"="Sample2"))

# Melt diet distances
diet_dist_melt <- diet_dist %>% as.matrix %>% as.data.frame %>% rownames_to_column("Sample1") %>%
                melt(idvar = "Sample1") %>% rename(Sample2 = variable, diet_distance = value)

# Merge with host distances
dist_df <- samples_dist_melt %>% full_join(mb_dist_melt, by=c("Sample1"="Sample1", "Sample2"="Sample2")) %>%
            full_join(diet_dist_melt, by=c("Sample1"="Sample1", "Sample2"="Sample2"))

# Add info on diet and habitat
eco_info <- phy_sp_f@sam_data %>% data.frame %>% select(diet.general, habitat.general) %>%
                rownames_to_column("Sample") %>% rename(host_diet = diet.general,
                                                        host_habitat = habitat.general)

# Add to dist df 
dist_df <- dist_df %>% left_join(eco_info, by=c("Sample1"="Sample")) %>%
            left_join(eco_info, by=c("Sample2"="Sample"), suffix = c("1","2"))  %>%
            mutate(diet_comparison = factor(case_when(host_diet1 == host_diet2 ~ "Same",
                                               host_diet1 %in% c("Herbivore", "Frugivore") & host_diet2 %in% c("Herbivore", "Frugivore") ~ "Similar",
                                               host_diet1  %in% c("Omnivore", "Animalivore") & host_diet2 %in% c("Omnivore", "Animalivore") ~ "Similar",
                                               host_diet1 != host_diet2 ~ "Different"), c("Same", "Similar", "Different")),
                   habitat_comparison = ifelse(host_habitat1 == host_habitat2, "Same", "Different")) %>%
            # Remove intraspecific comparisons
            filter(host_distance > 0)
write.csv(dist_df, file.path(subdir, "host_microbiome_diet_distances.csv"), row.names = FALSE, quote = FALSE)

# CLR
p_clr <- ggplot(dist_df, aes(x=host_distance, y=microbiome_distance_clr, colour = diet_distance)) +
  geom_jitter(width = 5, alpha = 0.1, size = 0.5) +
  scale_colour_gradient2(high = "#F04C3B", mid = "#FFF88B", low = "#106B85", midpoint = median(dist_df$diet_distance), name = "Diet distance") +
  geom_smooth(method="lm", aes(linetype = diet_comparison), colour = "black", linewidth = 1) +
  scale_linetype_manual(values = c("Same" = "solid", "Similar" = "dashed", "Different" = "dotted"),
                        labels = c("Same" = "Same diet", "Similar" = "Similar diet","Different" = "Different diet"), name = "") +
  labs(x="", y="Distance of\nCLR-transformed abundances") +
  theme(legend.position="none") +
  annotate("text", x = 250, y = 30, label = paste0("Partial Mantel test,\nSpearman r = ", round(mantel_mvp_clr$statistic, 2), ", p = ", mantel_mvp_clr$signif), size = 3)

# PhILR
p_philr <- ggplot(dist_df, aes(x=host_distance, y=microbiome_distance_philr, colour = diet_distance)) +
  geom_jitter(width = 5, alpha = 0.1, size = 0.5) +
  scale_colour_gradient2(high = "#F04C3B", mid = "#FFF88B", low = "#106B85", midpoint = median(dist_df$diet_distance), name = "Diet distance") +
  geom_smooth(method="lm", aes(linetype = diet_comparison), colour = "black", linewidth = 1) +
  scale_linetype_manual(values = c("Same" = "solid", "Similar" = "dashed", "Different" = "dotted"),
                        labels = c("Same" = "Same diet", "Similar" = "Similar diet","Different" = "Different diet"), name = "") +
  labs(x="Host phylogenetic distance", y="Distance of\nPhILR-transformed abundances") +
  theme(legend.position="bottom", legend.direction = "vertical") +
  annotate("text", x = 250, y = 250, label = paste0("Partial Mantel test,\nSpearman r = ", round(mantel_mvp_philr$statistic, 2), ", p = ", mantel_mvp_philr$signif), size = 3)

p <- plot_grid(p_clr, p_philr, nrow = 2, rel_heights = c(1, 1.5), axis = "rl", align = "v")

ggsave(file.path(subdir, "host_microbiome_distance_correlation.png"), p, width=5, height=8)

# The same but only within diets
dist_df_same_diet <- dist_df %>% filter(diet_comparison == "Same") %>% filter(host_diet1 != "Omnivore") # Remove omnivores because there is only one comparison

p_clr <- ggplot(dist_df_same_diet, aes(x=host_distance, y=microbiome_distance_clr, colour = host_diet1, shape = habitat_comparison)) +
  geom_jitter(width = 5, alpha = 0.2, size = 0.5) +
  scale_color_manual(values = diet_palette, name = "Diet category") +
  scale_shape_manual(values = c("Same" = 16, "Different" = 3), label = c("Same" = "Same habitat", "Different" = "Marine vs. Terrestrial"), name = "") +
  geom_smooth(method="lm", aes(group = host_diet1)) +
  labs(x="", y="Distance of\nCLR-transformed abundances") +
  theme(legend.position="none")

p_philr <- ggplot(dist_df_same_diet, aes(x=host_distance, y=microbiome_distance_philr, colour = host_diet1, shape = habitat_comparison)) +
  geom_jitter(width = 5, alpha = 0.2, size = 0.5) +
  scale_color_manual(values = diet_palette, name = "Diet category") +
  scale_shape_manual(values = c("Same" = 16, "Different" = 3), label = c("Same" = "Same habitat", "Different" = "Marine vs. Terrestrial"), name = "") +
  geom_smooth(method="lm", aes(group = host_diet1)) +
  labs(x="Host phylogenetic distance", y="Distance of\nPhILR-transformed abundances") +
  theme(legend.position="bottom")

p <- plot_grid(p_clr, p_philr, nrow = 2, rel_heights = c(1, 1.5), axis = "rl", align = "v")

ggsave(file.path(subdir, "host_microbiome_distance_correlation_same_diet.png"), p, width=5, height=8)

#### Tree congruence ####

# Create a list of samples per species
links <- data.frame(phy_sp_f@sam_data) %>% 
        select(Species) %>%
        mutate(Species = recode(Species, "Sus scrofa domesticus" ="Sus scrofa")) %>%
        rename("phy1" = "Species") %>%
        rownames_to_column("phy2") %>% select(phy1, phy2)

links$colour <- species_palette[as.character(links$phy1)]

# Remove tips that are not in the metadata
host_consensus <- drop.tip(host_consensus, tip = setdiff(host_consensus$tip.label, links$phy1))

# Microbiome tree from CLR distances
mb_tree_clr <- nj(as.dist(mb_dist_clr))

# Cophyloplot
coph <- cophylo(tr1=host_consensus, tr2=mb_tree_clr, assoc=links)

png(file.path(subdir, "cophyloplot_clr.png"), width = 2200, height = 2200)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.type="curved",link.lwd=4, link.lty="solid", link.col=links$colour)
# Add title
#text(x = 0.5, y = 0.95, labels = paste("p-value:", format(cod_results$p.value[i], digits = 3)), pos = 2, cex = 1.5, col = "black")
dev.off()

# Microbiome tree from philr distances
mb_tree_philr <- nj(as.dist(mb_dist_philr))

# Cophyloplot
coph <- cophylo(tr1=mb_tree_philr, tr2=host_consensus, assoc=links, rotate = TRUE)

png(file.path(subdir, "cophyloplot_philr.png"), width = 2200, height = 2200)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, tr1=mb_tree_philr, link.type="curved",link.lwd=4, link.lty="solid", link.col=links$colour)
# Add title
#text(x = 0.5, y = 0.95, labels = paste("p-value:", format(cod_results$p.value[i], digits = 3)), pos = 2, cex = 1.5, col = "black")
dev.off()
