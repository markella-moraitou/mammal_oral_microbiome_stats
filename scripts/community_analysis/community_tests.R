##### COMMUNITY-WIDE TESTS #####

#### Assess which factors drive microbial community composition using various tests ####

################
#### SET UP ####
################
library(dplyr)
library(tidyr)
library(reshape2)
library(tibble)
library(phyloseq)
library(ape)
library(picante)
#library(microbiome)
library(stringr)
library(vegan)
library(ecodist)

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis"))
subdir <- normalizePath(file.path(outdir, "community_tests"))
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

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

# Extract OTU table and sample data
otu_table <- as.data.frame(phy_sp_f@otu_table)
sample_data <- as.data.frame(phy_sp_f@sam_data)

# Transpose the OTU table
otu_table_t <- t(otu_table)

# Host phylogeny
host_consensus <- read.tree(file.path(outdir, "host_consensus.tre"))

# Bacterial phylogeny (GTDB)
bac_tree <- read.tree(file.path(phydir, "phy_tree.tree"))

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

#### Presence-absence
set.seed(123)
perm <- adonis2(t(otu_table(phy_sp_f)) ~ order + diet + habitat + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "jaccard")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_pa_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(t(otu_table(phy_sp_f)) ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "jaccard")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_pa_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

#### CLR-transformed data
set.seed(123)

perm <- adonis2(t(otu_table(phy_sp_f_clr)) ~ order + diet + habitat + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_clr_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(t(otu_table(phy_sp_f_clr)) ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_clr_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

#### PhILR-transformed data
set.seed(123)

perm <- adonis2(otu_table(phy_sp_philr) ~ order + diet + habitat + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_philr_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(otu_table(phy_sp_philr) ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_philr_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

#########################################
#### MULTIPLE REGRESSION ON MATRICES ####
#########################################

#### Get phylogenetic distances between samples ####

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
phy_dist_melt <- host_dist_melt %>%
                    right_join(samples_to_species, by=c("Item1"="Species"), relationship = "many-to-many") %>%
                    right_join(samples_to_species, by=c("Item2"="Species"), relationship = "many-to-many", suffix = c("1","2")) %>%
                    select(Sample1, Sample2, host_distance)

# Get phylogenetic distance matrix
phy_dist <- dcast(phy_dist_melt, Sample1 ~ Sample2, value = "host_distance") %>% column_to_rownames("Sample1") %>% as.dist

#### Get the rest of explanatory matrices ####

## Diet distances
diet <- sample_data %>% data.frame %>% select(cp, cf, nfe, ee, ash)
diet_dist <- vegdist(diet, method = "euclidean")

## Habitat distances: 0 if the same, 1 if different
habitat <- sample_data$habitat.general
habitat_dist <- as.dist(outer(habitat, habitat, FUN = Vectorize(function(x, y) as.numeric(x != y))))

## Ruminant or non-ruminant distances: 0 if the same, 1 if different
ruminant <- sample_data$digestion == "Ruminant"
ruminant_dist <- as.dist(outer(ruminant, ruminant, FUN = Vectorize(function(x, y) as.numeric(x != y))))

## Species distances: 0 if the same, 1 if different
species <- sample_data$Species
species_dist <- as.dist(outer(species, species, FUN = Vectorize(function(x, y) as.numeric(x != y))))

#### Get microbiome distances ####

# Fix tree labels
#tree <- bac_tree
#tree$tip.label <- gsub("_", " ", tree$tip.label) %>% gsub("*", "", .)

# Get OTU table only of Bacteria
#phy_bac <- phy_sp_f %>% subset_taxa(superkingdom == "Bacteria")
#taxa_names(phy_bac) <-  gsub("*", "", taxa_names(phy_bac), fixed=TRUE)
#otu_table_bac <- phy_sp_f %>% subset_taxa(superkingdom == "Bacteria") %>% otu_table

#length(rownames(otu_table_bac))
#length(tree$tip.label)
#length(intersect(rownames(otu_table_bac), tree$tip.label))
#setdiff(rownames(otu_table_bac), tree$tip.label)
#setdiff(tree$tip.label,rownames(otu_table_bac))

# Get phyloseq objects at different taxonomic levels
for (tax in c("species", "genus", "family", "order", "class", "phylum")) {
    dist_list <- list()
    # For species level, keep as is
    if (tax == "species") {
        phy <- phy_sp_f
        phy_clr <- phy_sp_f_clr
        #phy_philr <- phy_sp_philr # Skip philr for now until I understand how to change the taxrank
    } else {
        # For unnormalised and CLR-normalised, agglomerate at the given taxonomic level
        cat("Aggregating to", tax, "\n")
        phy <- tax_glom(phy_sp_f, taxrank = tax)
        phy_clr <- tax_glom(phy_sp_f_clr, taxrank = tax)
        # For philr, we need to get the phylogenetic tree at the given taxonomic level
        #tax_table <- as.data.frame(phy_sp_f@tax_table)[,c("species", tax)]
        #tree_new <- tree
        # Collapse all trees to the given taxonomic level
        #for (taxname in unique(tax_table[,tax])) {
        #    tips_to_merge <- tax_table$species[tax_table[,tax] == taxname]
        #    mrca <- getMRCA(tree_new, tips_to_merge)
        #    if (length(tips_to_merge) > 1) {
        #        tree_new <- collapse.singles(tree_new, tips_to_merge)
        #    } else {
        #        tree_new$tip.label[tree_new$tip.label == tips_to_merge] <- taxname
        #    }
        #}
    }
    # Calculate distances
    cat("Calculating distances for", tax, "\n")
    ## Jaccard distances (presence/absence)
    mb_dist_jaccard <- vegdist(t(phy@otu_table), method="jaccard")
    ## Aitchison distances (accounts for abundances and compositionality)
    mb_dist_clr <- vegdist(t(phy_clr@otu_table), method="euclidean")
    ## Euclidean philr distances (accounts for phylogeny, abundances, and compositionality)
    #mb_dist_philr <- vegdist(phy_sp_philr@otu_table, method="euclidean")
    ## Get richness differences
    richness <- estimate_richness(phy, measures = "Observed") %>% rownames_to_column("Sample")
    mb_dist_richness <- vegdist(richness$Observed, method="euclidean")
    ## Faith's PD (only for bacteria)
    #pd(otu_table, tree, include.root=TRUE)
    ### CONTINUE HERE .......###
    dist_list[["jaccard"]] <- mb_dist_jaccard
    dist_list[["clr"]] <- mb_dist_clr
    #dist_list[["philr"]] <- mb_dist_philr
    dist_list[["richness"]] <- mb_dist_richness
    assign(paste0("mb_dist_", tax), dist_list)
}

#### Run MRM ####
set.seed(123)

mrm_results <- list()
for (dist_list in c("mb_dist_species", "mb_dist_genus", "mb_dist_family", "mb_dist_order", "mb_dist_class", "mb_dist_phylum")) {
    tax = str_remove(dist_list, "mb_dist_")
    dist_list <- get(dist_list)
    for (dist_name in names(dist_list)) {
        mb_dist <- dist_list[[dist_name]]
        cat("Running MRM for", dist_name, "in", tax, "level", "\n")
        mrm <- MRM(mb_dist ~ phy_dist + diet_dist + habitat_dist + ruminant_dist, nperm = 1000)
        df <- mrm$coef %>% as.data.frame() %>% rownames_to_column("Variable")
        df$Taxonomic_level <- tax
        df$Distance <- dist_name
        mrm_results[[paste(tax, dist_name, sep="_")]] <- df
    }
}

# Merge all dataframes
mrm_results_df <- bind_rows(mrm_results) %>% filter(Variable != "Int") %>%
                    mutate(sig = case_when(pval < 0.001 ~ "***",
                                           pval < 0.01 ~ "**",
                                           pval < 0.05 ~ "*",
                                           TRUE ~ "")) %>%
                    mutate(Variable = str_remove(Variable, "_dist")) %>%
                    rename(coef = mb_dist) %>%
                    mutate(Taxonomic_level = factor(Taxonomic_level, levels = c("species", "genus", "family", "order", "class", "phylum"))) %>%
                    mutate(Variable = factor(Variable, levels = c("phy", "diet", "ruminant", "habitat")))

write.csv(mrm_results_df, file = file.path(subdir, "mrm_results_without_species.csv"), row.names = FALSE, quote = TRUE)

#### Plots ####
var_cols <- c("species" = "#A41B1B",
              "phy" = "#E35D5D",
              "diet" = "#2C9F2C",
              "ruminant" = "#C78837",
              "habitat" = "#295B7F")

# Plot result summary as boxplot
p <- ggplot(mrm_results_df,
            aes(x = Variable, y = coef, fill = Variable)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = var_cols) +
    geom_text(aes(label = sig), size = 6) +
    facet_grid(Distance ~ Taxonomic_level, scales = "free_y") +
    labs(y = "MRM coefficient", x = "") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p, filename = file.path(subdir, "mrm_results_barplot_without_species.png"), width = 10, height = 8)

# Plot distance correlations for species-level distances #
dist_df <-
    data.frame(jaccard = as.vector(mb_dist_jaccard),
               clr = as.vector(mb_dist_clr),
               philr = as.vector(mb_dist_philr),
               richness = as.vector(mb_dist_richness),
               phy = as.vector(phy_dist),
               diet = as.vector(diet_dist),
               habitat = as.vector(habitat_dist))

dist_df_melt <- melt(dist_df, measure.vars = c("jaccard", "clr", "philr", "richness"),
                     variable.name = "mb_distance", value.name = "mb_dist_value") %>%
                     mutate(habitat = as.factor(habitat))

p <- ggplot(dist_df_melt, aes(y = mb_dist_value, x = phy, colour = diet, shape = habitat)) +
    geom_jitter(size = 0.5, alpha = 0.1, width = 5, height = 0) +
    scale_color_viridis_c() +
    facet_wrap(~ mb_distance, scales = "free_y")

ggsave(p, filename = file.path(subdir, "mrm_distances.png"), width = 8, height = 6)
