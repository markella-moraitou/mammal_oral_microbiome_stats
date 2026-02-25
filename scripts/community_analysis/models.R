##### COMMUNITY-WIDE MODELS #####

#### Assess which factors drive microbial community composition using various models ####

################
#### SET UP ####
################

library(dplyr)
library(tidyr)
library(reshape2)
library(tibble)
library(phyloseq)
library(philr)
library(TreeTools)
library(ape)
library(picante)
library(phytools)
library(phangorn)
library(microbiome)
library(stringr)
library(vegan)
library(ecodist)
library(conflicted)

# Explicitly resolve the conflict
conflict_prefer("as.phylo", "ape")

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis"))
subdir <- normalizePath(file.path(outdir, "model_results"))
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

source(file.path("..", "phylo_functions.R"))

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
#### PREP DATA ####
###################

#### Get phylogenetic distances between samples ####

# Change tip labels to match metadata
host_consensus$tip.label <- host_consensus$tip.label %>%
                            gsub(pattern="Procolobus_badius", replacement="Piliocolobus_foai") %>%
                            gsub(pattern="Otaria_bryonia", replacement="Otaria_flavescens") %>%
                            gsub(pattern="_", replacement=" ", fixed=TRUE)

# Species distances
host_dist_melt <- cophenetic(host_consensus) %>% as.data.frame() %>% rownames_to_column("Item1") %>%
                melt(idvar = "Item1") %>% rename(Item2 = variable, host_distance = value) %>%
                # Scale distance
                mutate(host_distance = scale(host_distance))

# Samples to species
samples_to_species <- phy_sp_f@sam_data %>% data.frame %>% select(Species) %>% 
                      mutate(Species = recode(Species, "Sus domesticus" ="Sus scrofa")) %>%
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
diet <- sample_data %>% data.frame %>% select(Animal, PlantO, Fruit, Seed)
diet_dist <- scale(vegdist(diet, method = "euclidean"))

## Habitat distances: 0 if the same, 1 if different
habitat <- sample_data$habitat.general
habitat_dist <- outer(habitat, habitat, FUN = Vectorize(function(x, y) as.numeric(x != y)))
rownames(habitat_dist) <- rownames(sample_data)
colnames(habitat_dist) <- rownames(sample_data)
habitat_dist <- as.dist(habitat_dist)

## Ruminant or non-ruminant distances: 0 if the same, 1 if different
ruminant <- sample_data$digestion == "Ruminant"
ruminant_dist <- outer(ruminant, ruminant, FUN = Vectorize(function(x, y) as.numeric(x != y)))
rownames(ruminant_dist) <- rownames(sample_data)
colnames(ruminant_dist) <- rownames(sample_data)
ruminant_dist <- as.dist(ruminant_dist)

## Species distances: 0 if the same, 1 if different
species <- sample_data$Species
species_dist <- outer(species, species, FUN = Vectorize(function(x, y) as.numeric(x != y)))
rownames(species_dist) <- rownames(sample_data)
colnames(species_dist) <- rownames(sample_data)
species_dist <- as.dist(species_dist)

##########################
#### DEFINE FUNCTIONS ####
##########################

# Function to get microbiome distances at different taxonomic levels
# and for different distance metrics

collapse_tree_to_rank <- function(tree, tax_table, tax) {
    taxa <- unique(tax_table[,tax])
    tree_coll <- tree
    for (i in 1:length(taxa)) {
        taxname <- taxa[i]
        tips_to_merge <- tax_table$species[tax_table[,tax] == taxname] %>% intersect(tree_coll$tip.label)
        if (length(tips_to_merge) > 1) {
            # Identify nodes to collapse
            mrca <- getMRCA(tree_coll, tips_to_merge)
            cat(paste("Collapsing", length(tips_to_merge), "tips for", tax, taxname, "\n"))
            # Pick a random tip to keep and collapse others into it
            n <- sample(1:length(tips_to_merge), 1)
            tip_to_keep <- tips_to_merge[n]  
            tips_to_drop <- tips_to_merge[-n]
            tree_coll <- drop.tip(tree_coll, tips_to_drop)
            # Rename the kept tip to the taxonomic name
            tree_coll$tip.label[tree_coll$tip.label == tip_to_keep] <- taxname
        }  else if (length(tips_to_merge) == 1) {
            # If there is only one tip, just rename it
            tree_coll$tip.label[tree_coll$tip.label == tips_to_merge] <- taxname
        }
    }
    return(tree_coll)
}

# Write function to get microbiome distances at different taxonomic levels
get_mb_dist <- function(phy, tax, method = c("jaccard", "aitchison", "philr"), tree = NULL, scale = TRUE) {
    # Agglomerate to higher taxonomic level if needed
    if (tax == "species") {
        phy_raw <- phy
        collapsed_tree <- tree
    } else {
        cat("Agglomerating to", tax, "level\n")
        phy_raw <- tax_glom(phy, taxrank = tax)
        phy_raw@tax_table[, tax] <- make.unique(phy_raw@tax_table[, tax])
        taxa_names(phy_raw) <- phy_raw@tax_table[, tax]
        # If tree has been provided, collapse it to the given taxonomic level
        if (!is.null(tree)) {
            cat("Collapsing tree to", tax, "level\n")
            tax_table <- as.data.frame(phy@tax_table)[,c("species", tax)]
            collapsed_tree <- collapse_tree_to_rank(tree, tax_table, tax)
        }
    }
    if ("aitchison" %in% method) {
        # CLR normalisation
        phy_aitchison <- microbiome::transform(phy_raw, "clr")
    }
    if ("philr" %in% method) {
        # Keep only tips present in dataset
        tips_in_data <- intersect(phy_raw@tax_table[, tax], collapsed_tree$tip.label)
        tree_filt <- collapsed_tree %>% drop.tip(setdiff(collapsed_tree$tip.label, tips_in_data))
        # Philr normalisation
        phy_philr <- phy_raw %>% prune_taxa(tips_in_data, .)
        phy_philr@phy_tree <- tree_filt
        phy_philr <- philr(phy_philr, pseudocount=10^-5)
    }
    # Calculate distances
    dist_list <- list("raw" = "jaccard",
                      "aitchison" = "euclidean",
                      "philr" = "euclidean")
    mbdist_list <- list()
    if ("jaccard" %in% method) {
        mbdist <- vegdist(t(phy_raw@otu_table), method="jaccard")
        (scale == TRUE) & (mbdist <- as.dist(scale(mbdist))) # Scale microbiome distances
        mbdist_list[["jaccard"]] <- mbdist
    } 
    if ("aitchison" %in% method) {
        mbdist <- vegdist(t(phy_aitchison@otu_table), method="euclidean")
        (scale == TRUE) & (mbdist <- as.dist(scale(mbdist))) # Scale microbiome distances
        mbdist_list[["aitchison"]] <- mbdist
    }
    if ("philr" %in% method) {
        mbdist <- vegdist(phy_philr, method="euclidean")
        (scale == TRUE) & (mbdist <- as.dist(scale(mbdist))) # Scale microbiome distances
        mbdist_list[["philr"]] <- mbdist
    }
    return(mbdist_list)
}

# Write function to run MRM with random subsampling of samples per species
sampling_mrm <- function(mb_dist, predictor_dists = list(), niter = 100, sample_size = 5, sample_map) {
    mrm_res <- list()
    for (iter in 1:niter) {
        # Randomly subsamples 5 samples per species
        sampled_samples <- sample_map %>%
                            group_by(Species) %>%
                            slice_sample(n = sample_size) %>%
                            pull(Sample)
        # Subset microbiome distances to include only selected samples
        mb_dist_filt <- as.dist(as.matrix(mb_dist)[sampled_samples, sampled_samples])
        # Subset predictor distances to include only selected samples
        formula_parts <- c()
        for (matrix_name in names(predictor_dists)) {
            matrix <- get(matrix_name)
            matrix_filt <- as.dist(as.matrix(matrix)[sampled_samples, sampled_samples])
            assign(paste0(matrix_name, "_filt"), matrix_filt)
            if (all(as.matrix(matrix_filt) == 0)) {
                next
            } else {
                formula_parts <- c(formula_parts, paste0(matrix_name, "_filt"))
            }
        }
        formula_str <- as.formula(paste("mb_dist_filt ~", paste(formula_parts, collapse = " + ")))
        mrm <- MRM(formula_str, nperm = 1000)
        df <- mrm$coef %>% as.data.frame() %>% rownames_to_column("Variable") %>%
                mutate(Variable = str_remove(Variable, "_dist_filt"))
        df$iteration <- iter
        mrm_res[[iter]] <- df
    }
    mrm_df <- bind_rows(mrm_res)
    # Summarise to get min, q1, median, q3, max of the coefficients and median p-values
    mrm_summ <- mrm_df %>%
                group_by(Variable) %>%
                summarise(mb_dist_min = min(mb_dist_filt),
                        mb_dist_q1 = quantile(mb_dist_filt, 0.25),
                        mb_dist_median = median(mb_dist_filt),
                        mb_dist_q3 = quantile(mb_dist_filt, 0.75),
                        mb_dist_max = max(mb_dist_filt),
                        pval = median(pval), .groups = "drop")
    return(mrm_summ)
}

###################################
#### Get distances and run MRM ####
###################################

# Calculate microbiome distances at different taxonomic levels and for different distance metrics

tax_ranks <- c("species", "genus", "family", "order", "class", "phylum")

dist_list <- list()
for (tax in tax_ranks) {
   cat("Calculating distances at", tax, "level\n")
    method = c("jaccard", "aitchison", "philr")
    tree = bac_tree
    tree$tip.label <- gsub("_", " ", tree$tip.label)
    dist <- get_mb_dist(phy_sp_f, tax, method, tree)
    dist_list[[tax]] <- dist
}

saveRDS(dist_list, file = file.path(subdir, "mb_distances_all_levels.RDS"))

# Run MRM for each taxonomic level and each dataset, for each distance metric
# randomly sampling 5 samples per species

set.seed(123)

distances_wo_species <- list("phy_dist" = phy_dist,
                            "diet_dist" = diet_dist,
                            "habitat_dist" = habitat_dist,
                            "ruminant_dist" = ruminant_dist)

distances_with_species <- list(
                            "species_dist" = species_dist,
                            "phy_dist" = phy_dist,
                            "diet_dist" = diet_dist,
                            "habitat_dist" = habitat_dist,
                            "ruminant_dist" = ruminant_dist)

phy_list <- c("phy_sp_f", "phy_artio", "phy_prim", "phy_habitat")

mrm_results_nospecies <- list()
mrm_results_species <- list()

for(i in 1:length(dist_list)) {
    for (dataset in phy_list) {
        # For each dataset for each taxonomic level
        tax <- names(dist_list)[i] # Taxonomic level
        distances <- dist_list[[i]] # List of distance matrices (jaccard, aitchison, philr)
        phy <- get(dataset) # Phyloseq object
        sample_species_map <- data.frame(phy@sam_data) %>% rownames_to_column("Sample") %>%
                              select(Sample, Species) # Dataframe mapping Sample to Species
        for (method in names(distances)) {
            cat("Running MRM for", method, "in", tax, "level in", dataset, "without species distance matrix\n")
            mb_dist <- distances[[method]]
            df <- sampling_mrm(mb_dist, predictor_dists = distances_wo_species, niter = 100, sample_map = sample_species_map)
            df$Taxonomic_level <- tax
            df$Distance <- method
            df$Dataset <- dataset
            mrm_results_nospecies[[paste(dataset, tax, method, sep="-")]] <- df
            cat("Running MRM for", method, "on", tax, "level in", dataset, "with species distance matrix\n")
            df <- sampling_mrm(mb_dist, predictor_dists = distances_with_species, niter = 100, sample_map = sample_species_map)
            df$Taxonomic_level <- tax
            df$Distance <- method
            df$Dataset <- dataset
            mrm_results_species[[paste(dataset, tax, method, sep="-")]] <- df
        }
    }
}

# Display results without species matrix
mrm_results_nospecies_df <- bind_rows(mrm_results_nospecies) %>% dplyr::filter(Variable != "Int") %>%
                    mutate(sig = case_when(pval < 0.05 ~ "*",
                                           TRUE ~ "")) %>%
                    mutate(Distance = factor(recode(Distance,
                                            "jaccard" = "Jaccard",
                                            "aitchison" = "Aitchison",
                                            "philr" = "Philr"),
                                            levels = c("Jaccard", "Aitchison", "Philr"))) %>%
                    mutate(Taxonomic_level = factor(Taxonomic_level, levels = c("species", "genus", "family", "order", "class", "phylum"))) %>%
                    mutate(Dataset = factor(recode(Dataset,
                                            "phy_sp_f" = "Entire dataset",
                                            "phy_artio" = "Artiodactyla",
                                            "phy_prim" = "Primates",
                                            "phy_carni" = "Carnivora",
                                            "phy_habitat" = "Marine v. Terrest."),
                                            levels =c("Entire dataset", "Artiodactyla", "Carnivora", "Primates", "Marine v. Terrest."))) %>%
                    mutate(Variable = str_remove(Variable, "_dist_filt")) %>%
                    mutate(Variable = factor(recode(Variable,
                                            "phy" = "Host phylogeny",
                                            "diet" = "Diet",
                                            "ruminant" = "Ruminant status",
                                            "habitat" = "Habitat"),
                                            levels = c("Host phylogeny", "Diet",  "Ruminant status", "Habitat")))

colnames(mrm_results_nospecies_df) <- gsub('mb_dist', 'coef', colnames(mrm_results_nospecies_df), fixed=TRUE)

write.csv(mrm_results_nospecies_df, file = file.path(subdir, "mrm_results_without_species.csv"), row.names = FALSE, quote = TRUE)

mrm_results_nospecies_df <- mrm_results_nospecies_df %>% dplyr::filter(Taxonomic_level %in% c("species", "genus", "order", "phylum"))

var_cols <- c("Species identity" = "#A41B1B",
              "Host phylogeny" = "#E35D5D",
              "Diet" = "#2C9F2C",
              "Ruminant status" = "#C78837",
              "Habitat" = "#295B7F")

# Plot result summary as barplot
p <- ggplot(mrm_results_nospecies_df,
            aes(x = Variable, y = coef_median, fill = Variable, group = Taxonomic_level)) +
    geom_bar(aes(alpha = Taxonomic_level), stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = coef_q1, ymax = coef_q3), position = position_dodge(width = 0.9), width = 0.2) +
    scale_fill_manual(values = var_cols) +
    scale_alpha_discrete(range = c(0.4, 1), name = "Microbial taxonomic level") +
    geom_text(aes(label = sig, y = ifelse(coef_median > 0, coef_q3 + 0.1, coef_q1 - 0.1)),
                  position = position_dodge(width = 0.9), size = 4) +
    facet_grid(Distance ~ Dataset, scale = "free", space = "free") +
    labs(y = "MRM coefficient", x = "") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(hjust = 1),
          panel.grid = element_blank()) +
    guides(fill = "none")

ggsave(p, filename = file.path(subdir, "mrm_results_without_species_barplot.png"), width = 8, height = 10)

# Display results with species matrix
mrm_results_species_df <- bind_rows(mrm_results_species) %>% dplyr::filter(Variable != "Int") %>%
                    mutate(sig = case_when(pval < 0.05 ~ "*",
                                           TRUE ~ "")) %>%
                    mutate(Distance = factor(recode(Distance,
                                            "jaccard" = "Jaccard",
                                            "aitchison" = "Aitchison",
                                            "philr" = "Philr"),
                                            levels = c("Jaccard", "Aitchison", "Philr"))) %>%
                    mutate(Taxonomic_level = factor(Taxonomic_level, levels = c("species", "genus", "family", "order", "class", "phylum"))) %>%
                    mutate(Dataset = factor(recode(Dataset,
                                            "phy_sp_f" = "Entire dataset",
                                            "phy_artio" = "Artiodactyla",
                                            "phy_prim" = "Primates",
                                            "phy_carni" = "Carnivora",
                                            "phy_habitat" = "Marine v. Terrest."),
                                            levels =c("Entire dataset", "Artiodactyla", "Carnivora", "Primates", "Marine v. Terrest."))) %>%
                    mutate(Variable = str_remove(Variable, "_dist_filt")) %>% 
                    mutate(Variable = factor(recode(Variable,
                                            "species" = "Species identity",
                                            "phy" = "Host phylogeny",
                                            "diet" = "Diet",
                                            "ruminant" = "Ruminant status",
                                            "habitat" = "Habitat"),
                                            levels = c("Species identity", "Host phylogeny", "Diet",  "Ruminant status", "Habitat")))

colnames(mrm_results_species_df) <- gsub('mb_dist', 'coef', colnames(mrm_results_species_df), fixed=TRUE)

write.csv(mrm_results_species_df, file = file.path(subdir, "mrm_results_with_species.csv"), row.names = FALSE, quote = TRUE)

mrm_results_species_df <- mrm_results_species_df %>% dplyr::filter(Taxonomic_level %in% c("species", "genus", "order", "phylum"))

# Plot result summary as barplot
p <- ggplot(mrm_results_species_df,
            aes(x = Variable, y = coef_median, fill = Variable, group = Taxonomic_level)) +
    geom_bar(aes(alpha = Taxonomic_level), stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = coef_q1, ymax = coef_q3), position = position_dodge(width = 0.9), width = 0.2) +
    scale_fill_manual(values = var_cols) +
    scale_alpha_discrete(range = c(0.4, 1), name = "Microbial taxonomic level") +
    geom_text(aes(label = sig, y = ifelse(coef_median > 0, coef_q3 + 0.1, coef_q1 - 0.1)),
                  position = position_dodge(width = 0.9), size = 4) +
    facet_grid(Distance ~ Dataset, scale = "free", space = "free") +
    labs(y = "MRM coefficient", x = "") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(hjust = 1),
          panel.grid = element_blank()) +
    guides(fill = "none")

ggsave(p, filename = file.path(subdir, "mrm_results_with_species_barplot.png"), width = 8, height = 10)

#########################
#### Tree congruence ####
#########################

## Prep input
samples_to_species$Common.name <- sample_data$Common.name[match(samples_to_species$Species, sample_data$Species)]

# Create a list of samples per species
links <- samples_to_species %>%
        mutate(colour = species_palette[as.character(Species)]) %>%
        rename("phy1" = "Common.name", "phy2" = "Sample") %>%
        select(phy1, phy2, colour)

# Microbiome tree from Jaccard distances
mb_dist_jaccard <- vegdist(t(phy_sp_f@otu_table), method="jaccard")
mb_tree_jaccard <- nj(mb_dist_jaccard)
mb_dist_jaccard_matrix <- as.matrix(mb_dist_jaccard)

# Microbiome tree from CLR distances
mb_dist_clr <- vegdist(t(phy_sp_f_clr@otu_table), method="euclidean")
mb_tree_clr <- nj(mb_dist_clr)
mb_dist_clr_matrix <- as.matrix(mb_dist_clr)

# Microbiome tree from philr distances
mb_dist_philr <- vegdist(phy_sp_philr@otu_table, method="euclidean")
mb_tree_philr <- nj(mb_dist_philr)
mb_dist_philr_matrix <- as.matrix(mb_dist_philr)

# Get host distances per host species
host_dist <- dcast(host_dist_melt, Item1 ~ Item2, value = "host_distance") %>% column_to_rownames("Item1")

colnames(host_dist) <- sample_data$Common.name[match(colnames(host_dist), sample_data$Species)]
rownames(host_dist) <- sample_data$Common.name[match(rownames(host_dist), sample_data$Species)]

host_dist <- as.dist(host_dist)

# Get diet matrix again, but per host species
diet <- sample_data %>% data.frame %>% select(Common.name, Animal, PlantO, Fruit, Seed) %>%
        dplyr::filter(Common.name %in% links$phy1) %>%
        as_tibble %>% unique %>%
        column_to_rownames("Common.name")

diet_sp_dist <- vegdist(diet, method = "euclidean")

#### PARAFIT ####
## Compare host phylogeny and microbiome diversification
# Create association matrix (binary matrix: rows = hosts, columns = microbes)

parafit_results <- data.frame(mb_distance = character(),
                              host_distance = character(),
                              ParaFitGlobal = numeric(),
                              p.global = numeric())

assoc <- links %>% select(phy1, phy2) %>% mutate(assoc = 1) %>% pivot_wider(names_from = phy2, values_from = assoc, values_fill = 0) %>%
          column_to_rownames("phy1") %>% as.matrix()

# Jaccard
res <- parafit(host.D = as.matrix(host_dist), para.D = mb_dist_jaccard_matrix, HP = assoc, nperm = 1000, seed = 123, correction="cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "Jaccard",
                                    host_distance = "Phylogeny",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

# CLR
res <- parafit(host.D = as.matrix(host_dist), para.D = mb_dist_clr_matrix, HP = assoc, nperm = 1000, seed = 123, correction="cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "CLR",
                                    host_distance = "Phylogeny",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

# PhILR
res <- parafit(host.D = as.matrix(host_dist), para.D = mb_dist_philr_matrix, HP = assoc, nperm = 1000, seed = 123, correction="cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "PhILR",
                                    host_distance = "Phylogeny",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

## Compare host diet and microbiome diversification

# Jaccard
res <- parafit(host.D = as.matrix(diet_sp_dist), para.D = mb_dist_jaccard_matrix, HP = assoc, nperm = 1000, seed = 123, correction = "cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "Jaccard",
                                    host_distance = "Diet",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

# CLR
res <- parafit(host.D = as.matrix(diet_sp_dist), para.D = mb_dist_clr_matrix, HP = assoc, nperm = 1000, seed = 123, correction = "cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "CLR",
                                    host_distance = "Diet",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

# PhILR
res <- parafit(host.D = as.matrix(diet_sp_dist), para.D = mb_dist_philr_matrix, HP = assoc, nperm = 1000, seed = 123, correction = "cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "PhILR",
                                    host_distance = "Diet",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

write.csv(parafit_results, file = file.path(subdir, "parafit_results.csv"), row.names = FALSE, quote = TRUE)

#### Plot cophyloplots ####
# Cophyloplot Jaccard

host_consensus$tip.label <- sample_data$Common.name[match(host_consensus$tip.label, sample_data$Species)]

coph <- cophylo(tr1=host_consensus, tr2=mb_tree_jaccard, assoc=links)

png(file.path(subdir, "cophyloplot_jaccard.png"), width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()

coph <- cophylo(tr1=nj(diet_sp_dist), tr2=mb_tree_jaccard, assoc=links)

png(file.path(subdir, "cophyloplot_jaccard_diet.png"), width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()

# Cophyloplot CLR
coph <- cophylo(tr1=host_consensus, tr2=mb_tree_clr, assoc=links)

png(file.path(subdir, "cophyloplot_clr.png"), width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()

coph <- cophylo(tr1=nj(diet_sp_dist), tr2=mb_tree_clr, assoc=links)

png(file.path(subdir, "cophyloplot_clr_diet.png"), width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()

# Cophyloplot PhILR
coph <- cophylo(tr1=host_consensus, tr2=mb_tree_philr, assoc=links, rotate = TRUE)

png(file.path(subdir, "cophyloplot_philr.png"),width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, tr1=mb_tree_philr, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()

coph <- cophylo(tr1=nj(diet_sp_dist), tr2=mb_tree_philr, assoc=links)

png(file.path(subdir, "cophyloplot_philr_diet.png"), width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()

### Rerun tests with one microbiome per species

n_iter <- 100
tree_res <- data.frame(ParaFitGlobal = numeric(n_iter), pval = numeric(n_iter), rf_dist = numeric(n_iter))

set.seed(123)

# Also plot the sampled trees
pdf(file.path(subdir, "cophyloplot_clr_subsets.pdf"))
par(mar=c(5, 4, 4, 2) + 0.1)

for (i in 1:n_iter) {
  cat(i, "... ")
  # Sample one microbiome per species
  sampled_mbs <- links %>% group_by(phy1) %>% slice_sample(n = 1) %>% ungroup() %>% pull(phy2)
  # Create new association matrix
  links_sampled <- links %>% dplyr::filter(phy2 %in% sampled_mbs)
  assoc_sampled <- links_sampled %>%
                   select(phy1, phy2) %>% mutate(assoc = 1) %>%
                   pivot_wider(names_from = phy2, values_from = assoc, values_fill = 0) %>%
                   column_to_rownames("phy1") %>% as.matrix()
  # Subset microbiome dist
  mb_dist_clr_matrix_sampled <- mb_dist_clr_matrix[sampled_mbs, sampled_mbs]
  # Run parafit
  res <- parafit(host.D = host_dist, para.D = mb_dist_clr_matrix_sampled, HP = assoc_sampled, nperm = 1000, correction="cailliez")
  # Also get Robinson-Foulds distance
  mb_tree_clr_sampled <- drop.tip(mb_tree_clr, tip = setdiff(mb_tree_clr$tip.label, sampled_mbs))
  copy <- mb_tree_clr_sampled
  copy$tip.label <- links$phy1[match(copy$tip.label, links$phy2)]
  rf_dist <- RF.dist(host_consensus, copy, normalize = TRUE)
  tree_res[i, "ParaFitGlobal"] <- res$ParaFitGlobal
  tree_res[i, "pval"] <- res$p.global
  tree_res[i, "rf_dist"] <- rf_dist
  # Cophyloplot
  coph <- cophylo(tr1=host_consensus, tr2=mb_tree_clr_sampled, assoc=links_sampled, rotate = TRUE)
  plot(coph, tr1=mb_tree_philr, link.lwd=4, link.lty="solid", link.col=links_sampled$colour)
  title(main = paste("Iteration", i, "- p =", round(res$p.global, 3), "- RF dist =", round(rf_dist, 3)), 
        cex.main = 0.5, col.main = "black")
}
dev.off()
cat("done\n")

write.csv(tree_res, file = file.path(subdir, "parafit_results_subsampled.csv"), row.names = FALSE, quote = TRUE)