##### COMMUNITY-WIDE MODELS #####

#### Assess which factors drive microbial functional composition using various models ####

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
outdir <- normalizePath(file.path("..", "..", "output", "function"))
subdir <- normalizePath(file.path(outdir, "model_results"))
datadir <- normalizePath(file.path(outdir, "data")) # Directory with phyloseq objects

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
phy_gene_f <- readRDS(file.path(datadir, "phy_gene_f.RDS"))
phy_gene_f_clr <- readRDS(file.path(datadir, "phy_gene_f_clr.RDS"))

phy_pathway <- readRDS(file.path(outdir, "pathway_completeness", "phy_pathway.RDS"))
phy_pathway_clr <- readRDS(file.path(outdir, "pathway_completeness", "phy_pathway_clr.RDS"))

phy_gifts_el <- readRDS(file.path(datadir, "phy_gifts_el.RDS"))

# Host phylogeny
host_consensus <- read.tree(file.path(outdir, "..", "community_analysis", "host_consensus.tre"))

###################
#### PREP DATA ####
###################

sample_data <- data.frame(phy_gene_f@sam_data)

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
samples_to_species <- phy_gene_f@sam_data %>% data.frame %>% select(Species) %>% 
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
diet <- sample_data %>% data.frame %>% select(Animalivory, PlantO, Frugivory, Seed)
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

# Write function to get microbiome distances at different taxonomic levels
get_mb_dist <- function(phy, method = c("jaccard", "aitchison", "euclidean"), scale = TRUE) {
    if ("aitchison" %in% method) {
        # CLR normalisation
        phy_clr <- microbiome::transform(phy, "clr")
    }
    # Calculate distances
    dist_list <- list("raw" = "jaccard",
                      "aitchison" = "euclidean")
    mbdist_list <- list()
    otu_table <- phy@otu_table
    if (taxa_are_rows(phy)) {
        otu_table <- t(otu_table)
    }
    if ("jaccard" %in% method) {
        mbdist <- vegdist(otu_table, method="jaccard")
        (scale == TRUE) & (mbdist <- as.dist(scale(mbdist))) # Scale microbiome distances
        mbdist_list[["jaccard"]] <- mbdist
    } 
    if ("aitchison" %in% method) {
        mbdist <- vegdist(otu_table, method="euclidean")
        (scale == TRUE) & (mbdist <- as.dist(scale(mbdist))) # Scale microbiome distances
        mbdist_list[["aitchison"]] <- mbdist
    }
    if ("euclidean" %in% method) {
        mbdist <- vegdist(otu_table, method="euclidean")
        (scale == TRUE) & (mbdist <- as.dist(scale(mbdist))) # Scale microbiome distances
        mbdist_list[["euclidean"]] <- mbdist
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

# Calculate jaccard and aitchison distances for genes and pathways

dist_list <- list()

dist_list[["gene"]] <- get_mb_dist(phy_gene_f, method = c("jaccard", "aitchison"))
dist_list[["pathway"]] <- get_mb_dist(phy_pathway, method = c("jaccard", "aitchison"))
dist_list[["gifts"]] <- get_mb_dist(phy_gifts_el, method = c("jaccard", "euclidean"))

saveRDS(dist_list, file = file.path(subdir, "mb_func_distances.RDS"))

# Run MRM for each genes and pathways, for each distance metric
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

phy_list <- c("phy_gene_f", "phy_pathway", "phy_gifts_el")

mrm_results_nospecies <- list()
mrm_results_species <- list()

for (i in 1:length(phy_list)) {
    dataset <- phy_list[i]
    distances <- dist_list[[i]] # List of distance matrices (jaccard, aitchison, philr)
    phy <- get(dataset) # Phyloseq object
    sample_species_map <- data.frame(phy@sam_data) %>% rownames_to_column("Sample") %>%
                          select(Sample, Species) # Dataframe mapping Sample to Species
    for (method in names(distances)) {
        cat("Running MRM for", method, "on", dataset, "without species distance matrix\n")
        mb_dist <- distances[[method]]
        df <- sampling_mrm(mb_dist, predictor_dists = distances_wo_species, niter = 100, sample_map = sample_species_map)
        df$Distance <- method
        df$Dataset <- dataset
        mrm_results_nospecies[[paste(dataset, method, sep="-")]] <- df
        cat("Running MRM for", method, "on", dataset, "with species distance matrix\n")
        df <- sampling_mrm(mb_dist, predictor_dists = distances_with_species, niter = 100, sample_map = sample_species_map)
        df$Distance <- method
        df$Dataset <- dataset
        mrm_results_species[[paste(dataset, method, sep="-")]] <- df
    }
}

# Display results without species matrix
mrm_results_nospecies_df <- bind_rows(mrm_results_nospecies) %>% dplyr::filter(Variable != "Int") %>%
                    mutate(sig = case_when(pval < 0.05 ~ "*",
                                           TRUE ~ "")) %>%
                    # Show Euclidean and Aitchison together
                    mutate(Distance = case_when(Distance %in% c("euclidean", "aitchison") ~ "Euclidean/Aitchison",
                                                TRUE ~ Distance)) %>%
                    mutate(Distance = factor(recode(Distance,
                                            "jaccard" = "Jaccard",
                                            "aitchison" = "Euclidean/Aitchison"),
                                            levels = c("Euclidean/Aitchison", "Jaccard"))) %>%
                    mutate(Dataset = factor(recode(Dataset,
                                            "phy_gene_f" = "Genes",
                                            "phy_pathway" = "KEGG Pathways",
                                            "phy_gifts_el" = "Element-level GIFTs"),
                                            levels =c("Genes", "KEGG Pathways", "Element-level GIFTs"))) %>%
                    mutate(Variable = str_remove(Variable, "_dist_filt")) %>%
                    mutate(Variable = factor(recode(Variable,
                                            "phy" = "Host phylogeny",
                                            "diet" = "Diet",
                                            "ruminant" = "Ruminant status",
                                            "habitat" = "Habitat"),
                                            levels = c("Host phylogeny", "Diet",  "Ruminant status", "Habitat")))

colnames(mrm_results_nospecies_df) <- gsub('mb_dist', 'coef', colnames(mrm_results_nospecies_df), fixed=TRUE)

write.csv(mrm_results_nospecies_df, file = file.path(subdir, "func_mrm_results_without_species.csv"), row.names = FALSE, quote = TRUE)

var_cols <- c("Species identity" = "#A41B1B",
              "Host phylogeny" = "#E35D5D",
              "Diet" = "#2C9F2C",
              "Ruminant status" = "#C78837",
              "Habitat" = "#295B7F")

# Plot result summary as barplot
p <- ggplot(mrm_results_nospecies_df,
            aes(x = Variable, y = coef_median, fill = Variable)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = coef_q1, ymax = coef_q3), position = position_dodge(width = 0.9), width = 0.2) +
    scale_fill_manual(values = var_cols) +
    geom_text(aes(label = sig, y = ifelse(coef_median > 0, coef_q3 * 1.1, coef_q1 * 1.1)),
                  position = position_dodge(width = 0.9), size = 4) +
    facet_grid(Distance ~ Dataset, scale = "free", space = "free_x") +
    labs(y = "MRM coefficient", x = "") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(hjust = 1),
          panel.grid = element_blank()) +
    guides(fill = "none")

ggsave(p, filename = file.path(subdir, "func_mrm_results_without_species_barplot.png"), width = 8, height = 10)

# Display results with species matrix
mrm_results_species_df <- bind_rows(mrm_results_species) %>% dplyr::filter(Variable != "Int") %>%
                    mutate(sig = case_when(pval < 0.05 ~ "*",
                                           TRUE ~ "")) %>%
                    # Show Euclidean and Aitchison together
                    mutate(Distance = case_when(Distance %in% c("euclidean", "aitchison") ~ "Euclidean/Aitchison",
                                                TRUE ~ Distance)) %>%
                    mutate(Distance = factor(recode(Distance,
                                            "jaccard" = "Jaccard",
                                            "aitchison" = "Euclidean/Aitchison"),
                                            levels = c("Euclidean/Aitchison", "Jaccard"))) %>%
                    mutate(Dataset = factor(recode(Dataset,
                                            "phy_gene_f" = "Genes",
                                            "phy_pathway" = "KEGG Pathways",
                                            "phy_gifts_el" = "Element-level GIFTs"),
                                            levels =c("Genes", "KEGG Pathways", "Element-level GIFTs"))) %>%
                    mutate(Variable = str_remove(Variable, "_dist_filt")) %>% 
                    mutate(Variable = factor(recode(Variable,
                                            "species" = "Species identity",
                                            "phy" = "Host phylogeny",
                                            "diet" = "Diet",
                                            "ruminant" = "Ruminant status",
                                            "habitat" = "Habitat"),
                                            levels = c("Species identity", "Host phylogeny", "Diet",  "Ruminant status", "Habitat")))

colnames(mrm_results_species_df) <- gsub('mb_dist', 'coef', colnames(mrm_results_species_df), fixed=TRUE)

write.csv(mrm_results_species_df, file = file.path(subdir, "func_mrm_results_with_species.csv"), row.names = FALSE, quote = TRUE)

# Plot result summary as barplot
p <- ggplot(mrm_results_species_df,
            aes(x = Variable, y = coef_median, fill = Variable)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = coef_q1, ymax = coef_q3), position = position_dodge(width = 0.9), width = 0.2) +
    scale_fill_manual(values = var_cols) +
    geom_text(aes(label = sig, y = ifelse(coef_median > 0, coef_q3 * 1.1, coef_q1 * 1.1)),
                  position = position_dodge(width = 0.9), size = 4) +
    facet_grid(Distance ~ Dataset, scale = "free", space = "free_x") +
    labs(y = "MRM coefficient", x = "") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(hjust = 1),
          panel.grid = element_blank()) +
    guides(fill = "none")

ggsave(p, filename = file.path(subdir, "func_mrm_results_with_species_barplot.png"), width = 8, height = 10)

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

# Gene tree based on Aitchison distances
mb_dist_gene <- vegdist(t(phy_gene_f_clr@otu_table), method="euclidean")
mb_tree_gene <- nj(mb_dist_gene)
mb_dist_gene_matrix <- as.matrix(mb_dist_gene)

# Pathway tree based on Aitchison distances
mb_dist_pathway <- vegdist(t(phy_pathway_clr@otu_table), method="euclidean")
mb_tree_pathway <- nj(mb_dist_pathway)
mb_dist_pathway_matrix <- as.matrix(mb_dist_pathway)

# GIFT tree based on Euclidean distances
mb_dist_gift <- vegdist(phy_gifts_el@otu_table, method="euclidean")
mb_tree_gift <- nj(mb_dist_gift)
mb_dist_gift_matrix <- as.matrix(mb_dist_gift)

# Get host distances per host species
host_dist <- dcast(host_dist_melt, Item1 ~ Item2, value = "host_distance") %>% column_to_rownames("Item1")

colnames(host_dist) <- sample_data$Common.name[match(colnames(host_dist), sample_data$Species)]
rownames(host_dist) <- sample_data$Common.name[match(rownames(host_dist), sample_data$Species)]

host_dist <- as.dist(host_dist)

# Get diet matrix again, but per host species
diet <- sample_data %>% data.frame %>% select(Common.name, Animalivory, PlantO, Frugivory, Seed) %>%
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

## Genes vs host phylogeny
res <- parafit(host.D = as.matrix(host_dist), para.D = mb_dist_gene_matrix, HP = assoc, nperm = 1000, seed = 123, correction="cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "Genes",
                                    host_distance = "Phylogeny",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

## Genes vs host diet
res <- parafit(host.D = as.matrix(diet_sp_dist), para.D = mb_dist_gene_matrix, HP = assoc, nperm = 1000, seed = 123, correction = "cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "Genes",
                                    host_distance = "Diet",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

## Pathways vs host phylogeny
res <- parafit(host.D = as.matrix(host_dist), para.D = mb_dist_pathway_matrix, HP = assoc, nperm = 1000, seed = 123, correction="cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "Pathways",
                                    host_distance = "Phylogeny",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

## Pathways vs host diet
res <- parafit(host.D = as.matrix(diet_sp_dist), para.D = mb_dist_pathway_matrix, HP = assoc, nperm = 1000, seed = 123, correction = "cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "Pathways",
                                    host_distance = "Diet",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

## GIFTs vs host phylogeny
res <- parafit(host.D = as.matrix(host_dist), para.D = mb_dist_gift_matrix, HP = assoc, nperm = 1000, seed = 123, correction="cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "GIFTs",
                                    host_distance = "Phylogeny",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

## GIFTs vs host diet
res <- parafit(host.D = as.matrix(diet_sp_dist), para.D = mb_dist_gift_matrix, HP = assoc, nperm = 1000, seed = 123, correction = "cailliez")

parafit_results <- rbind(parafit_results,
                         data.frame(mb_distance = "GIFTs",
                                    host_distance = "Diet",
                                    ParaFitGlobal = res$ParaFitGlobal,
                                    p.global = res$p.global))

write.csv(parafit_results, file = file.path(subdir, "func_parafit_results.csv"), row.names = FALSE, quote = TRUE)

#### Plot cophyloplots ####

host_consensus$tip.label <- sample_data$Common.name[match(host_consensus$tip.label, sample_data$Species)]

# Cophyloplot Genes
coph <- cophylo(tr1=host_consensus, tr2=mb_tree_gene, assoc=links)

png(file.path(subdir, "func_cophyloplot_gene.png"), width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()

coph <- cophylo(tr1=nj(diet_sp_dist), tr2=mb_tree_gene, assoc=links)

png(file.path(subdir, "func_cophyloplot_gene_diet.png"), width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()

# Cophyloplot Pathways
coph <- cophylo(tr1=host_consensus, tr2=mb_tree_pathway, assoc=links)

png(file.path(subdir, "func_cophyloplot_pathway.png"), width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()

coph <- cophylo(tr1=nj(diet_sp_dist), tr2=mb_tree_pathway, assoc=links)

png(file.path(subdir, "func_cophyloplot_pathway_diet.png"), width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()

# Cophyloplot GIFTs

coph <- cophylo(tr1=host_consensus, tr2=mb_tree_gift, assoc=links)

png(file.path(subdir, "func_cophyloplot_gift.png"), width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()

coph <- cophylo(tr1=nj(diet_sp_dist), tr2=mb_tree_gift, assoc=links)

png(file.path(subdir, "func_cophyloplot_gift_diet.png"), width = 800, height = 500)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(coph, link.lwd=4, link.lty="solid", link.col=links$colour)
dev.off()
