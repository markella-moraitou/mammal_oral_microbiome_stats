##### FUNCTION VS TAXONOMY PROCRUSTES #####

#### Use Procrustes analysis to compare functional and taxonomic composition

################
#### SET UP ####
################

library(dplyr)
library(phyloseq)
library(microbiome)
library(microViz)
library(tidyr)
library(stringr)
library(tibble)
library(vegan)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata
outdir <- normalizePath(file.path("..", "..", "output", "function"))
taxdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Output of community analysis
subdir <- normalizePath(file.path(outdir, "func_tax_procrustes"))

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Phyloseq objects
phy_gene_f <- readRDS(file.path(outdir, "data", "phy_gene_f.RDS"))
phy_gene_f_clr <- readRDS(file.path(outdir, "data", "phy_gene_f_clr.RDS"))

phy_sp_f <- readRDS(file.path(taxdir, "phyloseq_objects", "phy_sp_f.RDS"))
phy_sp_f_clr <- readRDS(file.path(taxdir, "phyloseq_objects", "phy_sp_f_clr.RDS"))

##########################
#### DEFINE FUNCTION  ####
##########################

func_tax_procrustes <- function(phy_func, phy_tax, dist, out_suffix) {
    # Run PCA on both datasets
    cat("Running PCoA on functional data...\n")
    ord_func <- dist_calc(phy_func, dist) %>% ord_calc(method = "PCoA", centre = TRUE, scale. = TRUE)
    cat("Running PCoA on taxonomic data...\n")
    ord_tax <- dist_calc(phy_tax, dist) %>% ord_calc(method = "PCoA", centre = TRUE, scale. = TRUE)
    # Procrustes analysis
    pro <- procrustes(X = ord_tax@ord, Y = ord_func@ord, symmetric = FALSE)
    # Test significance
    set.seed(123)
    test <- protest(X = ord_tax@ord, Y = ord_func@ord, scores = "sites", permutations = 999)
    ss <- test$ss
    pval <- test$signif
    # Plot results
    png(file.path(subdir, paste0("procrustes_viz_", out_suffix, ".png")), width = 6, height = 6, units = "in", res = 300)
    plot(pro, kind = 1)
    title(main = paste0("Procrustes analysis (", dist, ")\n",
                        "SS = ", round(ss, 4), ", p = ", pval))
    dev.off()
    # Plot residuals per species
    res <- residuals(pro) %>% as.data.frame %>% rownames_to_column()
    colnames(res) <- c("Sample", "Residual")
    # Add metadata
    res <- res %>% left_join(select(data.frame(phy_tax@sam_data),
                                    c(new_name, diet.general, Species, Order, Order_grouped)),
                             by = c("Sample" = "new_name"))
    # Plot and save
    p <- ggplot(res, aes(x = Species, y = Residual, fill = diet.general)) +
        geom_boxplot() +
        facet_grid(cols = vars(Order_grouped), scales = "free_x", space = "free_x") +
        scale_fill_manual(values = diet_palette)
    ggsave(file.path(subdir, paste0("procrustes_residuals_", out_suffix, ".png")), p, width = 8, height = 6)
}

#############
#### RUN ####
#############

# Keep samples that are in both datasets
shared_samples <- intersect(sample_names(phy_gene_f_clr), sample_names(phy_sp_f_clr))

# Species level (Jaccard distances)
func_tax_procrustes(
            phy_tax = prune_samples(shared_samples, phy_sp_f),
            phy_func = prune_samples(shared_samples, phy_gene_f),
            dist = "jaccard", out_suffix = "jaccard_species")

# Species level (Aitchison distances)
func_tax_procrustes(
            phy_tax = prune_samples(shared_samples, phy_sp_f_clr),
            phy_func = prune_samples(shared_samples, phy_gene_f_clr),
            dist = "euclidean", out_suffix = "aitchison_species")

#### Get genus level data
phy_gen_f <- tax_glom(phy_sp_f, taxrank = "genus", NArm = FALSE)

phy_gen_f_clr <- microbiome::transform(phy_gen_f, "clr")

# Genus level (Jaccard distances)
func_tax_procrustes(
            phy_tax = prune_samples(shared_samples, phy_gen_f),
            phy_func = prune_samples(shared_samples, phy_gene_f),
            dist = "jaccard", out_suffix = "jaccard_genus")

# Genus level (Aitchison distances)
func_tax_procrustes(
            phy_tax = prune_samples(shared_samples, phy_gen_f_clr),
            phy_func = prune_samples(shared_samples, phy_gene_f_clr),
            dist = "euclidean", out_suffix = "aitchison_genus")
