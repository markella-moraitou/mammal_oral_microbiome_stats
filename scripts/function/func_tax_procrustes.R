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
library(ggnewscale)

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
phy_gene_f_clr <- readRDS(file.path(outdir, "data", "phy_gene_f_clr.RDS"))

phy_pathway_clr <- readRDS(file.path(outdir, "pathway_completeness", "phy_pathway_clr.RDS"))

phy_gifts_el <- readRDS(file.path(outdir, "data", "phy_gifts_el.RDS"))

phy_sp_f_clr <- readRDS(file.path(taxdir, "phyloseq_objects", "phy_sp_f_clr.RDS"))
phy_sp_f <- readRDS(file.path(taxdir, "phyloseq_objects", "phy_sp_f.RDS"))

#####################
####  PREP DATA  ####
#####################

# Recode order and habitat as TRUE and FALSE
phy_gene_f_clr <- phy_gene_f_clr %>%
         ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  Ruminant = (digestion == "Ruminant"),
                  Marine = (habitat.general == "Marine"))

phy_pathway_clr <- phy_pathway_clr %>%
         ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  Ruminant = (digestion == "Ruminant"),
                  Marine = (habitat.general == "Marine"))

phy_gifts_el <- phy_gifts_el %>%
         ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  Ruminant = (digestion == "Ruminant"),
                  Marine = (habitat.general == "Marine"))

phy_sp_f_clr <- phy_sp_f_clr %>%
         ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  Ruminant = (digestion == "Ruminant"),
                  Marine = (habitat.general == "Marine"))

# Species traits to use as constraints
species_traits <- c("Artiodactyla", "Perissodactyla", "Primates",
                    "Ruminant", "Marine", "Frugivory", "Animalivory")

##########################
#### DEFINE FUNCTION  ####
##########################

func_tax_procrustes <- function(phy_func, phy_tax, out_suffix) {
    # Scale and centre the data
    phy_func <- microbiome::transform(phy_func, "standardize")
    phy_tax <- microbiome::transform(phy_tax, "standardize")
    # Run PCA on both datasets
    cat("Running RDA on functional data...\n")
    ord_func <- phy_func %>% ord_calc(constraints = species_traits, method = "RDA", centre = TRUE)
    cat("Running RDA on taxonomic data...\n")
    ord_tax <- phy_tax %>% ord_calc(constraints = species_traits, method = "RDA", centre = TRUE)
    
    ## Procrustes analysis
    pro <- procrustes(X = ord_tax@ord, Y = ord_func@ord, symmetric = FALSE)
    
    # Test significance
    set.seed(123)
    test <- protest(X = ord_tax@ord, Y = ord_func@ord, scores = "sites", permutations = 999)
    ss <- test$ss
    pval <- test$signif
    
    ## Plot residuals per species
    res <- residuals(pro) %>% as.data.frame %>% rownames_to_column()
    colnames(res) <- c("Sample", "Residual")
    # Add metadata
    res <- res %>% left_join(select(data.frame(phy_tax@sam_data),
                                    c(new_name, diet.general, Species, Order, Order_grouped)),
                             by = c("Sample" = "new_name"))
    
    # Plot and save
    p <- ggplot(res, aes(y = Species, x = Residual, fill = diet.general)) +
        geom_boxplot() +
        facet_grid(rows = vars(Order_grouped), scales = "free_y", space = "free_y") +
        scale_fill_manual(values = diet_palette, name = "Host diet")
    ggsave(file.path(subdir, paste0("procrustes_residuals_", out_suffix, ".png")), p, width = 8, height = 6)
    
    ## Plot RDA comparison
    pro_df <- cbind(data.frame(pro$Yrot), data.frame(pro$X)) %>% rownames_to_column()
    colnames(pro_df) <- c("Sample", "Func1", "Func2", "Tax1", "Tax2")
    
    # Get species with higest residuals to colour
    highlight_species <- res %>% group_by(Species) %>% summarise(mean_res = mean(Residual)) %>%
        arrange(desc(mean_res)) %>% slice_head(n = 4) %>% pull(Species)
    
    pro_df <- pro_df %>% left_join(select(data.frame(phy_tax@sam_data),
                                    c(new_name, diet.general, Species, Order, Order_grouped)),
                             by = c("Sample" = "new_name")) %>%
        mutate(highlight = case_when(Species %in% highlight_species ~ Species, TRUE ~ "Other")) %>%
        mutate(highlight = factor(highlight, levels = c(levels(highlight_species), "Other")))
    
    species_palette_mod <- append(species_palette, setNames(c("black"), "Other"))
    
    p <- ggplot(data = pro_df) +
        # Add ellipses for each diet group
        stat_ellipse(aes(x = Tax1, y = Tax2, fill = diet.general), geom = "polygon", type = "t", alpha = 0.3) +
        scale_fill_manual(values = diet_palette, name = "Host diet") +
        # Add arrows
        geom_segment(aes(x = Tax1, y = Tax2, xend = Func1, yend = Func2, colour = highlight),
                    arrow = arrow(length=unit(0.10,"cm")), linewidth = 0.2) +
        scale_colour_manual(values = species_palette_mod, name = "Species") +
        labs(title = paste0("Procrustes analysis\nSS = ", round(ss, 4), ", p = ", pval),
             x = "RDA 1", y = "RDA2")
    
    ggsave(file.path(subdir, paste0("procrustes_viz_", out_suffix, ".png")), p, width = 6, height = 5)
}

#############
#### RUN ####
#############

# Keep samples that are in both datasets
shared_samples <- intersect(sample_names(phy_gene_f_clr), sample_names(phy_sp_f_clr))

# Species level - Genes
func_tax_procrustes(
            phy_tax = prune_samples(shared_samples, phy_sp_f_clr),
            phy_func = prune_samples(shared_samples, phy_gene_f_clr),
            out_suffix = "species_genes")

# Species level - Pathways
func_tax_procrustes(
            phy_tax = prune_samples(shared_samples, phy_sp_f_clr),
            phy_func = prune_samples(shared_samples, phy_pathway_clr),
            out_suffix = "species_pathways")

# Species level - GIFTS
func_tax_procrustes(
            phy_tax = prune_samples(shared_samples, phy_sp_f_clr),
            phy_func = prune_samples(shared_samples, phy_gifts_el),
            out_suffix = "species_gifts")

#### Get genus level data
phy_gen_f <- tax_glom(phy_sp_f, taxrank = "genus", NArm = FALSE)

phy_gen_f_clr <- microbiome::transform(phy_gen_f, "clr")
phy_gen_f_clr@sam_data <- phy_sp_f_clr@sam_data

# Genus level - Genes
func_tax_procrustes(
            phy_tax = prune_samples(shared_samples, phy_gen_f_clr),
            phy_func = prune_samples(shared_samples, phy_gene_f_clr),
            out_suffix = "genus_genes")

# Genus level - Pathways
func_tax_procrustes(
            phy_tax = prune_samples(shared_samples, phy_gen_f_clr),
            phy_func = prune_samples(shared_samples, phy_pathway_clr),
            out_suffix = "genus_pathways")

# Genus level - GIFTS
func_tax_procrustes(
            phy_tax = prune_samples(shared_samples, phy_gen_f_clr),
            phy_func = prune_samples(shared_samples, phy_gifts_el),
            out_suffix = "genus_gifts")
