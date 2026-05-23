##### DIFFERENTIAL ABUNDANCE TESTS #####

#### Run differential abundance analysis on filtered data

################
#### SET UP ####
################

library(dplyr)
library(phyloseq)
library(microbiome)
library(tidyr)
library(tibble)
library(stringr)
library(phyr)
library(phytools)
library(MCMCglmm)
library(parallel)
library(coda)
library(lme4)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis"))
subdir <- normalizePath(file.path(outdir, "differential_abundance"))
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

# Habitat OTU relations
habitats_table <- read.csv(file.path(outdir, "habitat_relations.csv"))

# Use OTU relations
uses_table <- read.csv(file.path(outdir, "use_relations.csv"))

# Phenotype OTU relations
phenotype_table <- read.csv(file.path(outdir, "phenotype_relations.csv"))

# Host phylogeny
host_consensus <- read.tree(file.path(outdir, "host_consensus.tre"))

#########################
#### PREP TEST INPUT ####
#########################

phy_genus <- phy_sp_f %>% tax_glom("genus")

phy_genus@tax_table[, "genus"] <- make.names(phy_genus@tax_table[, "genus"], unique = TRUE)
taxa_names(phy_genus) <- phy_genus@tax_table[, "genus"]

phy_genus_clr <- phy_genus %>% transform("clr")

# Combine Sirenia/Proboscidea in one clades
phy_genus_clr@sam_data$Order <- ifelse(phy_genus_clr@sam_data$Order %in% c("Sirenia", "Proboscidea"), "Sirenia/Proboscidea", as.character(phy_genus_clr@sam_data$Order))

# Turn habitat to factor
phy_genus_clr@sam_data$habitat.general <- factor(phy_genus_clr@sam_data$habitat.general, levels = c("Terrestrial", "Marine"))

# Test for ruminant differences
phy_genus_clr@sam_data$ruminant <- factor(ifelse(phy_genus_clr@sam_data$digestion == "Ruminant", "Ruminant", "Other"), levels = c("Other", "Ruminant"))

# Test for hypsodont differences
phy_genus_clr@sam_data$hypsodont <- factor(ifelse(grepl("hyps", phy_genus_clr@sam_data$molar_category), "Hypsodont", "Other"), levels = c("Other", "Hypsodont"))

# Melt
data <- psmelt(phy_genus_clr) %>%
        select(OTU, Abundance, Sample, Species, Order, diet.general, habitat.general, ruminant, Fruit, Animal, Fruit, Seed, unmapped_count, contig_reads_count)

data <- data %>%
    # Turn S. scrofa domesticus to S. scrofa to match tree
    mutate(Species = case_when(Species == "Sus domesticus" ~ "Sus scrofa",
                               TRUE ~ Species))

# Keep only most abundant taxa (first average by species to weight against uneven sampling, then average over species)
top <- data %>% group_by(OTU, Species) %>%
  summarise(av_abundance = mean(Abundance)) %>%
  group_by(OTU) %>% summarise(av_abundance = mean(av_abundance)) %>%
  arrange(desc(av_abundance)) %>%
  slice_head(n = 200) %>% pull(OTU)

data_top <- data %>% filter(OTU %in% top)

# Phylogeny
host_consensus$node.label <- paste0("node", c(1:length(host_consensus$node.label)))
host_consensus$tip.label <- gsub("_", " ", host_consensus$tip.label)
Ainv <- inverseA(host_consensus)$Ainv

##################################
#### APPROACH 1: PGLMM/lambda ####
##################################

# Approach inspired by Youngblut et al. 2019

# list of all unique genera
otus <- top

#### PGLMM: Test for ecology after accounting for phylogeny
pglmm_res <- data.frame(term = character(), 
                        coef = numeric(),
                        pval = numeric(),
                        OTU = character())

for (otu in otus) {
    cat("Running PGLMM for OTU", which(otus == otu), "of", length(otus), "--", otu, "\n")
    # Subset data for the current OTU
    data_filt <- data %>% filter(OTU == otu)
    # Run PGLMM
    model <- pglmm(Abundance ~ Animal + Fruit + habitat.general + ruminant +
                   (1 | Species__), 
                   data = data_filt, 
                   cov_ranef = list(Species = host_consensus),
                   family = "gaussian")
    # Extract results
    res <- cbind(model$B, model$B.pvalue) %>% as.data.frame %>%
            rownames_to_column %>% filter(rowname != "(Intercept)") %>%
            mutate(rowname = str_remove(str_remove(rowname, "habitat.general"), "ruminant"))
    colnames(res) <- c("term", "coef", "pval")
    # Add random phylogenetic effect
    #res <- rbind(res, data.frame(term = "phylogeny", coef = unname(model$s2r[2]), pval = NA))
    res$OTU <- otu
    # Combine results
    pglmm_res <- rbind(pglmm_res, res)
}

# Make wide
pglmm_wide <- pglmm_res %>%
    pivot_wider(names_from = term, values_from = c(coef, pval))

write.csv(pglmm_wide, file = file.path(subdir, "pglmm_results.csv"), quote = FALSE, row.names = FALSE)

#### Phylogenetic signal (Pagel's lambda) after regressing out ecology

phy_res <- data.frame(OTU = character(),
                       lambda = numeric(),
                       pval = numeric())

for (otu in otus) {
    cat("Calculating phylogenetic signal for OTU", which(otus == otu), "of", length(otus), "--", otu, "\n")
    # Subset data for the current OTU
    data_filt <- data %>% filter(OTU == otu)
    # Run PGLMM
    model <- lm(Abundance ~ Animal + Fruit + habitat.general + ruminant,
                   data = data_filt)
    # Extract residuals
    resids_df <- data.frame(Sample = data_filt$Sample,
                            Species = data_filt$Species,
                            residuals = residuals(model))
    resids_by_species <- resids_df %>%
                        group_by(Species) %>%
                        summarise(residuals = mean(residuals)) %>% as.data.frame %>%
                        column_to_rownames("Species")
    # Make vector
    resids_vec <- resids_by_species$residuals
    names(resids_vec) <- rownames(resids_by_species)
    # Calculate Pagel's lamda
    pagel <- phylosig(host_consensus, resids_vec, method="lambda", test = TRUE)
    res <- data.frame(OTU = otu,
                      lambda = pagel$lambda,
                      pval = pagel$P)
    phy_res <- rbind(phy_res, res)
}

write.csv(phy_res, file = file.path(subdir, "pagel_results.csv"), quote = FALSE, row.names = FALSE)

#### Combine results ####
combined_res <- phy_res %>% rename(coef_lambda = lambda,
                                   pval_lambda = pval) %>%
    full_join(pglmm_wide, by = "OTU") %>%
    # Add metadata
    left_join(unique(select(data.frame(phy_genus_clr@tax_table), -species)), by = c("OTU" = "genus"))

# Adjust p-values
combined_res <- combined_res %>%
    mutate(padj_Animal = p.adjust(pval_Animal, method = "holm"),
           padj_Fruit = p.adjust(pval_Fruit, method = "holm"),
           padj_Marine = p.adjust(pval_Marine, method = "holm"),
           padj_Ruminant = p.adjust(pval_Ruminant, method = "holm"),
           padj_lambda = p.adjust(pval_lambda, method = "holm"))

write.csv(combined_res, file = file.path(subdir, "pglmm_pagel_combined_results.csv"), quote = FALSE, row.names = FALSE)

combined_coef <- select(pivot_longer(combined_res, cols = starts_with("coef_"), names_to = "term", values_to = "coefficient"), "OTU", "term", "coefficient") %>%
    mutate(term = str_remove(term, "coef_"))

combined_pval <- select(pivot_longer(combined_res, cols = starts_with("pval_"), names_to = "term", values_to = "pval"), "OTU", "term", "pval") %>%
    mutate(term = str_remove(term, "pval_"))

combined_padj <- select(pivot_longer(combined_res, cols = starts_with("padj_"), names_to = "term", values_to = "padj"), "OTU", "term", "padj", "superkingdom":"family") %>%
    mutate(term = str_remove(term, "padj_"))

combined_long <- full_join(combined_coef, combined_pval, by = c("OTU", "term")) %>%
    full_join(combined_padj, by = c("OTU", "term")) %>%
    mutate(phylum_grouped = case_when(phylum %in% names(phylum_palette) ~ phylum,
                                      superkingdom == "Bacteria" ~ "Other Bacteria",
                                      superkingdom == "Archaea" ~ "Other Archaea")) %>%
    mutate(term = factor(term, levels = c("lambda", "Animal", "Fruit", "Ruminant", "Marine")),
           phylum_grouped = factor(phylum_grouped, levels = rev(names(phylum_palette))),
           significant = ifelse(pval < 0.05, ifelse(padj < 0.05, "yes", "pre-adjustment"), "no"))

# Plot
p <- ggplot(combined_long, aes(x = coefficient, y = phylum_grouped, colour = phylum_grouped, shape = significant)) +
    geom_point(alpha = 0.7, size = 1.5, position = position_jitter(height = 0.2)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    facet_grid(cols = vars(term), scales = "free_x") +
    scale_colour_manual(values = phylum_palette, name = "Phylum") +
    scale_shape_manual(values = c("yes" = 19, "pre-adjustment" = 6, "no" = 4), name = "Significant (padj < 0.05)") +
    theme(legend.position = "bottom", legend.direction = "vertical",
          axis.title.y = element_blank()) +
    guides(colour = guide_legend(ncol = 2, byrow = TRUE),
           shape = guide_legend(nrow = 2, byrow = TRUE))

ggsave(p, filename = file.path(subdir, "pglmm_pagel_combined.png"), width = 8, height = 6)

##############################
#### APPROACH 2: MCMCglmm ####
##############################

# Based on tutorial by Sweeny et al. 2023, mSystems

if (file.exists(file.path(subdir, "mcmcglmm_output.RDS"))) {
    cat("MCMCglmm output exists. Loading...\n")
    m <- readRDS(file.path(subdir, "mcmcglmm_output.RDS"))
} else {
    cat("Running MCMCglmm...\n")
    set.seed(14)
    # Prep formula
    response <- "Abundance"
    fixed <- c("Animal:OTU", "Fruit:OTU", "ruminant:OTU", "habitat.general:OTU")
    formula <- as.formula(paste(response, "~", paste(fixed, collapse = "+")))
    n_taxa <- length(unique(data_top$OTU))
    # Set priors
    prior = list(R = list(V = diag(1), nu = 0.002),
                 G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                          G2 = list(V = diag(n_taxa), nu = 0.002, alpha.mu = rep(0, n_taxa), alpha.V = diag(n_taxa)*100)))
    start.time <- Sys.time()
    cat("Starting MCMCglmm... ")
    cat(start.time, "\n")
    # Run MCMCglmm
    m <- mclapply(1:10, function(i) {
        MCMCglmm(fixed = formula,
           random = ~ OTU + idh(OTU):Species,
           ginverse = list(Species=Ainv),
           prior = prior,
           data = data_top,
           verbose = TRUE, # pr = TRUE, pl = TRUE,
           nitt = 40000,
           burnin = 20000,
           thin = 50)
    }, mc.cores = 10)
    end.time <- Sys.time()
    cat("Finished MCMCglmm. ")
    cat(end.time, "\n")
    time.lapsed <- end.time - start.time
    cat("Time lapsed for MCMCglmm: ")
    print(time.lapsed)
    # Save output
    saveRDS(m, file.path(subdir, "mcmcglmm_output.RDS"))
}

mlist <- lapply(m, function(model) model$Sol)
mlist <- do.call(mcmc.list, mlist)

# Diagnostics with gelman plot
pdf(file=file.path(subdir, "mcmcglmm_gelman_plots.pdf"))
par(mfrow=c(4,2), mar=c(2,2,1,2))
gelman.plot(mlist, auto.layout=F, bin.width = 1)
dev.off()

gelman.diag(mlist)

# Plot first chain
m1 = m[[1]]

pdf(file=file.path(subdir, "mcmcglmm_plots.pdf"))
par(mfrow = c(2,2))
plot(m1)
dev.off()

# Autocorrelation
diag(autocorr(m1$VCV)[2, , ])

# 95% Credible interval
HPDinterval(m1$VCV)

# Collect results into tables
fixed_results <- summary(m1)$solutions %>%
  data.frame %>% rownames_to_column("term") %>%
  mutate(OTU = str_extract(term, "OTU[^:]*")) %>% # Separate OTU
  mutate(term = str_remove(term, OTU) %>% str_remove(":")) %>%  # Separate term
  mutate(OTU = str_remove_all(OTU, "OTU|:")) %>% # Remove fluff from OTU name
  # Remove intercepts
  filter(!is.na(OTU))

random_results <- summary(m1)$Gcovariances %>%
  data.frame %>% rownames_to_column("term") %>%
  mutate(OTU = str_remove(term, ".Species")) %>% # Separate OTU
  mutate(term = str_remove(term, paste0(OTU, "."))) %>%  # Separate term
  mutate(OTU = str_remove_all(OTU, "OTU")) %>% # Remove fluff from OTU name
  mutate(pMCMC = NA) %>%
  filter(OTU != "")

# Combine resulrts
mcmc_res <- rbind(fixed_results, random_results) %>%
    group_by(term) %>%
    mutate(padj = p.adjust(pMCMC, "holm"))

write.csv(mcmc_res, file = file.path(subdir, "mcmcglmm_results.csv"), quote = FALSE, row.names = FALSE)

#### Extract diet and phylogeny effects per taxon ####

# Calculate phylogenetic lambda (λ) de Villemereuil and Nakagawa (2014).
# But it needs to be a different calculation for each taxon.

# Adapted from 
lambda <- m1$VCV %>% as.data.frame %>%
    pivot_longer(cols = contains("OTU"), names_to = "term", values_to = "Variance") %>%
    mutate(OTU = str_remove(term, ".Species")) %>% # Separate OTU
    mutate(term = str_remove(term, paste0(OTU, "."))) %>%  # Separate term
    mutate(OTU = str_remove_all(OTU, "OTU")) %>% # Remove fluff from OTU name
    filter(term != "OTU") %>%
    # Calculate lambda
    mutate(lambda = Variance/(Variance + units)) %>%
    # Average
    group_by(OTU) %>%
    summarise(lambda = mean(lambda, na.rm = TRUE))

write.csv(lambda, file = file.path(subdir, "mcmcglmm_lambda.csv"), quote = FALSE, row.names = FALSE)

# Calculate fixed effect R2 per Nakagawa & Schielzeth (2013)
# Get fitted values from fixed effects only
get_r2 <- function(n, model) {
    X_sub <- model$X[, grep(n, colnames(model$X))]
    sol_means <- colMeans(model$Sol[, grep(n, colnames(model$Sol))])
    fitted_fixed <- X_sub %*% sol_means
    # Variance of the fixed effects
    V_fixed <- var(as.numeric(fitted_fixed))
    #V_random <- sum(colMeans(model$VCV[, grepl(n, colnames(model$VCV)) | colnames(model$VCV) == "OTU"]))
    V_residual <- mean(model$VCV[, "units"])
    # Marginal R² (excluding random effects)
    R2 <- V_fixed / (V_fixed + V_residual)
    return(R2)
}

r2 <- data.frame(OTU = lambda$OTU, R2 = sapply(lambda$OTU, get_r2, model = m1))

write.csv(r2, file = file.path(subdir, "mcmcglmm_r2.csv"), quote = FALSE, row.names = FALSE)

# Combine lamda and R2
phylo_v_eco <-
    full_join(lambda, r2, by = "OTU") %>%
    # Add taxonomy
    left_join(phy_genus_clr@tax_table %>% as.data.frame %>% rownames_to_column("OTU"), by = "OTU") %>%
    mutate(phylum_grouped = case_when(phylum %in% names(phylum_palette) ~ phylum,
                                      superkingdom == "Bacteria" ~ "Other Bacteria",
                                      superkingdom == "Archaea" ~ "Other Archaea"))

p <- ggplot(phylo_v_eco, aes(x = lambda, y = phylum_grouped, colour = phylum_grouped)) +
    geom_jitter(size = 3, alpha = 0.8, height = 0.2, width = 0) +
    scale_colour_manual(values = phylum_palette, name = "Phylum")

ggsave(p, filename = file.path(subdir, "mcmcglmm_lambda.png"), width = 8, height = 6)

p <- ggplot(phylo_v_eco, aes(x = lambda, y = R2, colour = phylum_grouped)) +
    geom_point() +
    scale_colour_manual(values = phylum_palette, name = "Phylum") +
    labs(x = "Phylogenetic lambda (λ)", y = "Ecology effects (fixed R2)") +
    theme(legend.position = "bottom")

ggsave(p, filename = file.path(subdir, "mcmcglmm_lambda_v_r2.png"), width = 8, height = 6)

###########################################
#### COMBINE MCMCglmm and PGLMM/lambda ####
###########################################

# Use lambda instead of species for MCMCglmm results to match with PGLMM
lambda <- lambda %>% rename(coefficient = lambda) %>%
        mutate(term = "lambda",
               pval = NA)

mcmc_res_long <- mcmc_res %>% select(OTU, term, post.mean, pMCMC, padj) %>%
    mutate(term = str_remove(term, "habitat.general") %>% str_remove("ruminant")) %>%
    filter(term != "Species") %>%
    rename(coefficient = post.mean,
           pval = pMCMC) %>%
    rbind(lambda)

# t for traditional approach, b for Bayesian MCMCglmm
all_res <- full_join(combined_long, mcmc_res_long, by = c("OTU", "term"), suffix = c("_t", "_b"))

write.csv(all_res, file = file.path(subdir, "all_diffabund_results.csv"), quote = FALSE, row.names = FALSE)

# Do coefficients correlated?
p <- ggplot(all_res, aes(x = coefficient_t, y = coefficient_b, colour = term)) +
    geom_point() +
    facet_wrap(~ term, scales = "free") +
    labs(x = "PGLMM/Pagel's lambda coefficient", y = "MCMCglmm coefficient") +
    theme(legend.position = "none")

ggsave(p, filename = file.path(subdir, "all_diffabund_coeff_correlation.png"), width = 6, height = 6)

### Plot abundances ####

res_labels <- all_res %>%
            # Keep significant results only in either method (after p-value adjustment)
            filter(padj_t < 0.05 | pval_b < 0.05) %>%
            # label association (positive and negative)
            mutate(assoc = case_when(coefficient_t < 0 & coefficient_b < 0 ~ paste0(term, "-"),
                                     coefficient_t > 0 & coefficient_b > 0 ~ paste0(term, "+"),
                                     coefficient_t < 0 & coefficient_b > 0 ~ paste0(term, "mixed-+"),
                                     coefficient_t > 0 & coefficient_b < 0 ~ paste0(term, "mixed+-")),
                   signif = case_when(padj_t < 0.05 & pval_b < 0.05 ~ "both padj < 0.05",
                                      padj_t < 0.05 & pval_b < 0.05 ~ "padj_t < 0.05, pval_b < 0.05",
                                      padj_t < 0.05 ~ "padj_t < 0.05, pval_b ns",
                                      pval_t < 0.05 & pval_b < 0.05 ~ "pval_t < 0.05, pval_b < 0.05",
                                      pval_t >= 0.05 ~ "pval_t ns, pval_b < 0.05",
                                      # if an OTU doesn't come up as significant in either method (after adjusting), remove
                                      TRUE ~ "remove")) %>% 
            filter(signif != "remove") %>%
            # Then summarise all associations per taxon
            group_by(OTU) %>%
            summarise(label = paste(paste(assoc, signif, sep = ": "), collapse="\n")) %>%
            mutate(label = str_remove_all(label, "ruminant|habitat.general"))

# Get abundances per sample for the differentially abundant taxa
abundances <- phy_genus_clr@otu_table %>% t %>% data.frame %>% rownames_to_column("Sample") %>%
              pivot_longer(cols = -Sample, names_to = "OTU", values_to = "Abundance") %>%
              left_join(rownames_to_column(select(data.frame(phy_genus_clr@sam_data), Common.name, Order, diet.general), "Sample"), by = "Sample") %>%
              right_join(res_labels)

# Reorder host species
species_levels <- phy_genus_clr@sam_data %>% data.frame %>% arrange(as.character(Order), as.character(digestion), Common.name) %>% select(Order, Common.name) %>% unique
abundances$Common.name <- factor(abundances$Common.name, levels = species_levels$Common.name)

# Reorder genes
genus_levels <- abundances %>% arrange(label) %>% pull(OTU) %>% unique
abundances$OTU <- factor(abundances$OTU, levels = genus_levels)

# Plot
order_palette2 <- order_palette
order_palette2["Sirenia/Proboscidea"] <- order_palette2["Sirenia"]

p <- ggplot(abundances, aes(x = Common.name, y = Abundance, colour = Order, fill = diet.general)) +
    geom_boxplot(alpha = 0.8, size = 0.5) +
    scale_colour_manual(values = order_palette2, name = "Order") +
    scale_fill_manual(values = diet_palette, name = "Diet") +
    facet_wrap(~ paste(as.character(OTU), label, sep = "\n"), ncol = 6, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 8),
          legend.position = "bottom") + ylab("CLR-transformed abundances") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave(p, filename = file.path(subdir, "all_diffabund_abundances.png"), width = 20, height = 20)
