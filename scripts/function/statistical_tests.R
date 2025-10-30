##### STATISTICAL TESTS #####

#### Run PERMANOVA, differential abundance analysis and other tests on gene data

################
#### SET UP ####
################

library(dplyr)
library(phyloseq)
library(microbiome)
library(tidyr)
library(stringr)
library(vegan)
library(ANCOMBC)
library(MCMCglmm)
library(parallel)
library(coda)
library(tibble)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata
outdir <- normalizePath(file.path("..", "..", "output", "function"))
taxdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Output of community analysis
subdir <- normalizePath(file.path(outdir, "statistical_tests"))

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
phy_gene <- readRDS(file.path(outdir, "phy_gene.RDS"))
phy_gene_clr <- readRDS(file.path(outdir, "phy_gene_clr.RDS"))

# Low content sample list
low_content <- read.table(file.path(outdir, "low_content_samples.txt"), header = TRUE)

# Host phylogeny
host_consensus <- read.tree(file.path(taxdir, "host_consensus.tre"))

####################
#### PREP DATA  ####
####################

# Remove low content samples
phy_gene <- phy_gene %>% prune_samples(!sample_names(phy_gene) %in% low_content$x, .)
phy_gene_clr <- phy_gene_clr %>% prune_samples(!sample_names(phy_gene_clr) %in% low_content$x, .)

# Keep only genes with at least 10% prevalence
phy_gene_clr <- phy_gene_clr %>% subset_taxa(prevalence(phy_gene)>=0.1)

# Extract OTU table and sample data
otu_table <- as.data.frame(phy_gene_clr@otu_table)
sample_data <- as.data.frame(phy_gene_clr@sam_data)

# Transpose the OTU table
otu_table_t <- t(otu_table)

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

#### CLR-transformed data
set.seed(123)

perm <- adonis2(otu_table_t ~ order + diet + habitat + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_clr_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(otu_table_t ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_clr_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

#########################
#### PREP TEST INPUT ####
#########################

# Combine Sirenia/Proboscidea in one clades
phy_gene_clr@sam_data$Order <- ifelse(phy_gene_clr@sam_data$Order %in% c("Sirenia", "Proboscidea"), "Sirenia/Proboscidea", as.character(phy_gene_clr@sam_data$Order))

# Turn habitat to factor
phy_gene_clr@sam_data$habitat.general <- factor(phy_gene_clr@sam_data$habitat.general, levels = c("Terrestrial", "Marine"))

# Test for ruminant differences
phy_gene_clr@sam_data$ruminant <- factor(ifelse(phy_gene_clr@sam_data$digestion == "Ruminant", "Ruminant", "Other"), levels = c("Other", "Ruminant"))

# Test for hypsodont differences
phy_gene_clr@sam_data$hypsodont <- factor(ifelse(grepl("hyps", phy_gene_clr@sam_data$molar_category), "Hypsodont", "Other"), levels = c("Other", "Hypsodont"))

# OTU table
otu <- phy_gene_clr@otu_table %>% data.frame %>%
    rownames_to_column("OTU") %>% pivot_longer(cols = where(is.numeric), names_to = "Sample", values_to = "Abundance")

# Metadata
sam <- phy_gene_clr@sam_data %>% data.frame %>%
    select(Species, Order, diet.general, habitat.general, ruminant, Total_abundance, cf, cp, nfe, ee)%>% rownames_to_column("Sample")

# Combine data
data <- left_join(otu, sam, relationship = "many-to-one") %>%
  # Have sus scrofa as sus domesticus
  mutate(Species = case_when(Species == "Sus domesticus" ~ "Sus scrofa",
                             TRUE ~ Species))

data_wide <- data %>% pivot_wider(names_from = "OTU", values_from = "Abundance") %>% data.frame
colnames(data_wide) <- gsub(colnames(data_wide), pattern = ".", replacement = "_", fixed = TRUE)

# Get a metrics for diet
data_diet <- data_wide %>% select(Species, diet_general, cp, cf, ee, nfe) %>% unique %>%
  mutate(animalivory = log((cp + ee)/(cf+nfe)),
         frugivory = log((nfe + 5)/(cf + 5)))

p <- ggplot(data_diet, aes(x = animalivory, y = frugivory, colour = diet_general)) +
  geom_point(size = 2) + theme_bw() +
  geom_text(aes(label = Species)) +
  scale_color_manual(values = diet_palette)

ggsave(p, filename = file.path(subdir, "diet_variables.png"), width = 8, height = 8)

# Add to phyloseq
data_diet <- rbind(data_diet, data_diet[which(data_diet$Species == "Sus scrofa"),])
data_diet[nrow(data_diet), "Species"] <- "Sus domesticus"
phy_gene_clr@sam_data$animalivory <- data_diet$animalivory[match(phy_gene_clr@sam_data$Species, data_diet$Species)]
phy_gene_clr@sam_data$frugivory <- data_diet$frugivory[match(phy_gene_clr@sam_data$Species, data_diet$Species)]

# Add variables to big data table
data_wide <- data_diet %>% right_join(data_wide)

# Phylogeny
host_consensus$node.label <- paste0("node", c(1:length(host_consensus$node.label)))
host_consensus$tip.label <- gsub("_", " ", host_consensus$tip.label)
Ainv <- inverseA(host_consensus)$Ainv

######################
#### RUN MCMCglmm ####
######################

if (file.exists(file.path(subdir, "mcmcglmm_output.RDS"))) {
    cat("MCMCglmm output exists. Loading...\n")
    m <- readRDS(file.path(subdir, "mcmcglmm_output.RDS"))
} else {
    cat("Running MCMCglmm...\n")
    set.seed(14)
    responses <- intersect(colnames(data_wide), gsub(" ", "_", unique(otu$OTU)))
    formula <- as.formula(paste(
        "cbind(",
        paste(responses, collapse = ", "),
        ")  ~ -1 + trait + trait:animalivory + trait:frugivory + trait:habitat_general + trait:Order + trait:ruminant"
        ))
    m <- mclapply(1:10, function(i) {
    MCMCglmm(formula,
           random = ~idh(trait):Species,
           ginverse=list(Species=Ainv),
           rcov = ~idh(trait):units,
           data = data_wide, family = rep("gaussian", length(responses)),
           verbose = TRUE,
           nitt = 13000,
           burnin = 3000,
           thin = 10)
    }, mc.cores = 10)
    saveRDS(m, file.path(subdir, "mcmcglmm_output.RDS"))
}

mlist <- lapply(m, function(model) model$Sol)
mlist <- do.call(mcmc.list, mlist)

# Diagnostics with gelman plot
pdf(file=file.path(subdir, "mcmcglmm_gelman_plots.pdf"))
par(mfrow=c(4,2), mar=c(2,2,1,2))
gelman.plot(mlist, auto.layout=F)
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
fixed_results <- summary(m1)$solutions %>% data.frame %>% rownames_to_column("term") %>%
  mutate(term = str_remove(term, "trait")) %>%
  separate(term, into = c("gene", "term"), sep = ":", fill = "right") %>%
  # Remove intercepts
  filter(!is.na(term))

random_results <- summary(m1)$Gcovariances %>%
  data.frame %>% rownames_to_column("term") %>%
  mutate(term = str_remove(term, "trait")) %>%
  separate(term, into = c("gene", "term"), sep = "[.]") %>%
  mutate(pMCMC = NA)

# Combine resulrts
mcmc_res <- rbind(fixed_results, random_results) %>%
    group_by(term) %>%
    mutate(qval = p.adjust(pMCMC, "fdr"))

write.csv(mcmc_res, file = file.path(subdir, "mcmcglmm_results.csv"), quote = FALSE, row.names = FALSE)

#### Plot differentially abundant genes ####
# There is no significant results, so keep the 100 smallest qvalues
mcmc_signif <- mcmc_res %>%
            filter(!term %in% c("Species", "Total_abundance")) %>%
            filter(qval < 0.05)

mcmc_res_label <- mcmc_signif %>%
            # label association (positive and negative)
            mutate(assoc = case_when(post.mean < 0 ~ paste0(term, "-"),
                                     post.mean > 0 ~ paste0(term, "+"))) %>%
            mutate(assoc = str_remove(assoc, "Order|habitat_general|ruminant")) %>%
            # Then summarise all associations per taxon
            group_by(gene) %>%
            summarise(label = paste(assoc, collapse=" "))

# Get abundances per sample for the differentially abundant genes
abundances <- phy_gene_clr@otu_table %>% t %>% data.frame %>% rownames_to_column("Sample") %>%
              pivot_longer(cols = -Sample, names_to = "gene", values_to = "Abundance") %>%
              left_join(rownames_to_column(select(data.frame(phy_gene@sam_data), Common.name, Order, diet.general), "Sample"), by = "Sample")

mcmc_abund <- abundances %>% right_join(mcmc_res_label)
mcmc_abund$gene_description <- as.vector(phy_gene@tax_table[match(mcmc_abund$gene, phy_gene@tax_table[,"gene_id"]), "gene_description"])

# Reorder host species
species_levels <- phy_gene_clr@sam_data %>% data.frame %>% arrange(as.character(Order), as.character(digestion), Common.name) %>% select(Order, Common.name) %>% unique
mcmc_abund$Common.name <- factor(mcmc_abund$Common.name, levels = species_levels$Common.name)

# Reorder genes
gene_levels <- mcmc_abund %>% arrange(label) %>% pull(gene_description) %>% unique
mcmc_abund$gene_description <- factor(mcmc_abund$gene_description, levels = gene_levels)

# Plot
order_palette2 <- order_palette
order_palette2["Sirenia/Proboscidea"] <- order_palette2["Sirenia"]

p <- ggplot(mcmc_abund, aes(x = Common.name, y = Abundance, colour = Order, fill = diet.general)) +
    geom_boxplot(alpha = 0.8, size = 0.5) +
    scale_colour_manual(values = order_palette2, name = "Order") +
    scale_fill_manual(values = diet_palette, name = "Diet") +
    facet_wrap(~ paste(paste(str_trunc(as.character(gene_description), 49), gene), label, sep = "\n"), ncol = 3, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 8),
          legend.position = "bottom") + ylab("CLR-transformed abundances") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave(p, filename = file.path(subdir, "mcmc_abundances.png"), width = 12, height = 40)

#####################
#### RUN ANCOMBC ####
#####################

phy_gene@sam_data <- phy_gene_clr@sam_data

# Check if output is there, if not run
if(file.exists(file.path(subdir, "ancombc_results.rds"))) {
    cat("Output already there, so analysis will not be repeated! Loading existing ANCOMBC results...\n")
    out <- readRDS(file.path(subdir, "ancombc_results.rds"))
} else {
    out <- ancombc2(data = phy_gene,
                 fix_formula = "animalivory + frugivory + habitat.general + ruminant", #  + Order",
                 rand_formula = "(1|Species)",
                 p_adj_method = "fdr", prv_cut = 0.1, # Some parameters. prv_cut=0.10 means that taxa found in less than 10% of samples are ignored
                 struc_zero = FALSE,
                 global = FALSE,
                 lib_cut = 0,
                 verbose = TRUE)
    # Save output
    saveRDS(out, file.path(subdir, "ancombc_results.rds"))
    write.csv(out$res, file.path(subdir, "ancombc_res.csv"), row.names = FALSE, quote = TRUE)
    write.csv(out$ss_tab, file.path(subdir, "ancombc_ss_tab.csv"), row.names = FALSE, quote = TRUE)
    write.csv(out$samp_frac, file.path(subdir, "ancombc_samp_frac.csv"), row.names = TRUE, quote = TRUE)
}

#### Summary plots ####

# Get ancom results in a long format
lfc <-
    out$res %>%
    select(taxon, contains("lfc")) %>%
    pivot_longer(cols = contains("lfc"), names_to = "term", values_to = "lfc") %>%
    mutate(term = str_remove(term, "lfc_")) %>%
    mutate(term = str_remove(term, "Order|ruminant|diet.general|habitat.general|hypsodont"))

pvalues <-
    out$res %>%
    select(taxon, contains("p_")) %>%
    pivot_longer(cols = contains("p_"), names_to = "term", values_to = "pval") %>%
    mutate(term = str_remove(term, "p_")) %>%
    mutate(term = str_remove(term, "Order|ruminant|diet.general|habitat.general|hypsodont"))

qvalues <- 
    out$res %>%
    select(taxon, contains("q_")) %>%
    pivot_longer(cols = contains("q_"), names_to = "term", values_to = "qval") %>%
    mutate(term = str_remove(term, "q_")) %>%
    mutate(term = str_remove(term, "Order|ruminant|diet.general|habitat.general|hypsodont"))

sensitivity <- 
    out$res %>%
    select(taxon, contains("passed_ss_")) %>%
    pivot_longer(cols = contains("passed_ss_"), names_to = "term", values_to = "sensitivity") %>%
    mutate(term = str_remove(term, "passed_ss_")) %>%
    mutate(term = str_remove(term, "Order|ruminant|diet.general|habitat.general|hypsodont"))

ancombc_res_long <- lfc %>% full_join(pvalues) %>% full_join(qvalues) %>% full_join(sensitivity) %>%
                filter(term != "(Intercept)") %>%
                rename("gene" = "taxon")

write.csv(ancombc_res_long, file = file.path(subdir, "ancombc_long.csv"), quote = FALSE, row.names = FALSE)

# Volcano plot
p <- filter(ancombc_res_long) %>%
    ggplot(aes(x = lfc, y = -log10(pval), shape = sensitivity, colour = (qval <= 0.05))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    facet_wrap(~ term, ncol = 3) +
    scale_shape_manual(values = c("cross", "circle")) +
    scale_colour_manual(values = c("TRUE" = "#F6BD23", "FALSE" = "#0C1950"), name = "q-value <= 0.05") +
    labs(x = "Log-fold Change", y = "-log10 p-value") +
    theme(legend.position = "bottom")

ggsave(p, filename = file.path(subdir, "ancombc_volcano.png"), width = 8, height = 5)

#### Plot differentially abundant taxa ####

# Get significant results
ancombc_signif <- filter(ancombc_res_long) %>% filter(pval <= 0.05 & sensitivity == TRUE)

ancombc_res_label <- ancombc_signif %>%
            # label association (positive and negative)
            mutate(assoc = case_when(lfc < 0 ~ paste0(term, "-"),
                                     lfc > 0 ~ paste0(term, "+"))) %>%
            # Then summarise all associations per taxon
            group_by(gene) %>%
            summarise(label = paste(assoc, collapse=" "))

# Get abundances per sample for the differentially abundant taxa
abundances <- phy_gene_clr@otu_table %>% t %>% data.frame %>% rownames_to_column("Sample") %>%
              pivot_longer(cols = -Sample, names_to = "gene", values_to = "Abundance") %>%
              left_join(rownames_to_column(select(data.frame(phy_gene@sam_data), Common.name, Order, diet.general), "Sample"), by = "Sample")

ancombc_abund <- abundances %>% right_join(ancombc_res_label)
ancombc_abund$gene_description <- as.vector(phy_gene@tax_table[match(ancombc_abund$gene, phy_gene@tax_table[,"gene_id"]), "gene_description"])

# Reorder host species
species_levels <- phy_gene_clr@sam_data %>% data.frame %>% arrange(as.character(Order), as.character(digestion), Common.name) %>% select(Order, Common.name) %>% unique
ancombc_abund$Common.name <- factor(ancombc_abund$Common.name, levels = species_levels$Common.name)

# Reorder genes
gene_levels <- ancombc_abund %>% arrange(label) %>% pull(gene_description) %>% unique
ancombc_abund$gene_description <- factor(ancombc_abund$gene_description, levels = gene_levels)

# Plot
order_palette2 <- order_palette
order_palette2["Sirenia/Proboscidea"] <- order_palette2["Sirenia"]

p <- ggplot(ancombc_abund, aes(x = Common.name, y = Abundance, colour = Order, fill = diet.general)) +
    geom_boxplot(alpha = 0.8, size = 0.5) +
    scale_colour_manual(values = order_palette2, name = "Order") +
    scale_fill_manual(values = diet_palette, name = "Diet") +
    facet_wrap(~ paste(paste(str_trunc(as.character(gene_description), 49), gene), label, sep = "\n"), ncol = 3, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 8),
          legend.position = "bottom") + ylab("CLR-transformed abundances") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave(p, filename = file.path(subdir, "ancombc_abundances.png"), width = 12, height = 40)

##################################
#### COMPARE ANCOMBC-MCMCglmm ####
##################################

diffabund_comparison <- full_join(mcmc_res_label, ancombc_res_label, by = "gene", suffix = c("_mcmc", "_ancombc"))

write.csv(diffabund_comparison, file = file.path(subdir, "ancombc_mcmcglmm_comparison.csv"), quote = FALSE, row.names = FALSE)
