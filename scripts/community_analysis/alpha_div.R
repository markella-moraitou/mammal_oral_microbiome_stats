##### ALPHA DIVERSITY #####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(stringr)
library(vegan)
library(rphylopic)
library(ape)
library(picante)
library(phyr)
library(phytools)
library(MCMCglmm)
library(parallel)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis"))
subdir <- normalizePath(file.path(outdir, "alpha_diversity")) # subdirectory for the output of this script
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

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

# Bacterial phylogeny (GTDB)
bac_tree <- read.tree(file.path(phydir, "phy_tree.tree"))

# Host phylogeny
host_consensus <- read.tree(file.path(outdir, "host_consensus.tre"))

############################
#### RAREFACTION CURVES ####
############################

set.seed(1)

#### Filtered dataset ####

rare_results <- data.frame(Sample = character(),
                           S = numeric(),
                           se = numeric(),
                           subsample = numeric(),
                           stringsAsFactors = FALSE)

max <- 10^floor(log10(mean(sample_sums(phy_sp_f))))
step <- max / 20

for (s in seq(0, max, by=step)) {
  cat("Rarefying to", s, "sequences per sample...\n")
  rare <- rarefy(data.frame(phy_sp_f@otu_table), MARGIN = 2, sample = s, se = FALSE) %>% data.frame()
  colnames(rare) <- "S"
  rare <- rare %>% rownames_to_column(var = "Sample")
  rare$subsample <- s
  # Append to results
  rare_results <- rbind(rare_results, rare)
}

# Add metadata
metadata <- data.frame(Sample = sample_names(phy_sp_f),
                        lib_size = sample_sums(phy_sp_f),
                        Species = phy_sp_f@sam_data$Species,
                        Order_grouped = phy_sp_f@sam_data$Order_grouped,
                        stringsAsFactors = FALSE)

# Don't show subsamples larger than the achieved library size
rare_results_filt <- rare_results %>% left_join(metadata, by = "Sample") %>%
  filter(subsample <= lib_size)

write.csv(rare_results_filt, file.path(subdir, "rarefaction_results_filt.csv"), row.names = FALSE)

species_medians <-
  # Get the average species richness at max subsample
  rare_results_filt %>%
  group_by(Species, Order_grouped) %>%
  filter(subsample == max(rare_results_filt$subsample)) %>%
  summarise(median_S = mean(S)) %>%
  # add phylopic info
  left_join(phylopics)

# Plot rarefaction curves
p <- ggplot(rare_results_filt) +
  geom_line(aes(x = subsample, y = S, group = Sample, colour = Species)) +
  geom_phylopic(data = species_medians,
                 aes(x = max-step*2, y = median_S + 5, uuid = uid,
                     height = max(rare_results_filt$S)/2, color = Species),
                 alpha = 0.8, vjust = 0, hjust = 0) +
  scale_color_manual(values=species_palette, name = "Order") +
  facet_grid(~ Order_grouped) +
  theme(legend.position = "none") +
  xlab("Number of sequences sampled") +
  ylab("Observed species richness") +
  theme(legend.position = "none")

ggsave(file.path(subdir, "rarefaction_curves_filt.png"), p, width=8, height=8)

#### Raw dataset ####

rare_results <- data.frame(Sample = character(),
                           S = numeric(),
                           se = numeric(),
                           subsample = numeric(),
                           stringsAsFactors = FALSE)

max <- 10^floor(log10(mean(sample_sums(phy_sp))))
step <- max / 20

for (s in seq(0, max, by=step)) {
  cat("Rarefying to", s, "sequences per sample...\n")
  rare <- rarefy(data.frame(phy_sp@otu_table), MARGIN = 2, sample = s, se = FALSE) %>% data.frame()
  colnames(rare) <- "S"
  rare <- rare %>% rownames_to_column(var = "Sample")
  rare$subsample <- s
  # Append to results
  rare_results <- rbind(rare_results, rare)
}

# Add metadata
metadata <- data.frame(Sample = sample_names(phy_sp),
                        lib_size = sample_sums(phy_sp),
                        Species = phy_sp@sam_data$Species,
                        Order_grouped = phy_sp@sam_data$Order_grouped,
                        stringsAsFactors = FALSE)

# Don't show subsamples larger than the achieved library size
rare_results_filt <- rare_results %>% left_join(metadata, by = "Sample") %>%
  filter(subsample <= lib_size)

write.csv(rare_results_filt, file.path(subdir, "rarefaction_results_raw.csv"), row.names = FALSE)

species_medians <-
  # Get the average species richness at max subsample
  rare_results_filt %>%
  group_by(Species, Order_grouped) %>%
  filter(subsample == max(rare_results_filt$subsample)) %>%
  summarise(median_S = mean(S)) %>%
  # add phylopic info
  left_join(phylopics)

# Plot rarefaction curves
p <- ggplot(rare_results_filt) +
  geom_line(aes(x = subsample, y = S, group = Sample, colour = Species)) +
  geom_phylopic(data = species_medians,
                 aes(x = max-step*2, y = median_S + 5, uuid = uid,
                     width = max(rare_results_filt$S)/2, color = Species),
                 alpha = 0.8, vjust = 0, hjust = 0) +
  scale_color_manual(values=species_palette, name = "Order") +
  facet_grid(~ Order_grouped) +
  theme(legend.position = "none") +
  xlab("Number of sequences sampled") +
  ylab("Observed species richness") +
  theme(legend.position = "none")

ggsave(file.path(subdir, "rarefaction_curves_raw.png"), p, width=8, height=8)

#########################
#### ALPHA DIVERSITY ####
#########################

set.seed(1)

# Rarefy to the depth suggested by rarefaction curves
rar_level <- 200000
phy_sp_rarefied <- rarefy_even_depth(subset_samples(phy_sp_f, sample_sums(phy_sp_f) > rar_level), sample.size = rar_level, rngseed = 1)

alpha_div <- data.frame(estimate_richness(phy_sp_f, measures = c("Observed"))) %>%
  rownames_to_column(var = "Sample") %>% rename(filt = Observed) %>%
  left_join(data.frame(estimate_richness(phy_sp_rarefied, measures = c("Observed"))) %>%
              rownames_to_column(var = "Sample") %>% rename(filt_rarefied = Observed),
            by = "Sample")

# Also raw data
rar_level <- 500000
phy_sp_raw_rarefied <- rarefy_even_depth(subset_samples(phy_sp, sample_sums(phy_sp) > rar_level), sample.size = rar_level, rngseed = 1)

alpha_div <- left_join(alpha_div,
                      data.frame(estimate_richness(phy_sp, measures = c("Observed"))) %>%
                      rownames_to_column(var = "Sample") %>% rename(raw = Observed) %>%
              left_join(data.frame(estimate_richness(phy_sp_raw_rarefied, measures = c("Observed"))) %>%
                      rownames_to_column(var = "Sample") %>% rename(raw_rarefied = Observed),
                      by = "Sample"))

# Add metadata
alpha_div <- alpha_div %>%
  left_join(data.frame(phy_sp_f@sam_data) %>%
              rownames_to_column(var = "Sample") %>%
              select(Sample, Species, Common.name, Order, Order_grouped, diet.general, cp, cf, habitat.general, digestion),
            by = c("Sample"))

write.csv(alpha_div, file = file.path(subdir, "alpha_diversity.csv"), quote = FALSE, row.names = FALSE)

# Create a named vector for the relabeling
order_labels <- c("Rodentia" = "Rod.", "Carnivora" = "Carniv.", "Perissodactyla" = "Perissod.")

# Filtered
p <- ggplot(alpha_div, aes(x=Species, y=filt)) +
  geom_boxplot(aes(fill=diet.general)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=diet_palette, name = "Species") +
  scale_x_discrete(labels = setNames(phy_sp_f@sam_data$Common.name, phy_sp_f@sam_data$Species)) +
  facet_grid(Order_grouped ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(Order_grouped = as_labeller(order_labels, default = label_value))) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Observed species richness") +
  coord_flip()

ggsave(file.path(subdir, "alpha_diversity_filt.png"), p, width=8, height=10)

# Filtered & Rarefied
p <- ggplot(alpha_div, aes(x=Species, y=filt_rarefied)) +
  geom_boxplot(aes(fill=diet.general)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=diet_palette, name = "Species") +
  scale_x_discrete(labels = setNames(phy_sp_f@sam_data$Common.name, phy_sp_f@sam_data$Species)) +
  facet_grid(Order_grouped ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(Order_grouped = as_labeller(order_labels, default = label_value))) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Observed species richness (after rarefaction)") +
  coord_flip()

ggsave(file.path(subdir, "alpha_diversity_filt_rarefied.png"), p, width=8, height=10)

# Raw
p <- ggplot(alpha_div, aes(x=Species, y=raw)) +
  geom_boxplot(aes(fill=diet.general)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=diet_palette, name = "Species") +
  scale_x_discrete(labels = setNames(phy_sp@sam_data$Common.name, phy_sp@sam_data$Species)) +
  facet_grid(Order_grouped ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(Order_grouped = as_labeller(order_labels, default = label_value))) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Observed species richness") +
  coord_flip()

ggsave(file.path(subdir, "alpha_diversity_raw.png"), p, width=8, height=10)

# Raw & Rarefied
p <- ggplot(alpha_div, aes(x=Species, y=raw_rarefied)) +
  geom_boxplot(aes(fill=diet.general)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=diet_palette, name = "Species") +
  scale_x_discrete(labels = setNames(phy_sp@sam_data$Common.name, phy_sp@sam_data$Species)) +
  facet_grid(Order_grouped ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(Order_grouped = as_labeller(order_labels, default = label_value))) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Observed species richness (after rarefaction)") +
  coord_flip()

ggsave(file.path(subdir, "alpha_diversity_raw_rarefied.png"), p, width=8, height=10)

####################
#### FAITH'S PD ####
####################

bac_tree$tip.label <- gsub("_", " ", bac_tree$tip.label)

# Calculate Faith's Phylogenetic Diversity
otu_table <- t(as.matrix(phy_sp_f@otu_table))

# Get phylogenetic diversity
phy_div <- pd(otu_table, bac_tree) %>% rownames_to_column(var = "Sample") %>%
  left_join(data.frame(phy_sp_f@sam_data) %>%
              rownames_to_column(var = "Sample") %>%
              select(Sample, Species, Common.name, Order, Order_grouped, diet.general, cp, cf, habitat.general, digestion),
            by = c("Sample"))

write.csv(phy_div, file = file.path(subdir, "phylogenetic_diversity.csv"), quote = FALSE, row.names = FALSE)

# Plot phylogenetic diversity
p <- ggplot(phy_div, aes(x=Species, y=PD)) +
  geom_boxplot(aes(fill=diet.general)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=diet_palette, name = "Species") +
  scale_x_discrete(labels = setNames(phy_sp@sam_data$Common.name, phy_sp@sam_data$Species)) +
  facet_grid(Order_grouped ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(Order_grouped = as_labeller(order_labels, default = label_value))) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Faith's PD") +
  coord_flip()

ggsave(file.path(subdir, "phylogenetic_diversity.png"), p, width=8, height=10)

# Calculate ratio between PD and species richness
phy_div <- phy_div %>% mutate(ratio = PD/SR)

p <- ggplot(phy_div, aes(x=Species, y=ratio)) +
  geom_boxplot(aes(fill=diet.general)) +
  theme(legend.position = "none") +
  scale_fill_manual(values=diet_palette, name = "Species") +
  scale_x_discrete(labels = setNames(phy_sp@sam_data$Common.name, phy_sp@sam_data$Species)) +
  facet_grid(Order_grouped ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(Order_grouped = as_labeller(order_labels, default = label_value))) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  ylab("Faith's PD/Species richness ratio") +
  coord_flip()

ggsave(file.path(subdir, "alpha_vs_pd.png"))

###################
#### RUN TESTS ####
###################

# Fix tree labels
host_consensus$node.label <- paste0("node", c(1:length(host_consensus$node.label)))
host_consensus$tip.label <- gsub("_", " ", host_consensus$tip.label)

#### Run PGLMM ####
# to identify factors affecting alpha diversity and Faith's PD

# Alpha diversity
alpha_div <- alpha_div %>% mutate(Species = case_when(Species == "Sus scrofa domesticus" ~ "Sus scrofa",
                                                       TRUE ~ Species))

alpha_div$ruminant <- factor(ifelse(alpha_div$digestion == "Ruminant", "Ruminant", "Other"), levels = c("Other", "Ruminant"))

model <- pglmm(filt ~ cp + cf + habitat.general + ruminant + (1 | Species__), data = alpha_div, 
              cov_ranef = list(Species = host_consensus), family = "gaussian")

res_alpha <- cbind(model$B, model$B.pvalue) %>% as.data.frame %>%
            rownames_to_column %>% filter(rowname != "(Intercept)") %>%
            mutate(rowname = str_remove(str_remove(rowname, "habitat.general"), "ruminant"))

colnames(res_alpha) <- c("term", "coef", "pval")

res_alpha$response <- "alpha_diversity"

# Phylogenetic diversity

phy_div <- phy_div %>% mutate(Species = case_when(Species == "Sus scrofa domesticus" ~ "Sus scrofa",
                                                       TRUE ~ Species))

phy_div$ruminant <- factor(ifelse(phy_div$digestion == "Ruminant", "Ruminant", "Other"), levels = c("Other", "Ruminant"))

model <- pglmm(PD ~ cp + cf + habitat.general + ruminant + (1 | Species__), data = phy_div, 
              cov_ranef = list(Species = host_consensus), family = "gaussian")

res_phy <- cbind(model$B, model$B.pvalue) %>% as.data.frame %>%
            rownames_to_column %>% filter(rowname != "(Intercept)") %>%
            mutate(rowname = str_remove(str_remove(rowname, "habitat.general"), "ruminant"))

colnames(res_phy) <- c("term", "coef", "pval")

res_phy$response <- "phylogenetic_diversity"

res <- rbind(res_alpha, res_phy)

write.csv(res, file.path(subdir, "pglmm_results.csv"), row.names = FALSE, quote = FALSE)

#### Pagel's lambda on residuals ####

# Alpha diversity

model <- lm(filt ~ cp + cf + habitat.general + ruminant,
            data = alpha_div)

# Extract residuals
resids_df <- data.frame(Sample = alpha_div$Sample,
                        Species = alpha_div$Species,
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
res_alpha <- data.frame(response = "alpha_diversity",
                      lambda = pagel$lambda,
                      pval = pagel$P)

# Phylogenetic diversity
model <- lm(PD ~ cp + cf + habitat.general + ruminant,
            data = phy_div)

# Extract residuals
resids_df <- data.frame(Sample = phy_div$Sample,
                        Species = phy_div$Species,
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
res_phy <- data.frame(response = "phylogenetic_diversity",
                      lambda = pagel$lambda,
                      pval = pagel$P)

res_pagel <- rbind(res_alpha, res_phy)

write.csv(res_pagel, file.path(subdir, "pagel_results.csv"), row.names = FALSE, quote = FALSE)

#### MCMCglmm ####

cat("Running MCMCglmm...\n")
set.seed(14)

# Set priors
prior <- list(G = list(G1 = list(V = 1, nu = 0.002)),
              R = list(V = 1, nu = 0.002))

Ainv <- inverseA(host_consensus, nodes = "TIPS")$Ainv

m <- mclapply(1:10, function(i) {
      MCMCglmm(fixed = filt ~ cp + cf + habitat.general + ruminant,
         random = ~ Species,
         ginverse = list(Species=Ainv),
         prior = prior,
         data = alpha_div,
         verbose = TRUE, # pr = TRUE, pl = TRUE,
         nitt = 51000,
         burnin = 1000,
         thin = 100)
  }, mc.cores = 1)

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

# 95% Credible interval
HPDinterval(m1$VCV)

# Collect results into tables
fixed_results <- summary(m1)$solutions %>%
  data.frame %>% rownames_to_column("term") %>%
  mutate(term = str_remove(term, "habitat.general") %>% str_remove("ruminant")) %>%
  mutate(response = "alpha_diversity")

random_results <- summary(m1)$Gcovariances %>%
  data.frame %>% rownames_to_column("term") %>%
  mutate(pMCMC = NA) %>%
  mutate(response = "alpha_diversity")

# Combine resulrts
mcmc_res <- rbind(fixed_results, random_results) %>%
    group_by(term)

write.csv(mcmc_res, file = file.path(subdir, "mcmcglmm_alpha_results.csv"), quote = FALSE, row.names = FALSE)
