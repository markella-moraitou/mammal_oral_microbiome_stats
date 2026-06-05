#### PLOTS AND MODELS FOR NON KEGG ANNOTATIONS (MEROPS AND CAZY) ####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(lme4)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
datadir <- normalizePath(file.path("..", "..", "output", "function", "data"))
subdir <- normalizePath(file.path("..", "..", "output", "function", "non_kegg")) # subdirectory for the output of this script

dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Phyloseq objects
phy_gene_f <- readRDS(file.path(datadir, "phy_gene_f.RDS"))
phy_gene_f_clr <- readRDS(file.path(datadir, "phy_gene_f_clr.RDS"))

#############################
#### NON KEGG ANNOTATIOS ####
#############################

set.seed(123)

# Plot peptidase diversity
# Rarefy even depth before calculating richness
phy_merops <- phy_gene_f %>% rarefy_even_depth(2e+07) %>% subset_taxa(database == "MEROPS")

phy_merops_comp <- phy_gene_f %>% transform("compositional") %>% subset_taxa(database == "MEROPS")

merops_richness <- 
    estimate_richness(phy_merops, measures = "Observed") %>%
    rownames_to_column("Sample") %>%
    left_join(select(data.frame(sample_data(phy_merops)), Common.name, Animalivory, Frugivory, digestion, diet.general, Order, Order_grouped, habitat.general) %>%
                as.data.frame() %>% rownames_to_column("Sample")) %>%
    mutate(ruminant = (digestion == "Ruminant"),
           marine = (habitat.general == "Marine"))

write.csv(merops_richness, file.path(subdir, "peptidase_richness.csv"), row.names = FALSE)

model <- glmer(Observed ~ Animalivory + Frugivory + ruminant + marine + (1|Common.name), data = merops_richness, family = poisson)
res <- summary(model)$coefficients

label <- data.frame(res) %>% filter(Pr...z.. < 0.05) %>%
    rownames_to_column("Variable") %>%
    filter(Variable != "(Intercept)") %>%
    mutate(label = paste0(Variable, " coef = ", signif(Estimate, 2), ", p = ", signif(Pr...z.., 2))) %>%
    pull(label) %>%
    paste(collapse = "\n")

p <-
    ggplot(merops_richness, aes(x = as.factor(Animalivory), y = Observed, fill = Order)) +
    geom_boxplot() +
    scale_fill_manual(values = order_palette, name = "Taxonomic order") +
    labs(x = "Animal % in diet", y = "Peptidase richness (Observed)") +
    annotate("text", x = Inf, y = Inf,
              label = label, hjust = 1.1, vjust = 1.1, size = 3) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(file.path(subdir, "peptidase_richness.png"), width=5, height=4)

# Plot peptidase abundance

merops_abundance <- 
    data.frame(phy_merops_comp@otu_table) %>% rownames_to_column("OTU") %>%
    pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Abundance") %>%
    left_join(sample_data(phy_merops_comp) %>% data.frame() %>% rownames_to_column("Sample") %>%
                select(Sample, Common.name, Animalivory, diet.general, Order, Order_grouped)) %>%
    group_by(Common.name, Animalivory, diet.general, Order, Order_grouped) %>%
    # Get mean and quartile abundance of peptidases per species
    summarise(Abundance = mean(Abundance, na.rm = TRUE),
              q1 = quantile(Abundance, probs = 0.25),
              q3 = quantile(Abundance, probs = 0.75), .groups = "drop")

p <-
    ggplot(merops_abundance, aes(x = as.factor(Animalivory), y = Abundance, fill = Order)) +
    geom_boxplot() +
    scale_fill_manual(values = order_palette, name = "Taxonomic order") +
    labs(x = "Animal % in diet", y = "Peptidase relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(file.path(subdir, "peptidase_abundance.png"), width=5, height=4)

# Plot CAZy enzyme diversity

phy_cazy <- phy_gene_f %>% rarefy_even_depth(2e+07) %>% subset_taxa(database == "CAZY")

phy_cazy_comp <- phy_gene_f %>% transform("compositional") %>% subset_taxa(database == "CAZY")

cazy_richness <- 
    estimate_richness(phy_cazy, 1000, measures = "Observed") %>%
    rownames_to_column("Sample") %>%
    left_join(select(data.frame(sample_data(phy_merops)), Common.name, Animalivory, Frugivory, digestion, diet.general, Order, Order_grouped, habitat.general) %>%
                as.data.frame() %>% rownames_to_column("Sample")) %>%
    mutate(ruminant = (digestion == "Ruminant"),
           marine = (habitat.general == "Marine"))

write.csv(cazy_richness, file.path(subdir, "cazy_richness.csv"), row.names = FALSE)

model <- glmer(Observed ~ Animalivory + Frugivory + ruminant + marine + (1|Common.name), data = cazy_richness, family = poisson)
res <- summary(model)$coefficients

label <- data.frame(res) %>% filter(Pr...z.. < 0.05) %>%
    rownames_to_column("Variable") %>%
    filter(Variable != "(Intercept)") %>%
    mutate(label = paste0(Variable, " coef = ", signif(Estimate, 2), ", p = ", signif(Pr...z.., 2))) %>%
    pull(label) %>%
    paste(collapse = "\n")

p <-
    ggplot(cazy_richness, aes(x = as.factor(Animalivory), y = Observed, fill = Order)) +
    geom_boxplot() +
    scale_fill_manual(values = order_palette, name = "Taxonomic order") +
    annotate("text", x = Inf, y = Inf,
              label = label, hjust = 1.1, vjust = 1.1, size = 3) +
    labs(x = "Animal % in diet", y = "CAZy richness (Observed)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(file.path(subdir, "cazy_richness.png"), width=5, height=4)

# Plot CAZy abundance
cazy_abundance <- 
    data.frame(phy_cazy_comp@otu_table) %>% rownames_to_column("OTU") %>%
    pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Abundance") %>%
    left_join(sample_data(phy_cazy_comp) %>% data.frame() %>% rownames_to_column("Sample") %>%
                select(Sample, Common.name, Animalivory, diet.general, Order, Order_grouped)) %>%
    group_by(Common.name, Animalivory, diet.general, Order, Order_grouped) %>%
    # Get mean and quartile abundance of CAZy enzymes per species
    summarise(Abundance = mean(Abundance, na.rm = TRUE),
              q1 = quantile(Abundance, probs = 0.25),
              q3 = quantile(Abundance, probs = 0.75), .groups = "drop")

p <-
    ggplot(cazy_abundance, aes(x = factor(Animalivory), y = Abundance, fill = Order)) +
    geom_boxplot() +
    scale_fill_manual(values = order_palette, name = "Taxonomic order") +
    labs(x = "Animal % in diet", y = "CAZy relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(subdir, "cazy_abundance.png"), width=5, height=4)
