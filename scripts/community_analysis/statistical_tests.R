##### STATISTICAL TESTS #####

#### Run PERMANOVA, differential abundance analysis and other tests on filtered data using

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
library(car)
library(lme4)
library(tibble)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis"))
subdir <- normalizePath(file.path(outdir, "statistical_tests"))
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

#########################
#### PREP TEST INPUT ####
#########################

phy_genus <- phy_sp_f %>% tax_glom("genus")
phy_genus_clr <- phy_genus %>% transform("clr")

taxa_names(phy_genus_clr) <- phy_genus_clr@tax_table[, "genus"]

# Combine Sirenia/Proboscidea in one clades
phy_genus_clr@sam_data$Order <- ifelse(phy_genus_clr@sam_data$Order %in% c("Sirenia", "Proboscidea"), "Sirenia/Proboscidea", as.character(phy_genus_clr@sam_data$Order))

# Turn habitat to factor
phy_genus_clr@sam_data$habitat.general <- factor(phy_genus_clr@sam_data$habitat.general, levels = c("Terrestrial", "Marine"))

# Test for ruminant differences
phy_genus_clr@sam_data$ruminant <- factor(ifelse(phy_genus_clr@sam_data$digestion == "Ruminant", "Ruminant", "Other"), levels = c("Other", "Ruminant"))

# Test for hypsodont differences
phy_genus_clr@sam_data$hypsodont <- factor(ifelse(grepl("hyps", phy_genus_clr@sam_data$molar_category), "Hypsodont", "Other"), levels = c("Other", "Hypsodont"))

# OTU table
otu <- phy_genus_clr@otu_table %>% data.frame %>%
    rownames_to_column("OTU") %>% pivot_longer(cols = where(is.numeric), names_to = "Sample", values_to = "Abundance")

# Metadata
sam <- phy_genus_clr@sam_data %>% data.frame %>%
    select(Species, Order, diet.general, habitat.general, cf, cp, nfe, ee)%>% rownames_to_column("Sample")

# Combine data
data <- left_join(otu, sam, relationship = "many-to-one") %>%
  # Have sus scrofa as sus scrofa domesticus
  mutate(Species = case_when(Species == "Sus scrofa domesticus" ~ "Sus scrofa",
                             TRUE ~ Species))

data_wide <- data %>% pivot_wider(names_from = "OTU", values_from = "Abundance") %>% data.frame
colnames(data_wide) <- gsub(colnames(data_wide), pattern = ".", replacement = "_", fixed = TRUE)

# Get a metrics for diet
data_diet <- data_wide %>% select(Species, diet_general, cp, cf, ee, nfe) %>% unique %>%
  mutate(animalivory = log((cp + ee)/(cf+nfe)),
         frugivory = log(nfe/cf))

p <- ggplot(data_diet, aes(x = animalivory, y = frugivory, colour = diet_general)) +
  geom_point(size = 2) + theme_bw() +
  geom_text(aes(label = Species)) +
  scale_color_manual(values = diet_palette)

ggsave(p, filename = file.path(subdir, "diet_variables.png"), width = 8, height = 8)

# Add to phyloseq
data_diet <- rbind(data_diet, data_diet[which(data_diet$Species == "Sus scrofa"),])
data_diet[nrow(data_diet), "Species"] <- "Sus scrofa domesticus"
phy_genus_clr@sam_data$animalivory <- data_diet$animalivory[match(phy_genus_clr@sam_data$Species, data_diet$Species)]
phy_genus_clr@sam_data$frugivory <- data_diet$frugivory[match(phy_genus_clr@sam_data$Species, data_diet$Species)]

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
           nitt = 100000,
           burnin = 30000,
           thin = 100)
    }, mc.cores = 1)
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
  separate(term, into = c("OTU", "term"), sep = ":", fill = "right") %>%
  mutate(term = case_when(is.na(term) ~ "Intercept",
                          TRUE ~ term))

random_results <- summary(m1)$Gcovariances %>%
  data.frame %>% rownames_to_column("term") %>%
  mutate(term = str_remove(term, "trait")) %>%
  separate(term, into = c("OTU", "term"), sep = "[.]") %>%
  mutate(pMCMC = NA)

results <- rbind(fixed_results, random_results)
write.csv(results, file = file.path(subdir, "mcmcglmm_results.csv"), quote = FALSE, row.names = FALSE)

significant_results <- results %>% filter((l.95..CI < 0 & u.95..CI < 0) | (l.95..CI > 0 & u.95..CI > 0))

p <-
  ggplot(filter(significant_results, term %in% c("animalivory", "frugivory", "habitat_generalMarine")),
        aes(x = OTU, y = post.mean, fill = term, group = term)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  #geom_errorbar(aes(ymin = `l.95..CI`, ymax = `u.95..CI`), width = 0.1, position = "dodge") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("habitat_generalMarine" = "#553DCB", "animalivory" = "#FF6D5A", "frugivory" = "#B4F582"))

ggsave(p, filename = file.path(subdir, "mcmcglmm_significant_results.png"), width = 8, height = 5)

#####################
#### RUN ANCOMBC ####
#####################

phy_genus@sam_data <- phy_genus_clr@sam_data

# Check if output is there, if not run
if(file.exists(file.path(subdir, "ancombc_results.rds"))) {
    cat("Output already there, so analysis will not be repeated! Loading existing ANCOMBC results...\n")
        out <- readRDS(file.path(subdir, "ancombc_results.rds"))
} else {
    out <- ancombc2(data = phy_genus,
                 fix_formula = "animalivory + frugivory + Order + habitat.general + ruminant",
                 rand_formula = "(1|Species)",
                 tax_level = "species",
                 p_adj_method = "holm", prv_cut = 0.1, # Some parameters. prv_cut=0.10 means that taxa found in less than 10% of samples are ignored
                 group="Order",
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

res_long <- lfc %>% full_join(pvalues) %>% full_join(qvalues) %>% full_join(sensitivity) %>% filter(term != "(Intercept)")

write.csv(res_long, file = file.path(subdir, "ancomb_long.csv"), quote = FALSE, row.names = FALSE)

# Volcano plot
p <- filter(res_long) %>%
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

res_plot <- filter(res_long) %>% filter(qval <= 0.05 & sensitivity == TRUE) %>%
            # label association (positive and negative)
            mutate(assoc = case_when(lfc < 0 ~ paste0(term, "-"),
                                     lfc > 0 ~ paste0(term, "+"))) %>%
            # Then summarise all associations per taxon
            group_by(taxon) %>%
            summarise(label = paste(assoc, collapse=" "))

# Get abundances per sample for the differentially abundant taxa
abundances <- phy_sp_f_clr@otu_table %>% t %>% data.frame %>% rownames_to_column("Sample") %>%
              pivot_longer(cols = -Sample, names_to = "taxon", values_to = "Abundance") %>%
              left_join(select(data.frame(phy_ancombc@sam_data), new_name, Common.name, Order, diet.general), by = c("Sample" = "new_name")) %>%
              mutate(taxon = gsub(taxon, pattern = ".", replacement = " ", fixed = TRUE))

res_abundances <- abundances %>% right_join(res_plot) %>%
                    # Get genus
                    mutate(genus = str_remove(taxon, " .*"))

# Reorder host species
species_levels <- phy_ancombc@sam_data %>% data.frame %>% arrange(as.character(Order), as.character(digestion), Common.name) %>% select(Order, Common.name) %>% unique
res_abundances$Common.name <- factor(res_abundances$Common.name, levels = species_levels$Common.name)

# Reorder taxa
taxa_levels <- res_plot %>% arrange(label) %>% pull(taxon) %>% unique
res_abundances$taxon <- factor(res_abundances$taxon, levels = taxa_levels)

# Plot
order_palette2 <- order_palette
order_palette2["Sirenia/Proboscidea"] <- order_palette2["Sirenia"]

p <- ggplot(res_abundances, aes(x = Common.name, y = Abundance, colour = Order, fill = diet.general)) +
    geom_boxplot(alpha = 0.8, size = 0.5) +
    scale_colour_manual(values = order_palette2, name = "Order") +
    scale_fill_manual(values = diet_palette, name = "Diet") +
    facet_wrap(~ paste(taxon, label, sep = "\n"), ncol = 5, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 8),
          legend.position = "bottom") + ylab("CLR-transformed abundances") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave(p, filename = file.path(subdir, "ancombc_abundances.png"), width = 12, height = 30)

#### Plot diff abund taxa with habitat, uses, phenotypes ####
habitat_associations <-
            res_plot %>%
            left_join(select(habitats_table, c(taxon, OBT, obt_type, occurences))) %>%
            # Group OBTs into larger habitats
            mutate(OBT = case_when(OBT %in% c("dental plaque", "mouth") ~ "oral",
                                       OBT %in% c("marine water", "deep sea") ~ "marine",
                                       OBT %in% c("mammalian", "wild animal", "mammalian livestock") ~ "animal",
                                       TRUE ~ OBT)) %>%
            # Only select some of the most informative terms
            filter(OBT %in% c("oral", "animal", "soil", "marine", "rumen", "gut"))

use_associations <-
            res_plot %>%
            left_join(select(uses_table, c(taxon, OBT, obt_type, occurences)))

# Combine
associations <- full_join(habitat_associations, use_associations) %>% filter(!is.na(OBT)) %>% arrange(label, taxon)
associations$label <- factor(associations$label, levels = unique(associations$label))

write.csv(associations, file = file.path(subdir, "ancomb_diffabund_associations.csv"), quote = FALSE, row.names = FALSE)

p <- ggplot(associations, aes(y = taxon, x = OBT, size = occurences, colour = OBT)) +
        geom_point() +
        scale_color_manual(values = c("oral" = "#AE1E3D", "animal" = "#BD6E20", "rumen" = "#A4B81F", "gut" = "#BD9F20", "soil" = "#56A71C", "marine" = "#156B73",
                                     "antimicrobial" = "gray50", "fermentative" = "#3A459C", "proteolytic activity" = "#B41F34", "lipolytic activity" = "#BD9120")) +
        facet_grid(rows = vars(label), cols = vars(obt_type), scales = "free", space = "free") +
        theme(strip.text.y = element_text(angle = 0, size = 8), axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8, hjust = 1),
        legend.position = "none")

ggsave(p, filename = file.path(subdir, "ancombc_diffabund_associations.png"), width = 8, height = 8)

##########################
#### MULTIPLE ANOVAS #####
##########################

manyanovas <- function(physeq, formula) {
    phymelt <- psmelt(physeq)
    # Initialize an empty data frame to store results
    results <- data.frame()
    for (taxon in taxa_names(physeq)) {
        cat("Running ANOVA for taxon:", taxon, "\n")
        # Extract the abundance data for the current taxon
        df <- phymelt %>% filter(OTU == taxon)
        # Run the ANOVA
        model <- Anova(lm(formula, data = df), type=2)
        res <- data.frame(taxon = taxon, terms = rownames(model), F.value = model$`F value`, p.value = model$`Pr(>F)`)
        # Add the results to the results data frame
        results <- rbind(results, res)
    }
    return(results)
}

# Get metadata from ancombc
phy_sp_f_clr@sam_data <- phy_ancombc@sam_data

# Run multiple anovas
anovas_res <- manyanovas(phy_sp_f_clr, formula = "Abundance ~ Order + diet.general + habitat.general + ruminant")

# Plot results
anova_res_wide <- anovas_res %>% pivot_wider(names_from = terms, values_from = c(F.value, p.value))

# Add total abundance and genus information
anova_res_wide$total_abundance <- taxa_sums(phy_sp_f)[anova_res_wide$taxon]
anova_res_wide$genus <- as.data.frame(tax_table(phy_sp_f))[anova_res_wide$taxon, "genus"]

# Get grouped genus from OTUs that came up as differentially abundant
#abundant_genera <- anova_res_wide %>% group_by(genus) %>% summarise(total_abundance = sum(total_abundance)) %>%
#    arrange(desc(total_abundance)) %>% slice_head(n = 15) %>% pull(genus)

diff_abund_genera <- res_plot$taxon %>% str_remove(" .*") %>% unique

anova_res_wide$genus_grouped <- ifelse(anova_res_wide$genus %in% diff_abund_genera, anova_res_wide$genus, "Other")
anova_res_wide$genus_grouped <- factor(anova_res_wide$genus_grouped, levels = c(as.character(diff_abund_genera), "Other"))

write.csv(anova_res_wide, file = file.path(subdir, "taxon_anova_results.csv"), row.names = FALSE, quote = TRUE)

# Plot F-value for order and diet
max_axis <- max(c(anova_res_wide$F.value_Order, anova_res_wide$F.value_diet.general))

p <- ggplot(aes(x = F.value_Order, y = F.value_diet.general), data = anova_res_wide) +
    geom_hex(bins = 20) +
    scale_fill_viridis_c("H") +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", alpha = 0.7) +
    facet_wrap(~ genus_grouped, ncol = 4) +
    theme(legend.position = "none") +
    labs(x = "F-value Order", y = "F-value Diet") +
    xlim(NA, max_axis) + ylim(NA, max_axis)

ggsave(p, filename = file.path(subdir, "taxon_anova_order_diet.png"), width = 8, height = 8)

# Plot F-value for diet and habitat
max_axis <- max(c(anova_res_wide$F.value_habitat.general, anova_res_wide$F.value_diet.general))

p <- ggplot(aes(x = F.value_habitat.general, y = F.value_diet.general), data = anova_res_wide) +
    geom_hex(bins = 20) +
    scale_fill_viridis_c("H") +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", alpha = 0.7) +
    facet_wrap(~ genus_grouped, ncol = 4) +
    theme(legend.position = "none") +
    labs(x = "F-value Habitat", y = "F-value Diet") +
    xlim(NA, max_axis) + ylim(NA, max_axis)

ggsave(p, filename = file.path(subdir, "taxon_habitat_diet.png"), width = 8, height = 8)
