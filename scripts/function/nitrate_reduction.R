##### NITRATE REDUCTION
#### Specifically test if there are differences in nitrate reduction genes
#### Between high and low altitude species

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(phyloseq)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "function")) # subdirectory for the output of this script
datadir <- normalizePath(file.path(outdir, "data")) # Directory with data files
subdir <- normalizePath(file.path(outdir, "nitrate_reduction")) # Subdirectory for nitrate reduction analysis

here <- getwd() # Get current working directory

dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

# Get ordination functions
source(file.path("..","ordination_functions.R"))

#######################
#####  LOAD INPUT #####
#######################

# Phyloseq object
phy_gene_f <- readRDS(file.path(outdir, "data", "phy_gene_f.RDS"))
phy_gene_f_clr <- readRDS(file.path(outdir, "data", "phy_gene_f_clr.RDS"))

###############################
#### EXTRACT DATA AND PLOT ####
###############################

# Genes needed to convert nitrate to nitrite, and nitrite to nitric oxide

nitrate_reduc_kos <- c("narG"="K00370", "narH"="K00371", "narI"="K00374", "narX"="K07673", "narJ" = "K00373",
                     "napA"="K02567", "napB"="K02568", 
                     "nirK"="K00368", "nirS"="K15864")

nitrate_reduc_abund <- phy_gene_f_clr %>% prune_taxa(nitrate_reduc_kos, .) %>%
      psmelt %>% select(Sample, OTU, Abundance, Species, Common.name, diet.general, Locality, gene_name, module) %>%
      rename(KO = OTU)

nitrate_reduc_abund$gene_code <- factor(names(nitrate_reduc_kos)[match(nitrate_reduc_abund$KO, nitrate_reduc_kos)],
                                        level = names(nitrate_reduc_kos))

write.csv(nitrate_reduc_abund, file.path(subdir, "nitrate_reduction_gene_abundance.csv"), row.names = FALSE)

# Keep specific species pair for hypothesis testing
high_low_pairs <- c("Marmot", "Paca", "Reticulated giraffe", "Okapi", "Eastern gorilla", "Western gorilla")

nitrate_reduc_filt <-
    nitrate_reduc_abund %>% filter(Common.name %in% high_low_pairs) %>%
    # Out of Eastern gorillas, keep only mountain gorillas
    filter(Common.name != "Eastern gorilla" | Locality %in% c("Mount Mikeno", "Virunga Massif")) %>%
    mutate(category = case_when(Common.name %in% c("Marmot", "Eastern gorilla") ~ "Altitude-adapted",
                                Common.name == "Reticulated giraffe" ~ "Hypertension-adapted",
                                Common.name %in% c("Paca", "Okapi", "Western gorilla") ~ "Baseline"),
           clade = case_when(Common.name %in% c("Marmot", "Paca") ~ "Rodents",
                             Common.name %in% c("Reticulated giraffe", "Okapi") ~ "Giraffids",
                             Common.name %in% c("Eastern gorilla", "Western gorilla") ~ "Gorillas")) %>%
    mutate(category = factor(category, levels = c("Baseline", "Altitude-adapted", "Hypertension-adapted")),
           clade = factor(clade, levels = c("Gorillas", "Rodents", "Giraffids"))) %>%
    arrange(category, clade) %>%
    mutate(Common.name = factor(gsub(" ", "\n", Common.name), levels = unique(gsub(" ", "\n", Common.name))))

sample_n <- nitrate_reduc_filt %>% select(Sample, Common.name, category) %>% distinct() %>%
            group_by(Common.name, category) %>% summarise(n = n())

# Perform one-sided Wilcoxon test for each gene_code and clade
wilcox_results <- nitrate_reduc_filt %>%
    group_by(gene_code, clade) %>%
    arrange(category) %>%
    summarise(
        W = wilcox.test(Abundance ~ category, alternative = "less")$statistic,
        p_value = wilcox.test(Abundance ~ category, alternative = "less")$p.value,
        .groups = "drop"
    )

wilcox_results$adjusted_p_value <- p.adjust(wilcox_results$p_value, method = "holm")

wilcox_results$signif <- ifelse(wilcox_results$adjusted_p_value < 0.05, "***",
                                ifelse(wilcox_results$p_value < 0.05, "*", NA))

# Save Wilcoxon test results to a CSV file
write.csv(wilcox_results, file.path(subdir, "wilcox_test_results.csv"), row.names = FALSE)

p <- ggplot(nitrate_reduc_filt, aes(x = Common.name, y = Abundance, fill = category, colour = category)) +
    geom_boxplot() +
    geom_text(data = wilcox_results, aes(x = 1.5, y = 10, label = signif), size = 6, colour = "#FF9900", inherit.aes = FALSE) +
    scale_fill_manual(values = c("Altitude-adapted" = "#2A4D6E", "Hypertension-adapted" = "#AA4639", "Baseline" = "grey70")) +
    scale_color_manual(values = c("Altitude-adapted" = "#133453", "Hypertension-adapted" = "#802115", "Baseline" = "grey50")) +
    facet_grid(cols = vars(clade), rows = vars(gene_code), scales = "free") +
    scale_y_continuous(breaks = c(0, 10)) +
    theme(legend.position = "top", legend.direction = "vertical", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text.y = element_text(angle = 0)) +
    labs(y = "CLR abundance", x = "", fill = "", colour = "")

ggsave(p, filename = file.path(subdir, "nitrate_reduction_genes.png"), width = 4, height = 8)
