##### NITRATE REDUCTION
#### Specifically test if there are differences in nitrate reduction genes
#### Between high and low altitude species

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(microbiome)
library(ggnewscale)
library(stringr)

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

gene_str <- read.table(file.path(datadir, "gene_abundance_stratified_modified.tsv"),
                      quote = "", comment.char = "", header = TRUE, sep = "\t")

###############################
#### EXTRACT DATA AND PLOT ####
###############################

# Genes needed to convert nitrate to nitrite, and nitrite to nitric oxide

nitrate_reduc_kos <- c("narX" = "K07673", "narG, narZ, nxrA" = "K00370", "narH, narY, nxrB" = "K00371", "narI, narV"="K00374", "narJ, narW" = "K00373",
                     "napA"="K02567", "napB"="K02568", "napD" = "K02570", "napE" = "K02571",  "napF" = "K02572"
                     #"modA" = "K02020", "mogA" = "K03831", "moaB" = "K03638", "moaE" = "K03635", "moeA" = "K03750", "moeB" = "K21029", "moaA" = "K03639"
                     #"nirK"="K00368", "nirS"="K15864"
                     )

nitrate_reduc_abund <- phy_gene_f %>% transform("compositional") %>% prune_taxa(nitrate_reduc_kos, .) %>%
      psmelt %>% select(Sample, OTU, Abundance, Species, Common.name, diet.general, Locality, gene_name, module) %>%
      rename(KO = OTU)

nitrate_reduc_abund$gene_code <- factor(names(nitrate_reduc_kos)[match(nitrate_reduc_abund$KO, nitrate_reduc_kos)],
                                        level = names(nitrate_reduc_kos))

write.csv(nitrate_reduc_abund, file.path(subdir, "nitrate_reduction_gene_abundance.csv"), row.names = FALSE)

nitrate_reduc_abund$label <- paste(nitrate_reduc_abund$KO, nitrate_reduc_abund$gene_code, sep="\n")

# Keep specific species pair for hypothesis testing
high_low_pairs <- c("Alpine marmot", "Lowland paca", "Giraffe", "Okapi", "Eastern gorilla", "Western gorilla", "Argali", "Chamois", "Domestic sheep")

nitrate_reduc_filt <-
    nitrate_reduc_abund %>% filter(Common.name %in% high_low_pairs) %>%
    # Out of Eastern gorillas, keep only mountain gorillas
    filter(Common.name != "Eastern gorilla" | Locality %in% c("Mount Mikeno", "Virunga Massif")) %>%
    mutate(category = case_when(Common.name %in% c("Alpine marmot", "Eastern gorilla", "Argali", "Chamois") ~ "Altitude-adapted",
                                Common.name == "Giraffe" ~ "Hypertension-adapted",
                                Common.name %in% c("Lowland paca", "Okapi", "Western gorilla", "Domestic sheep") ~ "Baseline"),
           clade = case_when(Common.name %in% c("Alpine marmot", "Lowland paca") ~ "Rodents",
                             Common.name %in% c("Giraffe", "Okapi") ~ "Giraffids",
                             Common.name %in% c("Eastern gorilla", "Western gorilla") ~ "Gorillas",
                             Common.name %in% c("Argali", "Chamois", "Domestic sheep") ~ "Caprinae")) %>%
    mutate(category = factor(category, levels = c("Baseline", "Altitude-adapted", "Hypertension-adapted")),
           clade = factor(clade, levels = c("Gorillas", "Rodents", "Caprinae", "Giraffids"))) %>%
    arrange(category, clade) 

# Also get sum of nitrate reductases abundances
nitrate_reduc_sum <-
    nitrate_reduc_filt %>% group_by(Sample, Species, Common.name, diet.general, category, clade) %>%
    summarise(Abundance = sum(Abundance), label = "All nitrate\nreductases")

nitrate_reduc_all <- rbind(select(nitrate_reduc_filt, c(Sample, Species, Common.name, diet.general, category, clade, Abundance, label)),
                            nitrate_reduc_sum)

sample_n <- nitrate_reduc_all %>% select(Sample, Common.name, category) %>% distinct() %>%
            group_by(Common.name, category) %>% summarise(n = n())

# Perform one-sided Wilcoxon test for each gene_code and clade
wilcox_results <- nitrate_reduc_all %>%
    group_by(label, clade) %>%
    mutate(increase = ifelse(mean(Abundance[category != "Baseline"]) > mean(Abundance[category == "Baseline"]), TRUE, FALSE)) %>%
    group_by(label, clade) %>%
    arrange(category) %>%
    mutate(
        W = wilcox.test(Abundance ~ category)$statistic,
        p_value = wilcox.test(Abundance ~ category)$p.value,
        y = max(Abundance),
        x = n_distinct(Species)/2 + 0.4,
        .groups = "drop"
    ) %>%
    select(label, clade, increase, W, p_value, y, x) %>% unique

wilcox_results$adjusted_p_value <- p.adjust(wilcox_results$p_value, method = "holm")

wilcox_results$signif <- ifelse(wilcox_results$adjusted_p_value < 0.05, 2, ifelse(wilcox_results$p_value < 0.05, 1, 0))
wilcox_results$increase <- ifelse(wilcox_results$signif>0, wilcox_results$increase, "")

# Save Wilcoxon test results to a CSV file
write.csv(wilcox_results, file.path(subdir, "wilcox_test_results.csv"), row.names = FALSE)

p <- ggplot(nitrate_reduc_all, aes(x = Common.name, y = Abundance*100, fill = category, colour = category)) +
    geom_boxplot() +
    scale_fill_manual(values = c("Altitude-adapted" = "#2A4D6E", "Hypertension-adapted" = "#AA4639", "Baseline" = "grey70"), name = "") +
    scale_color_manual(values = c("Altitude-adapted" = "#133453", "Hypertension-adapted" = "#802115", "Baseline" = "grey50"), name = "") +
    new_scale_fill() +
    new_scale_colour() +
    geom_point(data = wilcox_results, aes(x = x, y = y*130, shape = increase, fill = increase, size = signif), inherit.aes = FALSE) +
    scale_fill_manual(values = c("TRUE" = "#33BF33", "FALSE" = "#B00707"), name = "", labels = c("Increase", "Decrease")) +
    scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 25), name = "", labels = c("Increase", "Decrease")) +
    scale_size(range = c(-1, 3), name = "significance", breaks = c(1,2), labels = c("p-val < 0.05", "p-adj < 0.05")) +
    ylim(0, NA) + scale_y_continuous(breaks = c(0.1, 0.5)) +
    facet_grid(cols = vars(clade), rows = vars(label), scales = "free") +
    theme(legend.position = "top", legend.direction = "vertical", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text.y = element_text(angle = 0, size = 8)) +
    labs(y = "Rel. abundance (%)", x = "", fill = "", colour = "") +
    guides(size = guide_legend(override.aes = list(shape = 16)),
           fill = guide_legend(override.aes = list(size = 2)))

ggsave(p, filename = file.path(subdir, "nitrate_reduction_genes.png"), width = 8, height = 9)

# Plot contributions of different taxa
nitrate_abund_str <- gene_str %>% filter(gene_id %in% nitrate_reduc_kos & Sample %in% phy_gene_f@sam_data$Ext.ID) %>%
    select(Sample, gene_id, genus, mapped_reads) %>% rename(KO = gene_id, Ext.ID = Sample) %>%
    left_join(select(rownames_to_column(data.frame(phy_gene_f@sam_data), "Sample"), c(Sample, Ext.ID))) %>%
    left_join(select(nitrate_reduc_filt, Species, Sample, KO, gene_code)) %>%
    filter(!is.na(genus))

write.csv(nitrate_abund_str, file.path(subdir, "nitrate_reduction_gene_stratified.csv"), row.names = FALSE)

# Keep only genes with significant comparison and the same species pairs as in the previous plot
signif_genes <- wilcox_results %>% filter(signif > 0) %>% pull(label) %>% str_remove("\n.*")
nitrate_abund_signif <- nitrate_abund_str %>% filter(Sample %in% nitrate_reduc_filt$Sample & KO %in% signif_genes) %>%
    mutate(Species = factor(Species, levels = unique(nitrate_reduc_filt$Species)))

# Identify top contributing taxa
top_taxa <- nitrate_abund_signif %>%
    group_by(genus) %>% filter(genus != "no support") %>%
    summarise(total_mapped_reads = sum(mapped_reads), .groups = "drop") %>%
    slice_max(order_by = total_mapped_reads, n = 8)

nitrate_abund_grouped <- nitrate_abund_signif %>%
    mutate(genus_grouped = factor(ifelse(genus %in% top_taxa$genus, genus, "Other"), levels = c(top_taxa$genus, "Other"))) %>%
    group_by(Sample, Species, gene_code, genus_grouped) %>%
    summarise(total_mapped_reads = sum(mapped_reads), .groups = "drop") %>%
    group_by(Species, gene_code) %>%
    # get relative abundances
    mutate(total_mapped_reads = total_mapped_reads / sum(total_mapped_reads) * 100)

p <- ggplot(nitrate_abund_grouped, aes(x = gene_code, y = total_mapped_reads, fill = genus_grouped)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#EC8904", "#9DDE04", "#0C5B98", "#B80375", "#ECDE04", "#EC2C04", "#560F9F", "#03AB42", "grey"), name = "") +
    facet_grid(rows = vars(Species)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            strip.text.y = element_text(angle = 0),
            legend.position = "top", legend.direction = "horizontal") +
    labs(x = "", y = "Total mapped reads", fill = "Genus")

ggsave(p, filename = file.path(subdir, "nitrate_reduction_gene_stratified.png"), width = 8, height = 9)
