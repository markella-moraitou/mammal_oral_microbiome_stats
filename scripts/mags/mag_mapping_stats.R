##### MAG presence-absence #####

#### Identify which MAGs are present in which samples using mapping information ####

#### LOAD PACKAGES ####
library(ape)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata
magdir <- normalizePath(file.path("..", "..", "output", "mags")) # Directory with MAG analysis output
subdir <- normalizePath(file.path("..", "..", "output", "mags", "mag_mapping_stats")) # subdirectory for the output of this script
taxdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Directory with output from community-level taxonomic analysis

if(!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))

#######################
#####  LOAD INPUT #####
#######################

# MAG metadata
bac_meta <- read.table(file.path(magdir, "bac_meta.tsv"), sep="\t", header=TRUE)
ar_meta <- read.table(file.path(magdir, "ar_meta.tsv"), sep="\t", header=TRUE)

# MAG trees
bac_tree <- read.tree(file = file.path(magdir, "bac_tree.tree"))
ar_tree <- read.tree(file = file.path(magdir, "ar_tree.tree"))

# Sample metadata
sample_meta <- read.table(file.path(taxdir, "metadata_all.tsv"), sep="\t", header=TRUE)

# List of dereplicated bins
drep_bins <- read.table(file.path(indir, "dereplicated_bins_list.txt"), header=FALSE, sep="\t") %>% pull(V1)

# Mapping stats
map_stats <- read.table(file.path(indir, "mapping_stats_bins.txt"), sep="\t", header=TRUE)

# Host phylogeny
host_trees <- read.nexus(file.path(indir, "mammal_vertlife.nex"))

##############################
#### GET HOST TREE LABELS ####
##############################

# Get consensus tree and fix tip labels
host_consensus <- consensus(host_trees)

# Change tip labels to match metadata
host_consensus$tip.label <- host_consensus$tip.label %>%
                            gsub(pattern="Tapirus_indicus", replacement="Acrocodia_indica") %>%
                            gsub(pattern="Procolobus_badius", replacement="Piliocolobus_foai") %>%
                            gsub(pattern="Otaria_bryonia", replacement="Otaria flavescens") %>%
                            gsub(pattern="_", replacement=" ", fixed=TRUE)

host_labels <- rev(host_consensus$tip.label)
sus <- which(host_labels == "Sus scrofa")
host_labels <- append(host_labels[1:sus], c("Sus domesticus")) %>%
                append(host_labels[(sus + 1):length(host_labels)]) %>%
                append(c("Extraction blank", "Library blank", "Environmental control"))

######################################
#### PLOT COMPLETENESS & IDENTITY ####
######################################

#### Collect metadata ####

# Combine bin metadata and keep only dereplicated bins
bin_meta <- rbind(bac_meta %>% mutate(domain="Bacteria"), ar_meta %>% mutate(domain="Archaea")) %>%
  filter(bin %in% drep_bins) %>%
  filter(Completeness >= 90 & Contamination <= 5)

write.csv(bin_meta, file = file.path(subdir, "hq_mag_metadata.csv"), row.names = FALSE)

# Add metadata to mapping stats
map_stats_meta <- map_stats %>%
    # Get host species info
    left_join(sample_meta[,c("Ext.ID", "Species", "Genus", "Order")],
              by = c("sample" = "Ext.ID"), relationship = "many-to-one") %>%
    # Order host species by order
    mutate(Species = factor(Species, levels = host_labels)) %>%
    # Get bin label and completeness
    inner_join(bin_meta[,c("bin", "label", "phylum", "Completeness", "Contamination")]) %>%
    # Label the sample where a bin was assembled from
    mutate(assembly_sample = case_when(sample == str_remove(str_remove(bin, "^.*-"), "\\..*$") ~ TRUE,
                                       TRUE ~ FALSE)) %>%
    rename(host_species = Species, host_genus = Genus, host_order = Order)

# Order bins by appearance in the tree
bin_order <- append(bac_tree$tip.label, ar_tree$tip.label)
bin_order <- bin_order[bin_order %in% map_stats_meta$label]

map_stats_meta$label <- factor(map_stats_meta$label, levels = bin_order)

#### Completeness per gene set and phylum

p <- ggplot(data = filter(bac_meta, is.tip), aes(y = Completeness, x = phylum)) +
        geom_boxplot()

ggsave(file.path(subdir, "completeness_per_phylum.png"), p, width = 6, height = 4)

# Keep only HQ MAGs as well as identities above 90%
map_stats_hq <- map_stats_meta %>% 
    filter(identity >= 90) 

#### Plot coverage - identity - number of reads ####

# Set identity and coverage thresholds based on mappings stats in assembly samples
min_id <- 98
min_cov <- 75
min_reads <- 10^4

p1 <- ggplot(data = map_stats_hq, aes(x = mapped_reads, y = coverage, colour = identity)) +
        geom_point(alpha = 0.5, size = 0.1) +
        geom_point(data = filter(map_stats_hq, assembly_sample), alpha = 1, size = 0.2, colour = "red") +
        scale_x_log10() +
        geom_vline(xintercept = min_reads, linewidth = 0.5, linetype = "dashed", color = "black") +
        geom_hline(yintercept = min_cov, linewidth = 0.5, linetype = "dashed", color = "black") +
        scale_colour_viridis_b(name = "Identity", option = "plasma") +
        theme(legend.position = "top", legend.text = element_text(size = 6, angle = 45, hjust = 1))

p2 <- ggplot(data = map_stats_hq, aes(x = identity, y = coverage, colour = mapped_reads)) +
        geom_point(alpha = 0.5, size = 0.1) +
        geom_point(data = filter(map_stats_hq, assembly_sample), alpha = 1, size = 0.2, colour = "red") +
        scale_x_log10() +
        geom_vline(xintercept = min_id, linewidth = 0.5, linetype = "dashed", color = "black") +
        scale_colour_viridis_b(name = "Mapped reads", option = "mako", trans = "log10") +
        theme(legend.position = "top", legend.text = element_text(size = 6, angle = 45, hjust = 1))

p3 <- ggplot(data = map_stats_hq, aes(x = coverage, fill = assembly_sample)) +
        geom_histogram(alpha=0.5) +
        geom_vline(xintercept = min_cov, linewidth = 0.5, linetype = "dashed", color = "black") +
        scale_fill_manual(name = "Assembly sample", values = c("TRUE" = "red", "FALSE" = "grey")) +
        theme(legend.position = "top")

p4 <- ggplot(data = map_stats_hq, aes(x = identity, fill = assembly_sample)) +
        geom_histogram(alpha=0.5) +
        geom_vline(xintercept = min_id, linewidth = 0.5, linetype = "dashed", color = "black") +
        scale_fill_manual(name = "Assembly sample", values = c("TRUE" = "red", "FALSE" = "grey")) +
        theme(legend.position = "top")

p <- plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv", axis = "tblr",
               rel_heights = c(1,1), rel_widths = c(1,1))

ggsave(file.path(subdir, "hq_mag_mapping_stats.png"), plot = p, width = 12, height = 12)

#### Summarise mapping stats per host species ####
map_stats_per_host <- map_stats_hq %>%
    group_by(host_species, host_order, label) %>%
    summarise(mean_identity = mean(identity),
              median_identity = median(identity),
              mean_coverage = mean(coverage),
              median_coverage = median(coverage),
              mean_reads = mean(mapped_reads),
              median_reads = median(mapped_reads),
              assembly_species = any(assembly_sample),
              max_identity = max(identity),
              max_coverage = max(coverage),
              max_reads = max(mapped_reads))

write.csv(map_stats_per_host, file = file.path(subdir, "hq_mag_mapping_stats_per_host.csv"), row.names = FALSE)

#### Plot mean and median stats per bin per host species ####
# Filter only good mappings
map_stats_per_host_filt <-
        map_stats_per_host %>% filter(mean_identity > 97 & mean_reads > 1000)

# Plot mean stats
p <- ggplot(data = map_stats_per_host_filt, aes(x = host_species, y = label)) +
    geom_point(aes(size = mean_reads, colour = mean_coverage)) +
    geom_point(data = filter(map_stats_per_host, assembly_species), shape = 8, size = 3, color = "black") +
    facet_grid(cols = vars(host_order), scales = "free_x", space = "free_x") +
    scale_size_continuous(name = "Mean mapped reads", trans = "log10") +
    scale_colour_viridis_b(name = "Mean coverage", option = "plasma") +
    theme(legend.position = "top", legend.direction = "vertical",
         axis.text.y = element_text(size = 6), axis.text.x = element_text(hjust = 1)) +
    xlab("Host species") + ylab("MAG (high-quality: >90% completeness, <5% contamination)")

ggsave(file.path(subdir, "hq_mag_cov_and_id_mean.png"), plot = p, width = 10, height = 40)

# Plot median stats
map_stats_per_host_filt <-
        map_stats_per_host %>% filter(median_identity > 97 & median_reads > 1000)

p <- ggplot(data = map_stats_per_host_filt, aes(x = host_species, y = label)) +
    geom_point(aes(size = median_reads, colour = median_coverage)) +
    geom_point(data = filter(map_stats_per_host, assembly_species), shape = 8, size = 3, color = "black")  +
    facet_grid(cols = vars(host_order), scales = "free_x", space = "free_x") +
    scale_size_continuous(name = "Median mapped reads", trans = "log10") +
    scale_colour_viridis_b(name = "Median coverage", option = "plasma") +
    theme(legend.position = "top", legend.direction = "vertical",
         axis.text.y = element_text(size = 6), axis.text.x = element_text(hjust = 1)) +
    xlab("Host species") + ylab("MAG (high-quality: >90% completeness, <5% contamination)")

ggsave(file.path(subdir, "hq_mag_cov_and_id_median.png"), plot = p, width = 10, height = 40)

# Plot max stats
map_stats_per_host_filt <-
        map_stats_per_host %>% filter(max_identity > 97 & max_reads > 1000)

p <- ggplot(data = map_stats_per_host_filt, aes(x = host_species, y = label)) +
    geom_point(aes(size = max_reads, colour = max_coverage)) +
    geom_point(data = filter(map_stats_per_host, assembly_species), shape = 8, size = 3, color = "black")  +
    facet_grid(cols = vars(host_order), scales = "free_x", space = "free_x") +
    scale_size_continuous(name = "Max mapped reads", trans = "log10") +
    scale_colour_viridis_b(name = "Max coverage", option = "plasma") +
    theme(legend.position = "top", legend.direction = "vertical",
         axis.text.y = element_text(size = 6), axis.text.x = element_text(hjust = 1)) +
    xlab("Host species") + ylab("MAG (high-quality: >90% completeness, <5% contamination)")

ggsave(file.path(subdir, "hq_mag_cov_and_id_max.png"), plot = p, width = 10, height = 40)

#### Define presence-absence of MAGs per host species ####

mag_pres_per_sp <- 
            # Keep only mappings passing the thresholds
            map_stats_hq %>% filter((identity >= min_id & coverage >= min_cov & mapped_reads >= min_reads) | assembly_sample) %>%
            group_by(host_species, host_genus, host_order, label, bin) %>%
            # Get maximum identity, coverage and mapped reads per species. Also indicate if the MAG was assembled from that species
            summarise(max_identity = max(identity),
                      max_coverage = max(coverage),
                      max_reads = max(mapped_reads),
                      assembly_species = any(assembly_sample))

write.csv(mag_pres_per_sp, file = file.path(subdir, "hq_mag_presence_per_host.csv"), row.names = FALSE)

mag_pres_summ <- mag_pres_per_sp %>% group_by(label) %>%
  summarise(species = n_distinct(host_species),
            genera = n_distinct(host_genus),
            orders = n_distinct(host_order)) %>% arrange(desc(species))

write.csv(mag_pres_summ, file = file.path(subdir, "hq_mag_presence_summary.csv"), row.names = FALSE)

# Plot
p <- ggplot(data = mag_pres_per_sp, aes(x = host_species, y = label)) +
    geom_point(aes(size = max_reads, colour = max_coverage)) +
    geom_point(data = filter(mag_pres_per_sp, assembly_species), shape = 8, size = 3, color = "black")  +
    facet_grid(cols = vars(host_order), scales = "free_x", space = "free_x") +
    scale_size_continuous(name = "Max mapped reads", trans = "log10") +
    scale_colour_viridis_b(name = "Max coverage", option = "plasma") +
    theme(legend.position = "top", legend.direction = "vertical",
         axis.text.y = element_text(size = 6), axis.text.x = element_text(hjust = 1)) +
    xlab("Host species") + ylab("MAG (high-quality: >90% completeness, <5% contamination)")

ggsave(file.path(subdir, "hq_mag_pres_abs.png"), plot = p, width = 10, height = 40)

