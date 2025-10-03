##### MAG presence-absence #####

#### Identify which MAGs are present in which samples using mapping information ####

#### LOAD PACKAGES ####
library(ape)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
#library(ggtree)
#library(ggtreeExtra)
#library(ggnewscale)
#library(RColorBrewer)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "mags", "mag_mapping_stats")) # subdirectory for the output of this script
taxdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Directory with output from community-level taxonomic analysis

if(!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))

#######################
#####  LOAD INPUT #####
#######################

# MAG metadata
bac_meta <- read.table(file.path(subdir, "bac_meta.tsv"), sep="\t", header=TRUE)
ar_meta <- read.table(file.path(subdir, "ar_meta.tsv"), sep="\t", header=TRUE)

# MAG trees
bac_tree <- read.tree(file = file.path(subdir, "bac_tree.tree"))
ar_tree <- read.tree(file = file.path(subdir, "ar_tree.tree"))

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
                            gsub(pattern="Equus_quagga", replacement="Equus_burchellii") %>%
                            gsub(pattern="Procolobus_badius", replacement="Piliocolobus_foai") %>%
                            gsub(pattern="Otaria_bryonia", replacement="Otaria_byronia") %>%
                            gsub(pattern="_", replacement=" ", fixed=TRUE)

host_labels <- rev(host_consensus$tip.label)
sus <- which(host_labels == "Sus scrofa")
host_labels <- append(host_labels[1:sus], c("Sus scrofa domesticus")) %>%
                append(host_labels[(sus + 1):length(host_labels)]) %>%
                append(c("Extraction blank", "Library blank", "Environmental control"))

######################################
#### PLOT COMPLETENESS & IDENTITY ####
######################################

#### Collect metadata ####

# Combine bin metadata and keep only dereplicated bins
bin_meta <- rbind(bac_meta %>% mutate(domain="Bacteria"), ar_meta %>% mutate(domain="Archaea")) %>%
  filter(bin %in% drep_bins)

# Add metadata to mapping stats
map_stats_meta <- map_stats %>%
    # Get host species info
    left_join(sample_meta[,c("Ext.ID", "Species", "Order")],
              by = c("sample" = "Ext.ID"), relationship = "many-to-one") %>%
    # Order host species by order
    mutate(Species = factor(Species, levels = host_labels)) %>%
    # Get bin label and completeness
    left_join(bin_meta[,c("bin", "label", "Completeness", "Contamination")]) %>%
    # Label the sample where a bin was assembled from
    mutate(assembly_sample = case_when(sample == str_remove(str_remove(bin, "^.*-"), "\\..*$") ~ TRUE,
                                       TRUE ~ FALSE)) %>%
    rename(host_species = Species)

# Order bins by appearance in the tree
bin_order <- append(bac_tree$tip.label, ar_tree$tip.label)
bin_order <- bin_order[bin_order %in% map_stats_meta$label]

map_stats_meta$label <- factor(map_stats_meta$label, levels = bin_order)

# Keep only HQ MAGs as well as identities above 90%
map_stats_hq <- map_stats_meta %>% 
    filter(Completeness >= 90 & Contamination <= 5 & identity >= 90) 

#### Plot coverage - identity - number of reads ####
p1 <- ggplot(data = map_stats_hq, aes(x = mapped_reads, y = coverage, colour = identity)) +
        geom_point(alpha = 0.5, size = 0.1) +
        geom_point(data = filter(map_stats_hq, assembly_sample), alpha = 1, size = 0.2, colour = "red") +
        scale_x_log10() +
        scale_colour_viridis_b(name = "Identity", option = "plasma") +
        theme(legend.position = "top", legend.text = element_text(size = 6, angle = 45, hjust = 1))

p2 <- ggplot(data = map_stats_hq, aes(x = identity, y = coverage, colour = mapped_reads)) +
        geom_point(alpha = 0.5, size = 0.1) +
        geom_point(data = filter(map_stats_hq, assembly_sample), alpha = 1, size = 0.2, colour = "red") +
        scale_x_log10() +
        geom_vline(xintercept = 98, linewidth = 0.5, linetype = "dashed", color = "black") +
        scale_colour_viridis_b(name = "Mapped reads", option = "mako", trans = "log10") +
        theme(legend.position = "top", legend.text = element_text(size = 6, angle = 45, hjust = 1))

p <- plot_grid(p1, p2, ncol = 1, align = "v", axis = "lr",
               rel_heights = c(1,1))

ggsave(file.path(subdir, "hq_mag_mapping_stats.png"), plot = p, width = 6, height = 12)

#### Summarise mapping stats per host species ####
map_stats_per_host <- map_stats_hq %>%
    group_by(host_species, label) %>%
    summarise(mean_identity = mean(identity),
              median_identity = median(identity),
              mean_coverage = mean(coverage),
              median_coverage = median(coverage),
              mean_reads = mean(mapped_reads),
              median_reads = median(mapped_reads),
              assembly_species = any(assembly_sample))

#### Plot mean and median stats per bin per host species ####
# Keep only id > 98% or mapping stats for the assembly host species
map_stats_per_host_filt <-
        map_stats_per_host %>% filter(mean_identity > 98)

# Plot mean stats
p <- ggplot(data = map_stats_per_host_filt, aes(x = host_species, y = label)) +
    geom_point(aes(size = mean_reads, colour = mean_coverage)) +
    geom_point(data = filter(map_stats_per_host, assembly_species), shape = 8, size = 3, color = "black") +
    scale_size_continuous(name = "Mean mapped reads", trans = "log10") +
    scale_colour_viridis_b(name = "Mean coverage", option = "plasma") +
    theme(legend.position = "top", legend.direction = "vertical")

ggsave(file.path(subdir, "hq_mag_cov_and_id_mean.png"), plot = p, width = 10, height = 40)

# Plot median stats
map_stats_per_host_filt <-
        map_stats_per_host %>% filter(median_identity > 98)

p <- ggplot(data = map_stats_per_host_filt, aes(x = host_species, y = label)) +
    geom_point(aes(size = median_reads, colour = median_coverage)) +
    geom_point(data = filter(map_stats_per_host, assembly_species), shape = 8, size = 3, color = "black") +
    scale_size_continuous(name = "Median mapped reads", trans = "log10") +
    scale_colour_viridis_b(name = "Median coverage", option = "plasma") +
    theme(legend.position = "top", legend.direction = "vertical")

ggsave(file.path(subdir, "hq_mag_cov_and_id_median.png"), plot = p, width = 10, height = 40)
