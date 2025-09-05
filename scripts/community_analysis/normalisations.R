##### NORMALISATIONS #####

#### Access OTU tree from Open Tree of Life, perform philr and clr normalisations
#### and create subsets for further analysis

################
#### SET UP ####
################

library(dplyr)
library(phyloseq)
library(tidyr)
library(readr)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(microbiome)
library(stringr)
library(philr)
library(microbiome)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis"))
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Load phy
phy_sp_f <- readRDS(file.path(phydir, "phy_sp_f.RDS"))

# Links to download GTDB tree
url <- "https://data.gtdb.ecogenomic.org/releases/release220/220.0/"

#######################
#### GET TAXA TREE ####
#######################

# Download tree and metadata
download.file(paste0(url, "bac120_r220.tree.gz"), file.path(outdir, "bac120_r220.tree.gz"))
options(timeout = 600)
download.file(paste0(url, "bac120_taxonomy_r220.tsv.gz"), file.path(outdir, "bac120_taxonomy_r220.tsv.gz"))

bac_tree <- read.tree(gzfile(file.path(outdir, "bac120_r220.tree.gz")))
bac_meta <- read_tsv(file.path(outdir, "bac120_taxonomy_r220.tsv.gz"), col_names = FALSE)

# Keep only metadata in tree
bac_meta_f <- bac_meta %>% filter(X1 %in% bac_tree$tip.label)

# Split taxonomy column
bac_meta_f <- bac_meta_f %>% separate(X2, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";") %>%
  mutate(s = str_remove(s, "s__")) %>% mutate(s = str_remove(s, "_[A-Z]+"))

# Filter metadata to only include taxa in phyloseq object
bac_taxa <- unique(str_remove(taxa_names(subset_taxa(phy_sp_f, superkingdom != "Archaea")), "\\*"))

table(!bac_taxa %in% bac_meta_f$s)
bac_taxa[!bac_taxa %in% bac_meta_f$s]

# Extract tip labels from the species in the dataset
bac_meta_f <- filter(bac_meta_f, s %in% bac_taxa) %>% select(X1, p, s) %>%
  # pick one representative if there are multiple (they should be monophyletic, so it shouldn't matter)
  group_by(s) %>% slice(1) %>% ungroup()

# Subset tree
tree <- drop.tip(bac_tree, setdiff(bac_tree$tip.label, bac_meta_f$X1))
tree$tip.label <- bac_meta_f$s[match(tree$tip.label, bac_meta_f$X1)]
tree$node.label <- paste0("N", 1:tree$Nnode) # Not very informative for now

write.tree(tree, file = file.path(phydir, "phy_tree.tree"))

#################################
#### CREATE PHYLOSEQ OBJECTS ####
#################################

# CLR transformation
phy_sp_f_clr <- phy_sp_f %>% microbiome::transform("clr")

# PhiLR transformation
phy_sp_philr <- phyloseq(otu_table(phy_sp_f), sample_data(phy_sp_f), tax_table(phy_sp_f), phy_tree(tree))
philr_otu <- philr(phy_sp_philr, pseudocount=10^-5)

write.table(philr_otu, file = file.path(phydir, "phy_sp_philr_OTU.tsv"), sep = "\t", row.names = TRUE, quote = FALSE)

phy_sp_philr <- phyloseq(otu_table(philr_otu, taxa_are_rows = FALSE), sample_data(phy_sp_f))

# Save phyloseq objects
saveRDS(phy_sp_f, file.path(phydir, "phy_sp_f.RDS"))
saveRDS(phy_sp_f_clr, file.path(phydir, "phy_sp_f_clr.RDS"))
saveRDS(phy_sp_philr, file.path(phydir, "phy_sp_philr.RDS"))

########################
#### CREATE SUBSETS ####
########################

# Create function to mirror subset of phyloseq object on other transformed objects
mirror_subset <- function(phy_subset, phy_norm, taxa = FALSE) {
  phy_norm_new <- phy_norm %>% prune_samples(sample_names(phy_norm) %in% sample_names(phy_subset), .)
  if (taxa) {
    cat("Also subsetting taxa\n")
    phy_norm_new <- phy_norm_new %>% prune_taxa(taxa_names(phy_norm) %in% taxa_names(phy_subset), .)
  }
  return(phy_norm_new)
}

## Wild animals only
phy_wild <- phy_sp_f %>% subset_samples(!(Species %in% c("Ovis aries", "Equus caballus", "Sus scrofa domesticus")))
phy_wild <- phy_wild %>% subset_taxa(taxa_sums(phy_wild) > 0)

phy_wild_clr <- mirror_subset(phy_subset=phy_wild, phy_norm=phy_sp_f_clr, TRUE)
phy_wild_philr <- mirror_subset(phy_subset=phy_wild, phy_norm=phy_sp_philr, FALSE)

## Marine terrestrial comparisons
phy_habitat <- phy_sp_f %>%
  subset_samples(Species %in% c("Otaria byronia", "Meles meles", "Ursus arctos", "Orcinus orca", "Hippopotamus amphibius", "Sus scrofa", "Dugong dugon", "Loxodonta africana"))
phy_habitat <- phy_habitat %>% subset_taxa(taxa_sums(phy_habitat) > 0)

phy_habitat_clr <- mirror_subset(phy_subset=phy_habitat, phy_norm=phy_sp_f_clr, TRUE)
phy_habitat_philr <- mirror_subset(phy_subset=phy_habitat, phy_norm=phy_sp_philr, FALSE)

# Artiodactyla
phy_artio <- phy_sp_f %>%
  subset_samples(Order == "Artiodactyla")
phy_artio <- phy_artio %>% subset_taxa(taxa_sums(phy_artio) > 0)

phy_artio_clr <- mirror_subset(phy_subset=phy_artio, phy_norm=phy_sp_f_clr, TRUE)
phy_artio_philr <- mirror_subset(phy_subset=phy_artio, phy_norm=phy_sp_philr, FALSE)

# Carnivora
phy_carni <- phy_sp_f %>%
  subset_samples(Order == "Carnivora")
phy_carni <- phy_carni %>% subset_taxa(taxa_sums(phy_carni) > 0)

phy_carni_clr <- mirror_subset(phy_subset=phy_carni, phy_norm=phy_sp_f_clr, TRUE)
phy_carni_philr <- mirror_subset(phy_subset=phy_carni, phy_norm=phy_sp_philr, FALSE)

# Primates
phy_prim <- phy_sp_f %>%
  subset_samples(Order == "Primates")
phy_prim <- phy_prim %>% subset_taxa(taxa_sums(phy_prim) > 0)

phy_prim_clr <- mirror_subset(phy_subset=phy_prim, phy_norm=phy_sp_f_clr, TRUE)
phy_prim_philr <- mirror_subset(phy_subset=phy_prim, phy_norm=phy_sp_philr, FALSE)

# Planned contrasts in more deeply sampled Order (Artiodactyla, Carnivora, Primates)
phy_deep <- phy_sp_f %>%
  subset_samples(Order %in% c("Artiodactyla", "Carnivora", "Primates", "Perissodactyla"))
phy_deep <- phy_deep %>% subset_taxa(taxa_sums(phy_deep) > 0)

phy_deep_clr <- mirror_subset(phy_subset=phy_deep, phy_norm=phy_sp_f_clr, TRUE)
phy_deep_philr <- mirror_subset(phy_subset=phy_deep, phy_norm=phy_sp_philr, FALSE)

# Save phyloseq objects
for (obj in c("phy_wild", "phy_habitat", "phy_artio", "phy_carni", "phy_prim", "phy_deep",
               "phy_wild_clr", "phy_habitat_clr", "phy_artio_clr", "phy_carni_clr", "phy_prim_clr", "phy_deep_clr",
               "phy_wild_philr", "phy_habitat_philr", "phy_artio_philr", "phy_carni_philr", "phy_prim_philr", "phy_deep_philr")) {
  saveRDS(get(obj), file.path(phydir, paste0(obj, ".RDS")))
}

###################
#### PLOT TREE ####
###################

# Plot microbial tree and colour tips by phylum
tree_meta <- bac_meta_f %>% filter(s %in% tree$tip.label) %>% rename(tip.label = s) %>%
  mutate(p = str_remove(p, "p__") %>% str_remove("_[A-Z]+")) %>% select(tip.label, p)

p <- ggtree(tree)  %<+% tree_meta +
  geom_tiplab(size = 1, aes(colour = p), name = "Phylum") +
  scale_color_manual(values = phylum_palette) +
  ggtitle("Bacterial phylogenetic tree for OTUs in dataset") +
  theme(legend.position = "top") +
  guides(colour = guide_legend(nrow = 3, override.aes = list(size=3)))

# Add average abundance as barplot 
abund_df <- as.data.frame(otu_table(phy_sp_f)) %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "Sample", values_to = "Abundance") %>%
  group_by(OTU) %>%
  summarise(Mean_Abundance = mean(Abundance)) %>%
  rename(tip.label = OTU)

p <- p + geom_fruit(data = abund_df, geom = geom_bar,
                    mapping = aes(x = Mean_Abundance, y = tip.label),
                    orientation = "y", stat = "identity",
                    pwidth = 0.3, size = 0.1, fill = "darkgrey")

ggsave(p, filename = file.path(outdir, "phy_tree.png"), width = 5, height = 30, units = "in", dpi = 300)

# Remove downloaded files to save space
file.remove(file.path(outdir, "bac120_r220.tree.gz"))
file.remove(file.path(outdir, "bac120_taxonomy_r220.tsv.gz"))
