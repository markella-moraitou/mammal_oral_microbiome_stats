##### CLADES FOR POPGEN #####

#### Get codiversifying clades of MAGs to select lineages for population genomics analyses.
#### including reference genomes from GTDBtk

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(stringr)
library(ape)
library(phytools)
library(tibble)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "mags")) 
subdir <- normalizePath(file.path(outdir, "clades_for_popgen")) # subdirectory for the output of this script

# Create output directory
dir.create(subdir, showWarnings = FALSE)

source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))

setwd(subdir)

#######################
#####  LOAD INPUT #####
#######################

# MAG bacteria tree
bac_tree <- read.tree(file = file.path(outdir, "bac_tree.tree"))

# Entire bacteria tree
bac_tree_full <- read.tree(file = file.path(indir, "gtdbtk.bac120.decorated_itol.tree"))

# Bacterial MAG metadata
bac_meta <- read.table(file.path(outdir, "bac_meta.tsv"), sep="\t", header=TRUE)

# Codiversifying clades
bac_cod <- read.csv(file.path(outdir, "codiversification", "bac_cod_results.csv"), header=TRUE)

##############################
#### GET EXTENDED CLADES  ####
##############################

# Extract tips from codiversifying clades
cod_nodes <- bac_cod %>% filter(p.adjust < 0.05) %>% pull(node)

cod_tips <- lapply(cod_nodes, function(node) {
  tips <- getDescendants(bac_tree, Ntip(bac_tree) + which(bac_tree$node.label==node))
  tips <- tips[tips <= length(bac_tree$tip.label)]
  tip_names <- bac_tree$tip.label[tips]
  return(tip_names)
})

names(cod_tips) <- cod_nodes

cod_tips_df <- enframe(cod_tips, name = "node", value = "tip") %>%
  unnest(cols = c(tip)) %>%
  # Get MAG file names
  left_join(unique(select(bac_meta, c(label, bin, Completeness, Contamination))), by = c("tip" = "label"))

# Now extract clades with reference genomes from GTDBtk tree
cod_tips_full <- data.frame(node = character(), tip = character(), is.mag = logical(), mag_label = character(), stringsAsFactors = FALSE)

for (n in unique(cod_tips_df$node)) {
  # Get bins from that node
  bins <- cod_tips_df %>% filter(node == n) %>% pull(bin) %>% str_remove(".gz")
  # Extract the corresponding clade from the full tree including reference genomes
  node <- getMRCA(bac_tree_full, bins)
  tips <- getDescendants(bac_tree_full, node)
  tips <- tips[tips <= length(bac_tree_full$tip.label)]
  tip_names <- bac_tree_full$tip.label[tips] %>% str_remove("GB_|RS_")
  df <- data.frame(node = n, tip = tip_names, is.mag = tip_names %in% bins)
  df$tip <- ifelse(df$is.mag, paste0(df$tip, ".gz"), df$tip)
  df$mag_label <- cod_tips_df$tip[match(df$tip, cod_tips_df$bin, ".gz")]
  # Make sure to keep only HQ bins
  df$Completeness <- cod_tips_df$Completeness[match(df$tip, cod_tips_df$bin, ".gz")]
  df$Contamination <- cod_tips_df$Contamination[match(df$tip, cod_tips_df$bin, ".gz")]
  df <- df %>% filter((Completeness >= 90 & Contamination <= 5) | !is.mag )
  df <- df %>% select(-c(Completeness, Contamination))
  write.csv(df, file = file.path(subdir, paste0("clade_", n, "_with_refs.csv")), row.names = FALSE, quote = FALSE)
  cod_tips_full <- rbind(cod_tips_full, df)
}

# Summarise
cod_tips_summary <- cod_tips_full %>%
  group_by(node) %>%
  summarise(num_tips = n(),
            num_mags = sum(is.mag),
            num_refs = num_tips - num_mags)

write.csv(cod_tips_summary, file = "clades_with_refs_summary.csv", row.names = FALSE, quote = FALSE)
