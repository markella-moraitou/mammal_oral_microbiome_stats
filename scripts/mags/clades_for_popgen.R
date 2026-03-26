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
library(rlang)

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
bac_tree <- read.tree(file = file.path(outdir, "bac_tree_drep.tree"))

# Entire bacteria tree
bac_tree_full <- read.tree(file = file.path(indir, "gtdbtk.bac120.decorated_itol.tree"))

# Bacterial MAG metadata
bac_meta <- read.table(file.path(outdir, "bac_meta_drep.tsv"), sep="\t", header=TRUE)

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
  left_join(unique(select(bac_meta, c(label, bin, Completeness, Contamination, host_order))), by = c("tip" = "label"))

#### Get info on reference genomes ####
ref_info <- select(bac_meta, closest_genome_reference.1, species) %>% rename(accession = closest_genome_reference.1) %>%
  rbind(select(bac_meta, closest_placement_reference, species) %>% rename(accession = closest_placement_reference) ) %>%
  filter(!is.na(species)) %>% unique

other_ref <- bac_meta$other_related_references.genome_id.species_name.radius.ANI.AF. %>% unique %>%
  str_split("; ") %>% lapply(function(x) x[1:2]) %>% lapply(function(x) data.frame(accession = x[1], species = x[2])) %>% bind_rows() %>%
  mutate(species = str_remove(species, "s__")) %>% filter(!is.na(species)) %>% unique

ref_info <- rbind(ref_info, other_ref) %>% unique

write.csv(ref_info, file = file.path(subdir, "reference_genomes_info.csv"), row.names = FALSE, quote = FALSE)

# Now extract clades with reference genomes from GTDBtk tree
cod_tips_full <- data.frame(node = character(), tip = character(), is.mag = logical(), label = character(), stringsAsFactors = FALSE)

forward_branches <- list("Propionibacterium" = "host_order=='Artiodactyla'",
                         "Actinomyces.6" = "host_order=='Artiodactyla'",
                         "Corynebacterium.8" = "host_order=='Carnivora'")

for (n in unique(cod_tips_df$node)) {
  print(n)
  # Get bins from that node
  bins <- cod_tips_df %>% filter(node == n) %>% pull(bin) %>% str_remove(".gz")
  # Extract the corresponding clade from the full tree including reference genomes
  node <- getMRCA(bac_tree_full, bins)
  tips <- getDescendants(bac_tree_full, node)
  tips <- tips[tips <= length(bac_tree_full$tip.label)]
  tip_names <- bac_tree_full$tip.label[tips] %>% str_remove("GB_|RS_")
  df <- data.frame(node = n, tip = tip_names, is.mag = grepl("MEGAHIT", tip_names))
  df$tip <- ifelse(df$is.mag, paste0(df$tip, ".gz"), df$tip)
  # Add metadata
  df <- df %>% left_join(cod_tips_df, by = c("node" = "node", "tip" = "bin")) %>%
    rename(label = tip.y) %>%
    # Keep only dereplicated MAGs (i.e. those that are in the tree)
    filter(!is.na(label) | !is.mag)
  # Make sure to keep only HQ bins
  df <- df %>% filter((Completeness >= 90 & Contamination <= 5) | !is.mag )
  # Add reference genome labels
  df$label <- ifelse(df$is.mag, df$label,
                     ref_info$species[match(df$tip, ref_info$accession)])
  df$forward_branch <- NA
  # Define forward branches when relevant
  if (n %in% names(forward_branches)) {
    fg_tips <- cod_tips_df %>% filter(node == n) %>%
          filter(!!parse_expr((forward_branches[[n]]))) %>% pull(bin) %>% str_remove(".gz")
    fg_node <- getMRCA(bac_tree_full, fg_tips)
    fg_tips_all <- getDescendants(bac_tree_full, fg_node)
    fg_tips_all <- fg_tips_all[fg_tips_all <= length(bac_tree_full$tip.label)]
    fg_tip_names <- bac_tree_full$tip.label[fg_tips_all] %>% str_remove("GB_|RS_")
    df <- df %>% mutate(forward_branch = case_when(
      str_remove(tip, ".gz") %in% fg_tip_names ~ "#1",
      TRUE ~ "#0"
    ))
  }
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
