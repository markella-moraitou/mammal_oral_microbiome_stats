##### EXPORT CLADES #####

#### Get codiversifying clades of MAGs including reference genomes from GTDBtk
#### population genomics analyses and further mappings

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
library(ggtree)
library(ggnewscale)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "mags")) 
subdir <- normalizePath(file.path(outdir, "export_clades")) # subdirectory for the output of this script

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

other_ref <- bac_meta$other_related_references.genome_id.species_name.radius.ANI.AF. %>%
  unique %>% lapply(function(x) str_extract_all(x, "GC.*?; s__.*?;")[[1]]) %>%
  lapply(function(x) data.frame(genome = x)) %>% bind_rows() %>% filter(!is.na(genome)) %>%
  separate(genome, sep = "; ", into = c("accession", "species")) %>%
  mutate(species = str_remove(str_remove(species, "s__"), ";")) %>% unique

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

###########################################
#### EXTRA PLOT FOR PROPIONIBACTERIUM #####
###########################################

# Get bins from that node
prop_bins <- cod_tips_df %>% filter(node == "Propionibacterium") %>%
              select(tip, bin, host_order) %>% mutate(bin = str_remove(bin, ".gz")) %>%
              mutate(is.mag = TRUE)

#This time I want the parent node of the MRCA to get the outgroup
prop_node <- getMRCA(bac_tree_full, prop_bins$bin)
parent_node <- getParent(bac_tree_full, prop_node)

prop_tree <- extract.clade(bac_tree_full, parent_node)

# Rename tips
prop_tree$tip.label <- ifelse(prop_tree$tip.label %in% prop_bins$bin, prop_bins$tip[match(prop_tree$tip.label, prop_bins$bin)], 
                              prop_tree$tip.label)
prop_tree$tip.label <- str_remove(prop_tree$tip.label, "GB_|RS_")

# Get info on references
ref_prop <- ref_info %>% filter(accession %in% prop_tree$tip.label) %>% select(accession, species) %>%
  # Add manually
  rbind(data.frame(accession = c("GCF_000940845.1", "GCF_900111365.1"), species = c("Propionibacterium freudenreichii", "Propionibacterium cyclohexanicum"))) %>%
  mutate(is.mag = FALSE, host_order = "Public data") %>%
  mutate(Env = case_when(
    species == "Propionibacterium ruminifibrarum" ~ "cattle rumen",
    species == "Propionibacterium australiense" ~ "bovine lesions",
    species == "Propionibacterium acidifaciens" ~ "human mouth",
    species == "Propionibacterium freudenreichii" ~ "cheese",
    species == "Propionibacterium cyclohexanicum" ~ "spoiled juice"
  )) %>%
  mutate(label = paste0(species, " (", Env, ")"))

prop_tree$tip.label <- ifelse(prop_tree$tip.label %in% ref_prop$accession, ref_prop$label[match(prop_tree$tip.label, ref_prop$accession)],
                              prop_tree$tip.label)

# Drop tip labels starting with MEGAHIT, these were excluded during dereplication
prop_tree <- drop.tip(prop_tree, prop_tree$tip.label[grepl("MEGAHIT", prop_tree$tip.label)])

prop_meta <- full_join(prop_bins, rename(select(ref_prop, -c(species)), tip = label, bin = accession)) %>%
  mutate(host_order = factor(host_order, levels = c(names(order_palette), "Public data")))

p <- ggtree(prop_tree) %<+% prop_meta + 
  geom_tiplab(aes(label = label, colour = is.mag), nudge_x = 0.01, size = 3) +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey40"), labels = c("TRUE" = "MAG", "FALSE" = "Ref"), name = "") +
  new_scale_colour() +
  geom_tippoint(aes(colour = host_order), size = 2) +
  scale_colour_manual(values = append(order_palette, c("Public data" = "black")),
    labels = c("Artiodactyla" = "Ruminant Artiodactyla"), name = "Source") +
  xlim(0,0.8) +
  theme(legend.position = "left", legend.direction = "vertical",
        legend.margin = margin(0,0,0,0), plot.margin = margin(0,0,0,0)) +
  guides(colour = guide_legend(ncol = 1)) +
  geom_treescale(x = 0, y = 20, width = 0.1)

ggsave("tree_propionibacterium.png", plot = p, width = 6, height = 4)

write.csv(select(prop_meta, tip, bin), file = "tree_propionibacterium_bins.csv", row.names = FALSE, quote = FALSE)

#####################################
#### EXTRA PLOT FOR ACTINOMYCES #####
#####################################

# Get bins from that node
act_bins <- cod_tips_df %>% filter(node == "Actinomyces.6") %>%
              select(tip, bin, host_order) %>% mutate(bin = str_remove(bin, ".gz")) %>%
              mutate(is.mag = TRUE)

#This time I want the parent node of the MRCA to get the outgroup
act_node <- getMRCA(bac_tree_full, act_bins$bin)

act_tree <- extract.clade(bac_tree_full, act_node)

# Rename tips
act_tree$tip.label <- ifelse(act_tree$tip.label %in% act_bins$bin, act_bins$tip[match(act_tree$tip.label, act_bins$bin)], 
                              act_tree$tip.label)
act_tree$tip.label <- str_remove(act_tree$tip.label, "GB_|RS_")

# Get info on references
ref_act <- ref_info %>% filter(accession %in% act_tree$tip.label) %>% select(accession, species) %>%
  mutate(Env = case_when(species == "Actinomyces glycerinitolerans" ~ "sheep rumen",
                         species == "Actinomyces succiniciruminis" ~ "cattle rumen",
                         species == "Actinomyces procaprae" ~ "gazelle feces",
                         species == "Actinomyces ruminicola" ~ "cattle rumen",
                         species == "Actinomyces qiguomingii" ~ "antelope feces",
                         species == "Actinomyces ruminicola" ~ "cattle rumen",
                         species == "Actinomyces sp009930875" ~ "antelope feces",
                         species == "Actinomyces sp030530635" ~ "sheep saliva",
                         species == "Actinomyces sp023369555" ~ "rodent orbital abscess",
                         species == "Actinomyces israelii" ~ "",
                         species == "Actinomyces gerencseriae" ~ "human mouth")) %>%
  # Add manually
  rbind(data.frame(accession = c("GCA_937926195.1", "GCF_000159035.1", "GCF_000220835.1",
                                  "GCF_000269805.1", "GCF_000429225.1", "GCF_001553565.1",
                                  "GCF_013184985.2", "GCF_014595995.2", "GCF_015355765.1",
                                  "GCF_029851875.1", "GCF_030524025.1", "GCF_900637165.1"),
                   species = c("Actinomyces sp937926195", "Actinomyces urogenitalis", "Actinomyces sp000220835",
                              "Actinomyces massiliensis", "Actinomyces dentalis", "Actinomyces radicidentis",
                              "Actinomyces faecalis", "Actinomyces respiraculi", "Actinomyces haliotis",
                              "Actinomyces sp029851875", "Actinomyces sp030524025", "Actinomyces howellii"),
                    Env = c("human gut", "human urogenital tract", "human mouth",
                              "human blood", "human mouth", "human mouth",
                              "marmot feces", "marmot respiratory tract", "abalone gut",
                              "ferret", "tapir saliva", "cattle mouth"))) %>%
  mutate(is.mag = FALSE, host_order = "Public data") %>%
  mutate(label = paste0(species, " (", Env, ")"))

act_tree$tip.label <- ifelse(act_tree$tip.label %in% ref_act$accession, ref_act$label[match(act_tree$tip.label, ref_act$accession)],
                              act_tree$tip.label)

# Drop tip labels starting with MEGAHIT, these were excluded during dereplication
act_tree <- drop.tip(act_tree, act_tree$tip.label[grepl("MEGAHIT", act_tree$tip.label)])

act_meta <- full_join(act_bins, rename(select(ref_act, -c(species)), tip = label, bin = accession)) %>%
  mutate(host_order = factor(host_order, levels = c(names(order_palette), "Public data")))

p <- ggtree(act_tree) %<+% act_meta + 
  geom_tiplab(aes(label = label, colour = is.mag), nudge_x = 0.01, size = 3) +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey40"), labels = c("TRUE" = "MAG", "FALSE" = "Ref"), name = "") +
  new_scale_colour() +
  geom_tippoint(aes(colour = host_order), size = 2) +
  scale_colour_manual(values = append(order_palette, c("Public data" = "black")),
    labels = c("Artiodactyla" = "Ruminant Artiodactyla"), name = "Source") +
  xlim(0,0.8) +
  theme(legend.position = "left", legend.direction = "vertical",
        legend.margin = margin(0,0,0,0), plot.margin = margin(0,0,0,0)) +
  guides(colour = guide_legend(ncol = 1)) +
  geom_treescale(x = 0.7, y = 20, width = 0.1)

ggsave("tree_actinomyces.png", plot = p, width = 6, height = 4)

write.csv(select(act_meta, tip, bin), file = "tree_actinomyces_bins.csv", row.names = FALSE, quote = FALSE)
