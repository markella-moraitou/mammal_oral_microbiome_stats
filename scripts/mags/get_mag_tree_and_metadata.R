##### Get MAG tree #####

#### Filter GTDB-Tk output to keep only bins, and combine with metadata

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(renv)
library(stringr)
library(tibble)
library(ape)
library(tidytree)
library(tibble)
library(treeio)
library(phytools)
library(rentrez)

#### VARIABLES AND WORKING DIRECTORY ####

indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "mags")) # subdirectory for the output of this script

# Create subdirectory for output
dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

# Get ENTREZ KEY for taxonomic information
options(ENTREZ_KEY = Sys.getenv("API_KEY"))

#######################
#####  LOAD INPUT #####
#######################

# GTDB-Tk de_novo output
bac_tree <- read.tree(file.path(indir, "gtdbtk.bac120.decorated_itol.tree"))
ar_tree <- read.tree(file.path(indir, "gtdbtk.ar53.decorated_itol.tree"))

# bin metadata
mag_meta <- read.csv(file.path(indir, "mq_hq_bin_metadata.csv"), header=TRUE)

# Replace more than one instance of . with a single .
colnames(mag_meta) <- gsub("..", ".", colnames(mag_meta), fixed=TRUE) %>% gsub("...", ".", ., fixed=TRUE)

# Host metadata
host_taxonomy <- read.csv(file.path(indir, "host_taxonomy.csv"), header=TRUE)
host_traits <- read.csv(file.path(indir, "host_traits.csv"), header=TRUE)
host_diet <- read.csv(file.path(indir, "Lintulaakso_diet_filtered.csv"), header=TRUE)

# contig info
contig_info <- read.table(file.path(indir, "contig_info_per_sample.txt"), header=TRUE)

#################
#### PROCESS ####
#################

### Combine host metadata ###
colnames(host_taxonomy) <- paste0("host_", tolower(colnames(host_taxonomy)))

host_meta <- host_taxonomy[, colnames(host_taxonomy)!="host_species"] %>%
    left_join(host_traits, by=c("host_museum.species"="Species")) %>%
    left_join(select(host_diet, c("Species", "calculated_cluster_main_diet")), by=c("host_museum.species"="Species")) %>%
    rename("diet.general"="calculated_cluster_main_diet")

# Get taxon label: the most specific taxonomic level that is not empty and is not a made up GTDB name
get_taxon_label <- function(row) {
    # For each rank, return the first non-empty value that does not contain numbers
    for (rank in c("species", "genus", "family", "order", "class", "phylum")) {
        if (!is.na(row[[rank]]) & !grepl("[0-9]", row[[rank]])) {
            return(row[[rank]])
        }
    }
}

# Combine and clean metadata
meta <-
    # Add host metadata
    mag_meta %>% rename("host_species"="Species") %>%
    left_join(host_meta, by=c("host_species"="host_museum.species"), relationship = "many-to-one") %>%
    left_join(contig_info, by=c("Sample" = "sample"), relationship = "many-to-one") %>%
    separate(classification, sep=";", into=c("domain", "phylum", "class", "order", "family", "genus", "species")) %>%
    mutate(across(domain:species, ~ case_when(str_detect(.x, "^[dpcofgs]__$") ~ NA, TRUE ~ .x))) %>%
    mutate(is.mag=TRUE) %>%
    # Get taxonomy label
    mutate(classification=apply(., 1, get_taxon_label)) %>%
    group_by(classification) %>% mutate(label=paste(gsub(" ", "_", classification), row_number(), sep="_")) %>%
    # Remove *__ prefix
    mutate(across(domain:species, ~ str_remove(.x, "^[dpcofgs]__"))) %>% ungroup() %>%
    rename(min_damage_model_p = min,
                 q1_damage_model_p = q1,
                 mean_damage_model_p = mean,
                 median_damage_model_p = median,
                 q3_damage_model_p = q3,
                 max_damage_model_p = max)

## Get taxonomic ids
get_taxon_id <- function(taxon_name) {
  search_results <- entrez_search(db = "taxonomy", term = paste(taxon_name, "AND (Bacteria[Subtree] OR Archaea[Subtree]) "))
  if (length(search_results$ids) == 1) {
    return(search_results$ids[1])
  } else if (length(search_results$ids) == 0) {
    cat("No results for", taxon_name, "\n")
    return(NA)
  } else {
    cat("Multiple results for", taxon_name, "\n")
    return(NA)
  }
 }

# Get unique classifications
name_to_taxid = data.frame(name = unique(meta$classification))

# Get simplied classifications to facilitate searching the NCBI DB
# e.g. remove underscores after genus and species names, as well as the *__ prefix
name_to_taxid$simple_name <- name_to_taxid$name %>%
                            str_remove("^[dpcofgs]__") %>%
                            str_replace("_.+ ", " ") %>% # Remove underscore plus character(s) before a space
                            str_replace("_.+$", "") # Remove underscore plus character at the end of the string

# Get taxids                         
taxids <- sapply(name_to_taxid$simple_name, get_taxon_id)
name_to_taxid$taxids <- unlist(taxids)[match(name_to_taxid$simple_name, names(unlist(taxids)))]

# Save
write.csv(name_to_taxid, file = file.path(subdir, "name_to_taxid.csv"), row.names=FALSE, quote=FALSE)

# Add ids to metadata
meta$ncbi_taxid <- name_to_taxid$taxids[match(meta$classification, name_to_taxid$name)]

# Filter MAG metadata
bac_meta <- meta %>% filter(domain == "Bacteria")
ar_meta <- meta %>% filter(domain == "Archaea")

# Function to subset tree and add metadata
subset_tree <- function(tree, metadata) {
    # Get tree with only MAGs
    tips.to.drop <- tree$tip.label[!(tree$tip.label %in%  str_remove(meta$bin, ".gz"))]
    pruned_tree <- drop.tip(tree, tips.to.drop)
    # Rename labels
    pruned_tree$tip.label <- meta$label[match(pruned_tree$tip.label, str_remove(meta$bin, ".gz"))]
    #pruned_tree <- unroot(pruned_tree)
    # When the node label is missing, use the node number
    node_labs <- paste0("N", (Ntip(pruned_tree)+1):(Ntip(pruned_tree) + Nnode(pruned_tree)))
    pruned_tree$node.label <- ifelse(pruned_tree$node.label=="", node_labs, pruned_tree$node.label)
    return(pruned_tree)
}

# Bacteria tree
bac_tree_new <- subset_tree(bac_tree, bac_meta)

# Archaea tree
ar_tree_new <- subset_tree(ar_tree, ar_meta)

# Add some traits should be added to the entire clade that shares a common ancestor, not just the tips
add_node_metadata <- function(tree, metadata, traits) {
    # Add metadata to tree
    tree_tibble <- tree %>% as_tibble() %>% mutate(node=1:nrow(.))
    tree_tibble <- tree_tibble %>% left_join(metadata, by="label")
    tree_tibble$is.tip <- ifelse(tree_tibble$node <= Ntip(tree), TRUE, FALSE)
    # Add metadata to nodes that share a common ancestor
    for (tr in traits) {
        node_traits <- data.frame(node=integer(), val=character())
        for (val in unique(tree_tibble[!is.na(tree_tibble[[tr]]), ][[tr]])) {
            # Get all tips with this value
            tips <- tree_tibble$label[!is.na(tree_tibble[[tr]]) & tree_tibble[[tr]] == val]
            # Get the node of the most recent common ancestor
            node <- ifelse(length(tips) > 1, getMRCA(tree, tips), tree_tibble$node[tree_tibble$label == tips])
            # Get the descendants of this node
            desc <- getDescendants(tree, node)
            node_traits <- rbind(node_traits, data.frame(node=desc, val=rep(val, length(desc))))
        }
        # Add to tree tibble
        tree_tibble[[tr]] <- node_traits$val[match(tree_tibble$node, node_traits$node)]
    }
    metadata_with_nodes <- tree_tibble %>% select(-parent, -node, -branch.length)
    return(metadata_with_nodes)
}

bac_meta_new <- add_node_metadata(bac_tree_new, bac_meta, c("phylum", "class", "order", "family", "genus")) %>%
    # remove _[A-Z] suffix from phylum to species
    mutate(across(phylum:species, ~ str_remove(.x, "_[A-Z]$")))

ar_meta_new <- add_node_metadata(ar_tree_new, ar_meta, c("phylum", "class", "order", "family", "genus")) %>%
    # remove _[A-Z] suffix from phylum to species
    mutate(across(phylum:species, ~ str_remove(.x, "_[A-Z]$")))

#####################
#### SAVE OUTPUT ####
#####################

write.tree(bac_tree_new, file = file.path(subdir, "bac_tree.tree"))
write.tree(ar_tree_new, file = file.path(subdir, "ar_tree.tree"))

write.table(bac_meta_new, file = file.path(subdir, "bac_meta.tsv"), sep = "\t", row.names=FALSE, quote=FALSE)
write.table(ar_meta_new, file = file.path(subdir, "ar_meta.tsv"), sep = "\t", row.names=FALSE, quote=FALSE)
