##### CODIVERSIFICATION #####

#### Running MAG-host codiversification tests

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ape)
library(phytools)
library(tibble)
library(ggnewscale)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "mags")) 
taxdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Directory with taxonomic composition analysis
subdir <- normalizePath(file.path(outdir, "codiversification")) # subdirectory for the output of this script

# Create output directory
dir.create(subdir, showWarnings = FALSE)

source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))

setwd(subdir)

#######################
#####  LOAD INPUT #####
#######################

# Bin metadata
bac_meta <- read.table(file.path(outdir, "bac_meta_drep.tsv"), sep="\t", header=TRUE)
ar_meta <- read.table(file.path(outdir, "ar_meta_drep.tsv"), sep="\t", header=TRUE)

# MAG trees
bac_tree <- read.tree(file = file.path(outdir, "bac_tree_drep.tree"))
ar_tree <- read.tree(file = file.path(outdir, "ar_tree_drep.tree"))

# Fix node labels
bac_tree$node.label <- str_remove_all(bac_tree$node.label, "'")
ar_tree$node.label <- str_remove_all(ar_tree$node.label, "'")

# Host phylogeny
host_consensus <- read.tree(file.path(taxdir, "host_consensus.tre"))

# Habitat info
bac_habitats <- read.table(file.path(outdir, "bac_habitats.csv"), sep=",", header=TRUE)
ar_habitats <- read.table(file.path(outdir, "ar_habitats.csv"), sep=",", header=TRUE)

# MAGs per host species
mag_host <- read.table(file.path(outdir, "mag_mapping_stats", "hq_mag_presence_per_host.csv"), sep=",", header=TRUE)

#######################
#### TIDY UP DATA  ####
#######################

#### TREES ####
# Get consensus tree and fix tip labels
host_consensus$tip.label <- gsub("_", " ", host_consensus$tip.label)

#### MAG - HOST LINKS ####
mag_host <- mag_host %>% select(host_species, label, assembly_species) %>%
  # Also add Patescibacteria, which were not used for mapping as they didn't pass completeness filters
    rbind(
        filter(bac_meta, phylum == "Patescibacteria" & is.tip) %>% select(label, host_species) %>% mutate(assembly_species = TRUE)
    ) %>%
    # Change Sus domesticus to Sus scrofa to match tree
    mutate(host_species = case_when(host_species == "Sus domesticus" ~ "Sus scrofa",
                                    TRUE ~ host_species)) %>%
    group_by(host_species, label) %>% summarise(assembly_species = any(assembly_species)) %>% ungroup

hosts_per_mag <- mag_host %>% group_by(label) %>%
             summarise(n_hosts = n_distinct(host_species))

p <- ggplot(data=hosts_per_mag) + geom_histogram(aes(x=n_hosts), bins=50)

ggsave(plot=p, filename=file.path(subdir, "hosts_per_mag_histogram.png"), width=6, height=4)

links_all <- mag_host %>% select(host_species, label) %>%
    filter(!host_species %in% c("Extraction blank", "Library blank", "Environmental control"))

write.csv(links_all, file=file.path(subdir, "mag_host_links.csv"), row.names=FALSE, quote=FALSE)

##########################
#### DEFINE FUNCTIONS ####
##########################

# Select MAG clades to test
# Outputs lists, startimg with the most recent clade with at least min_tips number of tips
# Then go back until the maximum depth is reached
select_mag_clades <- function(tree, min_tips, max_depth) {
  cat(paste("Finding clades with at least", min_tips, "tips and a maximum depth of", max_depth, "\n"))
  # Start by identifying clades with at least 5 tips
  # but exclude parent nodes if the child node is already in the list
  initial_nodes <- c()
  cat("Finding initial nodes\n")
  for (node in Ntip(tree) + (tree$Nnode:1)) {
    descendants <- getDescendants(tree, node)
    if (length(descendants) >= min_tips & length(intersect(descendants, initial_nodes)) == 0) {
      initial_nodes <- append(node, initial_nodes)
    }
  }
  # Empty list to store nodes
  nodes_lists <- list()
  cat("Going back to maximum depth\n")
  # Traverse tree backwards from each initial node to the maximum depth, adding nodes to a lineage list
  for (node in initial_nodes) {
    initial_node <- node
    lineage <- c()
    depth <- node.depth(tree)[node]
    while (depth <= max_depth & depth < max(node.depth(tree))) {
      depth <- node.depth(tree)[node]
      # Add to list
      lineage <- append(lineage, node)
      # Get the parent node
      node <- getParent(tree, node)
      depth <- node.depth(tree)[node]
    }
    # Add to list using node labels
    node_label <- tree$node.label[initial_node - Ntip(tree)]
    nodes_lists[[node_label]] <- tree$node.label[lineage - Ntip(tree)]
  }
  return(nodes_lists)
}

# Get input for codiversification analysis: the subset of the MAG tree, ther corresponding host tree and the links
prepare_input <- function(node, mag_tree, host_mag_links, host_tree) {
  # Get subtree
  mag_subtree <- extract.clade(mag_tree, node)
  # Get bins and their host species based on keywords
  links <- host_mag_links %>% select(label, host_species) %>% filter(label %in% mag_subtree$tip.label)
  # Subset_host_tree
  hosts.to.drop <- host_tree$tip.label[!(host_tree$tip.label %in% links$host_species)]
  host_tree_filtered <- drop.tip(host_tree, hosts.to.drop)
  return(list(mag_subtree, host_tree_filtered, links))
}

# Get correlation statistics for subtree
dist_correlation <- function(mag_tree_sub, host_tree_sub, links, plot=FALSE) {
   # Get distance matrices
  host_dist <- cophenetic(host_tree_sub)[links$host_species, links$host_species] %>% as.dist
  mag_dist <- cophenetic(mag_tree_sub)[links$label, links$label] %>% as.dist
  # Run Mantel test
  r <- cor.test(host_dist, mag_dist, method="pearson", alternative = "greater")
  if (plot) {
    plot <- ggplot(data.frame(x=host_dist, y=mag_dist), aes(x=x, y=y)) +
            geom_point() + scale_x_continuous() + scale_y_continuous() +
            labs(x="Host distance", y="MAG distance") +
            labs(title = mag_tree_sub$node.label[1], subtitle = paste("r:", format(r$statistic, digits=3), "p:", format(r$p.value, digits=3)), size=5)
    ggsave(plot=plot, filename=file.path(intermediatedir, paste0(mag_tree_sub$node.label[1], "_cor.png")), width=4, height=4)
  }
  return(list(r=r$statistic, p.value=r$p.value))
}

# Permute links
permute_tree_tips <- function(tree) {
  tree_perm <- tree
  tree_perm$tip.label <- sample(tree$tip.label)
  return(tree_perm)
}

# First correlation calculation on all selected nodes
first_correlation_test <- function(mag_tree, host_mag_links, host_tree, nodes_to_test) {
  results = data.frame(node = character(), r = numeric(), p.value = numeric())
  # Get unique nodes to test
  unique_nodes <- unique(unlist(nodes_to_test))
  for (node in unique_nodes) {
    # Get input
    input <- prepare_input(node, mag_tree, host_mag_links, host_tree)
    mag_tree_sub <- input[[1]]
    host_tree_sub <- input[[2]]
    links <- input[[3]]
    if (length(host_tree_sub$tip.label) <= 2) {
      cat("Skipping", node, "with less or equal to 2 host tips\n")
      next
    }
    # Calculate the correlation between the host and MAG distance matrices and permute to get p-value
    res <- permutation_test(mag_tree_sub, host_tree_sub, links, nperm=500)
    r <- res$r
    p <- res$p.value
    results <- rbind(results, data.frame(node=node, r=r, p.value=p))
    cat(node, "- r:", r, " - p:", p, ": #tips:", length(mag_tree_sub$tip.label), "\n")    
  }
  return(results)
}

best_per_lineage <- function(nodes_to_test, first_results) {
  best_nodes <- c()
  n <- 1
  for (lineage_name in names(nodes_to_test)) {
    cat(paste0("Lineage ", n, ":", length(nodes_to_test), " - ", lineage_name, "\n"))
    lineage_nodes <- nodes_to_test[[lineage_name]]
    # Scan lineage and get the node with the highest correlation
    lineage_res <- filter(first_results, node %in% lineage_nodes)
    best_node <- lineage_res[which.max(lineage_res$r), "node"]
    best_nodes <- append(best_nodes, best_node)
    if (is.null(best_node)) {
      next
    }
    n <- n + 1
  }
  best_nodes <- unique(best_nodes)
  return(best_nodes)
}

permutation_test <- function(mag_tree_sub, host_tree_sub, links, nperm) {
    # Calculate the correlation between the host and MAG distance matrices
    r <- dist_correlation(mag_tree_sub, host_tree_sub, links, plot=TRUE)$r[[1]]
    # Permute host species and calculate correlation
    permuted_r <- replicate(nperm, dist_correlation(permute_tree_tips(mag_tree_sub), host_tree_sub, links)$r[[1]]) %>% 
                  unlist
    # Calculate p-value
    p_value <- sum(permuted_r >= r) / length(permuted_r)
    # Plot
    plot <- ggplot(data.frame(r=permuted_r), aes(x=r)) +
          geom_histogram(bins = 30) +
          geom_vline(xintercept=r, color="red") +
          labs(title = mag_tree_sub$node.label[1], subtitle=paste("r:", format(r, digits=3), "p:", format(p_value, digits=3)), size=5)
    ggsave(plot=plot, filename=file.path(intermediatedir, paste0(mag_tree_sub$node.label[1], "_perm.png")), width=4, height=4)
    
    return(list(r=r, p.value=p_value))
}

#################################
#### CODIVERSIFICATION TESTS ####
#################################

analysis <- function(mag_tree, host_mag_links, host_tree, nodes_to_test, intermediatedir) {
  # Create directory for intermediate results
  intermediatedir <<- intermediatedir
  dir.create(file.path(subdir, intermediatedir), showWarnings = FALSE)
  # Do first correlation test in all nodes
  cat("First correlation test\n")
  first_results <- first_correlation_test(mag_tree, host_mag_links, host_tree, nodes_to_test)
  write.csv(first_results, file = file.path(intermediatedir, "first_correlation_results.csv"), quote = FALSE, row.names = FALSE)
  # Get best nodes per lineage
  cat("Deciding best nodes per lineage\n")
  best_nodes <- best_per_lineage(nodes_to_test, first_results)
  # Then test only best nodes
  cat("Will test", length(best_nodes), "nodes\n")
  n <- 1
  total <- length(best_nodes)
  test_results <- data.frame()
  for (node in best_nodes) {
    cat(paste0("Node ", n, ":", total, " - ", node, "\n"))
    # Prepare analysis input
    input <- prepare_input(node, mag_tree, host_mag_links, host_tree)
    mag_tree_sub <- input[[1]]
    host_tree_sub <- input[[2]]
    links <- input[[3]]
    # Get some info about the clade
    info <- data.frame(node = node,
              mags_tips = length(mag_tree_sub$tip.label),
              host_tips = length(host_tree_sub$tip.label))
    # Run permutation test
    perm = permutation_test(mag_tree_sub, host_tree_sub, links, nperm=500)
    r <- perm$r
    p_value <- perm$p.value
    # Combine with node info
    res <- data.frame(r=r, p.value=p_value)
    res <- cbind(info, res)
    # Append to table
    test_results <- rbind(test_results, data.frame(res))
    n <- n + 1
  }
  # Account for multiple testing: Add a pseudocount because we have some p-values of 0
  test_results$p.value <- test_results$p.value + 1e-5
  test_results$p.adjust <- p.adjust(test_results$p.value, method="holm")
  return(test_results)
}

### For Bacteria ####
set.seed(123)
nodes_to_test <- select_mag_clades(bac_tree, min_tips=3, max_depth=50)

bac_cod_results <- analysis(bac_tree, links_all, host_consensus, nodes_to_test, "bac_cod_intermediate")

write.csv(bac_cod_results, file=file.path(subdir, "bac_cod_results.csv"), row.names=FALSE, quote=FALSE)

### For Archaea ####
set.seed(123)
nodes_to_test <- select_mag_clades(ar_tree, min_tips=3, max_depth=10)

ar_cod_results <- analysis(ar_tree, links_all, host_consensus, nodes_to_test, "ar_cod_intermediate")

write.csv(ar_cod_results, file=file.path(subdir, "ar_cod_results.csv"), row.names=FALSE, quote=FALSE)

cod_results_combined <- rbind(
  bac_cod_results %>% mutate(domain="Bacteria"),
  ar_cod_results %>% mutate(domain="Archaea")
)

# Add phylum
cod_results_combined <- cod_results_combined %>% left_join(select(rbind(bac_meta, ar_meta), label, phylum), by=c("node"="label")) %>%
    mutate(phylum = case_when(node == "Synergistales" ~ "Synergistota",
                              node == "Patescibacteria" ~ "Patescibacteria",
                              node == "o__Erysipelotrichales" ~ "Bacillota", TRUE ~ phylum))

write.csv(cod_results_combined, file=file.path(subdir, "cod_results_combined.csv"), row.names=FALSE, quote=FALSE)

######################
#### PLOT RESULTS ####
######################

#### Plot cophylogenies ####
cophyloplot <- function(host_tree, mag_tree, links, host_colours, host_linetypes) {
  links_rename <- links %>% rename("phy1"="label", "phy2"="host_species")
  # Add colours
  links_rename$colour <- host_colours[match(links_rename$phy2, names(host_colours))]
  # Add line types
  links_rename$linetype <- host_linetypes[match(links_rename$phy2, names(host_linetypes))]
  links_rename <- links_rename %>% as.matrix
  # Plot comparison
  coph <- cophylo(tr1=mag_tree, tr2=host_tree, assoc=links_rename)
  # Plot
  return(coph)
}

create_plots <- function(mag_tree, host_mag_links, host_tree, cod_results, outdir) {
  n=1
  for (i in 1:nrow(cod_results)) {
    node <- cod_results$node[i]
    cat(paste0(n, ":", nrow(cod_results), " - node ", node, "\n"))
    # Get input data
    input <- prepare_input(node, mag_tree, host_mag_links, host_tree)
    mag_tree_sub <- input[[1]]
    host_tree_sub <- input[[2]]
    links <- input[[3]]
    # Use dotted lines for non-assembly species
    linetypes <- links %>% left_join(mag_host) %>%
        mutate(line_type = case_when(assembly_species == FALSE ~ "dotted", TRUE ~ "solid"))
    linetypes <- setNames(linetypes$line_type, linetypes$host_species)
    # Plot cophyloplo
    coph <- cophyloplot(host_tree=host_tree_sub, mag_tree=mag_tree_sub, links=links, host_colours=species_palette, host_linetypes = linetypes)
    # Save plot
    png(file.path(subdir, outdir, paste(node, "cophyloplot.png", sep="_")), width = 600, height = 400)
    par(mar=c(6, 8, 4, 2) + 0.1)
    plot(coph, link.type="curved",link.lwd=4, link.lty=coph$assoc[,"linetype"], link.col=coph$assoc[,"colour"])
    # Add title
    text(x = 0.5, y = 0.95, labels = paste("Node:", node), pos = 2, cex = 1.5, col = "black")
    text(x = 0.5, y = 0.90, labels = paste("p-value:", format(cod_results$p.value[i], digits = 3)), pos = 2, cex = 1.5, col = "black")
    text(x = 0.5, y = 0.85, labels = paste("q-value:", format(cod_results$p.adjust[i], digits = 3)), pos = 2, cex = 1.5, col = "black")
    dev.off()
    n <- n + 1
  }
}

# Save plots
dir.create(file.path(subdir, "bac_cophylo_plots"), showWarnings = FALSE)

create_plots(bac_tree, links_all, host_consensus, bac_cod_results, "bac_cophylo_plots")

# Save plots
dir.create(file.path(subdir, "ar_cophylo_plots"), showWarnings = FALSE)

create_plots(ar_tree, ar_meta, host_consensus, ar_cod_results, "ar_cophylo_plots")

#### Plot entire tree ####

### Bacteria tree

# Add to tree metadata
bac_meta_plot <- bac_meta %>% left_join(bac_cod_results, by=c("label"="node")) %>%
  # For non significant results, turn r to NA
  mutate(r = case_when(p.adjust <= 0.05 ~ r, TRUE ~ NA)) %>%
  # Get -log10(p.adjust) for colouring
  mutate(neg_log10_p = case_when(!is.na(p.adjust) ~ -log10(p.adjust+10^-5), TRUE ~ NA)) %>%
  select(label, r, neg_log10_p, host_order, habitat.general)

codiv_nodes <- 
    bac_cod_results %>% filter(p.adjust <= 0.05 & r > 0) %>%
  pull(node)

# Get r values of descending nodes and tips of codiversifying clades
descending_r <- data.frame()

for (n in codiv_nodes) {
  subclade <- extract.clade(bac_tree, n)
  all_desc <- c(subclade$tip.label, subclade$node.label)
  node_r <- bac_cod_results %>% filter(node == n) %>% pull(r)
  df <- data.frame(label = all_desc, desc_r = node_r)
  descending_r <- rbind(descending_r, df) %>% filter(label != n)
}

descending_r <- descending_r %>%
  # If a tip is in multiple codiversifying clades, take the maximum r value
  group_by(label) %>%
  summarise(desc_r = max(desc_r, na.rm=TRUE))

bac_meta_plot <- left_join(bac_meta_plot, descending_r, by="label")

# Colour by order and habitat
bac_p <- ggtree(bac_tree, layout = "circular", size = 1.5, aes(colour = desc_r)) %<+% bac_meta_plot +
  scale_colour_gradient2(low = "yellow", mid = "orange", high = "red", midpoint = median(bac_meta_plot$r, na.rm = TRUE), na.value = "grey", name = "r coefficient") +
  geom_nodepoint(aes(fill=r, size=neg_log10_p), shape = 21, colour = "black") +
  scale_fill_gradient2(low = "yellow", mid = "orange", high = "red", midpoint = median(bac_meta_plot$r, na.rm = TRUE), na.value = "white", name = "r coefficient") +
  scale_size_continuous(name = "-log10(p adjusted)", range = c(1,5)) +
  new_scale_colour() +
  geom_tiplab(size=2.5, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_colour() +
  new_scale_fill() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order, color=habitat.general)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_colour_manual(values = habitat_palette, name = "Host habitat", na.value = "black") +
  scale_x_continuous(expand = c(0, 0)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(-6, -6, 0, -6), "cm"), # Remove margins
    legend.position="bottom",
    legend.direction="vertical",
    legend.text = element_text(size=20),
    legend.title = element_text(size=20)) +
  guides(fill = guide_legend(override.aes = list(size = 5)), 
  color = guide_legend(override.aes = list(size = 5)))

# Identify the node for each phylum
phylum_nodes <- bac_meta %>% 
  select(phylum, label) %>%
  left_join(as_tibble(bac_tree)[c("label", "node")]) %>%
  group_by(phylum) %>%
  filter(node == max(node, na.rm=TRUE)) %>%
  mutate(colour = phylum_palette[phylum]) %>%
  filter(!is.na(colour)) %>%
  # Need to do this manually because ggtree node numbers to not line up with the original tree node numbers
  mutate(node = case_when(
    phylum == "Actinomycetota" ~ 678,
    phylum == "Bacillota" ~ 569,
    phylum == "Pseudomonadota" ~ 938,
    phylum == "Desulfobacterota" ~ 925,
    phylum == "Bacteroidota" ~ 1056,
    phylum == "Synergistota" ~ 900
  ))

bac_p_node <- bac_p

for (i in 1:nrow(phylum_nodes)) {
  node <- phylum_nodes$node[i]
  label <- phylum_nodes$phylum[i]
  colour <- phylum_nodes$colour[i]
  bac_p_node <- bac_p_node +
      geom_cladelabel(node = node, label = label, colour = colour,
                      offset = 0.9, barsize = 3, fontsize = 0)
}

ggsave(bac_p_node, file=file.path(subdir, "bac_codiv_tree.png"), width = 20, height = 22)

#######################################
#### COMPARE CODIVERSIFYING VS NOT ####
#######################################

# Get lists of MAGs in codiversifying clades vs not codiversifying clades
# Find all descending tips
get_desctips <- function(tree, node) {
  sub_tree <- extract.clade(tree, node)
  return(sub_tree$tip.label)
}

codiv_tips <- lapply(codiv_nodes, get_desctips, tree=bac_tree) %>% unlist %>% unique

# Also add all remaining tips, indicating they are not codiversifying
codiv_tips_df <- data.frame(label = codiv_tips, 
                            codiversifying = TRUE)

non_codiv_tips <- bac_tree$tip.label[!(bac_tree$tip.label %in% codiv_tips)]

codiv_tips_df <- rbind(codiv_tips_df, data.frame(label = non_codiv_tips, 
                                               codiversifying = FALSE))

# Add habitat info
codiv_habitats <- bac_habitats %>% filter(label %in% codiv_tips_df$label) %>%
  filter(habitat %in% c("gut", "oral", "rumen", "soil") & occurences > 0.3) %>%
  select(label, habitat, occurences)

codiv_tips_df <- codiv_tips_df %>% left_join(codiv_habitats, by="label") %>%
    # Where not habitat info, fill in with "no info"
    mutate(habitat = case_when(is.na(habitat) ~ "no info", TRUE ~ habitat)) %>%
    mutate(habitat = factor(habitat, levels=c("oral", "rumen", "gut", "soil", "no info"))) %>%
    # Indicate MAGs counted more than once
    group_by(label) %>%
    mutate(multiple_habitats = case_when(n_distinct(habitat) > 1 ~ TRUE, TRUE ~ FALSE)) %>%
    mutate(codiversifying = factor(codiversifying, levels=c(TRUE, FALSE)))

# Add taxonomic info
codiv_tax <- bac_meta %>% filter(label %in% codiv_tips_df$label) %>%
  select(label, phylum, class, order, family, genus)

codiv_tips_df <- codiv_tips_df %>% left_join(codiv_tax, by="label")

# Save table
write.csv(codiv_tips_df, file=file.path(subdir, "bac_codiv_tips_info.csv"), row.names=FALSE, quote=FALSE)

#### Plot how many codiversifying MAGs are associated with each habitat ####

p <- ggplot(data = codiv_tips_df, aes(x = codiversifying)) +
  geom_bar(aes(fill = multiple_habitats)) +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "darkblue"), labels = c("FALSE" = "single habitat", "TRUE" = "multiple habitats"), name = "") +
  facet_grid(cols = vars(habitat), scales = "free_y") +
  theme(strip.text = element_text(size = 10), legend.position = "top", legend.title = element_blank(),
     axis.title.y = element_blank())

ggsave(plot=p, filename=file.path(subdir, "codiv_per_habitat.png"), width=3.5, height=2.5)

#### Plot how many codiversifying MAGs are associated with each phylum ####
p <- ggplot(data = unique(select(codiv_tips_df, c(label, phylum, codiversifying))), aes(x = codiversifying)) +
  geom_bar(aes(fill = phylum)) +
  scale_fill_manual(values = phylum_palette, na.value = "grey", name = "Phylum") +
  facet_grid(cols = vars(phylum), scale = "free_y") +
  theme(strip.text = element_text(angle = 90, size = 10), legend.position = "none",
        axis.title.y = element_blank())

ggsave(plot=p, filename=file.path(subdir, "codiv_per_phylum.png"), width=3.5, height=3)
