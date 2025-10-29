#### FUNCTIONS TO USE WITH PHYLOGENETIC TREES ####
library(ape)
library(phytools)
library(tidyverse)

# Create a function that will assign the traits of the tips to all parent nodes until the MRCA
tips_to_nodes <- function(tree, trait_table, trait) {
  # Create a table to store the info
  node_table <- data.frame(matrix(nrow = 0, ncol = 2))
  colnames(node_table) <- c("node", sym(trait))
  # For each unique trait value
  for (value in unique(trait_table[, trait])) {
    # Get tips
    tips <- trait_table %>% filter(!!sym(trait) == value) %>% pull(tip)
    # If there is only one tip, the trait doesn't get assigned to any parent nodes
    if (length(tips)==1) {
      node_table <- node_table %>% rbind(data.frame(node = tip_index <- which(tree$tip.label == tips),
                                                    trait_value = value) %>% 
                                           rename(!!trait := trait_value))
    } else {
      # Get their mrca
      anc <- getMRCA(tree, tips)
      # Get descendant nodes of the mrca
      desc <- getDescendants(tree, anc)
      # Append the node numbers and corresponding trait values to node_table
      node_table <- node_table %>% rbind(data.frame(node = desc,
                                                    trait_value = rep(value, length(desc))) %>% 
                                           rename(!!trait := trait_value))
    }
  }
  return(node_table)
}
