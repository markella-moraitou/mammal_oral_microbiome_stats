##### Plot MAG tree #####

#### Plot MAG trees for Bacteria and Archaea and annotate using relevant metadata

#### LOAD PACKAGES ####
library(ape)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(RColorBrewer)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "mags")) # subdirectory for the output of this script

source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))

#######################
#####  LOAD INPUT #####
#######################

# Bin metadata
bac_meta <- read.table(file.path(subdir, "bac_meta.tsv"), sep="\t", header=TRUE)
ar_meta <- read.table(file.path(subdir, "ar_meta.tsv"), sep="\t", header=TRUE)

# MAG trees
bac_tree <- read.tree(file = file.path(subdir, "bac_tree.tree"))
ar_tree <- read.tree(file = file.path(subdir, "ar_tree.tree"))

# Host phylogeny
host_trees <- read.nexus(file.path(indir, "mammal_vertlife.nex"))

# Habitat relations
habitat_relations <- read.csv(file.path(subdir, "habitat_relations.csv"))

# List of dereplicated bins
drep_bins <- read.table(file.path(indir, "dereplicated_bins_list.txt"), header=FALSE, sep="\t") %>% pull(V1)

###########################
#### KEEP ONLY HQ MAGS ####
###########################

## Keep only HQ MAGs and Patescibacteria (which are likely to have they completeness underestimated)
hq_bacs <- bac_meta %>% filter(Completeness >= 90 & Contamination <= 5 | (phylum == "Patescibacteria" & Contamination <= 5)) %>%
      filter(bin %in% drep_bins) %>% pull(label) 

hq_ars <- ar_meta %>% filter(Completeness >= 90 & Contamination <= 5) %>%
      filter(bin %in% drep_bins) %>% pull(label)

## Subset trees to only include HQ MAGs
bac_tree <- drop.tip(bac_tree, setdiff(bac_tree$tip.label, hq_bacs))
ar_tree <- drop.tip(ar_tree, setdiff(ar_tree$tip.label, hq_ars))

# Fix node labels
bac_tree$node.label <- str_remove_all(bac_tree$node.label, "'")
ar_tree$node.label <- str_remove_all(ar_tree$node.label, "'")

## Subset metadata to only include HQ MAGs
bac_meta <- bac_meta %>% filter(label %in% c(bac_tree$tip.label, bac_tree$node.label))
ar_meta <- ar_meta %>% filter(label %in% c(ar_tree$tip.label, ar_tree$node.label))

# Save trees and tables
write.tree(bac_tree, file = file.path(subdir, "bac_tree_drep.tree"))
write.tree(ar_tree, file = file.path(subdir, "ar_tree_drep.tree"))

write.table(bac_meta, file = file.path(subdir, "bac_meta_drep.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(ar_meta, file = file.path(subdir, "ar_meta_drep.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

########################
#### PROCESS TABLES ####
########################

## Habitat relations
# For each taxon, get the number of references to each habitat
habitats <- rbind(bac_meta, ar_meta) %>% filter(!is.na(bin)) %>% select(label, classification, domain) %>%
            # Get taxon to match with habitat_relations
            mutate(taxon = str_remove(classification, ".*__")) %>% left_join(habitat_relations, relationship = "many-to-many") %>%
            # Fill NAs with 0
            mutate(occurences = replace(occurences, is.na(occurences), 0)) %>%
            # Get percentage of occurences
            group_by(label) %>% mutate(occurences = occurences/sum(occurences)) %>%
            # Group OBTs into larger habitats
            mutate(habitat = case_when(OBT %in% c("dental plaque", "mouth") ~ "oral",
                                       OBT %in% c("marine water", "deep sea") ~ "marine",
                                       OBT %in% c("mammalian", "wild animal", "mammalian livestock") ~ "animal",
                                       TRUE ~ OBT)) %>%
            # Make sure that each taxon has an occurence value for each habitat, even if 0, instead of row missing
            tidyr::complete(taxon, habitat, fill = list(occurences = 0)) %>%
            group_by(label, habitat, domain) %>%
            summarise(occurences = sum(occurences)) %>%
            mutate(habitat = gsub(" ", "_", habitat)) %>%
            # Only select some of the most informative terms
            mutate(habitat = case_when(habitat == "laboratory_equipment" ~ "lab", TRUE ~ habitat)) %>%
            filter(habitat %in% c("oral", "animal", "soil", "marine", "rumen", "gut", "skin")) %>%
            mutate(habitat = factor(habitat, levels = c("oral", "animal", "rumen", "gut", "marine", "soil", "skin")))

bac_habitats <- habitats %>% filter(domain == "Bacteria")
write.csv(bac_habitats, file.path(subdir, "bac_habitats.csv"), row.names = FALSE)

ar_habitats <- habitats %>% filter(domain == "Archaea")
write.csv(ar_habitats, file.path(subdir, "ar_habitats.csv"), row.names = FALSE)

# Get specimen ages and damage for plotting
ages_damage <- bac_meta %>% rbind(ar_meta) %>% select(label, bin, domain, Sample, host_order, host_species, Year, host_species, min_damage_model_pmax, q1_damage_model_pmax, mean_damage_model_pmax, median_damage_model_pmax, q3_damage_model_pmax, max_damage_model_pmax) %>% filter(!is.na(Year) & Year != "") %>%
  filter(!is.na(bin)) %>%
  # If the value in the year column is not numeric, remove the non numeric symbols and convert to number
  mutate(Year_most_recent = case_when(is.na(as.numeric(Year)) ~ as.numeric(str_remove_all(Year, "<|\\?| \\(received\\)|[0-9][0-9][0-9][0-9]-")), TRUE ~ as.numeric(Year))) %>%
  mutate(Approximated_year = case_when(is.na(as.numeric(Year)) ~ "YES", TRUE ~ "NO")) %>%
  filter(!is.na(Year_most_recent))

####################
#### PLOT TREES ####
####################

#### Bacteria tree ####
# Colour by order and habitat
bac_p1 <- ggtree(bac_tree, layout = "circular", aes(color=phylum), size = 1) %<+%
    select(bac_meta, c(label, phylum, host_order, habitat.general)) +
  scale_colour_manual(values = phylum_palette, name = "MAG phylum", na.value = "black") +
  new_scale_color() +
  geom_tiplab(size=3, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order, color=habitat.general)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_colour_manual(values = habitat_palette, name = "Host habitat", na.value = "black") +
  scale_x_continuous(expand = c(0, 0)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(-6, -6, 0, -6), "cm"), # Remove margins
        legend.position="bottom",
        legend.direction = "vertical",
        legend.text = element_text(size=20),
        legend.title = element_text(size=20)) +
  guides(fill = guide_legend(override.aes = list(size = 5)), 
        color = guide_legend(override.aes = list(size = 5)))

bac_py <- filter(select(bac_meta, c(label, bin, median_damage_model_pmax)), !is.na(bin))

# Add info
bac_p <- bac_p1 +
  new_scale_colour() +
  # Add bubbleplot with habitat occurences
  geom_fruit(data = bac_habitats, geom=geom_point, mapping = aes(y=label, x=habitat, colour=habitat, size = occurences),
             offset = 0.06, pwidth = 0.08) +
  scale_color_manual(values = c("oral" = "#AE1E3D", "animal" = "#BD6E20", "rumen" = "#A4B81F", "gut" = "#BD9F20", "soil" = "#56A71C", "marine" = "#156B73", "skin" = "#ad703a"),
                     name = "MAG reported habitat") +
  scale_size(name = "MAG reported occurences") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  # Add barplot with damage patterns
  geom_fruit(data=bac_py, geom=geom_bar, mapping = aes(y=label, x = -log10(median_damage_model_pmax + 0.001)),
             stat = "identity", axis.params=list(axis = "x", text.size = 6, hjust = 1, vjust = 0., nbreak = 3),
             offset = 0.01, pwidth = 0.15, alpha = 0.3) +
  new_scale_color() +
  # Add point with collection year
  geom_fruit(data=ages_damage, geom=geom_point, mapping = aes(y=label, colour = Year_most_recent, shape = Approximated_year), size = 3,
             offset = 0.0) +
  scale_colour_viridis_c(option = "magma", name = "Year of\ncollection", na.value = "transparent") +
  scale_shape(name = "Year\napproximated")

ggsave(bac_p, file=file.path(subdir, "bac_genome_tree.png"), width = 20, height = 25)

#### Archaea tree ####
# Colour by order and habitat
ar_p1 <- ggtree(ar_tree)%<+%
    select(ar_meta, c(label, phylum, host_order, habitat.general)) +
  # Colour by phylum
  geom_tree(aes(color=phylum)) +
  scale_colour_manual(values = c("red", "blue", "green"), name = "MAG phylum", na.value = "black") +
  new_scale_color() +
  geom_tiplab(size=3, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_x_continuous(expand = c(0.04, 0.04)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(0, 5, 0, 1), "cm"),
        legend.position="bottom",
        legend.direction = "vertical",
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)), 
         color = guide_legend(override.aes = list(size = 2.5))) +
  xlim(0, 2.5)

ar_py <- filter(select(ar_meta, c(label, bin, median_damage_model_pmax)), !is.na(bin))

# Add info
ar_p <- ar_p1 +
  new_scale_colour() +
  # Add bubbleplot with habitat occurences
  geom_fruit(data = ar_habitats, geom=geom_point, mapping = aes(y=label, x=habitat, colour=habitat, size = occurences),
             offset = 0.18, pwidth = 0.1) +
  scale_color_manual(values = c("oral" = "#AE1E3D", "animal" = "#BD6E20", "rumen" = "#A4B81F", "gut" = "#BD9F20", "soil" = "#56A71C", "marine" = "#156B73", "skin" = "#ad703a"),
                     name = "MAG reported habitat") +
  scale_size(name = "MAG reported occurences") +
  # Add boxplot with damage patterns
  geom_fruit(data=ar_py, geom=geom_bar, mapping = aes(y=label, x = -log10(median_damage_model_pmax + 0.001)),
             stat = "identity", axis.params=list(axis = "x", text.size = 6, hjust = 1, vjust = 0., nbreak = 3),
             offset = 0.02, pwidth = 0.2, alpha = 0.3) +
  new_scale_color() +
  # Add point with collection year
  geom_fruit(data=ages_damage, geom=geom_point, mapping = aes(y=label, colour = Year_most_recent, shape = Approximated_year), size = 3,
             offset = 0) +
  scale_colour_viridis_c(option = "magma", name = "Year of\ncollection", na.value = "transparent") +
  scale_shape(name = "Year\napproximated")

ggsave(ar_p, file=file.path(subdir, "ar_genome_tree.png"), width = 13, height = 6)

###########################
#### PLOT AGE V DAMAGE ####
###########################

p <- filter(ages_damage, Approximated_year=="NO") %>%
     group_by(host_order) %>% filter(n_distinct(Year_most_recent) > 1) %>%
      ggplot(aes(x = Year_most_recent, group = Year_most_recent, y = median_damage_model_pmax, shape = Approximated_year, colour = host_order)) +
        geom_point(alpha = 0.5) +
        scale_colour_manual(values = order_palette, name = "Host order") +
        geom_smooth(method = "glm", inherit.aes = FALSE, aes(x = Year_most_recent, y = median_damage_model_pmax), colour = "black", linetype = "dotted") +
        facet_grid(. ~ host_order) + theme(legend.position = "none")

ggsave(p, file=file.path(subdir, "age_v_damage.png"), width = 8, height = 6)
