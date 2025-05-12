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

# Contig metadata
contigs <- read.csv(file.path(indir, "contig_metadata.csv"))

# Habitat relations
habitat_relations <- read.csv(file.path(subdir, "habitat_relations.csv"))

########################
#### PROCESS TABLES ####
########################

## Contig metdata
bac_contigs <- bac_meta %>% select(label, bin, Sample) %>% left_join(select(contigs, c(bin, sample, damage_model_pmax, pvalue, qvalue, RMSE)), by = c("bin" = "bin", "Sample" = "sample")) %>%
                group_by(bin) %>% mutate(median_pmax = median(damage_model_pmax), median_pvalue = median(pvalue), median_qvalue = median(qvalue), median_RMSE = median(RMSE))

ar_contigs <- ar_meta %>% select(label, bin, Sample) %>% left_join(select(contigs, c(bin, sample, damage_model_pmax, pvalue, qvalue, RMSE)), by = c("bin" = "bin", "Sample" = "sample")) %>%
                group_by(bin) %>% mutate(median_pmax = median(damage_model_pmax), median_pvalue = median(pvalue), median_qvalue = median(qvalue), median_RMSE = median(RMSE))

## Habitat relations
# For each taxon, get the number of references to each habitat
habitats <- rbind(bac_meta, ar_meta) %>% select(label, classification, domain) %>%
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
            filter(habitat %in% c("oral", "animal", "soil", "marine", "rumen", "gut"))

bac_habitats <- habitats %>% filter(domain == "Bacteria")
write.csv(bac_habitats, file.path(subdir, "bac_habitats.csv"), row.names = FALSE)

ar_habitats <- habitats %>% filter(domain == "Archaea")
write.csv(ar_habitats, file.path(subdir, "ar_habitats.csv"), row.names = FALSE)

# Get specimen ages and damage for plotting
ages_damage <- bac_meta %>% rbind(ar_meta) %>% select(label, bin, domain, Sample, host_order, host_species, Year, host_species) %>% filter(!is.na(Year) & Year != "") %>%
  left_join(select(contigs, c(bin, sample, damage_model_pmax, pvalue, qvalue, RMSE)), by = c("bin" = "bin", "Sample" = "sample")) %>%
  # If the value in the year column is not numeric, remove the non numeric symbols and convert to number
  mutate(Year_most_recent = case_when(is.na(as.numeric(Year)) ~ as.numeric(str_remove_all(Year, "<|\\?| \\(received\\)|[0-9][0-9][0-9][0-9]-")), TRUE ~ as.numeric(Year))) %>%
  mutate(Approximated_year = case_when(is.na(as.numeric(Year)) ~ "YES", TRUE ~ "NO")) %>%
  filter(!is.na(Year_most_recent) & !is.na(damage_model_pmax)) %>%
  group_by(label) %>% mutate(median_pmax = median(damage_model_pmax), median_pvalue = median(pvalue), median_qvalue = median(qvalue), median_RMSE = median(RMSE))

####################
#### PLOT TREES ####
####################

#### Bacteria tree ####
# Colour by order and habitat
bac_p <- ggtree(bac_tree, layout = "circular", aes(color=phylum), size = 1.5) %<+% bac_meta +
  scale_colour_manual(values = phylum_palette, name = "MAG phylum", na.value = "black") +
  new_scale_color() +
  geom_tiplab(size=2, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=2, stroke=0.7, 
                aes(fill=host_order, color=habitat.general)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_colour_manual(values = habitat_palette, name = "Host habitat", na.value = "black") +
  scale_x_continuous(expand = c(0, 0)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(-6, -6, -6, 1), "cm"), # Remove margins
        legend.position=c(0.06, 0.5),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20)) +
  guides(fill = guide_legend(override.aes = list(size = 5)), 
  color = guide_legend(override.aes = list(size = 5)))

# Add info
bac_p <- bac_p +
  # Add boxplot with damage patterns
  geom_fruit(data=bac_contigs, geom=geom_boxplot, mapping = aes(y=label, x=log10(damage_model_pmax + 0.001), group=label),
             axis.params=list(axis = "x", text.size = 6, hjust = 1, vjust = 0., nbreak = 3),
             offset = 0.3, pwidth = 0.2, alpha = 0.3) +
  new_scale_color() +
  # Add point with collection year
  geom_fruit(data=ages_damage, geom=geom_point, mapping = aes(y=label, colour = Year_most_recent, shape = Approximated_year), size = 2,
             offset = -0.2) +
  scale_colour_viridis_c(option = "magma", name = "Collection year", na.value = "transparent") +
  scale_shape(name = "", labels = c("Yes" = "Approximated", "No" = "From records")) +
  new_scale_colour() +
  # Add bubbleplot with habitat occurences
  geom_fruit(data = bac_habitats, geom=geom_point, mapping = aes(y=label, x=habitat, colour=habitat, size = occurences),
             offset = -0.2) +
  scale_color_manual(values = c("oral" = "#AE1E3D", "animal" = "#BD6E20", "rumen" = "#A4B81F", "gut" = "#BD9F20", "soil" = "#56A71C", "marine" = "#156B73"),
                     name = "MAG reported habitat") +
  scale_size(name = "MAG reported occurences") +
  guides(colour = guide_legend(override.aes = list(size = 2.5)))

ggsave(bac_p, file=file.path(subdir, "bac_genome_tree.png"), width = 22, height = 20)

#### Archaea tree ####
# Colour by order and habitat
ar_p <- ggtree(ar_tree) %<+% ar_meta +
  # Colour by phylum
  geom_tree(aes(color=phylum)) +
  scale_colour_manual(values = phylum_palette, name = "MAG phylum", na.value = "black") +
  new_scale_color() +
  geom_tiplab(size=3, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order, color=habitat.general)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_colour_manual(values = habitat_palette, name = "Host habitat", na.value = "black") +
  scale_x_continuous(expand = c(0.04, 0.04)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(0, 0, 0, 3), "cm"),
        legend.position=c(-0.01, 0.5),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)), 
         color = guide_legend(override.aes = list(size = 2.5))) 

# Add info
ar_p <- ar_p +
  # Add boxplot with damage patterns
  geom_fruit(data=ar_contigs, geom=geom_boxplot, mapping = aes(y=label, x=log10(damage_model_pmax + 0.001), group=label),
             axis.params=list(axis = "x", text.size = 6, hjust = 1, vjust = 0., nbreak = 3),
             offset = 0.5, pwidth = 0.2, alpha = 0.3) +
  scale_color_gradient2(low = "#11B200", mid = "#FFBA00", high = "#A40073",
                        midpoint = -1.3, name = "median q-value", na.value = "black", breaks = c(0, 0.05, 1), trans = "log10") +
  new_scale_color() +
  # Add bubbleplot with habitat occurences
  geom_fruit(data = ar_habitats, geom=geom_point, mapping = aes(y=label, x=habitat, colour=habitat, size = occurences),
             offset = -0.1, pwidth = 0.1) +
  scale_color_manual(values = c("oral" = "#AE1E3D", "animal" = "#BD6E20", "rumen" = "#A4B81F", "gut" = "#BD9F20", "soil" = "#56A71C", "marine" = "#156B73"),
                     name = "MAG reported habitat") +
  scale_size(name = "MAG reported occurences") +
  guides(colour = guide_legend(override.aes = list(size = 2.5)))

ggsave(ar_p, file=file.path(subdir, "ar_genome_tree.png"), width = 12, height = 8)

###########################
#### PLOT AGE V DAMAGE ####
###########################

p <- filter(ages_damage, Approximated_year=="NO") %>%
     group_by(host_order) %>% filter(n_distinct(Year_most_recent) > 1) %>%
      ggplot(aes(x = Year_most_recent, group = Year_most_recent, y = damage_model_pmax, shape = Approximated_year, colour = host_order)) +
        geom_point(alpha = 0.1) +
        scale_colour_manual(values = order_palette, name = "Host order") +
        geom_smooth(method = "lm", span = 5, inherit.aes = FALSE, aes(x = Year_most_recent, y = damage_model_pmax), colour = "black", linetype = "dotted") +
        facet_grid(. ~ host_order) + theme(legend.position = "none")

ggsave(p, file=file.path(subdir, "age_v_damage.png"), width = 8, height = 6)
