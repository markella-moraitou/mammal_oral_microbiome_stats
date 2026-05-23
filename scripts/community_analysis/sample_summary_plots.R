##### PLOT SAMPLE METADATA #####

#### Summarise information about the samples and diet ####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggnewscale)
library(ape)
library(ggtree)
library(phytools)
library(rphylopic)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "sample_summary_plot")) # subdirectory for the output of this script
phydir <- normalizePath(file.path("..", "..", "output", "community_analysis", "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

phy_sp_f <- readRDS(file.path(phydir, "phy_sp_f.RDS"))
metadata <- data.frame(phy_sp_f@sam_data)

# Host phylogeny
host_trees <- read.nexus(file.path(indir, "mammal_vertlife.nex"))

# Phylopic data for plotting
phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"))

########################
#### SUMMARISE DIET ####
########################

meta <- metadata %>% filter(!is.neg)

#### Lintulaakso diet PCA ####
lint_diet_data <- unique(meta[,c("Common.name", "Species", "Order_grouped", "cp", "ee", "cf", "ash", "nfe", "diet.general")]) %>%
              filter(!grepl("blank", Species) & !grepl("control", Species))
              
rownames(lint_diet_data) <- NULL
lint_diet_data <- lint_diet_data %>% column_to_rownames("Species")

diet_var <- lint_diet_data[,c("cp", "ee", "cf", "ash", "nfe")]
ord <- prcomp(diet_var)

# Get some info for plotting
var_explained <- round(ord$sdev^2 * 100 / sum(ord$sdev^2), 1) # Variance explained

loadings_matrix <- data.frame(Variables = rownames(ord$rotation[,c(1,2)]), ord$rotation[,c("PC1", "PC2")]) %>%
    arrange(desc(PC1^2 + PC2^2)) %>% head(3)

pics <- phylopics$uid[match(rownames(ord$x), phylopics$Species)]

labels <- ord$x %>% as.data.frame %>% mutate(label = case_when(PC1 > 10 ~ rownames(.),
                                                               PC2 == max(PC2) | PC2 == min(PC2) ~ rownames(.))) %>%
  pull(label)

# Plot
pca <- ggplot(aes(x = PC1, y = PC2, colour=lint_diet_data$diet.general), data = data.frame(ord$x)) +
  geom_phylopic(aes(uuid = pics), width = 9, position = position_jitter(width = 2, height = 2), alpha = 0.8) +
  scale_colour_manual(values = diet_palette, name = "") +
  xlab(paste("PC1 -", var_explained[1], "%")) +
  ylab(paste("PC2 -", var_explained[2], "%")) +
  geom_segment(data = loadings_matrix, aes(x = 0, y = 0, xend = (PC1*18),
                                       yend = (PC2*18)), arrow = arrow(length = unit(0.5, "picas")),
               color = "black") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  annotate("text", x = (loadings_matrix$PC1*20), y = (loadings_matrix$PC2*20),
           label = loadings_matrix$Variables)

ggsave(filename  =  file.path(subdir, "dietary_PCA_lintulaakso.png"), pca, width  =  4, height = 4)

#### Elton traits diet PCA ####
elton_diet_data <- 
  unique(meta[,c("Common.name", "Species", "Order_grouped",
              "Animal", "Fruit", "Seed", "Nect", "PlantO", "diet.general")]) %>%
              filter(!grepl("blank", Species) & !grepl("control", Species))
              
rownames(elton_diet_data) <- NULL
elton_diet_data <- elton_diet_data %>% column_to_rownames("Species")

diet_var <- elton_diet_data[,c("Animal", "Fruit", "Nect", "Seed", "PlantO")]

# Remove variables that are all zero
diet_var <- diet_var[, colSums(diet_var) > 0]

ord <- prcomp(diet_var)

# Get some info for plotting
var_explained <- round(ord$sdev^2 * 100 / sum(ord$sdev^2), 1) # Variance explained

loadings_matrix <- data.frame(Variables = rownames(ord$rotation[,c(1,2)]), ord$rotation[,c("PC1", "PC2")]) %>%
    # Keep only the 4 variables that explain most variance
    arrange(desc(PC1^2 + PC2^2)) %>% head(4)

pics <- phylopics$uid[match(rownames(ord$x), phylopics$Species)]

# Plot
pca <- ggplot(aes(x = PC1, y = PC2, colour=elton_diet_data$diet.general), data = data.frame(ord$x)) +
  geom_phylopic(aes(uuid = pics), width = 9, position = position_jitter(width = 2, height = 2), alpha = 0.8) +
  scale_colour_manual(values = diet_palette, name = "") +
  xlab(paste("PC1 -", var_explained[1], "%")) +
  ylab(paste("PC2 -", var_explained[2], "%")) +
  geom_segment(data = loadings_matrix, aes(x = 0, y = 0, xend = (PC1*50),
                                       yend = (PC2*50)), arrow = arrow(length = unit(0.5, "picas")),
               color = "black") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2)) +
  annotate("text", x = (loadings_matrix$PC1*52), y = (loadings_matrix$PC2*52),
           label = loadings_matrix$Variables, hjust = 0)

ggsave(filename  =  file.path(subdir, "dietary_PCA_elton.png"), pca, width  =  4, height = 4)

#### Lintulaakso diet summary ####

lint_diet_long <- lint_diet_data %>% rownames_to_column("Species") %>%
  mutate(Common.name = meta$Common.name[match(Species, meta$Species)]) %>%
  pivot_longer(c(cp, ee, cf, ash, nfe), values_to = "proportion", names_to = "nutrient") %>%
  mutate(nutrient = factor(nutrient, levels = c("ee", "cp", "nfe", "cf", "ash"))) %>%
  mutate(Order_grouped = recode(Order_grouped, "Perissodactyla" = "Peris.", "Rodentia" = "Ro.", "Carnivora" = "Carn.")) %>%
  arrange(nutrient)

nutrient_palette <- c(ash = "grey",
                      cf = "#3BA01B",
                      nfe = "#3A459C",
                      cp = "#B41F34",
                      ee = "#BD9120")

diet_barplot <- ggplot(lint_diet_long, aes(x = Common.name, y = proportion, fill = nutrient, group = diet.general)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
       axis.title.y = element_blank(),
       axis.title.x = element_blank(),
       legend.position = "top",
       legend.text = element_text(size = 8)) +
  facet_grid(cols = vars(Order_grouped), scales = "free", space = "free") +
  scale_fill_manual(values = nutrient_palette,
                    labels = c("ash" = "ash\n(inorganic)",
                               "cf" = "cf\n(crude fiber)",
                               "nfe" = "nfe\n(carbohydrates)",
                               "cp" = "cp\n(crude protein)",
                               "ee" = "ee\n(fat)"),
                    name = "")

ggsave(filename  =  file.path(subdir, "diet_barplot_lintulaakso.png"), diet_barplot, width  =  6, height = 5)

#### Elton diet summary ####

elt_diet_long <- elton_diet_data %>% rownames_to_column("Species") %>%
  mutate(Common.name = meta$Common.name[match(Species, meta$Species)]) %>%
  pivot_longer(c(Animal, Fruit, Seed, PlantO), values_to = "proportion", names_to = "dietary_element") %>%
  mutate(dietary_element = factor(dietary_element, levels = c("Animal", "Fruit", "Seed", "PlantO"))) %>%
  mutate(Order_grouped = recode(Order_grouped, "Perissodactyla" = "Peris.", "Rodentia" = "Ro.", "Carnivora" = "Carn.")) %>%
  arrange(dietary_element)

item_palette <- c(Seed = "#FFD457",
                      PlantO = "#85C43C",
                      Fruit = "#4D71A7",
                      Animal = "#D95244")

diet_barplot <- ggplot(elt_diet_long, aes(x = Common.name, y = proportion, fill = dietary_element, group = diet.general)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
       axis.title.y = element_blank(),
       axis.title.x = element_blank(),
       legend.position = "top",
       legend.text = element_text(size = 8)) +
  facet_grid(cols = vars(Order_grouped), scales = "free", space = "free") +
  scale_fill_manual(values = item_palette,
                    labels = c("Seed" = "Seed",
                               "Fruit" = "Fruit",
                               "PlantO" = "Other plant material",
                               "Animal" = "Animal material"),
                    name = "")

ggsave(filename  =  file.path(subdir, "diet_barplot_elton.png"), diet_barplot, width  =  6, height = 5)

#### Combine and save diet data ####

diet_data <- full_join(rownames_to_column(elton_diet_data, "Species"),
                      rownames_to_column(lint_diet_data, "Species")) %>%
                      select(Species, Common.name, Animal, Fruit, Seed, PlantO,
                            cp, ee, cf, ash, nfe, diet.general)

write.csv(diet_data, file = file.path(subdir, "diet_data.csv"), quote  =  FALSE, row.names  =  FALSE)

#### Order and diet ####
summ_dat <- meta %>% group_by(Order, diet.general) %>% 
  summarise(n_samples = n_distinct(new_name), n_species = n_distinct(Species))

write.csv(summ_dat, file = file.path(subdir, "data_summary.csv"), quote  =  FALSE, row.names  =  FALSE)

sums <- summ_dat %>% ungroup %>% summarise(total_species = sum(n_species), total_samples = sum(n_samples))

p_summary <-
  ggplot(aes(x = diet.general, y = Order, size = n_species, fill = n_samples, colour = n_samples), data = summ_dat) +
  geom_point(shape = 21) +
  scale_fill_continuous(low = "#CDB139", high = "#C73842", name = "Sample number", trans = "log", breaks = c(5, 25, 125)) +
  scale_colour_continuous(low = "#CDB139", high = "#C73842", name = "Sample number", trans = "log", breaks = c(5, 25, 125)) +
  scale_size(range  =  c(5, 18), breaks = c(min(summ_dat$n_species), mean(unique(summ_dat$n_species)), max(summ_dat$n_species)),
             name = "Species number") +
  new_scale(new_aes = "size") +
  geom_text(aes(label = paste0(n_species, "(", n_samples, ")"), x = diet.general, y = Order, size = n_species), fontface = "bold", colour = "black") +
  scale_size_continuous(range = c(2, 3), guide = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(caption = "Number of species (number of samples)")

ggsave(filename  =  file.path(subdir, "dataset_summary_plot.png"), p_summary, width  =  4, height = 4)

########################
#### HOST PHYLOGENY ####
########################

# Get consensus tree and fix tip labels
host_consensus <- consensus.edges(host_trees, method = "least.squares", if.absent="zero", rooted = FALSE, check.labels = TRUE)
host_consensus <- root(host_consensus, outgroup = "Macropus_giganteus", resolve.root = TRUE)
host_consensus <- force.ultrametric(host_consensus)

# Change tip labels to match metadata
host_consensus$tip.label <- host_consensus$tip.label %>%
                            gsub(pattern="Tapirus_indicus", replacement="Acrocodia_indica") %>%
                            gsub(pattern="Procolobus_badius", replacement="Piliocolobus_foai") %>%
                            gsub(pattern="Otaria_bryonia", replacement="Otaria_flavescens") %>%
                            gsub(pattern="_", replacement=" ", fixed=TRUE)

# Drop tips not in dataset
host_consensus <- drop.tip(host_consensus, setdiff(host_consensus$tip.label, meta$Species))
host_consensus$tip.label <- gsub("_", " ", host_consensus$tip.label)
write.tree(host_consensus, file = file.path(subdir, "..", "host_consensus.tre"))

# Get traits for tree tips
host_traits <- data.frame(tip = host_consensus$tip.label,
                          node=nodeid(host_consensus, host_consensus$tip.label),
                          Common.name = meta$Common.name[match(host_consensus$tip.label, meta$Species)],
                          Order = meta$Order[match(host_consensus$tip.label, meta$Species)],
                          diet.general = meta$diet.general[match(host_consensus$tip.label, meta$Species)],
                          habitat.general = meta$habitat.general[match(host_consensus$tip.label, meta$Species)],
                          phylopic = phylopics$uid[match(host_consensus$tip.label, phylopics$Species)])

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

# Add order trait to nodes
node_traits <- tips_to_nodes(host_consensus, host_traits, "Order") %>% left_join(host_traits) %>% arrange(node)

# Add traits to tree
host_consensus_traits <- full_join(host_consensus, node_traits, by = "node")

# Plot tree
p <- 
  ggtree(host_consensus_traits, aes(colour=Order), layout = "circular", open.angle = 180, size=1) + 
  geom_tiplab(size=4, nudge_x = 4.5, aes(label = Common.name)) + 
  geom_phylopic(aes(uuid=phylopic, fill=Order), width = 0.7, position=position_nudge(x=3)) +
  # Colour branches by order
  scale_color_manual(values = order_palette, name = "Taxonomic order", na.value = "black") +
  scale_fill_manual(values = order_palette, name = "Taxonomic order") +
  # Add points indicating diet and habitat
  new_scale_colour() +
  new_scale_fill() +
  geom_tippoint(shape=21, size=4, stroke=2, 
                aes(fill=diet.general, colour=habitat.general), position = position_nudge(x = 1)) +
  scale_fill_manual(values = diet_palette, name = "Diet", na.value = "black") +
  scale_colour_manual(values = habitat_palette, name = "Habitat", na.value = "black") +
  # Adjust theme
  theme_tree() +
  theme(plot.margin = margin(t=0, r=3, b=0, l=3, "cm"),
        legend.position = "none")

ggsave(filename =  file.path(subdir, "host_phylogeny.png"), p, width  =  8, height = 8)

##########################
#### SEQUENCING DEPTH ####
##########################

metadata$this_study <- ifelse(metadata$project_name == "Moraitou2025 and this study", "This study", "Previous studies")

p <- ggplot(metadata, aes(y = Common.name, x = unmapped_count, fill = this_study)) +
    geom_boxplot() +
    scale_x_continuous(trans = "log10") +
    scale_fill_manual(values = c(`This study` = "#73B55B", `Previous studies` = "#C96474"), name = "") +
    facet_grid(rows = vars(Order), space = "free", scales = "free")

ggsave(filename =  file.path(subdir, "read_count_per_species.png"), p, width  =  8, height = 12)
