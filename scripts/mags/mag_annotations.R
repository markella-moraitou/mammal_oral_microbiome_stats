##### MAG FUNCTIONS #####

#### Plot and analyse functional capacities of MAGs ####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(distillR)
library(stringr)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(phytools)
library(tibble)
library(ggnewscale)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "mags")) # general output directory
subdir <- normalizePath(file.path("..", "..", "output", "mags", "annotations")) # subdirectory for the output of this script

if(!dir.exists(subdir)) {
  dir.create(subdir, recursive = TRUE)
}

source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))

#######################
#####  LOAD INPUT #####
#######################

# Bin metadata
bac_meta <- read.table(file.path(outdir, "bac_meta.tsv"), sep="\t", header=TRUE)
ar_meta <- read.table(file.path(outdir, "ar_meta.tsv"), sep="\t", header=TRUE)

# MAG trees
bac_tree <- read.tree(file = file.path(outdir, "bac_tree.tree"))
ar_tree <- read.tree(file = file.path(outdir, "ar_tree.tree"))

# DRAM output
mag_dram <- read.table(file.path(indir, "MAG_annotations.tsv"),
            sep = "\t", header = TRUE, check.names=FALSE, quote = "", comment = "", fill = TRUE)
            
colnames(mag_dram)[1] <- "contig"

#### Process metadata table ####
meta <- rbind(bac_meta, ar_meta)
meta_tips <- meta %>% filter(!is.na(bin))

meta_tips <- meta_tips %>% mutate(phylum_grouped = case_when(phylum %in% names(phylum_palette) ~ phylum,
                                                            domain == "Bacteria" ~ "Other Bacteria",
                                                            domain == "Archaea" ~ "Other Archaea"))

###################
#### DistillR #####
###################

mag_dram_filt <- mag_dram %>% filter(fasta %in% meta_tips$bin)

if(file.exists(file.path(subdir, "GIFTs.csv"))) {
  cat("Loading GIFTs from file\n")
  GIFTs <- read.csv(file.path(subdir, "GIFTs.csv"), row.names = 1)
} else {
  cat("Running distillR to get GIFTs\n")
  GIFTs <- distill(mag_dram_filt, GIFT_db, genomecol=2, annotcol=c(9, 10, 20, 21, 22))
  write.csv(GIFTs, file.path(subdir, "GIFTs.csv"), row.names = TRUE)
}

#Aggregate bundle-level GIFTs into the compound level
GIFTs_elements <- to.elements(GIFTs, GIFT_db)

GIFTs_elements_long <- 
  GIFTs_elements %>% data.frame %>% rownames_to_column("Sample") %>%
  pivot_longer(cols = -c(Sample), names_to = "Code_element", values_to = "Completeness") %>%
  inner_join(select(meta_tips, c(label, bin, domain:species, phylum_grouped, host_species, host_common.name, host_order, diet.general)), by=c("Sample"="bin")) %>%
  left_join(unique(select(GIFT_db, c("Code_element", "Element", "Function", "Domain")))) %>%
  # Turn function into factor
  arrange(Domain, Function) %>% mutate(Function = factor(Function, levels = unique(Function)))

p <- ggplot(GIFTs_elements_long, aes(y=label, x=Element)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Compound", fill="Completeness") +
  facet_grid(rows = vars(phylum_grouped), cols = vars(Function), scales = "free", space = "free")
  
ggsave(p, filename = file.path(subdir, "GIFTs_elements.png"), width = 25, height = 20)

#Aggregate element-level GIFTs into the function level
GIFTs_functions <- to.functions(GIFTs_elements, GIFT_db)

GIFTs_functions_long <- 
  GIFTs_functions %>% data.frame %>% rownames_to_column("Sample") %>%
  pivot_longer(cols = -c(Sample), names_to = "Code_function", values_to = "Completeness") %>%
  inner_join(select(meta_tips, c(label, bin, domain:species, phylum_grouped, host_species, host_common.name, host_order, diet.general)), by=c("Sample"="bin")) %>%
  left_join(unique(select(GIFT_db, c("Code_function", "Function", "Domain"))))

p <- ggplot(GIFTs_functions_long, aes(y=label, x=Function)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Compound", fill="Completeness") +
  facet_grid(rows = vars(phylum_grouped), cols = vars(Domain), scales = "free", space = "free")

ggsave(p, filename = file.path(subdir, "GIFTs_functions.png"), width = 10, height = 15)

# Keep only degradation and biosynthesis elements for plotting
gifts_func_plotting <- GIFTs_functions_long %>%
  filter(Domain %in% c("Degradation", "Biosynthesis")) %>%
  arrange(Domain, Function) %>%
  mutate(Function = factor(Function, levels = unique(Function)),
         Completeness = as.numeric(Completeness))

gift_func_degradation <- gifts_func_plotting %>% filter(Domain == "Degradation") %>% mutate(Function = droplevels(Function))
gift_func_biosynthesis <- gifts_func_plotting %>% filter(Domain == "Biosynthesis") %>% mutate(Function = droplevels(Function))

# Same for elements
# Keep only degradation and biosynthesis elements for plotting
gift_elem_plotting <- GIFTs_elements_long %>%
  filter(Domain %in% c("Degradation", "Biosynthesis")) %>%
  arrange(Domain, Function, Element) %>%
  mutate(Element = factor(Element, levels = unique(Element)),
         Completeness = as.numeric(Completeness))

gift_elem_degradation <- gift_elem_plotting %>% filter(Domain == "Degradation") %>% mutate(Element = droplevels(Element))
gift_elem_biosynthesis <- gift_elem_plotting %>% filter(Domain == "Biosynthesis") %>% mutate(Element = droplevels(Element))

###################
#### PLOT TREE ####
###################

#### Bacteria ####
# Colour by host order and diet
p <- ggtree(bac_tree, layout = "circular", size = 1.5) %<+% select(bac_meta, c(label, host_order, habitat.general, diet.general)) +
  geom_tiplab(size=2, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order, color=diet.general)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_colour_manual(values = diet_palette, name = "Host diet", na.value = "black") +
  scale_x_continuous(expand = c(0, 0)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(-6, -6, -6, 1), "cm"), # Remove margins
  legend.position=c(0.05, 0.5),
  legend.text = element_text(size=20),
  legend.title = element_text(size=20)) +
  guides(fill = guide_legend(override.aes = list(size = 5)), 
  color = guide_legend(override.aes = list(size = 5)))

# Add GIFTs as heatmap
p <- p +
    new_scale_fill() +
    geom_fruit(data = gift_func_degradation, geom=geom_tile, mapping = aes(y=label, x=Function, fill=Completeness),
             offset = 0.06, pwidth = 0.15, axis.params=list(axis = "x", text.angle = 90, text.size = 2.5, colour = "white")) +
    scale_fill_viridis_c(name = "Completeness", option = "C") +
    geom_fruit(data = gift_func_biosynthesis, geom=geom_tile, mapping = aes(y=label, x=Function, fill=Completeness),
             offset = 0.02, pwidth = 0.15, axis.params=list(axis = "x", text.angle = 90, text.size = 2.5, colour = "white"))

ggsave(p, file=file.path(subdir, "bac_function_tree.png"), width = 22, height = 20)

#### Archaea ####
# Colour by host order and diet
p <- ggtree(ar_tree) %<+% select(ar_meta, c(label, host_order, habitat.general, diet.general)) +
  geom_tiplab(size=3, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order, color=diet.general)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_colour_manual(values = diet_palette, name = "Host diet", na.value = "black") +
  scale_x_continuous(expand = c(0.04, 0.04)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(0, 0, 0, 3), "cm"),
        legend.position=c(-0.01, 0.5),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)), 
         color = guide_legend(override.aes = list(size = 2.5))) 

# Add GIFTs as heatmap
p <- p +
    new_scale_fill() +
    geom_fruit(data = gift_func_degradation, geom=geom_tile, mapping = aes(y=label, x=Function, fill=Completeness),
             offset = 0.15, pwidth = 0.2, axis.params=list(axis = "x", text.angle = 90, text.size = 2.5, colour = "black")) +
    scale_fill_viridis_c(name = "Completeness", option = "C") +
    geom_fruit(data = gift_func_biosynthesis, geom=geom_tile, mapping = aes(y=label, x=Function, fill=Completeness),
             offset = 0.05, pwidth = 0.2, axis.params=list(axis = "x", text.angle = 90, text.size = 2.5, colour = "black"))

ggsave(p, file=file.path(subdir, "ar_function_tree.png"), width = 15, height = 10)

############################
#### STREPTOCOCCUS TREE ####
############################

# Subset to Streptococcus
strep_tree <- bac_tree %>% 
  keep.tip(bac_tree$tip.label[grep("Streptococcus", bac_tree$tip.label)])

# Colour by order and habitat
p <- ggtree(strep_tree) %<+% select(bac_meta, c(label, host_order, diet.general)) +
  geom_tiplab(size=3, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order, color=diet.general)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_colour_manual(values = diet_palette, name = "Host diet", na.value = "black") +
  scale_x_continuous(expand = c(0.04, 0.04)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(0, 0, 0, 3), "cm"),
        legend.position=c(-0.01, 0.5),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)), 
         color = guide_legend(override.aes = list(size = 2.5))) 

# Add GIFTs as heatmap
p <- p +
    new_scale_fill() +
    geom_fruit(data = gift_elem_degradation, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.1, pwidth = 2, axis.params=list(axis = "x", text.angle = 90, text.size = 2, colour = "black")) +
    scale_fill_viridis_c(name = "Completeness", option = "C") +
    geom_fruit(data = gift_elem_biosynthesis, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.05, pwidth = 2, axis.params=list(axis = "x", text.angle = 90, text.size = 2, colour = "black"))

ggsave(p, file=file.path(subdir, "streptococcus_function_tree.png"), width = 20, height = 15)

##########################
#### ACTINOMYCES TREE ####
##########################

# Subset to Actinomyces
actino_tree <- bac_tree %>% 
  keep.tip(bac_tree$tip.label[grep("Actinomyces", bac_tree$tip.label)])

# Colour by order and habitat
p <- ggtree(actino_tree) %<+% select(bac_meta, c(label, host_order, diet.general)) +
  geom_tiplab(size=3, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order, color=diet.general)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_colour_manual(values = diet_palette, name = "Host diet", na.value = "black") +
  scale_x_continuous(expand = c(0.04, 0.04)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(0, 0, 0, 3), "cm"),
        legend.position=c(-0.01, 0.5),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)), 
         color = guide_legend(override.aes = list(size = 2.5)))

# Add GIFTs as heatmap
p <- p +
    new_scale_fill() +
    geom_fruit(data = gift_elem_degradation, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.2, pwidth = 2, axis.params=list(axis = "x", text.angle = 90, text.size = 2, colour = "black")) +
    scale_fill_viridis_c(name = "Completeness", option = "C") +
    geom_fruit(data = gift_elem_biosynthesis, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.05, pwidth = 2, axis.params=list(axis = "x", text.angle = 90, text.size = 2, colour = "black"))

ggsave(p, file=file.path(subdir, "actinomyces_function_tree.png"), width = 20, height = 15)

#####################
#### ROTHIA TREE ####
#####################

# Subset to Rothia
rothia_tree <- bac_tree %>% 
  keep.tip(bac_tree$tip.label[grep("Rothia", bac_tree$tip.label)])

# Colour by order and habitat
p <- ggtree(rothia_tree) %<+% select(bac_meta, c(label, host_order, diet.general)) +
  geom_tiplab(size=3, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order, color=diet.general)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_colour_manual(values = diet_palette, name = "Host diet", na.value = "black") +
  scale_x_continuous(expand = c(0.04, 0.04)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(0, 0, 0, 3), "cm"),
        legend.position=c(-0.01, 0.5),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)), 
         color = guide_legend(override.aes = list(size = 2.5)))

# Add GIFTs as heatmap
p <- p +
    new_scale_fill() +
    geom_fruit(data = gift_elem_degradation, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.2, pwidth = 2, axis.params=list(axis = "x", text.angle = 90, text.size = 2, colour = "black")) +
    scale_fill_viridis_c(name = "Completeness", option = "C") +
    geom_fruit(data = gift_elem_biosynthesis, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.05, pwidth = 2, axis.params=list(axis = "x", text.angle = 90, text.size = 2, colour = "black"))

ggsave(p, file=file.path(subdir, "rothia_function_tree.png"), width = 20, height = 15)

########################
#### NEISSERIA TREE ####
########################

# Subset to Neisseria
neisseria_tree <- bac_tree %>% 
  keep.tip(bac_tree$tip.label[grep("Neisseria", bac_tree$tip.label)])

# Colour by order and habitat
p <- ggtree(neisseria_tree) %<+% select(bac_meta, c(label, host_order, diet.general)) +
  geom_tiplab(size=3, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order, color=diet.general)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_colour_manual(values = diet_palette, name = "Host diet", na.value = "black") +
  scale_x_continuous(expand = c(0.04, 0.04)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(0, 0, 0, 3), "cm"),
        legend.position=c(-0.01, 0.5),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)), 
         color = guide_legend(override.aes = list(size = 2.5)))

# Add GIFTs as heatmap
p <- p +
    new_scale_fill() +
    geom_fruit(data = gift_elem_degradation, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.2, pwidth = 2, axis.params=list(axis = "x", text.angle = 90, text.size = 2, colour = "black")) +
    scale_fill_viridis_c(name = "Completeness", option = "C") +
    geom_fruit(data = gift_elem_biosynthesis, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.05, pwidth = 2, axis.params=list(axis = "x", text.angle = 90, text.size = 2, colour = "black"))

ggsave(p, file=file.path(subdir, "neisseria_function_tree.png"), width = 20, height = 15)

################################
#### PROPIONIBACTERIUM TREE ####
################################

# Subset to Propionibacterium
propio_tree <- bac_tree %>% 
  keep.tip(bac_tree$tip.label[grep("Propionibacterium", bac_tree$tip.label)])

# Colour by order and habitat
p <- ggtree(propio_tree) %<+% select(bac_meta, c(label, host_order, diet.general)) +
  geom_tiplab(size=3, aes(colour=host_order)) +
  scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
  new_scale_color() +
  geom_tippoint(shape=21, size=3, stroke=1, 
                aes(fill=host_order, color=diet.general)) +
  scale_fill_manual(values = order_palette, name = "Host order", na.value = "black") +
  scale_colour_manual(values = diet_palette, name = "Host diet", na.value = "black") +
  scale_x_continuous(expand = c(0.04, 0.04)) +  # Adjust the x-axis scaling 
  theme(plot.margin = unit(c(0, 0, 0, 3), "cm"),
        legend.position=c(-0.01, 0.5),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5)), 
         color = guide_legend(override.aes = list(size = 2.5)))

# Add GIFTs as heatmap
p <- p +
    new_scale_fill() +
    geom_fruit(data = gift_elem_degradation, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.2, pwidth = 2, axis.params=list(axis = "x", text.angle = 90, text.size = 2, colour = "black")) +
    scale_fill_viridis_c(name = "Completeness", option = "C") +
    geom_fruit(data = gift_elem_biosynthesis, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.05, pwidth = 2, axis.params=list(axis = "x", text.angle = 90, text.size = 2, colour = "black"))

ggsave(p, file=file.path(subdir, "propionibacterium_function_tree.png"), width = 20, height = 15)
