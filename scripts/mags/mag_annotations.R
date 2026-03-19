##### MAG FUNCTIONS #####

#### Plot and analyse functional capacities of MAGs ####

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(readr)
library(distillR)
library(stringr)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(rphylopic)
library(cowplot)
library(phytools)
library(tibble)
library(ggnewscale)
library(vegan)

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
mag_dram <- read_tsv(file.path(indir, "MAG_annotations.tsv.gz"), quote = "", comment = "")

colnames(mag_dram)[1] <- "contig"

# Phylopics
phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

#### Process metadata table ####
meta <- rbind(bac_meta, ar_meta)
meta_tips <- meta %>% filter(!is.na(bin))

meta_tips <- meta_tips %>% mutate(phylum_grouped = case_when(phylum %in% names(phylum_palette) ~ phylum,
                                                            domain == "Bacteria" ~ "Other Bacteria",
                                                            domain == "Archaea" ~ "Other Archaea"))

meta_tips$phylum_grouped <- factor(meta_tips$phylum_grouped, levels = names(phylum_palette))

codiv_results <- read.csv(file.path(outdir, "codiversification", "bac_cod_results.csv")) %>%
      rbind(read.csv(file.path(outdir, "codiversification", "ar_cod_results.csv")))

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
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90)) +
  labs(y="Sample", x="Compound", fill="Completeness") +
  facet_grid(rows = vars(phylum_grouped),
             cols = vars(gsub(
                gsub(Function, pattern = " biosynthesis", replacement = "\nbiosynthesis"),
                pattern = " degradation", replacement = "\ndegradation")),
             scales = "free", space = "free")
  
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

###################
#### PLOT TREE ####
###################

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

#####################
#### PCA OF BINS ####
#####################

# Only for HQ bins
hq_meta <- meta_tips %>% filter(Completeness >= 90 & Contamination <= 5)

GIFTs_elements <- GIFTs_elements[hq_meta$bin, ]

## Plot PCAs of bins based on GIFT completeness

order_shape_scale <- c("Carnivora" = 4, "Primates" = 19, "Artiodactyla" = 5, "Perissodactyla" = 2, "Rodentia" = 1, "Rest" = 12, "Proboscidea_Sirenia" = 12)
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

arrows_func <- function(ord, axes = c(1,2), scale = 1) {
   scores = data.frame(vegan::scores(ord, choices = axes, display = "species")*scale)
   scores$distance = sqrt(scale*scores[,1]^2 + scores[,2]^2)
   scores <- scores[order(scores$distance, decreasing = TRUE),]
}

#### ALL BINS ####
ord <- pca(GIFTs_elements)

ord_df <- scores(ord)$sites %>% data.frame %>% rownames_to_column("bin") %>% left_join(hq_meta, by=c("bin"="bin"))

var_explained <- ord$CA$eig / sum(ord$CA$eig) * 100

# Get centroids for genera with at least 3 members
family_centroids <- ord_df %>% group_by(family, phylum_grouped) %>% summarise(PC1 = mean(PC1), PC2 = mean(PC2), n = n()) %>% filter(n >= 3)

p <- ggplot(data = ord_df, aes(x = PC1, y = PC2, colour = phylum_grouped, shape = host_order)) +
  geom_point(size=2, alpha=0.7) +
  scale_colour_manual(values = phylum_palette, name = "Phylum", na.value = "black") +
  scale_shape_manual(values = order_shape_scale, name = "Host order", na.value = 12) +
  geom_label(data = family_centroids, aes(x = PC1, y = PC2, label = family, color = phylum_grouped, shape = NULL), size=3, fill="white", alpha=0.7) +
  xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2], 1), "%)"))

# Plot arrows for 10 GIFTs with highest loadings
arrows <- arrows_func(ord, scale = 0.5) %>% head(10)
arrows <- arrows %>% rownames_to_column("Code_element") %>%
  left_join(unique(select(GIFT_db, c("Code_element", "Element", "Domain"))), by=c("Code_element"="Code_element")) %>%
  mutate(label = paste(Element, tolower(Domain)))

p <- p +
  geom_segment(inherit.aes = FALSE, data = arrows, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", alpha = 0.8) +
  geom_text(inherit.aes = FALSE, data = arrows, aes(x = PC1, y = PC2 , label = label,
            vjust = ifelse(PC2 < 0, 1, 0), hjust = ifelse(PC1 < 0, 1, 0)),
            size = 3, color = "black", alpha = 0.8)

ggsave(p, file=file.path(subdir, "GIFTs_PCA_bins.png"), width = 10, height = 10)

#### TOP FAMILIES ####
# Separately for the families with most bins
top_families <- family_centroids %>% arrange(desc(n)) %>% pull(family) %>% head(10)

plot_list <- list()

# Subset to each family and run PCA
for (fam in top_families) {
  filt_meta <- hq_meta %>% filter(family == fam)
  filt_gifts <- GIFTs_elements[filt_meta$bin, ]
  ord_fam <- pca(filt_gifts)
  ord_fam_df <- scores(ord_fam)$sites %>% data.frame %>% rownames_to_column("bin") %>% left_join(filt_meta, by=c("bin"="bin"))
  var_explained_fam <- ord_fam$CA$eig / sum(ord_fam$CA$eig) * 100
  p_fam <- ggplot(data = ord_fam_df, aes(x = PC1, y = PC2, colour = diet.general, shape = host_order)) +
    geom_point(size=2, alpha=0.7) +
    scale_colour_manual(values = diet_palette, name = "Host diet", na.value = "black") +
    scale_shape_manual(values = order_shape_scale, name = "Host order", na.value = 12) +
    xlab(paste0("PC1 (", round(var_explained_fam[1], 1), "%)")) +
    ylab(paste0("PC2 (", round(var_explained_fam[2], 1), "%)")) +
    ggtitle(paste0("Family: ", fam)) + theme(legend.position = "none")
  # Plot arrows for 5 GIFTs with highest loadings
  arrows <- arrows_func(ord_fam, scale = 0.5) %>% head(5)
  arrows <- arrows %>% rownames_to_column("Code_element") %>%
    left_join(unique(select(GIFT_db, c("Code_element", "Element", "Domain"))), by=c("Code_element"="Code_element")) %>%
    mutate(label = paste(Element, tolower(Domain)))
  p_fam <- p_fam +
    geom_segment(inherit.aes = FALSE, data = arrows, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black", alpha = 0.8) +
    geom_text(inherit.aes = FALSE, data = arrows, aes(x = PC1, y = PC2, label = label,
              vjust = ifelse(PC2 < 0, 1, 0), hjust = ifelse(PC1 < 0, 1, 0)),
              size = 1.5, color = "black", alpha = 0.8)
    plot_list[[fam]] <- p_fam
}

p <- plot_grid(plotlist = plot_list, ncol = 2)

ggsave(p, file=file.path(subdir, "GIFTs_PCA_top_families.png"), width = 10, height = 15)

#### For codiversifying clades

cod_nodes <- codiv_results %>% filter(p.adjust < 0.05) %>% pull(node) %>% unique()

# Get tips of codiversifying clades
cod_clades <- list()
for (node in cod_nodes) {
  tips <- extract.clade(bac_tree, node)$tip.label
  cod_clades[[as.character(node)]] <- tips
}

plot_list <- list()

# Subset to each family and run PCA
for (clade in names(cod_clades)) {
  filt_meta <- hq_meta %>% filter(label %in% cod_clades[[clade]])
  filt_gifts <- GIFTs_elements[filt_meta$bin, ]
  ord_clade <- pca(filt_gifts)
  ord_clade_df <- scores(ord_clade)$sites %>% data.frame %>% rownames_to_column("bin") %>% left_join(filt_meta, by=c("bin"="bin"))
  var_explained <- ord_clade$CA$eig / sum(ord_clade$CA$eig) * 100
  p_clade <- ggplot(data = ord_clade_df, aes(x = PC1, y = PC2, colour = host_order, shape = diet.general)) +
    geom_point(size=2, alpha=0.7) +
    scale_shape_manual(values = diet_shape_scale, name = "Host diet", na.value = "black") +
    scale_colour_manual(values = order_palette, name = "Host order", na.value = 12) +
    xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
    ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) +
    ggtitle(paste0("Clade: ", clade)) + theme(legend.position = "none")
  arrows <- arrows_func(ord_clade, scale = 0.5) %>% head(5)
  arrows <- arrows %>% rownames_to_column("Code_element") %>%
    left_join(unique(select(GIFT_db, c("Code_element", "Element", "Domain"))), by=c("Code_element"="Code_element")) %>%
    mutate(label = paste(Element, tolower(Domain)))
  p_clade <- p_clade +
    geom_segment(inherit.aes = FALSE, data = arrows, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black", alpha = 0.8) +
    geom_text(inherit.aes = FALSE, data = arrows, aes(x = PC1, y = PC2, label = label,
             vjust = ifelse(PC2 < 0, 1, 0), hjust = ifelse(PC1 < 0, 1, 0)),
             size = 1.5, color = "black", alpha = 0.8)
  plot_list[[clade]] <- p_clade
}

p <- plot_grid(plotlist = plot_list, ncol = 3)

ggsave(p, file=file.path(subdir, "GIFTs_PCA_codiversifying.png"), width = 15, height = 15)

####################
#### PLOT TREES ####
####################

# Plot heatmaps of GIFTs on trees of codiversifying clades
hq_meta$uid <- phylopics$uid[match(hq_meta$host_species, phylopics$Species)]

# Subset to each family and run PCA
for (clade in names(cod_clades)) {
  filt_meta <- hq_meta %>% filter(label %in% cod_clades[[clade]])
  
  # Keep GIFTs for degradation or biosynthesis that vary across this clade
  filt_degr <- gift_elem_degradation %>% filter(label %in% filt_meta$label) %>% filter(Code_element %in% colnames(GIFTs_elements)) %>%
    group_by(Code_element) %>% filter(var(Completeness, na.rm=TRUE) > 0 & max(Completeness > 0.8)) %>% ungroup() %>%
    mutate(Element = droplevels(Element))
  
  filt_bios <- gift_elem_biosynthesis %>% filter(label %in% filt_meta$label) %>% filter(Code_element %in% colnames(GIFTs_elements)) %>%
    group_by(Code_element) %>% filter(var(Completeness, na.rm=TRUE) > 0 & max(Completeness > 0.8)) %>% ungroup() %>%
    mutate(Element = droplevels(Element))
  
  n_degr <- nlevels(filt_degr$Element)
  n_bios <- nlevels(filt_bios$Element)
  n_elements <- n_degr + n_bios
  
  # Subset tree
  sub_tree <- bac_tree %>% keep.tip(filt_meta$label)
  
  # Plot tree
  p <- ggtree(sub_tree) %<+% select(filt_meta, c(label, host_order, diet.general, uid)) +
    geom_tiplab(size=5, aes(colour=host_order)) +
    scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
    new_scale_color() +
    geom_phylopic(aes(uuid = uid, fill = host_order), size = 0.5) +
    geom_tippoint(aes(color=diet.general), size = 1) +
    scale_colour_manual(values = order_palette, name = "Host order", na.value = "black") +
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
    geom_fruit(data = filt_degr, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.6, pwidth = 4*(n_degr/n_elements), axis.params=list(axis = "x", text.angle = 90, text.size = 6, colour = "black", hjust=1)) +
    scale_fill_viridis_c(name = "Completeness", option = "C") +
    geom_fruit(data = filt_bios, geom=geom_tile, mapping = aes(y=label, x=Element, fill=Completeness),
             offset = 0.1, pwidth = 4*(n_bios/n_elements), axis.params=list(axis = "x", text.angle = 90, text.size = 6, colour = "black", hjust=1)) +
    theme(plot.margin = unit(c(0, -2, 6, 3), "cm")) +  # Ensure bottom margin is large enough
    coord_cartesian(clip = "off")

  ggsave(p, file=file.path(subdir, paste0(clade, "_function_tree.png")), width = ceiling(n_elements/3) + 6, height = ceiling(length(sub_tree$tip.label)/3 + 3))
}
