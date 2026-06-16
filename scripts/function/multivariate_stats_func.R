##### EXPLORATORY PLOTS #####

#### Plots exploratory plots such as ordinations and heatmaps

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(phyloseq)
library(microViz)
library(microbiome)
library(rphylopic)
library(ggplot2)
library(RColorBrewer)
library(ggnewscale)
library(scales)
library(ggExtra)
library(distillR)
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
datadir <- normalizePath(file.path("..", "..", "output", "function", "data"))
pathdir <- normalizePath(file.path("..", "..", "output", "function", "pathway_completeness")) # Directory with pathway analysis output
subdir <- normalizePath(file.path("..", "..", "output", "function", "multivariate_func")) # subdirectory for the output of this script

dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

# Get ordination functions
source(file.path("..","ordination_functions.R"))

#######################
#####  LOAD INPUT #####
#######################

# Phyloseq objects
phy_gene_f <- readRDS(file.path(datadir, "phy_gene_f.RDS"))
phy_gene_f_clr <- readRDS(file.path(datadir, "phy_gene_f_clr.RDS"))

phy_pathway <- readRDS(file.path(pathdir, "phy_pathway.RDS"))
phy_pathway_clr <- readRDS(file.path(pathdir, "phy_pathway_clr.RDS"))

phy_gifts_el <- readRDS(file.path(datadir, "phy_gifts_el.RDS"))

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

###################
#### RDA PLOTS ####
###################

#### Get shape scales for plotting ####
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)
order_shape_scale <- c("Carnivora" = 4, "Primates" = 19, "Artiodactyla" = 5, "Perissodactyla" = 2, "Rodentia" = 1, "Rest" = 12, "Proboscidea_Sirenia" = 12)

#### GENE ABUNDANCE (CLR) ####

# Recode order and habitat as TRUE and FALSE
phy_gene_f_clr <- phy_gene_f_clr %>%
        ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  Ruminant = (digestion == "Ruminant"),
                  Marine = (habitat.general == "Marine"))

# Species traits to use as constraints
species_traits <- c("Artiodactyla", "Perissodactyla", "Primates",
                    "Ruminant", "Marine", "Frugivory", "Animalivory")

# Ordinate using all data
ord <- ord_calc(phy_gene_f_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:length(species_traits)))

ggsave(file.path(subdir, "screeplot_genes.png"), p, width=2, height=2)

## SAMPLE PLOTS

# Color by diet
p <- custom_ord_plot(phy_gene_f_clr, ord, colour="diet.general", shape="Order_grouped", type = "RDA")

ggsave(p, filename = file.path(subdir, "gene_ordination_diet.png"), width=8, height=6)

ggsave(p, filename = file.path(subdir, "gene_ordination_diet.svg"), width=8, height=6)

# Colour by order
p <- custom_ord_plot(phy_gene_f_clr, ord, colour="Order_grouped", shape="diet.general", type = "RDA")

ggsave(p, filename = file.path(subdir, "gene_ordination_order.png"), width=8, height=6)

ggsave(p, filename = file.path(subdir, "gene_ordination_order.svg"), width=8, height=6)

## PLOT ARROWS

# Get loading arrows coordinaties
arrows <- arrow_coord(ord@ord, axes = c(1, 2))

# Get gene category
arrows$category <- as.character(phy_gene_f_clr@tax_table[match(rownames(arrows),  rownames(phy_gene_f_clr@tax_table)), "category"])

arrows$to_plot <- (rownames(arrows) %in% head(rownames(arrows), 500))

# Save
write.csv(rownames_to_column(arrows, "gene"), file.path(subdir, "gene_ordination_arrows.txt"), quote = FALSE, row.names = FALSE)

# Keep only strongest associations
arrows_filt <- arrows %>% filter(to_plot) %>%
              select(contains(c("1", "2")), category)

# Group uncommon categories
common_categories <- table(arrows_filt$category) %>% sort(decreasing = TRUE) %>% head(6) %>% names

arrows_filt <- arrows_filt %>%
    mutate(category_grouped = factor(case_when(category %in% common_categories ~ category,
                                            TRUE ~ "Other"), levels = c(common_categories, "Other")))

# Set colours for categories using colour brewer
arrow_colours <- brewer.pal(n = length(unique(arrows_filt$category_grouped)), name = "Dark2")
names(arrow_colours) <- unique(arrows_filt$category_grouped) # Remove "Other" from names
arrow_colours["Other"] <- "grey90" # Set "Other" to grey

p <- ggplot(data = arrows_filt) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = category_grouped), linewidth = 0.5, alpha = 0.5) +
  scale_color_manual(values = arrow_colours, name = "Category") +
  xlab("RDA1") + ylab("RDA2") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 2))

ggsave(p, filename = file.path(subdir, "gene_ordination_arrows.png"), width=8, height=6)

#### PATH ABUNDANCE (CLR) ####

phy_pathway_clr <- phy_pathway_clr %>% 
        ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  Ruminant = (digestion == "Ruminant"),
                  Marine = (habitat.general == "Marine"))

# Ordinate using all data
ord <- ord_calc(phy_pathway_clr, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:length(species_traits)))

ggsave(file.path(subdir, "screeplot_pathways.png"), p, width=2, height=2)

## SAMPLE PLOTS

# Color by diet
p <- custom_ord_plot(phy_pathway_clr, ord, colour="diet.general", shape="Order_grouped", type = "RDA")

ggsave(file.path(subdir, "pathway_ordination_diet.png"), p, width=8, height=6)

# Colour by order
p <- custom_ord_plot(phy_pathway_clr, ord, colour="Order_grouped", shape="diet.general", type = "RDA")

ggsave(file.path(subdir, "pathway_ordination_order.png"), p, width=8, height=6)

## Plot arrows
# Get scores and add info
taxa_rda <- data.frame(vegan::scores(ord@ord, display="species", choices=1:3)) %>%
            cbind(tax_table(phy_pathway_clr)[, c("path_name", "path_class")]) %>%
            rownames_to_column("OTU") %>%
            # Get grouped categories
            mutate(category = str_remove(path_class, ".*; "))

# Identify 10 pathways with the highest correlation
top_taxa <- taxa_rda %>% arrange(desc(sqrt(RDA1^2 + RDA2^2))) %>%
              slice_head(n = 10) %>% pull(path_name)

taxa_rda <- taxa_rda %>% mutate(label = ifelse(path_name %in% top_taxa, path_name, NA)) %>%
            mutate(label = factor(path_name, levels = top_taxa)) %>%
            mutate(linetype = ifelse(path_name %in% top_taxa, "solid", "dashed"))

# Group uncommon categories for labelled taxa
common_categories <- table(taxa_rda$category[!is.na(taxa_rda$label)]) %>% sort(decreasing = TRUE) %>% head(7) %>% names

taxa_rda <- taxa_rda %>%
    mutate(category_grouped = factor(case_when(category %in% common_categories ~ category,
                                            TRUE ~ "Other"), levels = c(rev(common_categories), "Other")))  

# Save
write.csv(taxa_rda, file.path(subdir, "pathway_ordination_arrows.txt"), quote = FALSE, row.names = FALSE)

# Set colours for categories using colour brewer
arrow_colours <- brewer.pal(n = length(unique(taxa_rda$category_grouped))-1, name = "Dark2")
names(arrow_colours) <- unique(taxa_rda$category_grouped)[-length(unique(taxa_rda$category_grouped))] # Remove "Other" from names
arrow_colours["Other"] <- "grey70" # Set "Other" to grey

# Plot
set.seed(245)

p <- ggplot(taxa_rda, aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = category_grouped)) +
  geom_segment(linewidth = 0.5, alpha = 0.8, aes(linetype = linetype)) +
  scale_linetype_identity() +
  scale_color_manual(values = arrow_colours, name = "Pathway class") +
  geom_label(aes(label = label, x = RDA1, y = RDA2),
            size = 2, alpha = 0.5, vjust = ifelse(taxa_rda$RDA2 < 0, 1, 0)) +
  custom_theme() + xlab("RDA1 scores") + ylab("RDA2 scores") +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 7)) +
  guides(colour = guide_legend(ncol = 1)) +
  xlim(min(taxa_rda$RDA1)*1.5, max(taxa_rda$RDA1)*1.2)

ggsave(p, filename = file.path(subdir, "pathway_ordination_arrows.png"), width=8, height=6)

#### GIFTs (distillR) ####

phy_gifts_el <- phy_gifts_el %>% 
        ps_mutate(Artiodactyla = (Order == "Artiodactyla"),
                  Carnivora = (Order == "Carnivora"),
                  Perissodactyla = (Order == "Perissodactyla"),
                  Primates = (Order == "Primates"),
                  Ruminant = (digestion == "Ruminant"),
                  Marine = (habitat.general == "Marine"))

# Ordinate using all data
ord <- ord_calc(phy_gifts_el, constraints = species_traits, method = "RDA")

# Select variables and check for collinearity
ord_step <- step(ord@ord, scope = formula(ord@ord), test = "perm")
vif.cca(ord_step)

# Scree plot
p <- ord %>% ord_get() %>% plot_scree() + custom_theme() +
            xlim(paste0("RDA", 1:length(species_traits)))

ggsave(file.path(subdir, "screeplot_gifts.png"), p, width=2, height=2)

## SAMPLE PLOTS

# Color by diet
p <- custom_ord_plot(phy_gifts_el, ord, colour="diet.general", shape="Order_grouped", type = "RDA")

ggsave(file.path(subdir, "gift_ordination_diet.png"), p, width=8, height=6)

# Colour by order
p <- custom_ord_plot(phy_gifts_el, ord, colour="Order_grouped", shape="diet.general", type = "RDA")

ggsave(file.path(subdir, "gift_ordination_order.png"), p, width=8, height=6)

## Plot arrows

# Get loading arrows coordinaties
arrows <- arrow_coord(ord@ord, axes = c(1, 2))

element_info <- GIFT_db %>% select(Code_element, Element, Function, Domain) %>%
                  filter(Code_element %in% rownames(arrows)) %>%
                  distinct() %>% column_to_rownames("Code_element")

# Get gene category
arrows <- arrows %>% cbind(element_info[rownames(arrows),])

arrows$plot_label <- (rownames(arrows) %in% head(rownames(arrows), 10))

# Save
write.csv(taxa_rda, file.path(subdir, "gift_ordination_arrows.txt"), quote = FALSE, row.names = FALSE)

# Keep only strongest associations
arrows_sum <- arrows %>% group_by(Function, Domain) %>%
              summarise(RDA1 = mean(RDA1), RDA2 = mean(RDA2)) %>%
              arrange(desc(sqrt(RDA1^2 + RDA2^2))) %>% head(10)

p <- ggplot(data = arrows) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = Domain), linewidth = 0.5, alpha = 0.3, linetype = "dashed") +
  geom_label(data = arrows[arrows$plot_label,], aes(x = RDA1*1.1, y = RDA2*1.1, label = Element, colour = Domain), size = 2, fill = "white", alpha = 0.7) +
  geom_segment(data = arrows_sum, aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = Domain), linewidth = 1, alpha = 0.7) +
  geom_label(data = arrows_sum, aes(x = RDA1*1.1, y = RDA2*1.1, label = str_remove(Function, tolower(Domain)), fill = Domain), size = 2, colour = "white", alpha = 0.7) +
  scale_colour_manual(values = c("Degradation" = "#D55E00", "Biosynthesis" = "#0072B2", "Structure" = "#6A6A6A"), name = "") +
  scale_fill_manual(values = c("Degradation" = "#D55E00", "Biosynthesis" = "#0072B2", "Structure" = "#6A6A6A"), name = "") +
  xlab("RDA1") + ylab("RDA2") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 1))

ggsave(p, filename = file.path(subdir, "gift_ordination_arrows.png"), width=8, height=6)

#### Heatmaps ####
# Heatmaps based on the pathways and gifts with the largest loadings along the first axis

# Pathways
paths_plot <- taxa_rda %>% arrange(desc(RDA1^2)) %>% slice_head(n = 6) %>% pull(OTU)

paths_heatmap <- prune_taxa(paths_plot, phy_pathway_clr) %>% psmelt %>%
                group_by(OTU, path_name, path_class, Animalivory, Common.name, Species, diet.general) %>%
                summarise(Abundance = mean(Abundance)) %>%
                # Order element by animalivory
                arrange(desc(Animalivory), Common.name, path_name) %>%
                mutate(path_name = paste(path_name, str_remove(OTU, "path:"))) %>%
                mutate(Common.name = factor(Common.name, levels = unique(Common.name)),
                       path_name = factor(path_name, levels = unique(path_name)),
                       category = str_remove(path_class, ".*; ")) %>%
                mutate(diet_short = recode(diet.general,
                                            "Animalivore" = "An.",
                                            "Omnivore" = "Omn.",
                                            "Herbivore" = "Herbivore",
                                            "Frugivore" = "Frugivore")) %>%
                mutate(diet_short = factor(diet_short, levels = rev(levels(diet_short)))) %>% ungroup()

p1 <- ggplot(paths_heatmap, aes(x = Common.name, y = path_name, fill = Abundance)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Pathway\nabundance(CLR) ", breaks = c(-9, -4, 1)) +
    facet_grid(rows = vars(category),
               cols = vars(diet_short), scale = "free", space = "free",
               labeller = as_labeller(function(x) str_wrap(x, width = 25))) + 
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
    theme(legend.position = "right", legend.title = element_text(angle = -90, size = 10), legend.text = element_text(size = 10), legend.title.position = "right",
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 11), axis.title = element_blank(),
          strip.background.y = element_rect(fill = "white", colour = "white"),
          strip.text.y = element_text(angle = 0, hjust = 0), plot.margin = margin(t = 0, r = 10, b = 0, l = 20))

# GIFTs
gifts_plot <- arrows %>% arrange(desc(RDA1^2)) %>% slice_head(n = 6) %>% rownames

gifts_heatmap <- # Keep only elements with strong associations with RDA axes
                prune_taxa(gifts_plot, phy_gifts_el) %>% psmelt %>%
                rename(Code_element = OTU, Completeness = Abundance) %>%
                left_join(rownames_to_column(element_info, "Code_element")) %>%
                group_by(Element, Function, Domain, Animalivory, Common.name, Species, diet.general) %>%
                summarise(Completeness = mean(Completeness)) %>%
                # Order element by animalivory
                arrange(desc(Animalivory), Common.name, Domain, Element) %>%
                mutate(Common.name = factor(Common.name, levels = unique(Common.name)),
                       Element = factor(Element, levels = unique(Element))) %>%
                mutate(diet_short = recode(diet.general,
                                            "Animalivore" = "An.",
                                            "Omnivore" = "Omn.",
                                            "Herbivore" = "Herbivore",
                                            "Frugivore" = "Frugivore")) %>%
                mutate(diet_short = factor(diet_short, levels = rev(levels(diet_short)))) %>% ungroup()

p2 <- ggplot(gifts_heatmap, aes(x = Common.name, y = Element, fill = Completeness)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#160583", name = "Functional trait\ncompleteness ", breaks = c(0, 0.5, 1)) +
    facet_grid(rows = vars(Function),
               cols = vars(diet_short), scale = "free", space = "free",
               labeller = as_labeller(function(x) str_wrap(x, width = 25))) + 
    theme(legend.position = c(1.7, 0.4), legend.title = element_text(angle = -90, size = 10), legend.text = element_text(size = 10), legend.title.position = "right",
          strip.background = element_rect(fill = "white", colour = "white"),
          strip.text.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(),
          strip.text.y = element_text(angle = 0, hjust = 0), plot.margin = margin(t = 0, r = 0, b = 0, l = 20))

# Depth and number of genes
depth <- data.frame(phy_gifts_el@sam_data) %>% select(Common.name, diet.general, contig_count, Gene_richness, Total_abundance) %>%
    group_by(Common.name, diet.general) %>%
    summarise(Gene_richness = mean(Gene_richness), Total_abundance = mean(Total_abundance)) %>%
    mutate(Common.name = factor(Common.name, levels = unique(gifts_heatmap$Common.name))) %>%
    left_join(unique(select(gifts_heatmap, c(Common.name, diet_short))))

p3 <- ggplot(depth, aes(x = Common.name, y = Total_abundance, group = diet_short)) +
    geom_line() +
    scale_y_continuous(sec.axis = sec_axis(~./10^5)) +
    geom_line(aes(y = Gene_richness*10^5), colour = "red", linetype = "dashed") +
    facet_grid(cols = vars(diet_short), scale = "free", space = "free_x") + 
    geom_label(data = data.frame(label = c("Total abundance", "Gene richness"),
                                 x = c(10, 10), y = c(mean(depth$Total_abundance), mean(depth$Gene_richness)*10^5),
                                 diet_short = depth$diet_short[1:2]),
              aes(label = label, x = x, y = y), colour = c("black", "red"), alpha = 0.6, inherit.aes = TRUE) +
    theme(legend.position = "bottom", strip.text.x = element_blank(),
          axis.text.x = element_text(size = 10, angle = 45, hjust=1), axis.title = element_blank(),
          strip.text.y = element_text(angle = 0, size = 7), axis.text.y.right = element_text(color = "red"),
          plot.margin = margin(t = 0, r = 10, b = 20, l = 0))

p <- plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 0.9, 1))

ggsave(p, filename = file.path(subdir, "loadings_heatmap.png"), width=10, height=6)

###################
#### PERMANOVA ####
###################

sample_data <- as.data.frame(phy_gene_f_clr@sam_data)

# Explanatory variables
order <- sample_data$Order
diet <- sample_data$diet.general
habitat <- sample_data$habitat.general
ruminant <- (sample_data$digestion == "Ruminant")
hypsodont <- grepl("hyps", sample_data$molar_category)
species <- sample_data$Species

#### GENE DATA ####
set.seed(123)

otu_table <- t(as.data.frame(phy_gene_f_clr@otu_table))

perm <- adonis2(otu_table ~ order + diet + habitat + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_gene_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(otu_table ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_gene_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

#### PATHWAY DATA ####
set.seed(123)

otu_table <- t(as.data.frame(phy_pathway@otu_table))[rownames(sample_data),]

perm <- adonis2(otu_table ~ order + diet + habitat + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_pathway_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(otu_table ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_pathway_onlyspecies.csv"), row.names = TRUE, quote = TRUE)

#### GIFT DATA ####

set.seed(123)

otu_table <- as.data.frame(phy_gifts_el@otu_table)[rownames(sample_data),]

perm <- adonis2(otu_table ~ order + diet + habitat + ruminant + hypsodont,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_gift_allfactors.csv"), row.names = TRUE, quote = TRUE)

perm <- adonis2(otu_table ~ sample_data$Species,
        permutations = 1000, by = "margin", method = "euclidean")

write.csv(as.data.frame(perm), file = file.path(subdir, "permanova_gift_onlyspecies.csv"), row.names = TRUE, quote = TRUE)
