##### INPUT FOR FUNCTIONAL ANALYSIS #####

#### Prepares tables and phyloseq objects for analysis of functional annotations

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
#library(microViz)
#library(cowplot)
#library(microbiome)
#library(rphylopic)
#library(vegan)
#library(stringr)
#library(ggplot2)
#library(RColorBrewer)
library(distillR)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "function")) # subdirectory for the output of this script

dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

## Set up for plotting
#source(file.path("..", "plot_setup.R"))
#plot_setup(file.path("..", "..", "input", "palettes"))
#theme_set(custom_theme())

# Get ordination functions
#source(file.path("..","ordination_functions.R"))

#######################
#####  LOAD INPUT #####
#######################

raw <- read.table(file.path(indir, "sample_annotations.tsv"),
                  sep = "\t", header = TRUE, check.names=FALSE, quote = "", comment = "", fill = TRUE)


##################
#### DISTILLR ####
##################

# Add gene identifier

if(file.exists(file.path(subdir, "GIFTs.csv"))) {
  cat("Loading GIFTs from file\n")
  GIFTs <- read.csv(file.path(subdir, "GIFTs.csv"), row.names = 1)
} else {
  cat("Running distillR to get GIFTs\n")
  GIFTs <- distill(raw, GIFT_db, genomecol=2, annotcol=c(9, 10, 20, 21, 22))
  write.csv(GIFTs, file.path(subdir, "GIFTs.csv"), row.names = TRUE)
}

GIFTs <- GIFTs[-which(rownames(GIFTs)=="fasta"),]
#Aggregate bundle-level GIFTs into the compound level
GIFTs_elements <- to.elements(GIFTs, GIFT_db)

GIFTs_elements_long <- 
  GIFTs_elements %>% data.frame %>% rownames_to_column("Sample") %>%
  pivot_longer(cols = -c(Sample), names_to = "Code_element", values_to = "Completeness") %>%
  mutate(Sample = str_remove(Sample, "_final_contigs")) %>%
  left_join(select(meta, c(Ext.ID, new_name, Common.name, Order, diet.general)), by=c("Sample"="Ext.ID")) %>%
  left_join(unique(select(GIFT_db, c("Code_element", "Element", "Function", "Domain")))) %>%
  # Turn function into factor
  arrange(Domain, Function) %>% mutate(Function = factor(Function, levels = unique(Function)))

# Remove samples for which nucleotide biosynthesis related paths are under 90% on average complete for any
incomplete_communities <- GIFTs_elements_long %>% filter(Function == "Nucleic acid biosynthesis") %>%
                  group_by(Sample) %>% summarise(mean_completeness = mean(Completeness)) %>% filter(mean_completeness < 0.9) %>%
                  pull(Sample)

GIFTs_elements_long <- GIFTs_elements_long %>% filter(!Sample %in% incomplete_communities)

p <- ggplot(GIFTs_elements_long, aes(y=new_name, x=Element)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Compound", fill="Completeness") +
  facet_grid(rows = vars(Common.name), cols = vars(Function), scales = "free", space = "free")
  
ggsave(p, filename = file.path(subdir, "GIFTs_elements.png"), width = 25, height = 20)

#Aggregate element-level GIFTs into the function level
GIFTs_functions <- to.functions(GIFTs_elements, GIFT_db)

GIFTs_functions_long <- 
  GIFTs_functions %>% data.frame %>% rownames_to_column("Sample") %>%
  pivot_longer(cols = -c(Sample), names_to = "Code_function", values_to = "Completeness") %>%
  mutate(Sample = str_remove(Sample, "_final_contigs")) %>%
  left_join(select(meta, c(Ext.ID, new_name, Common.name, Order, diet.general)), by=c("Sample"="Ext.ID")) %>%
  left_join(unique(select(GIFT_db, c("Code_function", "Function", "Domain")))) %>%
  # Remove samples in incomplete communities
  filter(!Sample %in% incomplete_communities)

p <- ggplot(GIFTs_functions_long, aes(y=new_name, x=Function)) +
  geom_tile(aes(fill=Completeness)) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Function", fill="Completeness") +
  facet_grid(rows = vars(Common.name), cols = vars(Domain), scales = "free", space = "free")

ggsave(p, filename = file.path(subdir, "GIFTs_functions.png"), width = 10, height = 15)

#Aggregate function-level GIFTs into overall Biosynthesis, Degradation and Structural GIFTs
GIFTs_domains <- to.domains(GIFTs_functions, GIFT_db)

#### GIFT - DIET ASSOCIATIONS ####

nutr_degradations <- GIFTs_elements_long %>% filter(Domain == "Degradation")

p <- ggplot(nutr_degradations, aes(y = Completeness, x = diet.general)) +
      geom_violin() +
      geom_jitter(width = 0.2, height = 0.05, alpha = 0.5, size = 0.5) +
      facet_wrap(~Element)

ggsave(p, filename = file.path(subdir, "nutr_degradations.png"), width = 20, height = 20)

nutr_biosynthesis <- GIFTs_elements_long %>% filter(Domain == "Biosynthesis")

p <- ggplot(nutr_biosynthesis, aes(y = Completeness, x = diet.general)) +
      geom_violin() +
      geom_jitter(width = 0.2, height = 0.05, alpha = 0.5, size = 0.5) +
      facet_wrap(~Element)

ggsave(p, filename = file.path(subdir, "nutr_biosynthesis.png"), width = 20, height = 20)


###############
#### PLOT  ####
###############

#### Gene abundance per header ####
p <- ggplot(filter(distill_long), aes(x=Percentage, y=new_name, fill=header)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set1"))(11)) +
  facet_grid(rows = vars(Common.name), scales = "free_y", space = "free_y") +
  theme(strip.text.y = element_text(angle = 0),
        legend.position = "bottom", legend.direction = "vertical") +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE))

ggsave(p, filename = file.path(subdir, "gene_count.png"), width = 8, height = 8)

#### Pathway completeness ####
p <- ggplot(pathway_completeness, aes(y=new_name, x=feature)) +
  geom_tile(aes(fill=completeness)) +
  scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey50") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Feature", fill="Completeness") +
  facet_grid(rows = vars(Common.name), cols = vars(category), scales = "free", space = "free")

ggsave(p, filename = file.path(subdir, "pathway_completeness.png"), width = 20, height = 20)

#### Product presence ####
p <- ggplot(product_presence, aes(y=new_name, x=feature)) +
  geom_tile(aes(fill=presence)) +
  scale_fill_manual(values = c(`TRUE` = "green", `FALSE` = "grey90")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.y = element_text(angle = 0)) +
  labs(y="Sample", x="Feature", fill="Presence") +
  facet_grid(rows = vars(Common.name), cols = vars(category), scales = "free", space = "free")

ggsave(p, filename = file.path(subdir, "product_presence.png"), width = 20, height = 20)