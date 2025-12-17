##### DATA FOR FUNCTIONAL ANALYSIS #####

#### Prepares tables and phyloseq objects for analysis of functional annotations

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(cowplot)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "function")) # subdirectory for the output of this script
subdir <- normalizePath(file.path(outdir, "data")) # subdirectory for the output of this script
taxdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Directory with taxonomy analyses

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# DRAM and CAT output (stratified annotations)
annot_str <- read_tsv(file.path(indir, "gene_abundance_stratified.tsv.gz"))

# Info on genes
dram_summary <- read_tsv(file.path(indir, "DRAM_genome_summary_form.tsv"))

# Contaminant taxa list
abund_negs <- read.csv(file.path(taxdir, "assess_contamination", "abundant_in_negs_taxa.csv"))
low_prev <- read.csv(file.path(taxdir, "assess_contamination", "low_prevalence_taxa.csv"))

contam <- unique(c(abund_negs$x, low_prev$x))

# DRAM genome summary form
summary_form <- read.table(file.path(indir, "DRAM_genome_summary_form.tsv"),
                            sep = "\t", header = TRUE, check.names=FALSE, quote = "\"", comment = "")

# Sample metadata
metadata <- read.csv(file.path(indir, "sample_metadata.csv")) # Sample metadata
sample_tax <- read.csv(file.path(indir, "host_taxonomy.csv")) %>%  # Sample taxonomy
  select(species, genus, subfamily, infraorder, suborder, order, superorder, Common.name)
species_traits <- read.csv(file.path(indir, "host_traits.csv")) # Species habitats
quant_diet <- read.csv(file.path(indir, "Lintulaakso_diet_filtered.csv")) # Diet quantification data from Lintulaakso et al. 2023 paper
elton_traits <- read.csv(file.path(indir, "elton_traits.csv")) # Elton traits data
contig_info <- read.csv(file.path(indir, "contig_info_per_sample.txt"), sep ="\t") # Info on contigs per sample

###########################
#### PREPARE METADATA #####
###########################

#### Collect and tidy metadata ####
# Normalise colnames
colnames(sample_tax) <- str_to_title(colnames(sample_tax))

# Combine dietary components of carnivory in Elton traits
colnames(elton_traits) <- colnames(elton_traits) %>%
                      gsub("Diet.", "", .)

elton_traits <- elton_traits %>%
  mutate(Animal = Inv + Vend + Vect + Vfish + Scav + Vunk) %>%
  select(-c(Inv, Vend, Vect, Vfish, Scav, Vunk))

# Combine all sample and host species metadata in one big table
meta <- metadata

meta <-
  meta %>%
  left_join(sample_tax, relationship = "many-to-many", by=c("Species"="Species")) %>%
  left_join(species_traits, relationship = "many-to-many") %>%
  left_join(quant_diet, relationship="many-to-many") %>%
  left_join(elton_traits, by="Species") %>%
  left_join(contig_info, by = c("Ext.ID" = "sample")) %>%
  # Combine different rows for pooled samples
  group_by(Ext.ID) %>%
  summarise_all(.funs = function(.cols) {
    if (n_distinct(.cols) > 1) {
      paste(.cols, collapse = "&")
    } else {
      .cols[1]
    }
  }) %>%
  rename("diet.general"="calculated_cluster_main_diet") %>%
  as.data.frame

# Adjust dietary categories (trying to reconcile Lintulaakso et al. and Elton traits)
meta <- meta %>%
  mutate(diet.general = case_when(
    Species == "Giraffa camelopardalis" ~ "Herbivore",
    Species %in% c("Sus scrofa", "Sus domesticus") ~ "Omnivore",
    Species == "Marmota marmota" ~ "Frugivore",
    TRUE ~ diet.general
  ))

# Turn order and species to a factor and create a separate order for blanks and controls
meta$Order <- ifelse(meta$Species %in% c("Environmental control", "Extraction blank", "Library blank"),
                    "Control/blank", meta$Order)

ord_levels <- unique(sort(meta$Order))
ord_levels <- ord_levels[ord_levels != c("Control/blank")] %>% append("Control/blank")
meta$Order <- factor(meta$Order, levels=ord_levels)

spe_levels <- meta %>% arrange(Order, Species) %>% pull(Species) %>% unique
meta$Species <- factor(meta$Species, levels=spe_levels)

# Create a new Order column with rare orders grouped together for easier plotting
meta <- meta %>% group_by(Order) %>%
  mutate(Order_grouped = case_when(n_distinct(Species) > 1 ~ Order,
                                   TRUE ~ "Rest")) %>% as.data.frame

levels <- unique(meta$Order_grouped)
levels <- levels[which(! levels %in% c("Rest", "Control/blank"))] %>% append(c("Rest", "Control/blank"))

meta$Order_grouped <- factor(meta$Order_grouped, levels = levels)

# Turn diet categories to factors
meta$diet.general <- factor(meta$diet.general, levels=c("Herbivore", "Frugivore", "Omnivore", "Animalivore"))

## Get better samples names
rename <- meta %>%
  # Give common name to blanks and controls
  mutate(Common.name = case_when(Species %in% c("Environmental control", "Extraction blank", "Library blank") ~ Species, TRUE ~ Common.name)) %>%
  select(Ext.ID, Common.name) %>% rename("old_name"="Ext.ID") %>%
  # new names will consist of the first letter of the genus, the first three of the species epithet and a number
  group_by(Common.name) %>% mutate(num=row_number() %>% str_pad(width = 2, pad = "0")) %>%
  separate(col=Common.name, into=c("part1", "part2", "part3", "part4"), fill="left", sep=" ") %>%
  # Turn NAs to ""
  mutate(part1=ifelse(is.na(part1), "", part1),
         part2=ifelse(is.na(part2), "", part2),
         part3=ifelse(is.na(part3), "", part3),
         part4=ifelse(is.na(part4), "", part4)) %>%
  # Use first 4 letters of last adjective (part3) and the entire last word (part4)
  mutate(new_name=paste(str_to_lower(ifelse(str_length(part3)> 3, str_sub(part3, 1, 1), part3)),
                        str_to_title(part4), "_", num, sep="")) %>%
  select(old_name, new_name)

# Temporary bit for this unknown library
rename$new_name[rename$old_name=="maybe_DM_017"] <- paste0(rename$new_name[rename$old_name == "DM_017"], ".maybe")

meta$new_name <- rename$new_name[match(meta$Ext.ID, rename$old_name)]

write.csv(meta, file.path(subdir, "meta.csv"), row.names = FALSE)

########################################
#### GET GENE ABUNDANCES PER SAMPLE ####
########################################

# Remove taxa that were identified as contaminants by community-wide analysis
annot_mod <- annot_str %>%
    filter(!species %in% contam)

annot_mod <- annot_mod %>%
        # Remove _[A-Z] from taxon names
        mutate(phylum = str_remove_all(phylum, "_[A-Z]+")) %>%
        # For CAZY annotations, remove the subcategories and keep only the two top level categories
        # If this leads to two different descriptions in the same category, revert to long name
        mutate(gene_id_temp = case_when(grepl(";.*;", gene_id) ~ str_extract(gene_id, "^([^;]+; [^;]+)"),
                                   TRUE ~ gene_id)) %>%
        mutate(gene_id_temp = case_when(n_distinct(gene_description) > 1 ~ gene_id,
                                   TRUE ~ gene_id_temp)) %>%
        mutate(gene_id = gene_id_temp,
                gene_id_temp = NULL) %>%
        # Also, shorten gene names
        mutate(gene_name = gene_description) %>%
        mutate(gene_name = case_when(database == "CAZY" & grepl("; ", gene_name) ~ str_remove(gene_name, "; .*"),
                                     TRUE ~ gene_name)) %>%
        mutate(gene_name = case_when(database == "CAZY" & str_length(gene_name) > 100 ~ str_remove(gene_name, " .*"),
                                     database == "MEROPS" & str_length(gene_name) > 100 ~ str_remove(gene_name, "#*#.()"),
                                     TRUE ~ gene_name))

# Save the modified annotations
write.table(annot_mod, file.path(subdir, "gene_abundance_stratified_modified.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Sum gene abundances per sample and make wider
annot <- annot_mod %>% group_by(Sample, gene_id) %>%
        summarise(abundance = sum(mapped_reads)) %>%
        pivot_wider(names_from = Sample, values_from = abundance, values_fill = 0) %>%
        column_to_rownames("gene_id")

write.table(annot, file.path(subdir, "gene_abundance.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)

# Summarise the DRAM summary form
dram_info <- dram_summary %>% select(gene_id, module, sheet, header) %>%
        mutate(sheet = str_remove(sheet, " (.*)")) %>%
        # Filter out MISC (not so informative)
        filter(sheet != "MISC") %>%
        # Summarise when a gene is in many sheets
        group_by(gene_id) %>%
        summarise(sheet = paste(unique(sheet), collapse = "; "),
                  module = paste(unique(module), collapse = "; "),
                  header = paste(unique(header), collapse = "; "))

# Get gene info
annot_info <- annot_mod %>%
        select(database, gene_id, gene_name, gene_description) %>% unique %>%
        # Add info from DRAM summary
        left_join(dram_info, by = "gene_id", relationship = "one-to-one") %>%
        select(gene_id, database, sheet, header, module, gene_name, gene_description) %>%
        # Indicate a gene is a peptidase or CAZy enzyme
        rename(category = sheet) %>%
        mutate(category = case_when(database == "MEROPS" ~ "Peptidase",
                                    database == "CAZY" ~ "CAZy enzyme",
                                    TRUE ~ category))

write.table(annot_info, file.path(subdir, "gene_info.csv"),
            sep = "\t", quote = FALSE, row.names = TRUE)

##################################
#### PREPARE PHYLOSEQ OBJECT  ####
##################################

# Fix names
colnames(annot) <-  meta$new_name[match(colnames(annot), meta$Ext.ID)]

tbl <- annot %>% mutate(across(where(is.numeric), floor))
otu <- otu_table(tbl, taxa_are_rows = TRUE)

sam <- sample_data(column_to_rownames(meta, "new_name"))

tax <- tax_table(as.matrix(column_to_rownames(annot_info, "gene_id")))

phy_gene <- phyloseq(otu, sam, tax)

phy_gene@sam_data$is.neg <- grepl("control|blank", phy_gene@sam_data$Species, ignore.case = TRUE)

# Keep only genes with at least 100 reads in total
phy_gene <- prune_taxa(taxa_sums(phy_gene) > 100, phy_gene)

# Count total abundance of features in otu table
phy_gene@sam_data$Total_abundance <- colSums(phy_gene@otu_table)

# Count unique gene richness
phy_gene@sam_data$Gene_richness <- estimate_richness(phy_gene, measure="Observed")[[1]]

#### Identify low content samples
# Identify when number of genes plateaus
contigs_to_genes <- data.frame(phy_gene@sam_data) %>% select(contig_count, len_median, Total_abundance, Gene_richness, is.neg)

thres <- 10^6

p <- ggplot(contigs_to_genes, aes(x = Total_abundance, y = Gene_richness, colour = is.neg)) +
    geom_point() +
    scale_x_log10() +
    geom_vline(xintercept = thres) +
    scale_y_log10() +
    scale_color_manual(values = c('TRUE' = 'red', 'FALSE' = 'black'), name = "Negative control?") +
    ylab("Gene richness (number of unique genes)") + xlab("Total gene abundance (number of classified reads)")

ggsave(p, filename = file.path(subdir, "contigs_to_genes.png"))

# low content samples
low_content_samples <- names(which(colSums(phy_gene@otu_table) < thres))

write.csv(low_content_samples, file.path(subdir, "low_content_samples.txt"), quote = FALSE, row.names = FALSE)

# CLR normalisation
phy_gene_clr <- transform(phy_gene, "clr")

saveRDS(phy_gene, file.path(subdir, "phy_gene.RDS"))
saveRDS(phy_gene_clr, file.path(subdir, "phy_gene_clr.RDS"))

# Remove low content samples and genes that do not pass the filters
phy_gene_f <- prune_samples(!(sample_names(phy_gene) %in% low_content_samples), phy_gene)
phy_gene_f <- subset_samples(phy_gene_f, !is.neg)
phy_gene_f <- prune_taxa(taxa_sums(phy_gene_f) > 0, phy_gene_f)

##############################
#### FILTER BY PREVALENCE ####
##############################

# Prevalence per species
species <- phy_gene_f@sam_data$Species %>% levels
prevalence <- data.frame(row.names = species)

for (spe in species) {
  subset <- phy_gene_f %>% subset_samples(Species == spe)
  temp <- prevalence(subset, detection = 10^-5, include.lowest = TRUE, sort = TRUE) %>% data.frame() %>% t
  rownames(temp) <- spe
  prevalence <- rbind(prevalence, temp)
}

prevalence <- t(prevalence) %>% data.frame %>% arrange(desc(rowSums(.))) %>%
                rownames_to_column("taxon")

write.csv(prevalence, file.path(subdir, "gene_prevalence_per_species.csv"), quote = FALSE)

# Identify taxa that have at least 20% prevalence in a single species
low_prevalence_genes <- prevalence %>%
  rowwise() %>%
  filter(all(c_across(where(is.numeric)) < 0.2)) %>%
  pull(taxon)

write.csv(low_prevalence_genes, file.path(subdir, "low_prevalence_genes.txt"), quote = FALSE, row.names = FALSE)

# Remove low prevalence taxa
phy_gene_f <- prune_taxa(setdiff(taxa_names(phy_gene_f), low_prevalence_genes), phy_gene_f)

# Make sure there are no empty samples
phy_gene_f <- prune_samples(sample_sums(phy_gene_f) > 0, phy_gene_f)

# CLR normalisation
phy_gene_f_clr <- transform(phy_gene_f, "clr")

saveRDS(phy_gene_f, file.path(subdir, "phy_gene_f.RDS"))
saveRDS(phy_gene_f_clr, file.path(subdir, "phy_gene_f_clr.RDS"))

dbinfo <-
    data.frame(database = phy_gene_f@tax_table[, "database"],
              gene_abundance = taxa_sums(phy_gene_f)) %>% 
    group_by(database) %>%
    summarise(n_genes = n(),
              p_genes = n()/nrow(phy_gene_f@tax_table),
              total_abundance = sum(gene_abundance),
              p_abundance = total_abundance/sum(taxa_sums(phy_gene_f)))

write.csv(dbinfo, file.path(subdir, "database_info.csv"), row.names = FALSE)
