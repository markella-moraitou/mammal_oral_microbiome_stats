##### DATA FOR FUNCTIONAL ANALYSIS #####

#### Prepares tables and phyloseq objects for analysis of functional annotations

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(tibble)
library(phyloseq)
library(microbiome)
library(ggplot2)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "function")) # subdirectory for the output of this script

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# DRAM distillation and product outputs

distill <- read_excel(file.path(indir, "sample_metabolism_summary.xlsx")) %>%
            unique

product <- read.table(file.path(indir, "sample_product.tsv"),
                        sep = "\t", header = TRUE, check.names=FALSE)

# Sample metadata
metadata <- read.csv(file.path(indir, "sample_metadata.csv")) # Sample metadata
sample_tax <- read.csv(file.path(indir, "host_taxonomy.csv")) %>%  # Sample taxonomy
  select(museum.species, genus, subfamily, infraorder, suborder, order, superorder, Common.name)
species_traits <- read.csv(file.path(indir, "species_traits.csv")) # Species habitats
quant_diet <- read.csv(file.path(indir, "Lintulaakso_diet_filtered.csv")) # Diet quantification data from Lintulaakso et al. 2023 paper
contig_info <- read.csv(file.path(indir, "contig_info_per_sample.txt"), sep ="\t") # Info on contigs per sample

###########################
#### PREPARE METADATA #####
###########################

#### Collect and tidy metadata ####
# Normalise colnames
colnames(sample_tax) <- str_to_title(colnames(sample_tax))

# Combine all sample and host species metadata in one big table
meta <- metadata

meta <-
  meta %>%
  left_join(sample_tax, relationship = "many-to-many", by=c("Species"="Museum.species")) %>%
  left_join(species_traits, relationship = "many-to-many") %>%
  left_join(quant_diet, relationship="many-to-many") %>%
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

# Turn order and species to a factor and create a separate order for blanks and controls
meta$Order <- ifelse(meta$Species %in% c("Environmental control", "Extraction blank", "Library blank"),
                    "Control/blank", meta$Order)

ord_levels <- unique(sort(meta$Order))
ord_levels <- ord_levels[ord_levels != c("Control/blank")] %>% append("Control/blank")
meta$Order <- factor(meta$Order, levels=ord_levels)

spe_levels <- meta %>% arrange(Order, Species) %>% pull(Species) %>% unique
meta$Species <- factor(meta$Species, levels=spe_levels)

comname_levels <- meta %>% arrange(Order, Common.name) %>% pull(Common.name) %>% unique
meta$Common.name <- factor(meta$Common.name, levels=comname_levels)

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

write.csv(meta, file.path(outdir, "meta.csv"), row.names = FALSE)

#######################################
#### PREP DRAM OUTPUT FOR PLOTTING ####
#######################################

#### Distill table ####

# Make longer for plotting
distill_long <- distill %>%
                pivot_longer(cols = contains("final_contigs"), names_to = "Sample", values_to = "Count") %>%
                mutate(Sample = str_remove(Sample, "_final_contigs")) %>%
                mutate(Count = as.numeric(Count)) %>%
                # Get percentages
                group_by(Sample) %>%
                mutate(Percentage = Count * 100 / sum(Count)) %>%
                filter(sum(Count) > 0)

# Add sample metadata
distill_long <- select(meta, c(Ext.ID, new_name, Common.name, Order, diet.general)) %>%
                right_join(distill_long, by=c("Ext.ID"="Sample"))

write.csv(distill_long, file.path(outdir, "distill_long.csv"), quote = FALSE, row.names = FALSE)

#### Product table ####

# Fix names
product$genome <- product$genome %>% str_remove("_final_contigs")

# Make sure logical columns are read as logical
product <- product %>% mutate(across(where(~ all(. %in% c("True", "False"))), ~ as.logical(str_to_upper(.)))) 

### TEMPORARY?
product <- product[-which(product$genome == "fasta"),]

rownames(product) <- NULL

# Make longer
pathway_completeness <- product %>% select(genome, where(is.numeric)) %>%
  pivot_longer(cols = -c(genome), names_to = "feature", values_to = "completeness") %>%
  # Separate columns into categorie and feature
  separate(col=feature, into=c("category", "feature"), sep=": ", fill="left") %>%
  # In the cases where there is no overall category (these pathways that are grouped under module in the product HTML)
  # just use the feature name, simplified
  mutate(category = case_when(is.na(category) ~ str_remove(feature, ", .*") %>% str_remove(" (.*)"),
                              TRUE ~ category))

process_presence <- product %>% select(genome, where(is.logical)) %>%
  pivot_longer(cols = -c(genome), names_to = "feature", values_to = "presence") %>%
  # Separate columns into categorie and feature
  separate(col=feature, into=c("category", "feature"), sep=": ", fill="left")

#### Add metadata ####
pathway_completeness <- pathway_completeness %>%
  left_join(select(meta, c(Ext.ID, new_name, Common.name, Order, diet.general)), by=c("genome"="Ext.ID")) %>%
  mutate(sample_category = paste(Order, diet.general, sep="\n"))

write.csv(pathway_completeness, file.path(outdir, "pathway_completeness.csv"), row.names = FALSE, quote = FALSE)

process_presence <- process_presence %>%
  left_join(select(meta, c(Ext.ID, new_name, Common.name, Order, diet.general)), by=c("genome"="Ext.ID")) %>%
  mutate(sample_category = paste(Order, diet.general, sep="\n"))

write.csv(process_presence, file.path(outdir, "process_presence.csv"), row.names = FALSE, quote = FALSE)

##################################
#### PREPARE PHYLOSEQ OBJECTS ####
##################################

#### GENE PHYLOSEQ ####
## Based on distillation

# Get OTU table
gene_table <- distill %>% select(gene_id, contains("final_contigs")) %>% unique %>%
  # Some genes appear in more than one line (due to different names, modules) -- get average (the abundance should be the same in every replicate line)
  group_by(gene_id) %>% summarise(across(everything(), ~ mean(as.numeric(.)))) %>%
  column_to_rownames(var="gene_id") %>%
  # make sure all columns are numeric
  mutate(across(everything(), ~ as.numeric(.)))

# Fix names
colnames(gene_table) <-  meta$new_name[match(str_remove(colnames(gene_table), "_final_contigs"), meta$Ext.ID)]

# Get get gene taxonomy e.g. header, subheader, module
gene_tax <- distill %>% select(header, subheader, module, gene_description, gene_id) %>% unique %>%
        # When a gene id shows up in more than one module, collapse together module names
        group_by(gene_id) %>%
        summarise_all(.funs = function(.cols) {if (n_distinct(.cols) > 1) {paste(.cols, collapse = " & ")} else {.cols[1]}}) %>%
        column_to_rownames("gene_id") %>% as.matrix

otu <- otu_table(gene_table, taxa_are_rows = TRUE)

sam <- sample_data(column_to_rownames(meta, "new_name"))

tax <- tax_table(gene_tax)

phy_gene <- phyloseq(otu, sam, tax)

# Count total abundance of features in otu table
phy_gene@sam_data$Total_abundance <- colSums(phy_gene@otu_table)

# Count unique gene richness
phy_gene@sam_data$Gene_richness <- estimate_richness(phy_gene, measure="Observed")[[1]]

#### Identify low content samples
# Identify when number of genes plateaus
contigs_to_genes <- data.frame(phy_gene@sam_data) %>% select(contig_count, len_median, Total_abundance, Gene_richness)

p <- ggplot(contigs_to_genes, aes(x = Total_abundance, y = Gene_richness, colour = contig_count)) +
    geom_point() +
    scale_color_viridis_c()

ggsave(p, filename = file.path(outdir, "contigs_to_genes.png"))

# low content samples
low_content_samples <- names(which(colSums(phy_gene@otu_table) < 1000))

write.csv(low_content_samples, file.path(outdir, "low_content_samples.txt"), quote = FALSE, row.names = FALSE)

# CLR normalisation
phy_gene_clr <- transform(phy_gene, "clr")

saveRDS(phy_gene, file.path(outdir, "phy_gene.RDS"))
saveRDS(phy_gene_clr, file.path(outdir, "phy_gene_clr.RDS"))

#### FUNCTION PHYLOSEQ ####
## Based on product, combining both coverage and presence/absence info as presence absence

cov_table <- product %>% select(genome, where(is.numeric))
cov_table$new_name <- meta$new_name[match(cov_table$genome, meta$Ext.ID)]
cov_table <- cov_table %>% select(-genome) %>%
  column_to_rownames(var="new_name")

# Fix colnames
colnames(cov_table) <- str_remove(colnames(cov_table), ".*: ")

# Turn presence absence table to numeric (1 if TRUE, 0 if FALSE)
pres_table <- product %>% select(genome, where(is.logical))
pres_table$new_name <- meta$new_name[match(pres_table$genome, meta$Ext.ID)]
pres_table <- pres_table %>% select(-genome) %>%
  column_to_rownames(var="new_name") %>%
  # make sure all columns are numeric
  mutate(across(everything(), ~ as.numeric(.)))

# Fix colnames
colnames(pres_table) <- str_remove(colnames(pres_table), ".*: ")

# Combine
function_table <- cbind(cov_table, pres_table)

# Remove duplicate columns
function_table <- function_table[,!duplicated(colnames(function_table))]

# Get groupings by category
func_group <- rbind(select(pathway_completeness, c(category, feature)),
                    select(process_presence, c(category, feature))) %>%
                    # When a feature shows up in more than one category, collapse together category names
                    group_by(feature) %>%
                    summarise(category = case_when(n_distinct(category) > 1 ~ paste(unique(category), collapse = " & "),
                                                       TRUE ~ category[1])) %>%
                    column_to_rownames("feature") %>% as.matrix

otu <- otu_table(function_table, taxa_are_rows = FALSE)

tax <- tax_table(func_group)

phy_function <- phyloseq(otu, sam, tax)

saveRDS(phy_function, file.path(outdir, "phy_function.RDS"))

