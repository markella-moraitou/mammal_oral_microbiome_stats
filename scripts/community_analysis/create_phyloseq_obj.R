##### GENERATE PHYLOSEQ #####

#### Use CAT taxonomy, abundance based on contig mapping, and sample metadata
#### to generate a phyloseq object for further analysis

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(renv)
library(phyloseq)
library(stringr)
library(tibble)
library(taxize)
library(microbiome)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # subdirectory for the output of this script
phydir <- normalizePath(file.path(subdir, "phyloseq_objects")) # subdirectory for the phyloseq objects

dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
dir.create(phydir, showWarnings = FALSE)

options(ENTREZ_KEY = Sys.getenv("API_KEY"))

#######################
#####  LOAD INPUT #####
#######################

#### Define file paths
# OTU table
otu_table_path <- file.path(indir, "abundance_table_CAT.tsv")
tax_table_path = file.path(indir, "taxonomy_table_CAT.tsv")

# Sample metadata
metadata_path <- file.path(indir, "sample_metadata.csv")
sample_tax_path <- file.path(indir, "host_taxonomy.csv")
host_traits_path <- file.path(indir, "host_traits.csv")
quant_diet_path <- file.path(indir, "Lintulaakso_diet_filtered.csv")
rc_path <- file.path(indir, "read_count.csv")
rl_path <- file.path(indir, "read_length.csv")
decom_path <- file.path(indir, "decOM_output.csv")

#### Load files
# Load OTU table
otu_table <- read.table(otu_table_path, sep="\t", comment.char="", header=TRUE)
tax_table <- read.table(tax_table_path, sep="\t", comment.char="", header=TRUE)

# Load metadata
metadata <- read.csv(metadata_path) # Sample metadata
sample_tax <- read.csv(sample_tax_path) %>%  # Sample taxonomy
  select(museum.species, genus, subfamily, infraorder, suborder, order, superorder, Common.name)
host_traits <- read.csv(host_traits_path) # Species habitats
quant_diet <- read.csv(quant_diet_path) # Diet quantification data from Lintulaakso et al. 2023 paper
rc <- read.csv(rc_path) %>% rename_with( ~ paste0(., "_count")) # read count per step
rl <- read.csv(rl_path) %>% rename_with( ~ paste0(., "_avlength")) # average read length per step
decom <- read.csv(decom_path, row.names = NULL) %>% # decOM output
  # remove entries where no kmers have been counted
  filter(rowSums(!is.na(select(., starts_with("p_")))) > 0)

###############################################
#### COMBINE AND TIDY UP SAMPLE METADATA  #####
###############################################

# Normalise colnames
colnames(sample_tax) <- str_to_title(colnames(sample_tax))

# Combine all sample and host species metadata in one big table
meta <- metadata

meta <-
  meta %>%
  left_join(sample_tax, relationship = "many-to-many", by=c("Species"="Museum.species")) %>%
  left_join(host_traits, relationship = "many-to-many") %>%
  left_join(quant_diet, relationship="many-to-many") %>%
  left_join(rc, by=c("Ext.ID"="sample_count")) %>%
  left_join(rl, by=c("Ext.ID"="sample_avlength")) %>% 
  left_join(decom, by=c("Ext.ID"="Sink")) %>%
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

# Keep a version with all lab metadata (not just for samples in OTU table)
meta_all <- meta

################################
#### PREPARE TABLES FOR PS #####
################################

#### OTU TABLE ####

# May remove this later: decimals in abundance are throwing errors when calculating richness
# therefore I will round everything down to the next lower integer
# this means anything with abundance less than 1 will appear as absent (fine with that!)

tbl <- otu_table %>% column_to_rownames("lineage") %>%
  mutate(across(everything(), floor))

# Get new names
tbl <- tbl %>%
  rename_with(~rename$new_name[match(., rename$old_name)], everything())

OTU = otu_table(tbl, taxa_are_rows = TRUE)

#### SAMPLE DATA ####
# Keep only samples that exist in tbl
meta <- meta %>% filter(new_name %in% colnames(tbl))

# Get rownames
rownames(meta) <- meta$new_name

SAM = sample_data(meta)

#### TAX TABLE ####

taxonomy <- tax_table %>% select(c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "lineage")) %>%
          filter(lineage %in% taxa_names(OTU)) %>%
          # Remove suffixes e.g. _A from Bacillota_A
          mutate(across(everything(), ~ str_remove_all(., "_[A-Z]+"))) %>%
          as.matrix()

# Get taxa names as row names
row.names(taxonomy) <- taxonomy[, "lineage"]

TAX = tax_table(taxonomy) %>% unique

###################################
####  CREATE PHYLOSEQ TABLES  #####
###################################

# Create phyloseq object and filter out empty samples
cat("Creating phyloseq table\n")
phy <- phyloseq(OTU, SAM, TAX)

# To reduce the size of this table and make it more manageable, remove taxa with less than 500 abundance
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
phy <- subset_samples(phy, sample_sums(phy) > 0)

# Agglomerate to species level and keep only prokaryotes and archaea
cat("Agglomerating to species level\n")
phy_sp <- phy %>% tax_glom(taxrank="species") %>% subset_taxa(species != "no support") %>%
 subset_taxa(superkingdom %in% c("Bacteria", "Archaea"))

# Extra lineage column causes issues with some scripts
phy_sp@tax_table <- phy_sp@tax_table[,c("superkingdom", "phylum", "class", "order", "family", "genus", "species")]

phy_sp <- subset_taxa(phy_sp, taxa_sums(phy_sp) > 0)
phy_sp <- subset_samples(phy_sp, sample_sums(phy_sp) > 0)

# Get taxa names instead of ids
taxa_names(phy_sp) <- make.unique(as.vector(phy_sp@tax_table[,"species"]))

#### Collect some metadata

# Collect number of OTUs per sample
phy_sp@sam_data$taxa_raw <- estimate_richness(phy_sp, measures="Observed")$Observed

#Add column indicating samples and controls
phy_sp@sam_data$is.neg <- grepl("blank|control", phy_sp@sam_data$Order_grouped)

# Calculate oral to soil ratio according to DecOM results
phy_sp@sam_data <- phy_sp@sam_data %>% data.frame %>% mutate(oral_to_soil_ratio=(p_mOral + p_aOral)/p_Sediment.Soil) %>% sample_data

# CLR-normalisation
phy_sp_clr <- phy_sp %>% transform('clr')

#####################
#### SAVE OUTPUT ####
#####################

# Phyloseq objects
saveRDS(phy_sp, file.path(phydir, "phy_sp.RDS"))
saveRDS(phy_sp_clr, file.path(phydir, "phy_sp_clr.RDS"))

# Constituent parts of phyloseq object
write.table(otu_table(phy_sp), file.path(subdir, "phyloseq_objects", "phy_sp_OTU.tsv"), sep = "\t", row.names=TRUE, quote=FALSE)
write.table(tax_table(phy_sp), file.path(subdir, "phyloseq_objects", "phy_sp_TAX.tsv"), sep = "\t", row.names=TRUE, quote=FALSE)
write.table(data.frame(sample_data(phy_sp)), file.path(subdir, "phyloseq_objects", "phy_sp_SAM.tsv"), sep = "\t", row.names=TRUE, quote=TRUE)
write.table(otu_table(phy_sp_clr), file.path(subdir, "phyloseq_objects", "phy_sp_clr_OTU.tsv"), sep = "\t", row.names=TRUE, quote=TRUE)

# Entire taxonomy table
write.table(taxonomy, file.path(subdir, "taxonomy_all.tsv"), sep = "\t", row.names=FALSE, quote=TRUE)

# Entire metadata table
# Get a big metadata table with all lab samples
meta_all <- meta_all %>% left_join(data.frame(phy_sp@sam_data)) %>%
  mutate(is.neg=grepl("blank|control", Species))

write.table(meta_all, file.path(subdir, "metadata_all.tsv"), sep = "\t", row.names=FALSE, quote=TRUE)

##################################
#### GET TAXIDS FOR OMNICROBE ####
##################################

if (!file.exists(file.path(subdir, "names_to_ids_filt.csv"))) {
  taxnames <- taxa_names(phy_sp)
  # Modify taxnames to match NCBI database (e.g. remove spxxxx type epithets)
  taxsimple <- taxnames %>% str_remove(" sp[0-9]+") %>% str_remove("\\*")
  # Search database to get NCBI codes
  taxids <- sapply(unique(taxsimple), function(x) {
    tryCatch({
      i <- which(unique(taxsimple) == x)[1]
      cat("Searching for ", x, " (", i, "/", length(unique(taxsimple)), ")\n", sep="")
      id = get_ids(x, db="ncbi", simplify=FALSE)$ncbi[[1]]
      return(id)
      Sys.sleep(0.5) # To avoid overloading the server
    },
      error=function(e) {
        message("Error:", conditionMessage(e), "\n")
        Sys.sleep(30) # To avoid overloading the server
        return(NA)
        })
    })
  # Collect results in dataframe 
  taxids_df <- data.frame(searchnames = names(taxids), ids = unlist(taxids))
  # Identify those not found and search only with genus name
  missing = taxids_df %>% filter(is.na(ids) & grepl(" ", searchnames)) %>% pull(searchnames)
  taxids_df = taxids_df %>% filter(!searchnames %in% missing)
  replacement = str_remove(missing, " .*")
  # Search again
  taxids_miss = sapply(unique(replacement), function(x) {
    tryCatch({
      i <- which(unique(replacement) == x)[1]
      cat("Searching for ", x, " (", i, "/", length(unique(replacement)), ")\n", sep="")
      id = get_ids(x, db="ncbi", simplify=FALSE)$ncbi[[1]]
      return(id)
      Sys.sleep(0.5) # To avoid overloading the server
    },
      error=function(e) {
        message("Error:", conditionMessage(e), "\n")
        Sys.sleep(30) # To avoid overloading the server
        return(NA)
        })
    })
  # Collect results in dataframe
  taxids_miss_df <- data.frame(searchnames = names(taxids_miss), ids = unlist(taxids_miss))
  # Combine two table and add original taxon name
  taxids_df <- rbind(taxids_df, taxids_miss_df) %>% filter(!is.na(ids)) %>% unique
  names_to_ids <- data.frame(names = taxnames, searchnames = taxsimple) %>%
      left_join(taxids_df) %>% select(names, ids, searchnames)
  names_to_ids_filt <- names_to_ids %>% filter(!is.na(ids))
  # Save names and ids
  write.table(names_to_ids_filt, file.path(subdir, "names_to_ids_filt.csv"), sep=",", row.names=FALSE, quote=FALSE)
} else {
  cat("Not creating names_to_ids.csv because it already exists. Delete or rename it and rerun script to update it.")
}
