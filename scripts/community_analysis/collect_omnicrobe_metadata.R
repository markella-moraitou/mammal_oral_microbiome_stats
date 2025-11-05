##### Collect metadata from omnicrobe database #####

#### Get taxids based on taxon names from phyloseq object
#### Use them to get info on known habitats, phenotypes and use for all assembled MAGs

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(stringr)
library(taxize)
library(phyloseq)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # subdirectory for the output of this script
phydir <- normalizePath(file.path(subdir, "phyloseq_objects")) # subdirectory for the phyloseq objects

options(ENTREZ_KEY = Sys.getenv("API_KEY"))

##################################
#### GET TAXIDS FOR OMNICROBE ####
##################################

phy_sp <- readRDS(file.path(phydir, "phy_sp.RDS"))

# Keep only taxa present in samples
phy_sp <- subset_samples(phy_sp, is.neg == FALSE)
phy_sp <- prune_taxa(taxa_sums(phy_sp) > 0, phy_sp)

taxnames <- taxa_names(phy_sp)

if (!file.exists(file.path(subdir, "names_to_ids_filt.csv"))) {
  # Modify taxnames to match NCBI database (e.g. remove spxxxx type epithets)
  taxsimple <- taxnames %>% str_remove(" sp[0-9]+") %>% str_remove("\\*") %>% str_remove("_[A-Z]+")
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
  replacement = str_remove(missing, " .*") %>% str_remove("_[A-Z]+")
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

###############################
#### RUN OMNICROBE SCRIPTS ####
###############################

#### Use taxids to get info on known habitats, phenotypes and use for all identified taxa

# Need to run commands through system
system('scriptdir="..";
        outdir="../../output/community_analysis";
        cut -d, -f 2 $outdir/names_to_ids_filt.csv | sort | uniq > $outdir/taxids.tmp;
        sed -i '/NA/d' $outdir/taxids.tmp;
        sed -i '/ids/d' $outdir/taxids.tmp;
        python $scriptdir/access_omnicrobe_db.py $outdir/taxids.tmp $outdir "habitat";
        python $scriptdir/access_omnicrobe_db.py $outdir/taxids.tmp $outdir "phenotype";
        python $scriptdir/access_omnicrobe_db.py $outdir/taxids.tmp $outdir "use";
        rm $outdir/taxids.tmp')

cat("Finished getting taxids and omnicrobe data.\n")

