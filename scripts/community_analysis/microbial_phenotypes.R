##### TEST HYPOTHESIS OF MICROBIOME EVOLUTION #####
## This script goes beyond exploratory plots to 
## explore our hypotheses regarding the evolution of the microbiome
## in the context of the host phylogeny and ecology

################
#### SET UP ####
################

#### LOAD PACKAGES ####
library(dplyr)
library(tidyr)
library(tibble)
library(reshape2)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(vegan)
library(RColorBrewer)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis")) # Directory with output from the community analysis
subdir <- normalizePath(file.path("..", "..", "output", "community_analysis", "microbial_phenotypes")) # subdirectory for the output of this script
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

# Load all phyloseq objects in phydir
for (phy_file in list.files(phydir, pattern = "*.RDS")) {
  assign(gsub(".RDS", "", phy_file), readRDS(file.path(phydir, phy_file)))
}

# Microbial taxa known uses and habitats
taxa_uses <- read.csv(file.path(outdir, "use_relations.csv"), header=TRUE)

taxa_habitats <- read.csv(file.path(outdir, "habitat_relations.csv"), header=TRUE)

###################################################
#### MICROBIAL PHENOTYPES ASSOCIATED WITH DIET ####
###################################################

#### First compare different types of metabolic activity with diet

# Uses table
taxa_uses_wide <- taxa_uses %>%
          # Keep only taxa in filtered phyloseq
          filter(taxon %in% taxa_names(phy_sp_f)) %>%
          mutate(OBT = gsub(OBT, pattern = " ", replacement = "_")) %>%
          pivot_wider(id_cols = taxon, names_from = OBT, values_from = occurences) %>%
          rename("fermentative_activity" = "fermentative") %>%
          # Turn all values >= to 3 to 1 and all NAs to 0
          # this means a taxon counts as having an activity if it is reported in at least 3 reports for that taxon
          mutate(across(ends_with("activity"), ~ifelse(is.na(.), 0, ifelse(. >= 3, 1, 0)))) %>% 
          mutate(across(ends_with("activity"), as.numeric)) %>%
          # Remove taxa with all 0s
          filter(if_any(ends_with("activity"), ~. > 0))

# Combine phyloseq info with taxa uses
taxa_uses_abund <- psmelt(transform(phy_sp_f, "compositional")) %>% select(Sample, OTU, Abundance, Species, Common.name, Order, diet.general, Fruit, PlantO, Animal, Seed) %>%
                      right_join(taxa_uses_wide, by=c("OTU"="taxon")) %>%
                      # Turn NA to 0 in all numeric columns
                      mutate(across(ends_with("activity"), ~ifelse(is.na(.), 0, .)))

# Get mean abundances per taxon and save table
taxa_uses_mean_abund <- taxa_uses_abund %>%
                        group_by(OTU, fermentative_activity, antimicrobial_activity, proteolytic_activity, lipolytic_activity) %>%
                        summarise(mean_abundance = mean(Abundance))

write.csv(taxa_uses_mean_abund, file.path(subdir, "taxa_uses.csv"), row.names = FALSE, quote = FALSE)

## Does the abundance of taxa with unknown use vary by species?
known_uses_abund <- psmelt(transform(phy_sp_f, "compositional")) %>% select(Sample, OTU, Abundance, Order, Order_grouped, Common.name, diet.general, Fruit, PlantO, Animal, Seed) %>%
                    filter(OTU %in% taxa_uses_wide$taxon) %>% group_by(Sample, Order, Order_grouped, Common.name, diet.general) %>% summarise(known_abund = sum(Abundance)) %>%
                    # Shorten Order names for plotting
                    mutate(Order_grouped = recode(Order_grouped, "Perissodactyla" = "Peris.", "Rodentia" = "R.", "Carnivora" = "Carn."))

# Get average known-use abundance across samples
mean_abund <- mean(known_uses_abund$known_abund)

p <- ggplot(known_uses_abund, aes(x = known_abund, y = Common.name, fill = diet.general)) +
  facet_grid(rows=vars(Order_grouped), space = "free", scales="free") +
  geom_violin(aes(colour = diet.general)) +
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.5) +
  geom_vline(xintercept = mean_abund, linetype = "dotted") +
  scale_fill_manual(values = diet_palette) +
  scale_colour_manual(values = diet_palette) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x="Abundance of taxa\nwith known functions") +
  theme(legend.position="none", axis.title.y = element_blank(), axis.text.x = element_text(angle = 0))

ggsave(file.path(subdir, "known_function_abundance.png"), p, width=5, height=7)

## Does the abundance and diversity of proteolytic microorganisms correlate with the amount of crude protein in the diet?
proteolytic_abund <- taxa_uses_abund %>% filter(proteolytic_activity > 0) %>% group_by(Sample, Species, Order, Animal) %>% summarise(proteolytic_abund = sum(Abundance))

p <- ggplot(proteolytic_abund, aes(y = proteolytic_abund*100, x = Animal, colour = Order)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", inherit.aes = FALSE, aes(y = proteolytic_abund*100, x = Animal), colour = "black", linetype = "dotted") +
  scale_color_manual(values = order_palette) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(x="Crude protein in diet (%)", y="Porphyromonas\nabundance (%)") +
  theme(legend.position = "none")

ggsave(file.path(subdir, "proteolytic_abundance_animal.png"), p, width=5, height=2)

## Does the abundance and diversity of fermentative microorganisms correlate with the amount of carbohydrates in the diet?
fermentative_abund <- taxa_uses_abund %>% filter(fermentative_activity > 0) %>% group_by(Sample, Species, Order, Fruit) %>% summarise(fermentative_abund = sum(Abundance))

p <- ggplot(fermentative_abund, aes(y = fermentative_abund*100, x = Fruit, colour = Order)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", inherit.aes = FALSE, aes(y = fermentative_abund*100, x = Fruit), colour = "black", linetype = "dotted") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = order_palette, name = "") +
    labs(x="Nitrogen free extract in diet (%)", y="fermentative\nabundance (%)") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))

ggsave(file.path(subdir, "fermentative_abundance_Fruit.png"), p, width=5, height=3)

## Control: antimicrobial activity and diet
antimicrobial_abund <- taxa_uses_abund %>% filter(antimicrobial_activity > 0) %>% group_by(Sample, Species, Order, Animal, Fruit, Seed) %>% summarise(antimicrobial_abund = sum(Abundance))

p <- ggplot(antimicrobial_abund, aes(y = antimicrobial_abund*100, x = Animal, colour = Order)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", inherit.aes = FALSE, aes(y = antimicrobial_abund*100, x = Animal), colour = "black", linetype = "dotted") +
  scale_color_manual(values = order_palette) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(x="Crude protein in diet (%)", y="antimicrobial\nabundance(%)") +
  theme(legend.position = "none")

ggsave(file.path(subdir, "antimicrobial_abundance_Animal.png"), p, width=5, height=2)

# Get overall plot of activity proportions per sample
ntaxa_uses <-
    taxa_uses_abund %>% group_by(Sample, Common.name, Order, diet.general) %>%
      summarise(p_antimicrobial = sum(antimicrobial_activity > 0 & Abundance > 0)/sum(Abundance > 0),
                a_antimicrobial = sum(Abundance[antimicrobial_activity > 0]),
                p_fermentative = sum(fermentative_activity > 0 & Abundance > 0)/sum(Abundance > 0),
                a_fermentative = sum(Abundance[fermentative_activity > 0]),
                p_proteolytic = sum(proteolytic_activity > 0 & Abundance > 0)/sum(Abundance > 0),
                a_proteolytic = sum(Abundance[proteolytic_activity > 0]),
                p_lipolytic = sum(lipolytic_activity > 0 & Abundance > 0)/sum(Abundance > 0),
                a_lipolytic = sum(Abundance[lipolytic_activity > 0])) %>% 
                pivot_longer(cols = starts_with("p_") | starts_with("a_"), names_to = "activity", values_to = "value") %>%
                separate(activity, into = c("measure", "activity"), sep = "_", remove = FALSE) %>%
                mutate(measure = ifelse(measure == "p", "Proportion of taxa", "Proportion of abundance"))

# Summarise by species
ntaxa_uses_species <- 
    ntaxa_uses %>% group_by(Common.name, activity, measure, diet.general, Order) %>%
      summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE)) %>%
      ungroup() %>%
      # Shorten "Omnivore" for plotting
      mutate(diet.general = recode(diet.general, "Omnivore" = "Omn.", "Animalivore" = "Anim."))

p <- ggplot(ntaxa_uses_species, aes(x = Common.name, y = mean, fill = activity)) +
        geom_bar(stat = "identity", position = "dodge")  +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), colour = "black", position = position_dodge(width = 0.9), width = 0.25) +
        scale_fill_manual(values = c("antimicrobial" = "gray50",
                                     "fermentative" = "#3A459C",
                                     "proteolytic" = "#B41F34",
                                     "lipolytic" = "#BD9120")) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        facet_grid(cols = vars(diet.general),
                   rows = vars(measure),
                   scales = "free", space = "free_x") +
        theme(legend.position = "bottom", axis.title = element_blank(),
              axis.text.x = element_text(hjust = 1)) 

ggsave(file.path(subdir, "activity_proportions_per_species.png"), p, width=12, height=8)

#### Second, compare known taxon habitats with host habitat

# Uses table
taxa_habitats_wide <- taxa_habitats %>% mutate(OBT = gsub(OBT, pattern = " ", replacement = "_")) %>%
          # Keep only taxa in filtered phyloseq
          filter(taxon %in% taxa_names(phy_sp_f)) %>%
          # Combine marine_water and deep_sea as sea_water, and dental_plaque and mouth as oral
          mutate(OBT = case_when(OBT == "marine_water" | OBT == "deep_sea" ~ "sea_water",
                                 OBT == "dental_plaque" | OBT == "mouth" ~ "oral",
                                 TRUE ~ OBT)) %>%
          group_by(taxon, OBT) %>% summarise(occurences = sum(occurences)) %>%
          pivot_wider(id_cols = taxon, names_from = OBT, values_from = occurences) %>%
          # Turn all values >= 3 to 1 and all NAs to 0
          # this means a taxon counts as coming from a habitat if it is reported in at least 3 studies
          mutate(across(everything(), ~ ifelse(is.na(.), 0, ifelse(. >= 3, 1, 0)))) %>%
          # Remove taxa with all 0s
          filter(if_any(where(is.numeric), ~. > 0))

# Combine phyloseq info with taxa uses
taxa_habitats_abund <- psmelt(transform(phy_sp_f, "compositional")) %>% select(Sample, OTU, Abundance, Species, Common.name, Order, digestion, diet.general, habitat.general) %>%
                      right_join(taxa_habitats_wide, by=c("OTU"="taxon"))

# Get mean abundances per taxon and save table
taxa_habitats_mean_abund <- taxa_habitats_abund %>%
                        group_by(OTU) %>%
                        mutate(mean_abundance = mean(Abundance)) %>%
                        select(OTU, where(is.numeric), -Abundance) %>% unique

write.csv(taxa_habitats_mean_abund, file.path(subdir, "taxa_habitats.csv"), row.names = FALSE, quote = FALSE)

## Do marine hosts have more marine taxa and terrestrial hosts more soil taxa?
habitat_abund <- taxa_habitats_abund %>% filter((sea_water > 0 | soil > 0) & Species %in% unique(phy_habitat@sam_data$Species)) %>%
        group_by(Sample, Species, Common.name, Order, habitat.general) %>%
        summarise(marine_abund = sum(Abundance[sea_water > 0 & soil == 0]),
                  soil_abund = sum(Abundance[soil > 0 & sea_water == 0])) %>%
        pivot_longer(cols = ends_with("_abund"), names_to = "taxon_habitat", values_to = "abund") %>%
        mutate(taxon_habitat = gsub("_abund", "-only taxa", taxon_habitat)) %>%                       
        # Combine Probiscidea and Sirenia for plotting
        mutate(Order = recode(Order, "Proboscidea" = "Probosc./Sir.", "Sirenia" = "Probosc./Sir."))

p <- ggplot(habitat_abund, aes(y = abund, x = Common.name, fill = habitat.general)) +
  geom_boxplot(alpha = 0.8) +
  facet_grid(cols=vars(Order), rows=vars(taxon_habitat), scales="free", space="free_x") +
  scale_fill_manual(values = habitat_palette, name = "Habitat") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x="", y="Total rel. abundance of taxa") +
  theme(axis.text.x = element_text(hjust = 1))

ggsave(file.path(subdir, "habitat_abundance.png"), p, width=8, height=8)

## Do ruminants have more rumen taxa?
rumen_abund <- taxa_habitats_abund %>% filter(rumen > 0) %>%
                # Get genera identified in rumen
                mutate(genus = gsub(" .*", "", OTU)) %>% group_by(Sample, Common.name, Order, digestion, diet.general, genus) %>%
                summarise(Abundance = sum(Abundance)) %>% group_by(Common.name, Order, digestion, diet.general, genus) %>%
                summarise(Abundance = mean(Abundance)) %>%
                mutate(digestion = case_when(digestion == "" ~ "Non-ruminant", TRUE ~ digestion)) %>%
                # If Artiodactyla, group by digestion, otherwise group by diet
                mutate(group = case_when(Order == "Artiodactyla" ~ digestion,
                                         TRUE ~ diet.general)) %>%
                mutate(group = factor(group, levels = c("Ruminant", "Pseudoruminant", "Non-ruminant", "Herbivore", "Frugivore", "Omnivore", "Animalivore"))) %>%
                # shorten names for plotting
                mutate(group = recode(group, "Pseudoruminant" = "Ps.", "Non-ruminant" = "Non.", "Omnivore" = "Om.", "Animalivore" = "An.")) %>%
                select(Common.name, genus, Abundance, group)
 
# Group top 5 genera, group the rest as "Other"
top_genera <- rumen_abund %>% group_by(genus) %>% summarise(total = sum(Abundance, na.rm = TRUE)) %>% arrange(desc(total)) %>% top_n(5, total) %>% pull(genus)               

rumen_abund <- rumen_abund %>% mutate(genus = ifelse(genus %in% top_genera, genus, "Other")) %>%
                mutate(genus = factor(genus, levels = c("Other", rev(top_genera))))

genera_palette <- setNames(brewer.pal(length(top_genera), "Set3"), top_genera)
genera_palette["Other"] <- "grey70"

p <- ggplot(rumen_abund, aes(x = Abundance*100, y = Common.name, fill = genus)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = genera_palette, name = "") +
      facet_grid(rows = vars(group), scales = "free", space = "free") +
      theme(legend.position = "bottom") +
      labs(y="", x="Rel. abundance of\nrumen microorganisms (%)") +
      # Reverse guide order
      guides(fill = guide_legend(reverse = TRUE, nrow = 3))

ggsave(file.path(subdir, "rumen_abundance.png"), p, width=5, height=8)

#######################################
#### COLONISER STAGE - HYPSODONTY  ####
#######################################

#### Early colonisers: Streptococcus & Actinomyces ####

early_colonisers <- phy_sp_f %>% tax_glom("genus") %>% transform("compositional") %>% psmelt() %>%
                  filter(genus %in% c("Streptococcus", "Actinomyces")) %>%
                  group_by(Sample, Species, genus, molar_category, hypsodonty_index, Order) %>%
                  summarise(early_abundance = sum(Abundance))

# Plot by hypsodonty index (for ungulates only)
p <- ggplot(filter(early_colonisers, !is.na(hypsodonty_index)), aes(y = early_abundance*100, x = hypsodonty_index, group = Species, colour = Order)) +
  geom_boxplot() +
  geom_smooth(method = "lm", inherit.aes = FALSE, aes(y = early_abundance*100, x = hypsodonty_index), colour = "black", linetype = "dotted") +
  facet_grid(rows=vars(genus), scales = "free") +
  #scale_x_continuous(trans = "log10") +
  #scale_y_continuous(trans = "log10") +
  scale_color_manual(values = order_palette, name = "") +
  labs(x="Hypsodonty index", y="Early coloniser abundance (%)") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))

ggsave(file.path(subdir, "early_abund_by_hypsodonty.png"), p, width=5, height=5)

# Plot by molar category
p <- ggplot(filter(early_colonisers, molar_category!=""), aes(y = early_abundance*100, x = Species, fill = Order)) +
  geom_boxplot() +
  facet_grid(cols=vars(molar_category), rows=vars(genus), scales = "free", space = "free_x") +
  #scale_x_continuous(trans = "log10") +
  #scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = order_palette, name = "") +
  labs(x="Molar category", y="Early coloniser abundance (%)") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))

ggsave(file.path(subdir, "early_abund_by_molar_cat.png"), p, width=5, height=8)

#### Late colonisers: Treponema, Porphyromonas, Fusobacterium, Corynebacterium ####

late_colonisers <- phy_sp_f %>% tax_glom("genus") %>% transform("compositional") %>% psmelt() %>%
                  filter(genus %in% c("Prophyromonas", "Treponema", "Fusobacterium", "Corynebacterium")) %>%
                  group_by(Sample, Species, genus, molar_category, hypsodonty_index, Order) %>%
                  summarise(late_abundance = sum(Abundance))

# Plot by hypsodonty index (for ungulates only)
p <- ggplot(filter(late_colonisers, !is.na(hypsodonty_index)), aes(y = late_abundance*100, x = hypsodonty_index, group = Species, colour = Order)) +
  geom_boxplot() +
  facet_grid(rows=vars(genus), scales = "free") +
  geom_smooth(method = "lm", inherit.aes = FALSE, aes(y = late_abundance*100, x = hypsodonty_index), colour = "black", linetype = "dotted") +
  #scale_x_continuous(trans = "log10") +
  #scale_y_continuous(trans = "log10") +
  scale_color_manual(values = order_palette, name = "") +
  labs(x="Hypsodonty index", y="Late coloniser abundance (%)") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))

ggsave(file.path(subdir, "late_abund_by_hypsodonty.png"), p, width=5, height=8)

# Plot by molar category
p <- ggplot(filter(late_colonisers, molar_category!=""), aes(y = late_abundance*100, x = Species, fill = Order)) +
  geom_boxplot() +
  facet_grid(cols=vars(molar_category), rows=vars(genus), scales = "free", space = "free_x") +
  #scale_x_continuous(trans = "log10") +
  #scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = order_palette, name = "") +
  labs(x="Molar category", y="Late coloniser abundance (%)") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))

ggsave(file.path(subdir, "late_abund_by_molar_cat.png"), p, width=5, height=8)

