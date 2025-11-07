#### FUNCTIONS TO USE WITH PCA ORDINATIONS ####
library(wesanderson)
library(dplyr)
library(vegan)
library(tidyr)

## Add species centroid as phylopics
centroids <- function(ordination, phyloseq) {
  # Get ordination vectors as data frame and species info
  ord_df <- data.frame(vegan::scores(ordination, choices = c(1:4), display = "sites"))
  ord_df <- ord_df %>%
    cbind(select(data.frame(phyloseq@sam_data), c(Common.name, Species, Order_grouped, Order, diet.general, habitat.general)))
  # Group and calculate means
  centroids <- ord_df %>% group_by(Common.name, Species, Order_grouped, Order, diet.general, habitat.general) %>%
          summarise_all(.funs = mean, na.rm = TRUE)
  centroids <- left_join(centroids, phylopics, by = c("Species" = "Species"))
  return(centroids)
}

## Get loadings plot to display alongside PCA
loadings_plot <- function(ordination, axes, top_taxa = 20, shorten_names = TRUE) {
  # Get loadings
  loadings <- data.frame(vegan::scores(ordination, choices = axes, display = "species"))
  axes_names=paste0("PC", axes)
  eigval <- ordination$CA$eig[axes]
  # Keep only taxa with highest loadings
  # Get half of top taxa for each axis
  loadings_filt <- rbind(slice_max(loadings, order_by = abs(!!sym(axes_names[1])), n = ceiling(top_taxa/2)),
                         slice_max(loadings, order_by = abs(!!sym(axes_names[2])), n = floor(top_taxa/2) )) %>% 
    rownames_to_column("taxon") %>%
    # Get labels for plotting
    mutate(genus = gsub(" .*", "", taxon))
  if (shorten_names) {
    loadings_filt <- loadings_filt %>%
      # Create labels for taxa by keeping only first letter of genus
      mutate(label = gsub("([A-Z])[a-z]+", "\\1.", taxon))
  } else {
    loadings_filt$label <- loadings_filt$taxon
  }
  loadings_filt <- loadings_filt %>%
    # Arrange by first axis
    arrange(desc(!!sym(axes_names[1]))) %>% mutate(label = factor(label, levels = label))
  # Keep top 8 genera and group the rest as "Other"
  top_genera <- table(loadings_filt$genus) %>% sort(decreasing = TRUE) %>% names() %>% head(8)
  loadings_filt <- loadings_filt %>% mutate(genus = ifelse(genus %in% top_genera, genus, "Other"))
  # Get genus palette
  genera <- loadings_filt %>% filter(genus != "Other") %>% pull(genus) %>% unique
  genus_palette <- setNames(wes_palette("Darjeeling1", length(genera), type = "continuous"), genera)
  genus_palette["Other"] <- "grey"
  # Set genus as factor
  loadings_filt$genus <- factor(loadings_filt$genus, levels = names(genus_palette))
  # Melt and plot
  loadings_melt <- pivot_longer(loadings_filt, cols = starts_with("PC"), names_to = "PC", values_to = "loading")
  p <- ggplot(loadings_melt, aes(y = label, x = loading, fill = genus)) +
        geom_bar(stat = "identity", colour = "black") +
        facet_grid(~PC, scales = "free_y", space = "free_y", switch = "y") +
        geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
        scale_fill_manual(values = genus_palette, name = "Genus") +
        theme(axis.text.y = element_text(face = "italic", family = "serif", size = 8),
              axis.ticks = element_blank(), axis.title = element_blank(),
              panel.grid = element_blank(),
              strip.background = element_rect(colour = "white"),
              legend.position = "left")
  return(p)
}

## Get loadings plot to display alongside RDA
loadings_plot_rda <- function(ordination, axes, top_taxa = 20) {
  # Get loadings
  loadings <- data.frame(vegan::scores(ordination, choices = axes, display = "species"))
  axes_names=paste0("RDA", axes)
  eigval <- ordination$CA$eig[axes]
  # Keep only taxa with highest loadings
  # Get half of top taxa for each axis
  loadings_filt <- rbind(slice_max(loadings, order_by = !!sym(axes_names[1]), n = ceiling(top_taxa/2)),
                         slice_max(loadings, order_by = !!sym(axes_names[2]), n = floor(top_taxa/2) )) %>% 
    rownames_to_column("taxon") %>%
    # Get labels for plotting
    mutate(genus = gsub(" .*", "", taxon)) %>%
    # Create labels for taxa by keeping only first letter of genus
    mutate(label = gsub("([A-Z])[a-z]+", "\\1.", taxon)) %>%
    # Arrange by first axis
    arrange(desc(!!sym(axes_names[1]))) %>% mutate(label = factor(label, levels = label))
  # Keep top 8 genera and group the rest as "Other"
  top_genera <- table(loadings_filt$genus) %>% sort(decreasing = TRUE) %>% names() %>% head(8)
  loadings_filt <- loadings_filt %>% mutate(genus = ifelse(genus %in% top_genera, genus, "Other"))
  # Get genus palette
  genera <- loadings_filt %>% filter(genus != "Other") %>% pull(genus) %>% unique
  genus_palette <- setNames(wes_palette("Darjeeling1", length(genera), type = "continuous"), genera)
  genus_palette["Other"] <- "grey"
  # Set genus as factor
  loadings_filt$genus <- factor(loadings_filt$genus, levels = names(genus_palette))
  # Melt and plot
  loadings_melt <- pivot_longer(loadings_filt, cols = starts_with("RDA"), names_to = "RDA", values_to = "loading")
  p <- ggplot(loadings_melt, aes(y = label, x = loading, fill = genus)) +
        geom_bar(stat = "identity", colour = "black") +
        facet_grid(~RDA, scales = "free_y", space = "free_y", switch = "y") +
        geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
        scale_fill_manual(values = genus_palette, name = "Genus") +
        theme(axis.text.y = element_text(face = "italic", family = "serif", size = 8),
              axis.ticks = element_blank(), axis.title = element_blank(),
              panel.grid = element_blank(),
              strip.background = element_rect(colour = "white"),
              legend.position = "left")
  return(p)
}

arrow_coord <- function(ordination, axes) {
   # Get taxon pca scores
   scores = data.frame(vegan::scores(ord@ord, choices = axes, display = "species"))
   scores$distance = sqrt(scores[,1]^2 + scores[,2]^2)
   scores <- scores[order(scores$distance, decreasing = TRUE),]
   return(scores) 
}

library(RColorBrewer)
library(colorspace)

expand_palette <- function(subcategories_df, base_colours) {
  # Initialize an empty list to store the expanded palette
  expanded_palette <- c()
  categories <- unique(subcategories_df[1])[[1]]
  # Iterate over each category
  for (category in categories) {
    # Filter the subcategories for the current category
    subcats <- subcategories_df[subcategories_df[, 1] == category, 2][[1]]
    n_subcats <- length(subcats)
    # Calculate the number of lighter and darker shades
    n_lighter <- floor(n_subcats / 2)
    n_darker <- n_subcats - n_lighter -1
    
    # Generate lighter and darker shades
    lighter_shades <- lighten(base_colours[category], seq(0.1, 0.5, length.out = n_lighter))
    darker_shades <- rev(darken(base_colours[category], seq(0.1, 0.5, length.out = n_darker)))
    
    # Combine the shades with the base color in the middle
    palette <- c(darker_shades, base_colours[category], lighter_shades)
    palette <- palette[!is.na(palette)]
    names(palette) <- subcats
    # Assign colors to subcategories
    for (i in seq_along(subcats)) {
      expanded_palette <- append(expanded_palette, palette)
    }
  }
  return(expanded_palette)
}

custom_ord_plot <- function(phyloseq, ordination, colour_var, shape_var, arrows_scaling, type) {
  #### Get arrows ####
  # Get loading arrows coordinaties
  arrows <- arrow_coord(ord@ord, axes = c(1, 2, 3))
  # Get genus and phylum
  arrows$genus <- as.character(phyloseq@tax_table[match(rownames(arrows),  phyloseq@tax_table[, "species"]), "genus"])
  arrows$phylum <- as.character(phyloseq@tax_table[match(rownames(arrows),  phyloseq@tax_table[, "species"]), "phylum"])
  arrows$superkingdom <- as.character(phyloseq@tax_table[match(rownames(arrows),  phyloseq@tax_table[, "species"]), "superkingdom"])
  arrows$to_plot <- (rownames(arrows) %in% head(rownames(arrows), 100))
  # Group phyla for better plotting
  arrows <- arrows %>% mutate(phylum_grouped = factor(case_when(phylum %in% names(phylum_palette) ~ phylum,
                                                          superkingdom == "Bacteria" ~ "Other Bacteria",
                                                          superkingdom == "Archaea" ~ "Other Archaea"), levels = names(phylum_palette)))
  # Keep only strongest associations
  arrows_filt <- arrows %>% filter(to_plot) %>%
                select(contains(c("1", "2")), phylum_grouped)
  # Get species centroids
  centroids <- centroids(ord@ord, phyloseq)
  # Plot
  p <- ord_plot(ord, colour=colour_var, shape=shape_var, alpha = 0.8) +
    custom_theme() +
    geom_phylopic(data = centroids, aes_string(colour = colour_var), uuid = centroids$uid, fill = "transparent", width = 0.3)
  # Add the correct scales
  if (colour_var == "Order_grouped") {
    p <- p +
        scale_color_manual(values=order_palette, name = "Host order")
  } else if (colour_var == "diet.general") {
    p <- p +
        scale_color_manual(values=diet_palette, name = "Estimated diet")
  } else if (colour_var == "habitat.general") {
    p <- p +
        scale_colour_manual(values=habitat_palette, name = "Habitat")
  }
  if (shape_var == "diet.general") {
    p <- p +
        scale_shape_manual(values=diet_shape_scale, name = "Estimated diet")
  } else if (shape_var == "Order_grouped") {
    p <- p +
        scale_shape_manual(values=order_shape_scale, name = "Order")
  } else if (shape_var == "Common.name") {
    p <- p +
        scale_shape_manual(values=species_shape_scale, name = "Species")
  }
  # Add more layers
  p <- p +
    theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8)) +
    guides(shape = guide_legend(ncol = 2), colour = guide_legend(ncol = 2))
  # If PCA, add taxon arrows
  if (type == "PCA") {
    p <- p +
      new_scale_colour() +
      geom_segment(data = arrows_filt, aes(x = 0, y = 0, xend = PC1*arrows_scaling, yend = PC2*arrows_scaling, colour = phylum_grouped), linewidth = 0.5, alpha = 0.5) +
      scale_color_manual(values = phylum_palette, name = "Microbial phylum")
  }
  # If RDA add marginals
  if (type == "RDA") {
    p <- ggMarginal(p, type="violin", groupColour = TRUE, groupFill = TRUE, size=5)
  }
  return(p)
}

taxa_plot <- function(ord, phyloseq, ntaxa = 20) {
  # Get taxa scores and add taxonomic info
  taxa_rda <- data.frame(vegan::scores(ord@ord, display="species", choices=1:3)) %>%
              cbind(tax_table(phyloseq)[, c("superkingdom", "phylum", "genus", "species")]) %>%
              rownames_to_column("OTU") %>%
              # Get grouped phylum
              mutate(phylum_grouped = factor(case_when(phylum %in% names(phylum_palette) ~ phylum,
                                        superkingdom == "Bacteria" ~ "Other Bacteria",
                                        superkingdom == "Archaea" ~ "Other Archaea"), levels = names(phylum_palette)))
  # Identify 20 (or other defined number) genera with the highest correlation
  top_taxa <- taxa_rda %>% arrange(desc(sqrt(RDA1^2 + RDA2^2))) %>%
                slice_head(n = ntaxa) %>% pull(species)
  taxa_rda <- taxa_rda %>% mutate(label = ifelse(species %in% top_taxa, genus, NA)) %>%
              mutate(label = factor(species, levels = top_taxa)) %>%
              mutate(linetype = ifelse(species %in% top_taxa, "solid", "dashed"))  
  set.seed(245)
  p <- ggplot(taxa_rda, aes(x = 0, y = 0, xend = RDA1, yend = RDA2, colour = phylum_grouped)) +
    geom_segment(linewidth = 0.5, alpha = 0.8, aes(linetype = linetype)) +
    scale_linetype_identity() +
    scale_color_manual(values = phylum_palette, name = "Phylum") +
    geom_label(aes(label = label, x = RDA1, y = RDA2),
              size = 1.7, alpha = 0.5, vjust = ifelse(taxa_rda$RDA2 < 0, 1, 0)) +
    scale_size_continuous(range = c(0.01, 2), name = "Mean CLR-abundance") +
    custom_theme() + xlab("RDA1 scores") + ylab("RDA2 scores") +
    theme(legend.position = "bottom", legend.title = element_blank()) + guides(colour = guide_legend(nrow = 3)) +
    xlim(min(taxa_rda$RDA1)*1.1, max(taxa_rda$RDA1)*1.2)
  return(list(plot = p, data = taxa_rda))
}
