#### FUNCTIONS TO USE WITH PCA ORDINATIONS ####
library(wesanderson)
library(dplyr)
library(vegan)
library(tidyr)

## Add species centroid as phylopics
centroids <- function(ordination, phyloseq) {
  # Get ordination vectors as data frame and species info
  ord_df <- data.frame(scores(ordination, choices = c(1:4), display = "sites"))
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
  loadings <- data.frame(scores(ordination, choices = axes, display = "species"))
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
  loadings <- data.frame(scores(ordination, choices = axes, display = "species"))
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

arrow_coord <- function(ordination, phyloseq) {
   # Get fit vectors of taxa on ordination
   envfit_obj <- envfit(ordination, as.data.frame(t(otu_table(phyloseq))), permute = F, choices=c(1,2,3,4))
   # Extract relevant information for plotting (adapted from Jackie Zorz: https://jkzorz.github.io/2020/04/04/NMDS-extras.html)
   en_coord = as.data.frame(scores(envfit_obj, "vectors")) * ordiArrowMul(envfit_obj)
   # Add effect size and order and p-value
   en_coord$r <- envfit_obj$vectors$r
   en_coord$pval <- envfit_obj$vectors$pval
   # Get only significant fits and order by r
   sig <- names(which(envfit_obj$vectors$pval < 0.05))
   en_coord <- en_coord[sig,]
   en_coord <- en_coord[rev(order(en_coord$r)),]
   return(en_coord)
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
