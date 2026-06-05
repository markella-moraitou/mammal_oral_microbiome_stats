library(ggplot2)

###############
#### SETUP ####
###############

## SET THEME
custom_theme <- function() {
  theme_bw() +
    theme(
      # Customize other theme elements here
      axis.text=element_text(size = 12),
      axis.text.x=element_text(angle=90),
      axis.title=element_text(size = 14, face = "bold"),
      legend.title=element_text(size = 14, face = "bold"),
      legend.text=element_text(size=12),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "white", colour = "black")
    )
}

## SET THEME
rda_theme <- function() {
  theme_classic() +
    theme(
      # Customize other theme elements here
      axis.text=element_text(size = 12),
      axis.text.x=element_text(angle=90),
      axis.title=element_text(size = 14, face = "bold"),
      legend.title=element_text(size = 14, face = "bold"),
      legend.text=element_text(size=12),
      plot.background=element_rect(fill = "white", colour = "white")
    )
}

large_font <- function() {
  theme_bw() +
    theme(
      # Customize other theme elements here
      axis.text=element_text(size = 14),
      axis.text.x=element_text(angle=90),
      axis.title=element_text(size = 18, face = "bold"),
      legend.title=element_text(size = 14, face = "bold"),
      legend.text=element_text(size=12),
      strip.background = element_rect(fill = "white", colour = "black")
    )
}

plot_setup <- function(palette_dir) {
  # Set the working directory to the location of this script
  theme_set(custom_theme())

  ## LOAD PALETTES
  if (dir.exists(file.path(palette_dir))) {
    cat("Loading custom palettes\n")

    species_palette <<- read.csv(file=file.path(palette_dir, "species_palette.csv"))$x
    names(species_palette) <<- read.csv(file=file.path(palette_dir, "species_palette.csv"))$X

    order_palette <<- read.csv(file=file.path(palette_dir, "order_palette.csv"))$x
    names(order_palette) <<- read.csv(file=file.path(palette_dir, "order_palette.csv"))$X

    diet_palette <<- read.csv(file=file.path(palette_dir, "diet_palette.csv"))$x
    names(diet_palette) <<- read.csv(file=file.path(palette_dir, "diet_palette.csv"))$X

    habitat_palette <<- read.csv(file=file.path(palette_dir, "habitat_palette.csv"))$x
    names(habitat_palette) <<- read.csv(file=file.path(palette_dir, "habitat_palette.csv"))$X
    
    phylum_palette <<- read.csv(file=file.path(palette_dir, "phylum_palette.csv"))$x
    names(phylum_palette) <<- read.csv(file=file.path(palette_dir, "phylum_palette.csv"))$X
    
  } else {
    cat("No custom palettes found\n")
  }
}