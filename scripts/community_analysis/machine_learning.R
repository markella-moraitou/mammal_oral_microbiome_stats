##### MACHINE LEARNING #####

#### Use machine learning approaches like tsne and random forest

################
#### SET UP ####
################

library(dplyr)
library(phyloseq)
library(microbiome)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(tsne)
library(rphylopic)
library(mlr3)
library(mlr3viz)
library(mlr3tuning)
library(mlr3mbo)
library(mlr3learners)
library(mlr3extralearners)
library(mlr3tuningspaces)
library(patchwork)
library(paradox)
library(rlist)

#### VARIABLES AND WORKING DIRECTORY ####

# Directory and file paths paths
indir <- normalizePath(file.path("..", "..", "input")) # Directory with phyloseq output and sample metadata 
outdir <- normalizePath(file.path("..", "..", "output", "community_analysis"))
subdir <- normalizePath(file.path(outdir, "machine_learning"))
phydir <- normalizePath(file.path(outdir, "phyloseq_objects")) # Directory with phyloseq objects

# Create output directory if it doesn't exist
if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

## Set up for plotting
source(file.path("..", "plot_setup.R"))
plot_setup(file.path("..", "..", "input", "palettes"))
theme_set(custom_theme())

#######################
#####  LOAD INPUT #####
#######################

# Load phy
phy_sp_f <- readRDS(file.path(phydir, "phy_sp_f.RDS"))
phy_sp_f_clr <- readRDS(file.path(phydir, "phy_sp_f_clr.RDS"))

phylopics <- read.csv(file.path(indir, "palettes", "phylopics.csv"), stringsAsFactors = FALSE)

# Aggregate to genus level
phy_gen <- tax_glom(phy_sp_f, taxrank = "genus")
phy_gen_clr <- microbiome::transform(phy_gen, "clr")

# Extract relevant species metadata
species_meta <- data.frame(phy_gen@sam_data) %>% select(Species, Common.name, diet.general, Order) %>% unique

##################
#### RUN TSNE ####
##################

# Extract OTU table and sample data
otu_table <- as.data.frame(phy_gen_clr@otu_table)
sample_data <- as.data.frame(phy_gen_clr@sam_data)

# Transpose the OTU table
otu_table_t <- t(otu_table)

# Perform t-SNE
set.seed(123)

tsne_result <- tsne(otu_table_t)

#### Plot t-SNE results

tsne_df <- as.data.frame(tsne_result) %>% 
  rename(tsne1 = V1, tsne2 = V2)

tsne_df$Sample <- sample_names(phy_sp_f)
tsne_df$Species <- sample_data$Species
tsne_df$Order <- sample_data$Order_grouped
tsne_df$Diet <- sample_data$diet.general

# Get centroids
centroids <- tsne_df %>%
  group_by(Species, Order, Diet) %>%
  summarise(tsne1 = mean(tsne1), tsne2 = mean(tsne2)) %>%
  left_join(phylopics, by = "Species")

# Get shape scales for plotting
diet_shape_scale <- c("Animalivore" = 8, "Omnivore" = 9 , "Frugivore" = 2, "Herbivore" = 16)

uniq_species <- unique(subset_samples(phy_sp_f_clr, Order %in% c("Primates", "Carnivora", "Artiodactyla") | habitat.general == "Aquatic")@sam_data$Common.name)
species_shape_scale <- c(1:25, 35:35+26-length(uniq_species))
names(species_shape_scale) <- uniq_species

order_shape_scale <- c("Carnivora" = 4, "Primates" = 19, "Artiodactyla" = 5, "Perissodactyla" = 2, "Rodentia" = 1, "Rest" = 12)

# Plot
p <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = Order, shape = Diet)) +
  geom_point(size = 2, alpha = 0.5) +
  labs(title = "t-SNE of OTU Table", x = "t-SNE 1", y = "t-SNE 2") +
  scale_color_manual(values = order_palette) +
  scale_shape_manual(values=diet_shape_scale, name = "Estimated diet") +
  geom_phylopic(data = centroids, aes(colour = Order, uuid = uid), width = 3, alpha = 0.8)

ggsave(p, filename = file.path(subdir, "tsne_plot_order.png"), width = 8, height = 6)

p <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = Diet, shape = Order)) +
  geom_point(size = 2, alpha = 0.5) +
  labs(title = "t-SNE of OTU Table", x = "t-SNE 1", y = "t-SNE 2") +
  scale_color_manual(values = diet_palette) +
  scale_shape_manual(values=order_shape_scale, name = "Order") +
  geom_phylopic(data = centroids, aes(colour = Diet, uuid = uid), width = 3, alpha = 0.8)

ggsave(p, filename = file.path(subdir, "tsne_plot_diet.png"), width = 8, height = 6)

##########################
#### MACHINE LEARNING ####
##########################

#### Function for tuning ####
# Write function for tuning and saving results: this will be used to tune several candidate learners before benchmarking

tune_and_save <- function(learner, label, task, dir) {
  # If best params exist, load
  if(file.exists(file.path(dir, "best_params.yml"))) {
    best_params <<- list.load(file.path(dir, "best_params.yml"))
  } else {
    best_params <<- list()
  }
  # Define resampling strategy
  set.seed(123)
  rcv = rsmp("repeated_cv", repeats = 10, folds = 10)
  #rcv = rsmp("repeated_cv", repeats = 2, folds = 5)
  # Define terminator
  term = trm("evals", n_evals = 100)
  #term = trm("evals", n_evals = 5)
  # Define tuner: Bayesian optimisation
  tuner = tnr("mbo")
  
  # Tune learner on task
  instance = TuningInstanceBatchSingleCrit$new(task, learner, rcv, measure, term)
  tuner$optimize(instance)
  instance$result
  
  # Plot tuning results
  p <- wrap_plots(autoplot(instance))
  ggsave(p, filename  =  file.path(dir, paste0(label, "_tuning.png")), width  =  13, height = 8)
  
  # Add tuning to best_params list and save list
  best_params[[label]] <<- instance$result_learner_param_vals
  list.save(best_params, file.path(dir, "best_params.yml"))
  # Return best parameters
  return(instance$result_learner_param_vals)
}

#### CLASSIFICATION TASK ####

# Get response and predictors in one table
ml_data <- data.frame(phy_gen_clr@otu_table) %>% t %>%
          cbind(select(data.frame(phy_sp_f_clr@sam_data), c(Species))) %>% #, Order, diet.general, habitat.general))) %>%
          rownames_to_column("Sample") %>%
          # Keep only species with > 3 samples
          filter(Species %in% names(which(table(Species) > 3)))

colnames(ml_data) <- str_replace_all(colnames(ml_data), pattern = " ", replacement = ".") %>% str_replace_all(pattern = "-", replacement = "_")

# Turn habitat into factor
#ml_data$habitat.general <- factor(ml_data$habitat.general, levels = c("Terrestrial", "Marine"))
# Drop empty levels from response variable
ml_data$Species <- droplevels(ml_data$Species)

## Run to predict species
response <- "Species"
predictors <- setdiff(colnames(ml_data), c(response, "Sample")) %>% paste(., collapse = "+")

# Create task with indexing copies log as the target
task <- as_task_classif(as.formula(paste(response, "~", predictors)), data = ml_data)

#### Split training - testing ####
# Set species as a stratum (so both train and test datasets get a similar distribution of these factors)
task$set_col_roles("Species", c("target", "stratum"))

set.seed(123)
splits = partition(task, ratio = 0.67)

# Plot datasets to verify stratification
task_split <- task$data()
task_split$Species <- ml_data$Species
task_split$dataset <- ifelse(task$row_ids %in% splits$train, "train",
                             ifelse(task$row_ids %in% splits$test, "test", NA))

p <- ggplot(aes(x = dataset, fill = Species), data = task_split) +
    geom_bar() +
    scale_fill_manual(values = species_palette)

ggsave(p, filename = file.path(subdir, "ml_data_split.png"), width = 8, height = 6)

# Save
write.csv(task_split, file = file.path(subdir, "task_split.csv"), quote = FALSE, row.names = FALSE)

# Use only the training data for tuning
task_train = task$clone()$filter(splits$train)
task_test = task$clone()$filter(splits$test)

# Classification measure
measure = msr("classif.acc")

# Output directory
dir <- file.path(subdir, "ml_species")
if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

#### Tune kknn ####
lrn_kknn = lrn("classif.kknn",
               k = to_tune(1, 30),
               distance = to_tune(0, 5),
               kernel = to_tune(c("rectangular", "triangular", "epanechnikov")))

lrn_kknn$param_set$values <- tune_and_save(lrn_kknn, "classif.kknn", task_train, dir)

saveRDS(object = lrn_kknn, file = file.path(dir, "lrn_kknn.RDS"))

#### Tune ranger ####
lrn_ranger = lrn("classif.ranger",
                 num.trees = to_tune(100, 1000), 
                 mtry = to_tune(1, task_train$n_features), 
                 min.node.size = to_tune(p_int(3, 25)), 
                 sample.fraction = to_tune(0.5, 1), 
                 max.depth = to_tune(4, 25), 
                 splitrule = to_tune(c("gini", "extratrees")),
                 importance = "impurity")

# Tune on training dataset
lrn_ranger$param_set$values <- tune_and_save(lrn_ranger, "classif.ranger", task_train, dir)

saveRDS(object = lrn_ranger, file = file.path(dir, "lrn_ranger.RDS"))

#### Tune rpart ####
lrn_rpart = lrn("classif.rpart",
                cp = to_tune(0.001, 0.1),  # Complexity parameter
                minbucket = to_tune(5, 30),  # Minimum number of observations in any terminal node
                maxdepth = to_tune(1, 30)  # Maximum depth of any node of the final tree
                )

lrn_rpart$param_set$values <- tune_and_save(lrn_rpart, "classif.rpart", task_train, dir)

saveRDS(object = lrn_rpart, file = file.path(dir, "lrn_rpart.RDS"))

#### Tune GLM with regularisation ####
lrn_glmnet = lrn("classif.glmnet",
                lambda = to_tune(p_dbl(0.1, 1)), # Regularization parameter
                nlambda = to_tune(p_int(10, 100)), # Number of lambda values to generate
                standardize = TRUE
                )

lrn_glmnet$param_set$values <- tune_and_save(lrn_glmnet, "classif.glmnet", task_train, dir)

saveRDS(object = lrn_glmnet, file = file.path(dir, "lrn_glmnet.RDS"))

#### Benchmark ####
rcv = rsmp("repeated_cv", repeats = 10, folds = 10)

# Get benchmark design: Benchmark using the training set

# Great list of learners (load from RDS)
learner_files <- dir(dir, pattern = "lrn.*.RDS")
learners <- list()

# Add featureless learner
learners[["classif.featureless"]] = lrn("classif.featureless")

# Add learners to list
for (file in learner_files) {
  obj <- readRDS(file.path(dir, file))
  learners[[str_remove(file, ".RDS")]] <- obj
}

design = benchmark_grid(
  tasks = task_train,
  learners = learners,
  resamplings = rcv
)

# Save plot
set.seed(12)
bmr = benchmark(design)
p <- autoplot(bmr, measure = measure) +
  custom_theme() +
  scale_fill_viridis_d() +
  theme(plot.background = element_rect(fill = "white")) +
  theme(legend.location = "none")

ggsave(p, filename  =  file.path(dir, "benchmarking_classif.png"), width  =  13, height = 8)

# Save aggregate scores
bmr_scores <- bmr$aggregate(measure) %>% select(learner_id, measure$id)
write.csv(bmr_scores, file = file.path(dir, "bmr_classif_aggregate_scores.csv"), quote = FALSE, row.names = FALSE)

#### Evaluate using test dataset ####
evaluations <- list()

for (learner in learners) {
  set.seed(12)
  learner$train(task_train)
  evaluations[paste(learner$id, "test", sep = "_")] <- learner$predict(task_test)$score(measure)
  evaluations[paste(learner$id, "train", sep = "_")] <- learner$predict(task_train)$score(measure)
  # Convert the confusion matrix to a data frame for ggplot
  confusion_df <- as.data.frame(learner$predict(task_test)$confusion)
  # Add diet and order info
  confusion_meta <- confusion_df %>%
    left_join(species_meta, by = c("truth" = "Species"), relationship = "many-to-one") %>%
    left_join(species_meta, by = c("response" = "Species"), suffix = c(".truth", ".response"), relationship = "many-to-one") %>%
    # Is it same order and/or same diet
    mutate(class = case_when(truth == response ~ "Correct",
                             Order.truth == Order.response & diet.general.truth == diet.general.response ~ "Within order & diet",
                             Order.truth == Order.response ~ "Within order",
                             diet.general.truth == diet.general.response ~ "Within diet",
                             TRUE ~ NA)) %>%
    # Arrange by order and then diet
    arrange(Order.truth, diet.general.truth, Common.name.truth) %>%
    mutate(Common.name.truth = factor(Common.name.truth, levels = unique(Common.name.truth)),
           Common.name.response = factor(Common.name.response, levels = unique(Common.name.truth))) %>%
    # Calculate relative frequency of predictions (per truth class)
    group_by(truth) %>% mutate(Rel_Freq = Freq/sum(Freq)) %>%
    mutate(Rel_Freq_rounded = case_when(Rel_Freq > 0 ~round(Rel_Freq, 1)))
  
  # Plotting the confusion matrix
  p <- 
    ggplot(confusion_meta, aes(x = Common.name.response, y = Common.name.truth, fill = Rel_Freq, colour = class)) +
    geom_tile(linewidth = 0.5) +
    geom_text(aes(label = Rel_Freq_rounded), size = 3, color = "black") +
    scale_fill_gradient2(low = "white", mid = "#FCAA4C", high = "#FC714C", midpoint = 0.5) +
    scale_colour_manual(values = c("Correct" = "#84E936", "Within order & diet" = "#4FB700", "Within order" = "#0366A1", "Within diet" = "#FCEE00"),
                        name = "Classification", na.value = "transparent") +
    labs(title = "Confusion Matrix", x = "Predicted", y = "Actual") +
    theme(axis.text.x = element_text(hjust = 1))

  ggsave(p, filename  =  file.path(dir, paste0(learner$id, "_predictions.png")), width  =  13, height = 8)

}

# Plot
evaluations_df <- unlist(evaluations) %>% as.data.frame()
colnames(evaluations_df) <- measure$label
evaluations_df$learner <- rownames(evaluations_df) %>% str_remove(., "lrn_") %>% paste0("classif.", .) %>% str_remove("_.*")
evaluations_df$dataset <- rownames(evaluations_df) %>% str_remove(".*_")

write.csv(evaluations_df, file = file.path(dir, "model_evaluation_classif.csv"), quote = FALSE, row.names = FALSE)

p <- ggplot(aes(y = !!sym(measure$label), x = dataset, fill = learner), data = evaluations_df) +
  facet_grid(cols = vars(learner)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() + theme(legend.position = "none") 

ggsave(p, filename  =  file.path(dir, "model_evaluation_classif.png"), width  =  13, height = 8)
