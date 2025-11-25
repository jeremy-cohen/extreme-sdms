#### ============================================================
#### 0. Libraries & Paths
#### ============================================================
rm(list = ls())
library(ggplot2)
library(ape)
library(nlme)
# library(phytools) # optional
# library(geiger)   # not required

root <- "~/extremes"
data_directory <- file.path(root, "data")

#### ============================================================
#### 1. Data Loading & Preprocessing
#### ============================================================

# Read in main and auxiliary data
x <- read.csv(file.path(data_directory, "birds_VE8_comp.csv"))
x$sciname <- gsub(" ", "_", x$sciname)

x2 <- read.csv(file.path(data_directory, "birds_VE5b_comp.csv"))
x2$sciname <- gsub(" ", "_", x2$sciname)

# Join mass and habitat data
x$Habitat <- x2$Habitat[match(x$sciname, x2$sciname)]
x$Mass <- x2$Mass[match(x$sciname, x2$sciname)]
x$Trophic.Level <- x2$Trophic.Level[match(x$sciname, x2$sciname)]
x$Trophic.Niche <- x2$Trophic.Niche[match(x$sciname, x2$sciname)]

# Combine "Scavenger" into "Carnivore" in Trophic.Level
x$Trophic.Level[x$Trophic.Level == "Scavenger"] <- "Carnivore"

# Create a new habitat variable
x$Habitat_Group <- dplyr::case_when(
    x$Habitat %in% c("Forest", "Woodland") ~ "Forest",
    x$Habitat %in% c("Human modified", "Human Modified") ~ "Urban",
    x$Habitat %in% c("Marine", "Coastal", "Wetland", "Riverine") ~ "Water",
    x$Habitat %in% c("Grassland", "Rock", "Shrubland") ~ "Open",
    TRUE ~ "Other"
)

# Read and prune tree
tree <- ape::read.nexus(file.path(data_directory, "phylo.consensus.sumtrees.nex"))
tree_sub <- ape::drop.tip(tree, setdiff(tree$tip.label, x$sciname))

# source functions
source("~/extremes/scripts/00-functions.R")

#### ============================================================
#### 3. Analysis Pipeline
#### ============================================================

# User inputs
RESPONSES <- c("delta_range_area", "delta_prop_in_range")
BASE_PRED <- c("Mass", "Hand.Wing.Index")
EXTRA_PRED <- c("Habitat_Group", "Trophic.Level")
SEASONS <- c("all", "summer", "winter")
OUTLIER_ROW <- NA_integer_
CORR_KIND <- "pagel"

# Data prep
spdata <- prep_extremes(x, outlier_row = OUTLIER_ROW, VARIABLE = "extremes")
spdata <- standardize_sciname(spdata)
spdata$Mass <- log(spdata$Mass)
aligned <- align_tree_data(tree_sub, spdata)
tree_pruned <- aligned$tree
spdata <- aligned$data

#### ============================================================
#### 4. Model Fitting & Stepwise Selection
#### ============================================================

start_terms <- c(BASE_PRED, EXTRA_PRED, "season")
res_step <- backward_stepwise_pgls(
    df = spdata, tree = tree_pruned,
    response = "delta_prop_in_range",
    start_terms = start_terms,
    always_keep = "season",
    ssn = "all",
    corr_kind = "pagel",
    AIC_IMPROVE_MIN = 0,
    refit_REML = TRUE
)
cat("\n=== Backward stepwise (delta_prop_in_range; all seasons; Pagel λ) ===\n")
print(res_step$path)
print(summary(res_step$final))

res_step2 <- backward_stepwise_pgls(
    df = spdata, tree = tree_pruned,
    response = "delta_range_area",
    start_terms = start_terms,
    always_keep = "season",
    ssn = "all",
    corr_kind = "pagel",
    AIC_IMPROVE_MIN = 0,
    refit_REML = TRUE
)
cat("\n=== Backward stepwise (delta_range_area; all seasons; Pagel λ) ===\n")
print(res_step2$path)
print(summary(res_step2$final))

#### ============================================================
#### 5. Output Results
#### ============================================================

delta_area_table <- as.data.frame(summary(res_step2$final)$tTable)
delta_prop_table <- as.data.frame(summary(res_step$final)$tTable)
anova_delta_area <- as.data.frame(anova(res_step2$final))
anova_delta_prop <- as.data.frame(anova(res_step$final))

# Write to csv
write.csv(delta_area_table, "~/Downloads/E_delta_area_table.csv", row.names = TRUE)
write.csv(delta_prop_table, "~/Downloads/E_delta_prop_table.csv", row.names = TRUE)
write.csv(anova_delta_area, "~/Downloads/E_anova_delta_area.csv", row.names = TRUE)
write.csv(anova_delta_prop, "~/Downloads/E_anova_delta_prop.csv", row.names = TRUE)
