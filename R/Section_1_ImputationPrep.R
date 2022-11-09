# Section 1. Script for preparing data for imputation process including checking distributions and categorical frequencies.

# Required input for this script:
# dfRaw = Trait dataset with missing values.
# l_trees = Phylogenetic trees in Newick format (if using phylogenetic imputation).

# Output:
# CleanedCompleteCaseDataset.csv = Cleaned complete-case dataset.
# CleanedOriginalCaseDataset.csv = Cleaned original dataset with missing values.

# Acknowledgments. ----

# Trait data used in this study obtained from:
# Meiri S. Traits of lizards of the world: Variation around a successful evolutionary design. Glob Ecol Biogeogr. 2018;27(10):1168–72. 
# Meiri S. Data from: Traits of lizards of the world: Variation around a successful evolutionary design. Dryad Dataset [Internet]. 2019; Available from: https://doi.org/10.5061/dryad.f6t39kj

# Multiple sequence alignment obtained from:
# Pyron RA, Burbrink FT, Wiens JJ. A phylogeny and revised classification of Squamata, including 4161 species of lizards and snakes. BMC Evol Biol. 2013 Apr 29;13(1):93. 
# Pyron RA, Burbrink FT, Wiens JJ. Data from: A phylogeny and revised classification of Squamata, including 4161 species of lizards and snakes. Dryad Dataset [Internet]. 2013; Available from: https://doi.org/10.5061/dryad.82h0m

# Cleaned COI dataset provided by Dr. Cameron Nugent. Originally downloaded from the Barcode of Life Data (BOLD) System on March 4th, 2020 (associated with below citation):
# Citation: Nugent CA, Elliott TA, Ratnasingham S, Adamowicz SJ. coil: an R package for cytochrome c oxidase I (COI) DNA barcode data cleaning, translation, and error evaluation. Genome. 2020; 63:291-305.

### 1. Load libraries and functions. ----

library(ape)
library(data.table)
library(tidyverse)
source("R/Functions/DataHandling_Functions.R")
source("R/Functions/Phylo_Functions.R")

### 2. Dataset loading and variable assignment. ----

# Read in trait dataset you wish to impute. For an example, you can download the Meiri (2018) dataset available at the DOI mentioned above and subset to the columns you want to consider for imputation.
fileName <- file.choose()
dfRaw <- fread(fileName, data.table = F)
# Ensure blanks are NAs.
dfRaw[dfRaw == ""] <- NA
# Create vector of column names that contain taxonomic information. For example:
taxCols <- "species_name"
# Another example if using Meiri 2018 dataset:
colnames(dfRaw)
# Renaming column that contains  to "species_name (this is required for downstream analysis).
colnames(dfRaw)[1] <- "species_name"
# Cleaning up column names.
dfRaw <- CleanColumns(dfRaw)
# Subset by column name (removing misc column/columns we don't need but should check these in your final datasets prior to imputation e.g. for comments about data quality etc). For example:
dfRaw <- subset(dfRaw, select = -c(epithet, genus, valid_reptile_database_february_2018, year_of_description, country_described_from, known_only_from_the_only_type, mass_equation_feldman_et_al_2016_unless_stated, intercept, slope, activity_time_comments, substrate_comments, diet_comments, foraging_mode_comments, family, phylogeny, phylogenetic_data, remarks, references_biology_all_columns_except_m_n_and_o, references_svl_of_unsexed_individuals_neonates_and_hatchlings, references_svl_of_females, references_svl_of_males))

# Complete-case dataset. ---
# If already created, read in complete-case trait dataset.
fileName <- file.choose()
dfCC <- fread(fileName, data.table = F)
# Ensure blanks are NAs.
dfCC[dfCC == ""] <- NA

# OPTIONAL: Complete-case data-set creation. ---
# # To create a complete-case dataset, uncomment the following lines:
# ns <- seq(100, 5000, by = 100) ## These are the sample size thresholds to consider for your complete-case datasets. NOTE: In this example, the species name column is called "species_name".
# # Extract candidate trait names for this taxon (everything other than taxCols). I.e. the traits you want to *consider* for imputation.
# candTraits <- setdiff(colnames(dfRaw), taxCols)
# # Apply the CritNumberLoop() function to dfRaw to create a list of potential complete-case dataframes that contain data for at least three traits. This function writes the candidate complete-case dataframes to file.
# CritNumberLoop(data = dfRaw, traitCols = candTraits, taxCols = "species_name", critNumbers = ns)
# # Once selected, read in the complete-case dataset by uncommenting the following line:
# fileName <- file.choose()
# dfCC <- fread(fileName, data.table = F, header = T)

# Extract trait names for this taxon (everything other than taxCols). I.e. the traits you want to impute in dfRaw.
traits <- setdiff(colnames(dfCC), taxCols)
traits ## Double check that these are the correct trait names!
# Ensure species name column format in each dataframe matches. For example, to ensure species names have underscores and not spaces:
dfCC$species_name <- gsub(" ", "_", dfCC$species_name)
dfRaw$species_name <- gsub(" ", "_", dfRaw$species_name)
# Ensure columns in dfRaw match columns in dfCC.
dfCC <- dfCC[, c(taxCols, traits)]
dfRaw <- dfRaw[, c(taxCols, traits)]
# Remove columns without any trait information in dfRaw (after subsetting to candidate traits we want to test).
rmThese <- apply(dfRaw[, traits], 1, function(x) all(is.na(x)))
dfRaw <- dfRaw[!rmThese, ]

# NOTE: Check all classes of your traits now. Ensure they are of the correct class that we need (e.g. double or integer for numerical traits/character for categorical traits).
lapply(dfRaw[, traits], class)
lapply(dfCC[, traits], class)
# E.g. to convert to numerical:
#dfRaw$maximum_svl <- as.numeric(dfRaw$maximum_svl)
#dfCC$maximum_svl <- as.numeric(dfCC$maximum_svl)
# Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
l_traits <- BreakIntoTypes(dfRaw, traits)
# Extract numerical traits.
contTraits <- l_traits[[1]]
# Extract categorical traits.
catTraits <- l_traits[[2]]

# Combine complete-case and original data into list.
l_dfTraits <- list(Complete = dfCC, Original = dfRaw)

# OPTIONAL: If using phylogenetic imputation, uncomment the following lines.
# # Read in tree. For example:
# fileName <- file.choose()
# tree <- read.tree(fileName)
# # Ensure tip labels match species name columns.
# tree$tip.label <- gsub(" ", "_", tree$tip.label)
# # Match species names in tree to those in datasets. If sample size is less than 100, consider dropping traits here.
# l_matched <- mapply(DropAndMatch, data = l_dfTraits, MoreArgs = list(tree = tree), SIMPLIFY = F)
# # Extract updated trees (species not found in datasets were dropped).
# l_trees <- lapply(l_matched, function(x) x[[1]])
# # Extract updated datasets.
# l_dfTraits <- lapply(l_matched, function(x) x[[2]])

# Check sample sizes before proceeding!
lapply(l_dfTraits, nrow)

### 3. Check trait data. ----

# Descriptive info for numerical data:
l_contRes <- lapply(l_dfTraits, function(x) lapply(x[, contTraits], GetNumericalInfo))
# Descriptive info for categorical data:
l_catRes <- lapply(l_dfTraits, function(x) lapply(x[, catTraits], GetCategoricalInfo))

# Let's view the results.
for(i in 1:length(l_dfTraits)){
  # Take the name of the dataframe.
  dataName <- names((l_dfTraits))[[i]]
  # Take the ith contRes.
  contRes <- l_contRes[[i]]
  # For every continuous trait..
  for(cont in 1:length(contRes)){
    # Print the name of the trait.
    cat("Trait:", names((contRes))[[cont]])
    # Extract the result.
    result <- contRes[[cont]]
    # Print the results.
    print(result[1:6])
    # Plot the histogram.
    plot(result$Plot, col = "skyblue", main = paste(dataName, names((contRes))[[cont]], sep = " "))
  }
  
  # Take the ith catRes.
  catRes <- l_catRes[[i]]
  # For every categorical trait..
  for(cat in 1:length(catRes)){
    # Print the name of the trait.
    cat("Trait:", names((catRes))[[cat]])
    # Extract the result.
    result <- catRes[[cat]]
    # Print the results.
    print(result[1:4])
    # Convert frequency count into format amenable to plotting.
    dfCount  <- data.frame(result[[3]])
    # Plot barplot.
    barplot(height = dfCount$Freq, names.arg = dfCount$vector, col = "skyblue", main = paste(dataName, names((catRes))[[cat]], sep = " "))
  }
}

# Check for outliers in the numerical data (there will be two checks here, one for complete-case and one for original dataset).
l_dfOutliers <- lapply(contTraits, function(x) mapply(OutlierCheck, data = l_dfTraits, col = list(x), SIMPLIFY = F))
# Name according to contTraits.
names(l_dfOutliers) <- contTraits
# At this point, can check the names in l_dfOutliers to determine whether any of these are result of human error etc.
# For example:
# View(l_dfOutliers$smallest_clutch$Complete)

# Check for class imbalances in the complete-case data.
catFreqsCC <- mapply(ScreenCategories, variable = na.omit(l_dfTraits$Complete[, catTraits]), varName = catTraits)
# Check for class imbalances in the original data.
catFreqsO <- mapply(ScreenCategories, variable = na.omit(l_dfTraits$Original[, catTraits]), varName = catTraits)


# Alternatively, if you decide to remove a trait here, uncomment the following lines and remove in traits the original dataframes and create a new vector of traits for imputation. For example:
# dfRaw$clutch_size <- NULL
# dfCC$clutch_size <- NULL
# # Then rerun the following lines:
# # Get vector of traits.
# traits <- setdiff(colnames(dfCC), taxCols)
# # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
# l_traits <- BreakIntoTypes(dfRaw, traits)
# # Extract numerical traits.
# contTraits <- l_traits[[1]]
# # Extract categorical traits.
# catTraits <- l_traits[[2]]
# # AGAIN, remember to ensure the data are the correct classes!
# # Ensure columns in dfRaw match columns in dfCC.
# dfCC <- dfCC[, c(taxCols, traits)]
# dfRaw <- dfRaw[, c(taxCols, traits)]
# # Combine complete-case and original data into list.
# l_dfTraits <- list(Complete = dfCC, Original = dfRaw)

# At this point, you would remove categories that fall below the indicated threshold (e.g. comprising less than 10% of the data). Just make sure that when you're removing categories, it is in l_dfTraits as these are the datasets we will be working with now.
# Example:
#l_dfTraits$Original <- l_dfTraits$Original[!l_dfTraits$Original$insular_endemic %in% "unknown", ]

# If you remove observations here, it is recommended to run the ScreenCategories function again to ensure classes are balanced.
# Check for class imbalances in the complete-case data.
catFreqsCC <- mapply(ScreenCategories, variable = na.omit(l_dfTraits$Complete[, catTraits]), varName = catTraits)
catFreqsO <- mapply(ScreenCategories, variable = na.omit(l_dfTraits$Original[, catTraits]), varName = catTraits)

# If any data were removed, write the cleaned data to file for use in imputation simulations.
#fwrite(l_dfTraits$Complete, "Data/TraitData/CleanedCompleteCaseDataset.csv")
#fwrite(l_dfTraits$Original, "Data/TraitData/CleanedOriginalCaseDataset.csv")
