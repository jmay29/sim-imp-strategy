# Section 1. Script for preparing data for imputation process including checking distributions and categorical frequencies.

### 1. Load libraries and functions. ----

library(ape)
library(data.table)
library(tidyverse)
source("R/Functions/DataHandling_Functions.R")
source("R/Functions/Phylo_Functions.R")

### 2. Dataset loading and variable assignment. ----

# Read in the original trait dataset with missing data.
dfRaw <- fread("Data/TraitData/OriginalDataset.csv", data.table = F)
# Create vector of column names that contain taxonomic/misc information.
taxCols <- "species_name"
# Extract trait names for this taxon.
traits <- setdiff(colnames(dfRaw), taxCols)

# Complete-case dataset creation. ---
# If already created, read in the complete-case trait dataset.
dfCC <- fread("Data/TraitData/CompleteCaseDataset.csv", data.table = F)
# OPTIONAL:
# Let's create candidate complete-case datasets.
# Apply the CritNumberLoop() function to dfRaw to create a list of potential complete-case dataframes that contain data for at least three traits. This function writes the candidate complete-case dataframes to file. Uncomment the following line if you wish to use this function:
#CritNumberLoop(l_df = list(dfRaw), traitCols = traits, taxCols = "species_name")
# Once selected, read in the complete-case dataset by uncommenting the following line:
#dfCC <- file.choose()

# If using phylogenetic imputation, read in tree(s).
treeFiles <- list.files(path = "Data/RAxMLTrees/", pattern = "_ultra.tre")
l_trees <- lapply(treeFiles, function(x) read.tree(paste("Data/RAxMLTrees/", x, sep = "")))
# Name l_trees.
names(l_trees) <- treeFiles
# Match species names in trees to those in complete-case dataset. If less than 100, consider dropping traits here.
l_matched <- lapply(l_trees, DropAndMatch, data = dfCC)
# Extract updated trees (species not found in complete-case dataset were dropped from trees).
l_trees <- lapply(l_matched, function(x) x[[1]])
# Check number of species in the trees.
sapply(l_trees, function(x) length(x$tip.label)) ## Greater than 100? If not, consider excluding tree from imputation simulations.

# Combine complete-case and original data into list.
l_dfTraits <- list(Complete = dfCC, Original = dfRaw)

### 3. Check trait data. ----

# Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
l_traits <- BreakIntoTypes(dfCC, traits)
# Extract numerical traits.
contTraits <- l_traits[[1]]
# Extract categorical traits.
catTraits <- l_traits[[2]]

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
View(l_dfOutliers$smallest_clutch$Complete)

# Check for class imbalances in the complete-case data.
catFreqsCC <- lapply(dfCC[, catTraits], ScreenCategories)
print(catFreqsCC)
# Check for class imbalances in the original data.
catFreqsO <- lapply(dfRaw[, catTraits], ScreenCategories)
print(catFreqsO)

# At this point, you would remove categories that fall below the indicated threshold (e.g. comprising less than 10% of the data).
# Example:
dfRaw <- dfRaw[!dfRaw$insular_endemic == "unknown", ]

# If any data were removed, write the cleaned data to file for use in imputation simulations. No alterations to complete-case dataset so we can leave it as is.
fwrite(dfCC, "Data/TraitData/CleanedCompleteCaseDataset.csv")
fwrite(dfRaw, "Data/TraitData/CleanedOriginalCaseDataset.csv")
