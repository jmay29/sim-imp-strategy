# Section 1. Script for preparing data for imputation process including checking distributions and categorical frequencies.

# Required input for this script:
# dfRaw = Trait dataset with missing values.
# l_trees = Phylogenetic trees in Newick format (if using phylogenetic imputation).

# Output:
# CleanedCompleteCaseDataset.csv = Cleaned complete-case dataset.
# CleanedOriginalCaseDataset.csv = Cleaned original dataset with missing values.

# Acknowledgments. ----

# Trait data used in this study obtained from:
# Meiri, S. Traits of lizards of the world: Variation around a successful evolutionary design. Global Ecol Biogeogr. 2018; 27: 1168– 1172. https://doi.org/10.1111/geb.12773
# Meiri, Shai (2019), Data from: Traits of lizards of the world: variation around a successful evolutionary design, Dryad, Dataset, https://doi.org/10.5061/dryad.f6t39kj

# Multiple sequence alignment obtained from:
# Pyron, R.A., Burbrink, F.T. & Wiens, J.J. A phylogeny and revised classification of Squamata, including 4161 species of lizards and snakes. BMC Evol Biol 13, 93 (2013). https://doi.org/10.1186/1471-2148-13-93
# Pyron RA, Burbrink FT, Wiens JJ. Data from: A phylogeny and revised classification of Squamata, including 4161 species of lizards and snakes. Dryad Dataset [Internet]. 2013; Available from: https://doi.org/10.5061/dryad.82h0m

# Cleaned COI dataset provided by Dr. Cameron Nugent. Originally downloaded from the Barcode of Life Data (BOLD) System on March 4th, 2020.

### 1. Load libraries and functions. ----

library(ape)
library(data.table)
library(tidyverse)
source("R/Functions/DataHandling_Functions.R")
source("R/Functions/Phylo_Functions.R")

### 2. Dataset loading and variable assignment. ----

# Read in trait dataset you wish to impute. For an example, you can download the Meiri (2018) dataset available at the DOI mentioned above.
fileName <- file.choose()
dfRaw <- fread(fileName, data.table = F)
# Ensure blanks are NAs.
dfRaw[dfRaw == ""] <- NA
# Create vector of column names that contain taxonomic/misc information. For example:
taxCols <- "species_name"
# Extract trait names for this taxon (everything other than taxCols). I.e. the traits you want to consider for imputation.
traits <- setdiff(colnames(dfRaw), taxCols)

# Complete-case dataset. ---
# If already created, read in complete-case trait dataset.
fileName <- file.choose()
dfCC <- fread(fileName, data.table = F)
# Ensure blanks are NAs.
dfCC[dfCC == ""] <- NA

# Ensure species name column in each dataframe matches. For example, to ensure species names have underscores and not spaces:
dfCC$species_name <- gsub(" ", "_", dfCC$species_name)
dfRaw$species_name <- gsub(" ", "_", dfRaw$species_name)

# OPTIONAL: Complete-case data-set creation. ---
# If you don't want to allow for any missingness in your complete-case dataset:
# Apply the CritNumberLoop() function to dfRaw to create a list of potential complete-case dataframes that contain data for at least three traits. This function writes the candidate complete-case dataframes to file.
# Uncomment the following lines if you wish to use this function:
#ns <- seq(100, 5000, by = 100) ## These are the sample size thresholds to consider for your complete-case datasets.
#CritNumberLoop(data = dfRaw, traitCols = traits, taxCols = "species_name", critNumbers = ns)
# Once selected, read in the complete-case dataset by uncommenting the following line:
#fileName <- file.choose()
#dfCC <- fread(fileName)

# Combine complete-case and original data into list.
l_dfTraits <- list(Complete = dfCC, Original = dfRaw)

# OPTIONAL: If using phylogenetic imputation, uncomment the following lines.
# Read in tree. For example:
#fileName <- file.choose()
#tree <- read.tree(fileName)
# Ensure tip labels match species name columns.
#tree$tip.label <- gsub(" ", "_", tree$tip.label)
# Match species names in tree to those in datasets. If sample size is less than 100, consider dropping traits here. 
#l_matched <- mapply(DropAndMatch, data = l_dfTraits, MoreArgs = list(tree = tree), SIMPLIFY = F) 
# Extract updated trees (species not found in datasets were dropped).
#l_trees <- lapply(l_matched, function(x) x[[1]])
# Extract updated datasets.
#l_dfTraits <- lapply(l_matched, function(x) x[[2]])

# Check sample sizes before proceeding!
lapply(l_dfTraits, nrow)

# Another option here is to allow for some missingness in your complete-case dataset (i.e. a "nearly" complete-case dataset - see Optional_NearlyCompleteCaseCreation.R script if you wish to do this).

### 3. Check trait data. ----

# Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
l_traits <- BreakIntoTypes(l_dfTraits[[1]], traits)
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
catFreqsCC <- lapply(na.omit(l_dfTraits$Complete[, catTraits]), ScreenCategories)
print(catFreqsCC)
# Check for class imbalances in the original data.
catFreqsO <- lapply(na.omit(l_dfTraits$Original[, catTraits]), ScreenCategories)
print(catFreqsO)

# At this point, you would remove categories that fall below the indicated threshold (e.g. comprising less than 10% of the data). 
# Example:
l_dfTraits$Original <- l_dfTraits$Original[!l_dfTraits$Original$insular_endemic == "unknown", ]

# If you remove observations here, it is recommended to run the ScreenCategories function again to ensure classes are balanced.
#catFreqsCC <- lapply(na.omit(l_dfTraits$Complete[, catTraits]), ScreenCategories)
#print(catFreqsCC)
#catFreqsO <- lapply(na.omit(l_dfTraits$Original[, catTraits]), ScreenCategories)
#print(catFreqsO)

# If any data were removed, write the cleaned data to file for use in imputation simulations.
#fwrite(l_dfTraits$Complete, "Data/TraitData/CleanedCompleteCaseDataset.csv")
#fwrite(l_dfTraits$Original, "Data/TraitData/CleanedOriginalCaseDataset.csv")