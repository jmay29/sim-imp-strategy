# Section 2. This script imputes trait data using mean/mode replacement, k-nearest neighbour (KNN), missForest (RF), and multivariate imputation using chained equations (MICE) with and without phylogeny (in the form of phylogenetic eigenvectors derived from a phylogenetic tree). Given a complete-case dataset, it simulates data missing completely at random (MCAR), missing at random (MAR), and missingness not at random (MNAR). It then uses different imputation methods to fill in the missing values using the known observations. The outputs for this script are the imputation error rates for each trait in the dataset. These are written to the user's current working directory.

# Required input for this script:
# dfCC = Cleaned complete-case dataset from section 1 (CleanedCompleteCaseDataset.csv).
# dfRaw = Cleaned original dataset (with missing values) (CleanedOriginalCaseDataset.csv) from section 1.
# l_trees = Phylogenetic trees in Newick format (if using phylogenetic imputation).

# Output:
# Error rate csv files for each trait/missing pattern.
# Best parameter csv files for each trait/missing pattern.

# Acknowledgments. ----

# This script makes use of the following imputation functions:

# kNN() function in the "VIM" package for imputation.
# Citations: Kowarik A, Templ M. Imputation with the R Package VIM. J Stat Softw. 2016;74(7):1–16. 
# Templ M, Kowarik A, Alfons A, de Cillia G, Prantner B, Rannetbauer W. R package “VIM”: Visualization and Imputation of Missing Values [Internet]. 2021. Available from: https://cran.r-project.org/web/packages/VIM/VIM.pdf. R package v. 6.1.0.

# missForest() function in the "missForest" package for imputation.
# Citations: 40.	Stekhoven DJ. missForest: Nonparametric Missing Value Imputation using Random Forest. 2013.R package v. 1.4.
# Stekhoven DJ, Bühlmann P. MissForest—non-parametric missing value imputation for mixed-type data. Bioinformatics. 2012 Jan 1;28(1):112–8. 

# mice() function in the "mice" package for imputation.
# Citations: van Buuren S, Groothuis-Oudshoorn K. MICE: Multivariate Imputation by Chained Equations in R. J Stat Softw. 2011;45(3):1–67. R package v. 3.13.0.
# https://cran.r-project.org/web/packages/mice/mice.pdf.

### 1. Load libraries and functions. ----

library(ape)
library(data.table)
library(Information)
library(lmtest)
library(doParallel)
library(doRNG)
library(Information)
library(lmtest)
library(matrixStats)
library(mclust)
library(matrixStats)
library(Metrics)
library(mice)
library(missForest)
library(mltools)
library(naniar)
library(nnet)
library(phytools)
library(plotrix)
library(MPSEM)
library(rsample)
library(simputation)
library(tidyverse)
library(VIM)
source("R/Functions/DataHandling_Functions.R")
source("R/Functions/Imputation_Functions.R")
source("R/Functions/Phylo_Functions.R")
source("R/Functions/Simpute_Functions.R")
source("R/Functions/ParlMICEWrapper_Functions.R")

### 2. Data loading and variable assignment. ----

# Read in the cleaned original trait dataset.
fileName <- file.choose()
dfRaw <- fread(fileName, data.table = F)
# Ensure blanks are NAs.
dfRaw[dfRaw == ""] <- NA
# Read in cleaned complete-case trait dataset.
fileName <- file.choose()
dfCC <- fread(fileName, data.table = F)
# Ensure blanks are NAs.
dfCC[dfCC == ""] <- NA

# If using phylogenetic imputation, use the following lines. Customize according to your tree file and what you want to call your tree(s).
# Read in tree(s). For example:
treeFiles <- choose.files()
l_trees <- lapply(treeFiles, read.tree)
# Name l_trees according to file name.
names(l_trees) <- treeFiles
treeNames <- names(l_trees)
# Create vector of tree names (used to name error rate files later on). For example:
treeNames <- sapply(treeFiles, function(x) word(x, start = 2, sep = "_"))

# Which columns contain the taxonomic (or misc) information in dfCC?
taxCols <- "species_name"
# Which column contains species name information?
speciesCol <- "species_name"
# Extract trait names.
traits <- setdiff(colnames(dfCC), taxCols)
# Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
l_traits <- BreakIntoTypes(dfCC, traits)
# Extract numerical traits.
contTraits <- l_traits[[1]]
# Extract categorical traits.
catTraits <- l_traits[[2]]
# Ensure columns in dfRaw and dfCC match.
dfCC <- dfCC[, c(speciesCol, traits)]
dfRaw <- dfRaw[, c(speciesCol, traits)]
# Ensure matching species columns. We are adding underscores to the species names.
dfCC[, speciesCol] <- gsub(" ", "_", dfCC[, speciesCol])
dfRaw[, speciesCol] <- gsub(" ", "_", dfRaw[, speciesCol])
# Also ensure underscores in the tip labels of your tree if using phylogenetic imputation!

# Create an integer for the number of replicates to use for missingness simulations and imputation. Default here is 10 replicates.
r <- 2
# If simulating data missing completely at random (MCAR), create an integer for the maximum proportion of missingness to simulate (e.g. 0.4 would include simulations at 0.1, 0.2, 0.3. and 0.4 missingness). Default here is 0.4 missingness.
ml <- 0.3

# Missing not at random (MNAR) variable assignment. ---
# If simulating data missing NOT at random (MNAR), uncomment the following lines and adjust according to your needs:
# # Make list of quantile thresholds for removing data for each numerical trait (e.g. removing observations below the 10th quantile). Default here is removing data below the 10th quantile.
# contTraits
# l_quantiles <- as.list(rep(0.1, length(contTraits)))
# names(l_quantiles) <- contTraits
# # Make list of directions (e.g. greater or less than the quantile threshold) for which to introduce NAs for numerical traits. If "less" is chosen, data will be removed below the corresponding quantile. If "greater" is chosen, data will be removed above the corresponding quantile.
# l_directions <- as.list(rep("less", length(contTraits)))
# names(l_directions) <- contTraits
# # Make list indicating whether to use absolute values or not when removing data for each numerical trait (TRUE or FALSE). Default here is to not use absolute values.
# l_abs <- as.list(rep(F, length(contTraits)))
# names(l_abs) <- contTraits
# # For example, if we want to use absolute values for the latitude trait:
# l_abs$latitude_centroid_from_roll_et_al_2017 <- T
# # Make list of categories for which to introduce NAs for categorical traits. For example:
# catTraits
# lapply(dfCC[, catTraits], table)
# l_categories <- list(insular_endemic = "yes", activity_time = "Nocturnal")

# Lastly, get a list of the functions loaded into the global environment.
l_functions <- mget(lsf.str())

### 3. Simulation/imputation. ----

# Set a seed to ensure results can be replicated.
set.seed(547)
# Create vector of missingness patterns you wish to simulate. For example, we can simulate missingness MCAR and MAR:
missPatts <- c("MCAR", "MAR")
# Create empty list to hold results for each missingness pattern.
l_l_simputeResults <- CreateNamedList(listLength = length(missPatts), elementNames = names(missPatts))

# Are you using phylogenetic imputation?
phyImp <- T

# For every missingness pattern..
for(m in 1:length(missPatts)){
  # Take mth missingness pattern.
  MP <- missPatts[[m]]
  # Index to get the corresponding simulation/imputation (simpute) functions.
  index <- grep(paste("Simpute", MP, sep = ""), names(l_functions))
  l_simputeFunctions <- l_functions[index]
  # Create empty list to hold results for each function.
  l_simputeResults <- CreateNamedList(listLength = length(l_simputeFunctions), elementNames = names(l_simputeFunctions))
  # For every simpute function..
  for(f in 1:length(l_simputeResults)){
    # Take fth function.
    simputeFunc <- l_simputeFunctions[[f]]
    # Get name of function.
    funcName <- names(l_simputeFunctions)[[f]]
    # Get name of method.
    impMethod <- gsub(x = funcName, pattern = paste("Simpute", MP, sep = ""), replacement = "")
    # If the missingness pattern is MCAR..
    if(MP == "MCAR"){
      # Simpute without phylogeny.
      results <- simputeFunc(data = dfCC, vars = traits, int = r, missLevel = ml)
      # If the imputation method is MeanMode..
      if(impMethod == "MeanMode"){
        # Apply AverageErrors() function to write average error rates to file.
        AverageErrors(results = results$errorRates, data = dfCC, vars = traits, method = impMethod, missLevel = ml)
        # Append results to l_simputeResults.
        l_simputeResults[[f]] <- results
      } else {
        # Apply AverageErrors() function to write average error rates to file. Set paramTrack = T so we can track parameter values as well.
        AverageErrors(results = results$errorRates, data = dfCC, vars = traits, method = impMethod, missLevel = ml, paramTrack = T)
        # If phylogenetic imputation was chosen..
        if(phyImp == T){
          # Simpute with phylogeny.
          phyResult <- mapply(simputeFunc, tree = l_trees, MoreArgs = list(data = dfCC, vars = traits, int = r, missLevel = ml, phyImp = T), SIMPLIFY = F)
          # Extract error rates from result object.
          errors <- lapply(phyResult, function(x) x$errorRates)
          # Write average error rates to file.
          mapply(AverageErrors, results = errors, treeName = treeNames, MoreArgs = list(data = dfCC, vars = traits, method = impMethod, missLevel = ml, paramTrack = T), SIMPLIFY = F)
          # Combine results.
          results <- list(withoutPhy = results, withPhy = phyResult) 
        }
        # Append results to l_simputeResults.
        l_simputeResults[[f]] <- results
      }
      # If the missingness pattern is MAR..
    } else if(MP == "MAR") {
      # Simpute without phylogeny.
      results <- simputeFunc(data = dfCC, raw = dfRaw, vars = traits, int = r)
      # If the imputation method is MeanMode..
      if(impMethod == "MeanMode"){
        # Apply AverageErrorsBias() function to write average error rates to file.
        AverageErrorsBias(results = results$errorRates, data = dfCC, cont = contTraits, cat = catTraits, method = impMethod, type = MP)
        # Append results to l_simputeResults.
        l_simputeResults[[f]] <- results
      } else {
        # Apply AverageErrorsBias() function to write average error rates to file. Set paramTrack = T so we can track parameter values as well.
        AverageErrorsBias(results = results$errorRates, data = dfCC, cont = contTraits, cat = catTraits, method = impMethod, type = MP, paramTrack = T)
        # If phylogenetic imputation was chosen..
        if(phyImp == T){
          # Simpute with phylogeny.
          phyResult <- mapply(simputeFunc, tree = l_trees, MoreArgs = list(data = dfCC, raw = dfRaw, vars = traits, int = r, phyImp = T), SIMPLIFY = F)
          # Extract error rates from results object.
          errors <- lapply(phyResult, function(x) x$errorRates)
          # Write average error rates to file.
          mapply(AverageErrorsBias, results = errors, treeName = treeNames, MoreArgs = list(data = dfCC, cont = contTraits, cat = catTraits, method = impMethod, type = MP, paramTrack = T), SIMPLIFY = F)
          # Combine results.
          results <- list(withoutPhy = results, withPhy = phyResult) 
        }
        # Append l_results to l_simputeResults.
        l_simputeResults[[f]] <- results
      }
      # If the missingness pattern is MNAR..
    } else if(MP == "MNAR"){
      # Simpute without phylogeny.
      results <- simputeFunc(data = dfCC, vars = traits, int = r, quantiles = l_quantiles, directions = l_directions, absolutes = l_abs, categories = l_categories)
      # If the imputation method is MeanMode..
      if(impMethod == "MeanMode"){
        # Apply AverageErrorsBias() function to write average error rates to file.
        AverageErrorsBias(results = results$errorRates, data = dfCC, cont = contTraits, cat = catTraits, method = impMethod, type = MP)
        # Append results to l_simputeResults.
        l_simputeResults[[f]] <- results
      } else {
        # Apply AverageErrorsBias() function to write average error rates to file. Set paramTrack = T so we can track parameter values as well.
        AverageErrorsBias(results = results$errorRates, data = dfCC, cont = contTraits, cat = catTraits, method = impMethod, type = MP, paramTrack = T)
        # If phylogenetic imputation was chosen..
        if(phyImp == T){
          
          # Simpute with phylogeny.
          phyResult <- mapply(simputeFunc, tree = l_trees, MoreArgs = list(data = dfCC, vars = traits, int = r, quantiles = l_quantiles, directions = l_directions, absolutes = l_abs, categories = l_categories, phyImp = T), SIMPLIFY = F)
          # Extract error rates from results object.
          errors <- lapply(phyResult, function(x) x$errorRates)
          # Write average error rates to file.
          mapply(AverageErrorsBias, results = errors, treeName = treeNames, MoreArgs = list(data = dfCC, cont = contTraits, cat = catTraits, method = impMethod, type = MP, paramTrack = T), SIMPLIFY = F)
          # Combine results.
          results <- list(withoutPhy = results, withPhy = phyResult) 
        }
        # Append l_results to l_simputeResults.
        l_simputeResults[[f]] <- results
      }
    }
  }
  # Append results to l_l_simputeResults.
  l_l_simputeResults[[m]] <- l_simputeResults
}
