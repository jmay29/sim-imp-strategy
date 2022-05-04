# Section 2. This script imputes trait data using mean/mode replacement, k-nearest neighbour (KNN), missForest (RF), and multivariate imputation using chained equations (MICE) with and without phylogeny (in the form of phylogenetic eigenvectors derived from a phylogenetic tree). Given a complete-case dataset, it simulates data missing completely at random (MCAR), missing at random (MAR), and missingness not at random (MNAR). It then uses different imputation methods to fill in the missing values using the known observations. The outputs for this script are the imputation error rates for each trait in the dataset.

# Acknowledgments. ----

# This script makes use of the following imputation functions:

# kNN() function in the "VIM" package for imputation.
# Citations: Kowarik A, Templ M (2016). “Imputation with the R Package VIM.” Journal of Statistical Software, 74(7), 1–16. doi: 10.18637/jss.v074.i07.
# https://cran.r-project.org/web/packages/VIM/VIM.pdf

# missForest() function in the "missForest" package for imputation.
# Citations: Daniel J. Stekhoven (2013). missForest: Nonparametric Missing Value Imputation using Random Forest. R package version 1.4.
# Stekhoven D. J., & Buehlmann, P. (2012). MissForest - non-parametric missing value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.

# mice() function in the "mice" package for imputation.
# Citations: van Buuren S, Groothuis-Oudshoorn K (2011). “mice: Multivariate Imputation by Chained Equations in R.” Journal of Statistical Software, 45(3), 1-67. https://www.jstatsoft.org/v45/i03/.
# R package version 3.13.0. https://cran.r-project.org/web/packages/mice/mice.pdf

### 1. Load libraries and functions. ----

library(ape)
library(data.table)
library(Information)
library(lmtest)
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

### 2. Data loading and variable assignment. ----

# Read in cleaned complete-case trait dataset.
fileName <- file.choose()
dfCC <- read.csv(fileName)
# Ensure blanks are NAs.
dfCC[dfCC == ""] <- NA
# Read in the cleaned original trait dataset.
fileName <- file.choose()
dfRaw <- read.csv(fileName)
# Ensure blanks are NAs.
dfRaw[dfRaw == ""] <- NA

# If using phylogenetic imputation, uncomment the following lines.
# Read in trees. For example (these trees have all been matched to the trait data previously):
# treeFiles <- list.files(path = "Data/RAxMLTrees/", pattern = "_ultra.tre")
# l_trees <- lapply(treeFiles, function(x) read.tree(paste("Data/RAxMLTrees/", x, sep = "")))
# # Name l_trees.
# names(l_trees) <- treeFiles
# # Create vector of tree names (used to name error rate files later on).
# treeNames <- sapply(treeFiles, function(x) word(x, start = 2, sep = "_"))

# Extract trait names.
traits <- setdiff(colnames(dfCC), "species_name")
# Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
l_traits <- BreakIntoTypes(dfCC, traits)
# Extract numerical traits.
contTraits <- l_traits[[1]]
# Extract categorical traits.
catTraits <- l_traits[[2]]

# Create an integer for the number of replicates to use for missingness simulations and imputation. Default here is 10 replicates.
r <- 10
# If simulating data missing completely at random (MCAR), create an integer for the maximum proportion of missingness to simulate (e.g. 0.4 would include simulations at 0.1, 0.2, 0.3. and 0.4 missingness). Default here is 0.4 missingness.
ml <- 0.4

# Missing not at random (MNAR) variable assignment. ---
# If simulating data missing NOT at random (MNAR)..
# Make list of quantile thresholds for removing data for each numerical trait (e.g. removing observations below the 10th quantile). Default here is removing data below the 10th quantile.
l_quantiles <- as.list(rep(0.1, length(contTraits)))
names(l_quantiles) <- contTraits
# Make list of directions (e.g. greater or less than the quantile threshold) for which to introduce NAs for numerical traits. If "less" is chosen, data will be removed below the corresponding quantile. If "greater" is chosen, data will be removed above the corresponding quantile.
l_directions <- as.list(rep("less", length(contTraits)))
names(l_directions) <- contTraits
# Make list indicating whether to use absolute values or not when removing data for each numerical trait (TRUE or FALSE). Default here is to not use absolute values.
l_abs <- as.list(rep(F, length(contTraits)))
names(l_abs) <- contTraits
# For example, if we want to use absolute values for the latitude trait:
l_abs$latitude_centroid_from_roll_et_al_2017 <- T
# Make list of categories for which to introduce NAs for categorical traits. For example:
lapply(dfCC[, catTraits], table)
l_categories <- list(insular_endemic = "yes", activity_time = "Nocturnal")

# Lastly, get a list of the functions loaded into the global environment.
l_functions <- mget(lsf.str())

### 3. Simulation/imputation. ----

# Set a seed to ensure results can be replicated.
set.seed(547)
# Create vector of missingness patterns.
missPatts <- c("MCAR", "MAR", "MNAR")
# Create empty list to hold results for each missingness pattern.
l_l_simputeResults <- CreateNamedList(listLength = length(missPatts), elementNames = names(missPatts))

# Are you using phylogenetic imputation? Default is false.
phyImp <- F

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

