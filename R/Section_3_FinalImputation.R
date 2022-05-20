# Section 3. This script imputes raw dataset based on results of simulations in section 2.

# Required input for this script:
# dfCC = Cleaned complete-case dataset from section 1 (CleanedCompleteCaseDataset.csv).
# dfRaw = Cleaned original dataset (with missing values) (CleanedOriginalCaseDataset.csv) from section 1.
# l_trees = Phylogenetic trees in Newick format (if using phylogenetic imputation).
# Error rate csv files for each trait/missing pattern from section 2.
# Best parameter csv files for each trait/missing pattern from section 2.

# Output:
# Box plots comparing distributions of complete, original, and imputed numerical trait data.
# Bar plots comparing categorical frequencies of complete, original, and imputed categorical trait data.
# "DescriptiveStats.csv" = descriptive statistics of complete, original, and imputed numerical trait data.
# "ImputedData.csv" = final imputed dataset using best method

### Acknowledgments. ----

# Referenced these tutorials for customizing barplots and boxplots in R: 
# https://www.r-graph-gallery.com/boxplot.html
# https://www.r-graph-gallery.com/89-box-and-scatter-plot-with-ggplot2.html
# https://towardsdatascience.com/how-to-create-and-customize-bar-plot-using-ggplot2-package-in-r-4872004878a7
# https://stackoverflow.com/questions/52507719/how-to-change-ggplot2-boxplot-color-with-points
# https://stackoverflow.com/questions/37008705/ggplot-bar-chart-of-percentages-over-groups

# Referenced these tutorials for plots and descriptive stats in R: 
# https://statsandr.com/blog/descriptive-statistics-in-r/
# https://statsandr.com/blog/correlation-coefficient-and-correlation-test-in-r/
# http://www.sthda.com/english/wiki/scatter-plot-matrices-r-base-graphs

### 1. Load libraries and functions. ----

library(ape)
library(data.table)
library(lmtest)
library(Information)
library(matrixStats)
library(Metrics)
library(mice)
library(missForest)
library(mltools)
library(naniar)
library(nnet)
library(pastecs)
library(phytools)
library(plotrix)
library(MPSEM)
library(rsample)
library(simputation)
library(tidyverse)
library(VIM)
library(viridis)
source("R/Functions/DataHandling_Functions.R")
source("R/Functions/Imputation_Functions.R")
source("R/Functions/Phylo_Functions.R")
source("R/Functions/Simpute_Functions.R")

### 2. Data loading and variable assignment. ----
# Set the seed.
set.seed(547)

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

# Check dataframes to ensure they only have trait and taxonomy information!
colnames(dfCC)
colnames(dfRaw)

# Variable assignment. ---
# Assign order name.
order <- "Squamata"
# Look at column names of dfCC.
colnames(dfCC)
# Create vector of column names that contain taxonomic information.
taxCols <- c("species_name")
# Extract trait names.
traits <- setdiff(colnames(dfCC), taxCols)
# Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
l_traits <- BreakIntoTypes(dfCC, traits)
# Extract numerical traits.
contTraits <- l_traits[[1]]
# Extract categorical traits.
catTraits <- l_traits[[2]]
# Identify integer (count) traits, if any.
intTraits <- GetTraitNames(data = dfCC[, traits], class = "integer")

# Which missingness patterns were simulated?
missPatts <- c("MCAR", "MAR", "MNAR")

### 3. Handling error rate and parameter files. ----

# Error rate file handling. ---
# Section 2 wrote error files for each trait and imputation method. Now we must read those files into R. I have relocated the error files to a folder called Results/ErrorRates/All/.
# Read in the error rate files as dataframes.
errorFiles <- list.files(path = "Results/ErrorRates/All/", pattern = "_ErrorRates.csv")
l_dfErrors <- lapply(paste("Results/ErrorRates/All/", errorFiles, sep = ""), fread, data.table = F)
# Name according to file.
names(l_dfErrors) <- errorFiles
# Get names of all unique trait/imputation method combos (cleaning up the names here).
allMethodCombos <- unique(gsub(pattern = "_MCAR_|_MAR_|_MNAR_", replacement = "_", errorFiles))
# Remove .csv extension from names.
allMethodCombos <- unique(gsub(pattern = "_ErrorRates.csv", replacement = "", allMethodCombos))
# Create a list to hold all of the error rate dataframes for each method.
l_dfErrorRatesAll <- CreateNamedList(listLength = length(allMethodCombos), elementNames = allMethodCombos)
# For every method..
for(m in 1:length(allMethodCombos)){
  # Take mth method.
  method <- allMethodCombos[[m]]
  # Index to get the corresponding dfErrors for that method.
  index <- grep(pattern = paste(method, "_[A-Z]+_ErrorRates.csv", sep = ""), x = names(l_dfErrors))
  l_dfTraitErrors <- l_dfErrors[index]
  names(l_dfTraitErrors)
  # Rbind dataframes using rbindlist.
  dfMethod <- rbindlist(l_dfTraitErrors)
  # Append to l_dfErrorRatesAll.
  l_dfErrorRatesAll[[m]] <- dfMethod
}

# Parameter file handling. ---
# Repeat same steps to handle the parameter files.
# Read in the parameter files as dataframes.
paramFiles <- list.files(path = "Results/ErrorRates/All/", pattern = "_Parameters.csv")
l_dfParams <- lapply(paste("Results/ErrorRates/All/", paramFiles, sep = ""), fread, data.table = F)
# Name according to file.
names(l_dfParams) <- paramFiles
# Create a list to hold parameter dataframes for each method.
l_dfParamsAll <- CreateNamedList(listLength = length(allMethodCombos), elementNames = allMethodCombos)
# For every method..
for(m in 1:length(allMethodCombos)){
  # Take mth method.
  method <- allMethodCombos[[m]]
  # Index to get the corresponding dfParams for that method.
  index <- grep(pattern = paste(method, "_[A-Z]+_Parameters.csv", sep = ""), x = names(l_dfParams))
  l_dfTraitParams <- l_dfParams[index]
  names(l_dfTraitParams)
  # Rbind dataframes using rbindlist.
  dfMethod <- rbindlist(l_dfTraitParams)
  # Append to l_dfErrorRatesAll.
  l_dfParamsAll[[m]] <- dfMethod
}

### 4. Determining best imputation method and parameter values for the dataset. ----

# Convert categorical traits to factors in both dfCC and dfRaw.
dfRaw[, catTraits] <- lapply(dfRaw[, catTraits], as.factor)
dfCC[, catTraits] <- lapply(dfCC[, catTraits], as.factor)

# Let's determine the best imputation method for each missingness pattern.
l_missPatt <- as.list(l_dfErrorRatesAll[[1]]$missingness_level)
# Create list to hold winner votes for imputation method for each trait and missingness pattern.
l_l_winningImp <- CreateNamedList(listLength = length(l_missPatt), elementNames = l_missPatt)

# For each missingness pattern..
for(m in 1:length(l_missPatt)){
  # Take mth missPatt.
  missPatt <- l_missPatt[[m]]
  # Create list to hold error dataframes based on trait/method combo.
  l_dfTraitError <- CreateNamedList(listLength = length(allMethodCombos), elementNames = allMethodCombos)
  # For every method combo..
  for(c in 1:length(allMethodCombos)){
    # Extract error dataframe.
    dfTraitError <- l_dfErrorRatesAll[[c]]
    # Subset to missingness pattern results.
    dfTraitError <- dfTraitError[dfTraitError$missingness_level == missPatt, ]
    # Rename columns so we can merge with other dfTraitErrors easily.
    names(dfTraitError)[3:4] <- c("error", "SE")
    # Add a trait_method column.
    dfTraitError$trait_method <- names(l_dfErrorRatesAll)[c]
    # Append to lists.
    l_dfTraitError[[c]] <- dfTraitError
  }
  # Create lists to hold winner votes.
  l_winningImp <- CreateNamedList(listLength = length(traits), elementNames = traits)
  # For every trait..
  for(t in 1:length(traits)){
    # Take tth trait.
    trait <- traits[[t]]
    # Index to get the corresponding error dataframes for the trait.
    l_dfSubsetError <- l_dfTraitError[grep(trait, names(l_dfTraitError))]
    # Bind the lists of dataframes.
    dfTraitErrorAll <- bind_rows(l_dfSubsetError)
    # Create trait columns.
    dfTraitErrorAll$trait <- trait
    # Create method columns.
    dfTraitErrorAll$method <- gsub(pattern = paste(trait, "_", sep = ""), "", dfTraitErrorAll$trait_method)
    # Determine the method that resulted in the lowest error rate.
    l_winningImp[[t]] <- dfTraitErrorAll$method[which.min(dfTraitErrorAll$error)]
  }
  # Append to winner lists.
  l_l_winningImp[[m]] <- l_winningImp
}

# Determine the winner for each missingness pattern.
l_BEST <- lapply(l_l_winningImp, function(x) names(sort(table(unlist(x)), decreasing = T))[1])
# The best method overall is...drum roll..!:
optMethod <- names(sort(table(unlist(l_BEST)), decreasing = T))[1]
optMethod
# The best method for MAR (most biologically realistic scenario is)...more drum roll..!:
optMAR <- l_BEST[[which(names(l_BEST) == "MAR")]]
optMAR

# Let's find the parameter values that performed best for our optimal method.

# For every dataframe in l_dfParamsAll..
for(d in 1:length(l_dfParamsAll)){
  # Extract parameter dataframe.
  dfParamsMethod <- l_dfParamsAll[[d]]
  # Add a trait_method column.
  dfParamsMethod$trait_method <- names(l_dfParamsAll)[d]
  # Replace in l_dfParamsAll.
  l_dfParamsAll[[d]] <- dfParamsMethod
}

# Bind all dataframes in l_dfParamsAll.
dfParamsAll <- bind_rows(l_dfParamsAll)

# Subset to only those parameters used for our optimal method.
dfParamsOpt <- dfParamsAll[grepl(optMAR, x = dfParamsAll$trait_method) & dfParamsAll$missingness_level == "MAR", ]
# Create a trait column.
dfParamsOpt$trait <- gsub(pattern = paste("_", optMAR, sep = ""), replacement = "", dfParamsOpt$trait_method)
# Extract parameter values.
bestParams <- dfParamsOpt$param
# Name according to traits.
names(bestParams) <- dfParamsOpt$trait

# For those traits that could not be simulated MAR, select best parameter based on 0.1 MCAR results (e.g. values are less than 10% missing in these cases).
traitsNOMAR <- setdiff(traits, names(bestParams))
# Subset to MCAR and MNAR results.
dfParamsNOMAR <- dfParamsAll[grepl(optMAR, x = dfParamsAll$trait_method) & dfParamsAll$missingness_level == "0.1", ]
# Create a trait column.
dfParamsNOMAR$trait <- gsub(pattern = paste("_", optMAR, sep = ""), replacement = "", dfParamsNOMAR$trait_method)
# Subset to traitsNOMAR.
dfParamsNOMAR <- dfParamsNOMAR[dfParamsNOMAR$trait %in% traitsNOMAR, ]
# Extract parameter values.
bestParamsNOMAR <- dfParamsNOMAR$param
# Name according to traits.
names(bestParamsNOMAR) <- dfParamsNOMAR$trait
# Append to bestParams.
bestParams <- c(bestParams, bestParamsNOMAR)


### 5. Imputation prep. ----

# Here, if you want to proceed with phylogenetic imputation method, uncomment the following lines.
# Read in the tree you wish to use. For example:
#ultraTree <- read.tree("Data/RAxMLTrees/FinalTrees/Squamata_CMOS_ultraTree_finalImp.tre")
# Make sure species names in original dataset and tree match.
#ultraTree <- drop.tip(phy = ultraTree, tip = ultraTree$tip.label[!ultraTree$tip.label %in% dfRaw$species_name])
# Make sure the species data in dfRaw match the tip labels in ultraTree.
#dfRaw <- dfRaw[match(ultraTree$tip.label, dfRaw$species_name), ]
#all.equal(ultraTree$tip.label, dfRaw$species_name)

# Let's check the distributions of the data and for class imbalances in categorical data.
# Descriptive info for numerical data:
contRes <- lapply(dfRaw[, contTraits], GetNumericalInfo)
# Descriptive info for categorical data:
catRes <- lapply(dfRaw[, catTraits], GetCategoricalInfo)

# For every continuous trait..
for(cont in 1:length(contRes)){
  # Print the name of the trait.
  cat("Trait:", names((contRes))[[cont]])
  # Extract the result.
  result <- contRes[[cont]]
  # Print the results.
  print(result[1:6])
  # Plot the histogram.
  plot(result$Plot, col = "skyblue", main = paste(names((contRes))[[cont]], sep = " "))
}

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
  barplot(height = dfCount$Freq, names.arg = dfCount$vector, col = "skyblue", main = paste(names((catRes))[[cat]], sep = " "))
}

# Selecting predictors. ---
# Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait).
# Apply the SelectPredictors() function to obtain a list of predictors for each trait.
l_predictors <- SelectPredictors(dfRaw[, traits])

# Phylogenetic imputation prep. ---
# Set phylogenetic imputation to true if proceeding (default is F):
phyImp <- F
# If using phylogenetic imputation, uncomment the following lines:
# Append eigenvectors to dfRaw and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors. Note: may take a while depending on size of dataset.
l_evs <- AppendEigenvectors(data = dfRaw, vars = traits, tree = ultraTree, predictors = l_predictors)
# Extract list of dataframes (one dataset for each trait with appended eigenvectors).
l_dfRawEV <- l_evs[[1]]
# Extract updated list of predictors (with appended eigenvectors).
l_EVPredictors <- l_evs[[2]]


### 6. Final imputation. ----

# Replicate dfRaw so we can impute the missing values in a new dataset.
dfRawImp <- dfRaw

# For each trait..
for(t in 1:length(traits)) {
  # Take the name of the tth trait.
  trait <- traits[[t]]
  # If phylogenetic imputation was selected...
  if(phyImp == T){
    # Take the corresponding predictors with appended eigenvectors for that trait.
    preds <- l_EVPredictors[[grep(trait, names(l_EVPredictors))]]
  } else {
    # Take the corresponding predictors for that trait.
    preds <- l_predictors[[grep(trait, names(l_predictors))]]
  }
  # Take the best parameter for the trait.
  param <- bestParams[[grep(trait, names(bestParams))]]
  # If phylogenetic imputation was selected...
  if(phyImp == T){
    # Take the tth dfRawEV.
    dfRawMiss <- l_dfRawEV[[grep(trait, names(l_EVPredictors))]] 
  } else {
    # Make a copy of dfRaw to impute.
    dfRawMiss <- dfRaw
  }
  # If optimal method is KNN..
  if(grepl("KNN", optMAR)){
    # Impute using kNN and optimal value of k.
    dfImputed <- kNN(dfRawMiss, variable = trait, dist_var = preds, k = param)
    # If optimal method is RF..
  } else if(grepl("RF", optMAR)){
    # Impute the dataset (which contains the trait in question and its predictors) using missForest and optimal value of ntree.
    imputedRF <- missForest(as.data.frame(dfRawMiss[, c(trait, preds)]), ntree = param)
    # Access the imputed dataframe.
    dfImputed <- imputedRF$ximp
    # Round the imputed count data.
    dfImputed[, intTraits] <- lapply(dfImputed[, intTraits], function(x) as.integer(round(x)))
    # If optimal method is MICE..
  } else if(grepl("MICE", optMAR)){
    
    # If phylogenetic imputation was selected..
    if(phyImp == T){
      # MICE requires a predictor matrix for imputation. So, we will first merge all of the dataframe in l_dfRawEV into one dataframe. This is so we can identify the eigenvector columns needed for imputation of each trait and update this in the predictor matrix. This should also reduce computation time since we are just imputing one dataframe and not a list of dataframes. 
      dfMissEV <- l_dfRawEV[[1]]
      # Identify the eigenvector columns and remove those from dfMissEV.
      evIndex <- grep(pattern = "V_", colnames(dfMissEV))
      dfMissEV <- dfMissEV[, -evIndex]
      # Now, let's merge dfMissEV with the corresponding eigenvector columns in each dataframe in l_dfRawEV.
      for(e in 1:length(l_dfRawEV)){
        # Take the eth l_dfRawEV.
        dfMisseth <- l_dfRawEV[[e]]
        # Identify eigenvector columns in eth dataframe.
        index <- grep("V_", colnames(dfMisseth))
        evCols <- colnames(dfMisseth)[index]
        # Merge these columns with dfMisseth (including species_name).
        dfMissEV <- merge(dfMissEV, dfMisseth[, c("species_name", evCols)], by = "species_name")
      }
      # Create predictor matrix for use in MICE imputation with appended eigenvectors. We can indicate the column names of dfMissEV (excluding species_name) as names of variables to consider in imputation process. May take a while depending on size of dataset.
      matPredictors <- CreatePredictorMatrix(dfMissing = dfMissEV, cols = colnames(dfMissEV)[-1], predictors = l_EVPredictors)
      # Ensure all of the rows for the eigenvectors to 0 as they themselves will not be imputed.
      eigens <- grep(pattern = "V_", rownames(matPredictors))
      matPredictors[eigens, ] <- 0
      # Impute the datasets using optimal m and matPredictors (can take a while depending on size of dataset). May take a while depending on size of dataset.
      imputedMICE <- mice(dfMissEV[, c(colnames(matPredictors))], predictorMatrix = matPredictors, m = param, maxit = 10, print = FALSE)
      # Combine the imputed data into one dataframe.
      dfImputed <- CombineMIDataframes(midata = imputedMICE, method = "MICE", m = param, contVars = contTraits, catVars = catTraits)
    } else {
      # Create predictor matrix for use in MICE imputation (based on results of our trait predictor screening).
      matPredictors <- CreatePredictorMatrix(cols = traits, predictors = l_predictors, dfMissing = dfRawMiss)
      # Impute the datasets using optimal m and matPredictors (can take a while depending on size of dataset).
      imputedMICE <- mice(dfRawMiss[, c(colnames(matPredictors))], predictorMatrix = matPredictors, m = param, maxit = 10, print = FALSE)
      # Combine the imputed data into one dataframe.
      dfImputed <- CombineMIDataframes(midata = imputedMICE, method = "MICE", m = param, contVars = contTraits, catVars = catTraits)
    }
  }
  # Replace the data for the trait in dfRawImp with the newly imputed data in dfImputed.
  dfRawImp[, trait] <- dfImputed[, trait]
}

# Ensure dataframes are in the same order.
dfCC <- dfCC[c("species_name", traits)]
dfRaw <- dfRaw[c("species_name", traits)]
dfRawImp <- dfRawImp[c("species_name", traits)]
# Combine dataframes into list.
l_dfAll <- list(dfCC, dfRaw, dfRawImp)
# Name the list according to df.
names(l_dfAll) <- c("dfCC", "dfRaw", "dfRawImp")
# Now we want to edit the column names in each dataframe so we can easily identify which dataset the trait data are from. First, create list of strings to append to each col name (in order of l_dfAll).
colStrings <- list("cc", "raw", "raw_imp")
# Now apply PasteColNames() using mapply. Because we wrapped traits in a list() it is recycled for each iteration (our constants). Specified SIMPLIFY = F so it doesn't reduce results to a matrix.
l_dfAll <- mapply(PasteColNames, df = l_dfAll, colNames = list(traits), string = colStrings, sep = "_", SIMPLIFY = F)
# Merge all the dataframes by species name.
dfAll <- Reduce(function(...) merge(..., by = "species_name", all = T), l_dfAll)

### 7. Plots. ----

# In this section, we are plotting box plots for each continuous trait and bar plots for each categorical trait.
# Select the traits we want to plot (traits that actually had data imputed in dfRaw).
impTraits <- c("activity_time", "smallest_clutch", "largest_clutch", "female_svl")

# For every trait...
for(t in 1:length(impTraits)){
  # Take the tth trait.
  trait <- impTraits[[t]]
  # Get the corresponding column indices for the trait.
  index <- grep(pattern = trait, x = colnames(dfAll))
  # Subset out the columns.
  dfSubset <- dfAll[, index]
  # Create plot title.
  plotTitle <- paste("Final_Imputation", trait, ".tiff", sep = "_")
  # If the trait is continuous...
  if(trait %in% contTraits) {
    # Log transform data.
    dfSubset[, c(1:3)] <- lapply(dfSubset[, c(1:3)], log)
    # Get sample size counts for each column.
    sampleSizes <- sapply(dfSubset[, 1:3], function(x) sum(!is.na(x)))
    # Pivot dataframe to long form so we can more easily plot variables by group.
    dfPivot <- pivot_longer(dfSubset[, c(1:3)], cols = colnames(dfSubset[, c(1:3)]))
    # Convert to factor.
    dfPivot$name <- as.factor(dfPivot$name)
    # X-axis label containing sample size information.
    dataType <- c("Complete case", "Original data", "Imputed original data")
    # Paste sample size onto dataType vector.
    xlabel <- paste(dataType, "\n", "(n = ", sampleSizes, ")" , sep = "")
    
    # Boxplot.
    tiff(plotTitle, units = "in", width = 14, height = 8, res = 600)
    plot <- ggplot(dfPivot, aes(x = name, y = value, fill = name)) +
      geom_boxplot() +
      scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
      geom_point(position = position_jitterdodge(jitter.width = 2, dodge.width = 0), 
                 pch = 21, aes(fill = name, show.legend = F)) +
      labs(title = trait,
           y = "", x = "")  +
      scale_x_discrete(labels = xlabel) + theme_minimal() +
      theme(legend.position = "none") +
      theme(axis.text = element_text(size = 20, face = "bold"))
    print(plot)
    dev.off()
  } else if(trait %in% catTraits) {
    # Get sample size counts for each column.
    sampleSizes <- sapply(dfSubset, function(x) sum(!is.na(x)))
    # Create dataframe containing category counts for each group.
    dfPivot <- pivot_longer(dfSubset, cols = colnames(dfSubset))
    # Convert to factor.
    dfPivot$name <- as.factor(dfPivot$name)
    # Group by name and value.
    dfCount <- group_by(na.omit(dfPivot), name, value)
    # Get counts.
    dfCount <- summarise(dfCount, count = n())
    # Create proportion column.
    dfCount <- mutate(dfCount, prop = count/sum(count))
    # X-axis label containing sample size information.
    dataType <- c("Complete case", "Original data", "Imputed original data")
    # Paste sample size onto dataType vector.
    xlabel <- paste(dataType, "\n", "(n = ", sampleSizes, ")" , sep = "")
    
    # Plot.
    tiff(plotTitle, units = "in", width = 14, height = 8, res = 600)
    plot <- ggplot(data = dfCount, mapping = aes(x = name, y = prop, fill = value)) + 
      geom_bar(stat="identity", position = "dodge") + theme_minimal() +
      scale_fill_viridis(discrete = TRUE) +
      scale_x_discrete(labels = xlabel) +
      geom_text(aes(label = scales::percent(prop, accuracy = 0.1)), vjust = -.5, position = position_dodge(0.9), size = 6) +
      theme(axis.text=element_text(size = 20, face = "bold")) +
      labs(title = trait,
           y = "", x = "")
    print(plot)
    dev.off()
  }
}

### 8. Descriptive statistics. ----
contCols <- sapply(dfAll, is.numeric)
# Remove categorical as we are obtaining descriptive stats for numerical traits.
dfContSubset <- dfAll[, contCols]
# Apply stat.desc to dfContSubset to get some descriptive stats for each trait and dataset.
dfDescCont <- stat.desc(dfContSubset)
# Order alphabetically.
dfDescCont <- dfDescCont[, order(colnames(dfDescCont))]
# Write to file.
write.csv(dfDescCont, "DescriptiveStats.csv")

# Finally, write imputed dataframe to file.
#write.csv(dfRawImp, "ImputedData.csv")
