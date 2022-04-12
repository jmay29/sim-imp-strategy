# Section 4. This script creates the plots used in figures 3-5.

### Acknowledgments. ----

# The following tutorials were consulted:
# http://5.9.10.113/67164758/adding-hatches-or-patterns-to-ggplot-bars
# https://coolbutuseless.github.io/package/ggpattern/articles/geom-gallery-geometry.html#colour-example-1
# http://www.sthda.com/english/wiki/ggplot2-point-shapes
# https://towardsdatascience.com/how-to-create-and-customize-bar-plot-using-ggplot2-package-in-r-4872004878a7

### 1. Load libraries and functions. ----

library(ape)
library(caper)
library(car)
library(colorspace)
library(correlation)
library(corrplot)
library(data.table)
library(fastDummies)
library(ggpattern)
library(ggpol)
library(Hmisc)
library(Information)
library(lattice)
library(lmtest)
library(matrixStats)
library(memisc)
library(Metrics)
library(nnet)
library(pastecs)
library(phytools)
library(plotly)
library(plotrix)
library(psych)
library(rsample)
library(summarytools)
library(tidyverse)
library(VIM)
library(viridis)
source("R/Functions/DataHandling_Functions.R")
source("R/Functions/Imputation_Functions.R")
source("R/Functions/Phylo_Functions.R")
source("R/Functions/Plot_Functions.R")

### 2. Data loading and variable assignment. ----

# Read in cleaned complete-case trait dataset.
dfCC <- fread("Data/TraitData/CleanedCompleteCaseDataset.csv", data.table = F)
# Ensure blanks are NAs.
dfCC[dfCC == ""] <- NA
# Read in the cleaned original trait dataset.
dfRaw <- fread("Data/TraitData/CleanedOriginalCaseDataset.csv", data.table = F)
# Ensure blanks are NAs.
dfRaw[dfRaw == ""] <- NA

# Read in ultrametric gene trees.
treeFiles <- list.files(path = "Data/RAxMLTrees/", pattern = "_ultra.tre")
l_trees <- lapply(treeFiles, function(x) read.tree(paste("Data/RAxMLTrees/", x, sep = "")))
# Name l_trees.
names(l_trees) <- treeFiles
# Create vector of tree names (used to name error rate files later on).
treeNames <- sapply(treeFiles, function(x) word(x, start = 2, sep = "_"))

# Variable assignment. ---

# Assign order name.
order <- "Squamata"

# Create vector of column names that contain taxonomic information.
taxCols <- c("species_name", "class", "order", "family", "genus", "species")
# Extract trait names.
traits <- setdiff(colnames(dfCC), taxCols)
# Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
l_traits <- BreakIntoTypes(dfCC, traits)
# Extract numerical traits.
contTraits <- l_traits[[1]]
# Extract categorical traits.
catTraits <- l_traits[[2]]
# Identify integer (count) traits, if any.
intTraits <- GetTraitNames(df = dfCC[, traits], class = "integer")

# Assign phylogenetic signal measure to use (one of "K" or "lambda")
phyloSignal <- "lambda"

### 3. MCAR line plots without phylogeny. ----

# Section 2 wrote error files for each trait and imputation method. I have relocated the error files to a folder called Results/ErrorRates/All.

# Read in the error rate files as dataframes.
# Identify file names.
errorFiles <- list.files(path = "Results/ErrorRates/All/", pattern = "_ErrorRates.csv")
# Get MCAR results.
index <- grep("MCAR", x = errorFiles)
errorFiles <- errorFiles[index]
# Read in error files.
l_dfErrors <- lapply(paste("Results/ErrorRates/All/", errorFiles, sep = ""), fread, data.table = F)
# Name according to file.
names(l_dfErrors) <- errorFiles
  
# For each trait...
for(t in 1:length(traits)) {
  # Take the tth trait.
  trait <- traits[[t]]
  # Index to get the corresponding dfErrors for the trait.
  index <- grep(trait, names(l_dfErrors))
  l_dfTraitErrors <- l_dfErrors[index]
  names(l_dfTraitErrors)
  # Merge all of the dataframes using Reduce() (will get warning messsages here but we can ignore those).
  dfErrorTrait <- Reduce(function(...) merge(..., by = "missingness_level", all = T), l_dfTraitErrors)
  # Remove extraneous cols.
  index <- grep("V1.*", names(dfErrorTrait))
  dfErrorTrait <- dfErrorTrait[, -index]
  # Create a list of strings that we want to remove from the column names.
  l_rm <- list("_MCAR_ErrorRates.csv", trait, "^_")
  # Create a vector to hold the original column names.
  colNames <- names(l_dfTraitErrors)
  # Iterate over the lists, updating the column names with each pass.
  for(n in 1:length(l_rm)) {
    # Apply gsub function using the nth elements of the lists.
    colNames <- gsub(pattern = l_rm[[n]], replacement = "", x = colNames)
  }
  # As we also need to add columns to hold SE values, let's replicate the column names.
  newCols <- rep(colNames, each = 2)
  # Take every 2nd element.
  indices <- seq(2, length(newCols), 2)
  # Paste "_SE" onto every 2nd element.
  newCols[indices] <- paste(newCols[indices], "_SE", sep = "")
  # Replace in column names with newCols, starting at the 2nd column (because the first column is missingness_level).
  names(dfErrorTrait)[2:ncol(dfErrorTrait)] <- newCols
  names(dfErrorTrait)
    
  # As we want all of the plots to be on the same scales, extract highest error rates for continuous and categorical traits.
  maxError <- max(dfErrorTrait[, colNames])
  # Which row and column contains the maxError rate? To get the associated SE that we need for plotting, get row index and the column next to the one with maxError.
  rowIndex <- which(dfErrorTrait == maxError, arr.ind = T)[1]
  colIndex <- which(dfErrorTrait == maxError, arr.ind = T)[2] + 1
  # Extract SE for maxError.
  maxSE <- dfErrorTrait[rowIndex, colIndex]
  # Extract lowest error rates for continuous and categorical traits.
  minError <- min(dfErrorTrait[, colNames])
  # Find associated SE for minError.
  rowIndex <- which(dfErrorTrait == minError, arr.ind = T)[1]
  colIndex <- which(dfErrorTrait == minError, arr.ind = T)[2] + 1
  # Extract SE for minError.
  minSE <- dfErrorTrait[rowIndex, colIndex]
  # Add SE to maxError and minError.
  upper <- maxError + maxSE
  lower <- minError - minSE
 
  # Create title for the plots.
  plotTitle <- paste(order, trait, sep = "_")
  # If the trait is continuous...
  if(trait %in% contTraits) {
    # Apply the corresponding function for the taxon, specifying traitType = "numeric".
    PlotLine(df = dfErrorTrait, traitType = "numeric", plotType = "between", minError = lower, maxError = upper, title = plotTitle)
  } else if(trait %in% catTraits) {
    # Apply the corresponding function for the taxon, specifying traitType = "categorical".
    PlotLine(dfErrorTrait, traitType = "categorical", plotType = "between", minError = lower, maxError = upper, title = plotTitle)
  }
}

### 4. Bar plots. ----

# Read in the error rate files as dataframes.
# Identify file names.
errorFiles <- list.files(path = "Results/ErrorRates/All/", pattern = "_ErrorRates.csv")
# Read in error files.
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

# Create a list to hold plots.
l_l_plots <- vector(mode = "list", length = length(traits))

# For each trait...
for(t in 1:length(traits)) {
  # Take the tth trait.
  trait <- traits[[t]]
  # Index to get the corresponding dfErrors for the trait.
  index <- grep(trait, names(l_dfErrorRatesAll))
  l_dfTraitErrors <- l_dfErrorRatesAll[index]
  names(l_dfTraitErrors)
  # Convert to dataframe.
  l_dfTraitErrors <- lapply(l_dfTraitErrors, as.data.frame)
  
  # Merge all of the dataframes using Reduce(). suppressWarnings is used because column names are duplicated but we will fix this problem later on.
  dfErrorTrait <- suppressWarnings(Reduce(function(...) merge(..., by = "missingness_level", all = T), l_dfTraitErrors))
  # Remove extraneous cols.
  index <- grep("V1.*", names(dfErrorTrait))
  dfErrorTrait <- dfErrorTrait[, -index]
  # Create a list of strings that we want to remove from the column names.
  l_rm <- list(trait, "^_")
  # Create a vector to hold the original column names.
  colNames <- names(l_dfTraitErrors)
  # Iterate over the lists, updating the column names with each pass.
  for(n in 1:length(l_rm)) {
    # Apply gsub function using the nth elements of the lists.
    colNames <- gsub(pattern = l_rm[[n]], replacement = "", x = colNames)
  }
  # As we also need to add columns to hold SE values, let's replicate the column names.
  newCols <- rep(colNames, each = 2)
  # Take every 2nd element.
  indices <- seq(2, length(newCols), 2)
  # Paste "_SE" onto every 2nd element.
  newCols[indices] <- paste(newCols[indices], "_SE", sep = "")
  # Replace in column names with newCols, starting at the 2nd column (because the first column is missingness_level).
  names(dfErrorTrait)[2:ncol(dfErrorTrait)] <- newCols
  names(dfErrorTrait)
  
  # As we want all of the plots to be on the same scales, extract highest error rates for continuous and categorical traits.
  maxError <- max(dfErrorTrait[, colNames], na.rm = T)
  # Which row and column contains the maxError rate? To get the associated SE that we need for plotting, get row index and the column next to the one with maxError.
  rowIndex <- which(dfErrorTrait == maxError, arr.ind = T)[1]
  colIndex <- which(dfErrorTrait == maxError, arr.ind = T)[2] + 1
  # Extract SE for maxError.
  maxSE <- dfErrorTrait[rowIndex, colIndex]
  # Add SE to maxError.
  upper <- maxError + maxSE
  
  # Create list to hold plots.
  l_plots <- vector(mode = "list", length = length(dfErrorTrait$missingness_level))
  
  # For every missingness level...
  for(m in 1:length(dfErrorTrait$missingness_level)) {
    # If the trait is continuous...
    if(trait %in% contTraits) {
      # Apply PlotBar function, specifying traitType = "numeric".
      plot <- PlotBar(df = dfErrorTrait, traitType = "numeric", missLevel = m, title = trait, upperThresh = upper)
      # If the trait is categorical...
    } else if(trait %in% catTraits) {
      # Apply PlotBar function, specifying traitType = "categorical".
      plot <- PlotBar(dfErrorTrait, traitType = "categorical", missLevel = m, title = trait, upperThresh = upper)
    }
    # Append to l_plots.
    l_plots[[m]] <- plot
  }
  # Append to l_l_plots.
  l_l_plots[[t]] <- l_plots
}

### 5. Phylogenetic signal measurement & visualizations. ----

# Create a list to hold the phylogenetic signal estimates for each gene tree.
l_l_sig <- vector(mode = "list", length = length(l_trees))
# Name according to gene tree.
names(l_l_sig) <- names(l_trees)

# For every gene tree..
for(g in 1:length(l_trees)) {
  # Take the gth gene tree.
  tree <- l_trees[[g]]
  # Make sure the trait data and tree tips match.
  tree <- drop.tip(phy = tree, tip = tree$tip.label[!tree$tip.label %in% dfCC$species_name])
  # Make sure the order matches between the phylogeny and dataframe.
  dfPhySig <- dfCC[match(tree$tip.label, dfCC$species_name), ]
  # Convert categorical traits to factors.
  dfPhySig[, catTraits] <- lapply(dfPhySig[, catTraits], as.factor)
  # Identify multi-categorical traits as we need to make dummy variables for these traits (since D metric can only be measured on binary traits).
  multicat <- sapply(dfPhySig[, catTraits], function(x) length(levels(x)))
  # Name the vector according to catTraits.
  names(multicat) <- catTraits
  # Identify the names of the binary traits.
  binTraits <- names(which(multicat == 2))
  # Identify the names of the multi-categorical traits.
  multicat <- names(which(multicat > 2))
  # If there are multi-categorical traits...
  if(length(multicat) > 0) {
    # Use the function dummy_cols to create dummy columns to replace these traits.
    dfPhySig <- dummy_cols(dfPhySig, select_columns = multicat, remove_selected_columns = T, ignore_na = T)
    # Update traits column.
    traitsPS <- setdiff(colnames(dfPhySig), "species_name")
    # Identify factor traits.
    factorTraits <- sapply(dfPhySig[, traitsPS], is.factor)
    # Get the names of the columns that are TRUE.
    factorTraits <- names(which(factorTraits == TRUE))
    # Identify binary traits (omitting NAs) and append them to factorTraits as they will also be treated as factors.
    binTraits <- GetTraitNames(na.omit(dfPhySig[, traitsPS]), class = "binary")
    newCatTraits <- c(factorTraits, binTraits)
    # Make sure all binary traits are factor type.
    dfPhySig[, newCatTraits] <- lapply(dfPhySig[, newCatTraits], as.factor)
  }
  # Identify continuous traits.
  contTraits <- GetTraitNames(dfPhySig[, traitsPS], class = "numeric")
  # Create a list to hold the measures of phylogenetic signal.
  l_sig <- vector(mode = "list", length = length(traitsPS))
  # Name according to traitsPS.
  names(l_sig) <- traitsPS
  # For every trait...
  for(t in 1:length(traitsPS)) {
    # Get the name of the tth trait.
    trait <- traitsPS[[t]]
    # Apply TestPhyloSig function.
    sig <- TestPhyloSig(df = dfPhySig, colName = trait, phylo = tree, sig = phyloSignal)
    # Append to l_sig.
    l_sig[[t]] <- sig
  }
  # Append to l_l_sig.
  l_l_sig[[g]] <- l_sig
}

# Reformatting into dataframe for easier plotting.
# Create empty dataframe to store results.
dfMetrics <- data.frame(gene = character(), trait = character(), phylogenetic_signal= numeric(), p_value = numeric(), stringsAsFactors = FALSE) 

# For every element in l_l_sig..
for(g in 1:length(l_l_sig)){
  # Take name of gene.
  gene <- word(names(l_l_sig)[[g]], 2, sep = "_")
  # Take gth element in l_l_sig.
  l_sig <- l_l_sig[[g]]
  # For every trait in traitsPS..
  for(t in 1:length(traitsPS)){
    # Take the name of the trait.
    trait <- traitsPS[[t]]
    # Take the corresponding element of l_sig.
    sig <- l_sig[[grep(trait, names(l_sig))]]
    # If the trait is numeric..
    if(trait %in% contTraits){
      # Assign phyloSignal to metric.
      metric <- phyloSignal
      # If lambda was selected..
      if(metric == "lambda"){
        # Extract phylogenetic signal.
        value <- sig$lambda
        # If K was selected..
      } else if(metric == "K"){
        # Extract phylogenetic signal.
        value <- sig$K
      }
      # Extract p-value.
      pVal <- sig$P
      # Else if trait is in newCatTraits (multi cats broken down into binary)..
    } else if(trait %in% newCatTraits){
      # Assign D to metric.
      metric <- "D"
      # Extract phylogenetic signal.
      value <- sig$DEstimate
      # Extract p-value (Probability of E(D) resulting from no (random) phylogenetic structure).
      pVal <- sig$Pval1
    }
    # Combine results.
    result <- list(gene = gene, trait = trait, metric = metric, phylogenetic_signal = value, p_value = pVal)
    # Bind to dfMetrics.
    dfMetrics <- rbind(dfMetrics, result)
  }
}


# Plotting phylogenetic signal. ---
# Convert grouping variables to factor type. Here, I am also ensuring correct plotting order for the different variables. 
dfMetrics$gene <- factor(dfMetrics$gene, levels = c("COI", "CMOS", "RAG1"))
# Split dataframe into continuous and categorical dataframes.
# Continuous:
dfCont <- dfMetrics[dfMetrics$metric == phyloSignal, ]
# Replace latitude name.
dfCont$trait[dfCont$trait == "latitude_centroid_from_roll_et_al_2017"] <- "latitude"
# Categorical:
dfCat <- dfMetrics[dfMetrics$metric == "D", ]
# Convert trait variables to factor type and enter correct plotting order.
dfCont$trait <- factor(dfCont$trait, levels = c("largest_clutch", "smallest_clutch", "female_svl", "maximum_svl", "latitude"))
dfCat$trait <- factor(dfCat$trait, levels = c("activity_time_Diurnal", "activity_time_Nocturnal", "insular_endemic"))
# Remove NAs (resulting from extraneous category for multilevel traits).
dfCat <- na.omit(dfCat)
# Write results to file.
#write_csv(dfCont, file = "ContPhySig.csv")
#write_csv(dfCat, file = "CatPhySig.csv")

# Combine dataframes into list.
l_dfMetrics <- list(dfCont, dfCat)
# Create vector of plot titles.
plotTitles <- list("continuous", "categorical")
# Apply PlotPhySig function to list of dataframes.
mapply(PlotPhySig, df = l_dfMetrics, title = plotTitles)


# Plotting phylogenetic signal vs. error rate. ---

# Split dfMetrics by trait.
l_dfMetricTrait <- split(dfMetrics, dfMetrics$trait)
# Create empty list to store the dataframe that contains error ratios and phylogenetic signal for each trait.
l_dfTraitPhySig <- vector(mode = "list", length = length(traits))
# Name according to traits.
names(l_dfTraitPhySig) <- traits

# For each trait...
for(t in 1:length(traits)) {
  # Take the tth trait.
  trait <- traits[[t]]
  # Index to get the corresponding ErrorRatesAll for the trait (as these also include MAR and MNAR results).
  index <- grep(trait, names(l_dfErrorRatesAll))
  subset <- l_dfErrorRatesAll[index]
  names(subset)
  # Convert to data.frame. format.
  subset <- lapply(subset, as.data.frame)
  # Merge all of the dataframes using Reduce(). suppressWarnings is used because column names are duplicated but we will fix this problem later on.
  dfErrorTrait <- suppressWarnings(Reduce(function(...) merge(..., by = "missingness_level", all = T), subset))
  # Remove extraneous cols.
  index <- grep("V1.*", names(dfErrorTrait))
  dfErrorTrait <- dfErrorTrait[, -index]
  # Create a list of strings that we want to remove from the column names.
  l_rm <- list(trait, "^_")
  # Create a vector to hold the original column names.
  colNames <- names(subset)
  # Iterate over the lists, updating the column names with each pass.
  for(n in 1:length(l_rm)) {
    # Apply gsub function using the nth elements of the lists.
    colNames <- gsub(pattern = l_rm[[n]], replacement = "", x = colNames)
  }
  # As we also need to add columns to hold SE values, let's replicate the column names.
  newCols <- rep(colNames, each = 2)
  # Take every 2nd element.
  indices <- seq(2, length(newCols), 2)
  # Paste "_SE" onto every 2nd element.
  newCols[indices] <- paste(newCols[indices], "_SE", sep = "")
  # Replace in column names with newCols, starting at the 2nd column (because the first column is missingness_level).
  names(dfErrorTrait)[2:ncol(dfErrorTrait)] <- newCols
  names(dfErrorTrait)
  # Remove SE cols for now as we don't need them anymore.
  dfErrorTrait <- dfErrorTrait[, -grep("_SE", colnames(dfErrorTrait))]
  # Pivot the dataframe to get it ready for plotting.
  dfErrorTrait <- pivot_longer(dfErrorTrait, !missingness_level, names_to = "method", values_to = "error")
  # Add trait column.
  dfErrorTrait$trait <- trait
  # Initiate ratio column.
  dfErrorTrait$ratio <- NA ## Initiating column.
  # Split by missingness level.
  l_dfMiss <- split(dfErrorTrait, dfErrorTrait$missingness_level)
  
  # For every missingness pattern..
  for(m in 1:length(l_dfMiss)) {
    # Take dfMiss.
    dfMiss <- l_dfMiss[[m]] 
    # Make sure dfMiss is dataframe.
    dfMiss <- as.data.frame(dfMiss)
    # Ratio the errors for each method.
    dfMiss <- RatioError(dfMiss, method = "KNN")
    dfMiss <- RatioError(dfMiss, method = "RF")
    dfMiss <- RatioError(dfMiss, method = "MICE")
    # Remove NAs.
    dfMiss <- na.omit(dfMiss)
    # Add gene column.
    dfMiss$gene <- gsub(".+_", "", dfMiss$method)
    # Modify method column.
    dfMiss$method <- word(dfMiss$method, start = 1, sep = "_")
    # Replace in l_dfMiss.
    l_dfMiss[[m]] <- dfMiss
  }
  
  # Bind l_dfMiss.
  dfMissAll <- bind_rows(l_dfMiss, .id = "missingness_level")
  # If the trait is activity_time...
  if(trait == "activity_time") {
    # Index to get the phylogenetic signal metrics for activity_time_Diurnal (trait with most variability).
    index <- grep("activity_time_Diurnal", names(l_dfMetricTrait))
    dfMetricTrait <- l_dfMetricTrait[[index]]
    # Replace with "activity_time" for matching purposes.
    dfMetricTrait$trait <- "activity_time" 
  } else { 
    # Index to get the corresponding dfMetricTrait for the trait.
    index <- grep(trait, names(l_dfMetricTrait))
    dfMetricTrait <- l_dfMetricTrait[[index]]
  }
  # Merge dfTraitPhySig with dfMiss to get the phylogenetic signal metrics.
  dfTraitPhySig <- merge(dfMissAll, dfMetricTrait, by = "gene")
  # Append to l_dfTraitPhySig.
  l_dfTraitPhySig[[t]] <- dfTraitPhySig
}

# Bind l_dfTraitPhySig.
dfPlotAll <- bind_rows(l_dfTraitPhySig)
# Split into continuous and categorical dataframes.
dfPlotCont <- dfPlotAll[dfPlotAll$trait.x %in% contTraits, ]
dfPlotCat <- setdiff(dfPlotAll, dfPlotCont)
# Replace latitude trait name in dfPlotCont.
dfPlotCont$trait.x[dfPlotCont$trait.x == "latitude_centroid_from_roll_et_al_2017"] <- "latitude"

# Continuous. ---
# Split dfPlot by missingness_level.
l_dfPlotCont <- split(dfPlotCont, dfPlotCont$missingness_level)
# For every missingness level...
for(m in 1:length(l_dfPlotCont)) {
  # Take mth dataframe.
  dfPlot <- l_dfPlotCont[[m]]
  # Remove extra col.
  dfPlot$trait.y <- NULL
  # Convert character variables to factor.
  dfPlot[, c("gene", "method", "trait.x")] <- lapply(dfPlot[, c("gene", "method", "trait.x")], as.factor)
  
  # Split by gene.
  l_dfGene <- split(dfPlot, f = dfPlot$gene)
  # For each gene..
  for(g in 1:length(l_dfGene)){
    # Gene name of gene.
    gene <- names(l_dfGene)[[g]]
    # Get dataframe.
    dfGene <- l_dfGene[[g]]
    
    # Dot plot. ---
    # Create title.
    plotTitle <- paste(m, "cont", gene, ".tiff", sep = "")
    # Write to .tiff file.
    tiff(plotTitle, units = "in", width = 14, height = 8, res = 600)
    theme_set(theme_light())
    plot <- ggplot(dfGene, aes(x = phylogenetic_signal, y = ratio)) +
      geom_point(aes(shape = method, color = method), size = 20, position = "jitter") +
      scale_color_manual(values=c("#56B4E9", "#CC79A7", "#009E73")) +
      scale_shape_manual(values=c(17, 15, 16)) +
      labs(title = plotTitle, y = "", x = "") +
      geom_hline(yintercept = 1, color = "gray", size = 2) +
      theme(axis.text.y = element_text(vjust = 0.6, size = 28, face = "bold"), axis.text.x = element_text(vjust = 0.6, size = 28, face = "bold")) + theme(legend.text = element_text(face = "bold")) + guides(shape = guide_legend(override.aes = list(size = 5))) 
    # Print to tiff.
    print(plot)
    # Turn off dev.
    dev.off()
  }
}
  

# Categorical. ---
# Split dfPlotCat by missingness_level.
l_dfPlotCat <- split(dfPlotCat, dfPlotCat$missingness_level)
# For every missingness level...
for(m in 1:length(l_dfPlotCat)) {
  # Take mth dataframe.
  dfPlot <- l_dfPlotCat[[m]]
  # Remove extra col.
  dfPlot$trait.y <- NULL
  # Convert character variables to factor.
  dfPlot[, c("gene", "method", "trait.x")] <- lapply(dfPlot[, c("gene", "method", "trait.x")], as.factor)
  # Split dfTraitType by gene.
  l_dfGene <- split(dfPlot, f = dfPlot$gene)
  # For each gene..
  for(g in 1:length(l_dfGene)){
    
    # Gene name of gene.
    gene <- names(l_dfGene)[[g]]
    # Get dataframe.
    dfGene <- l_dfGene[[g]]

    # Dot plot. ---
    # Create title.
    plotTitle <- paste(m, "cat", gene, ".tiff", sep = "")
    # Write to .tiff file.
    tiff(plotTitle, units = "in", width = 14, height = 8, res = 600)
    theme_set(theme_light())
    plot <- ggplot(dfGene, aes(x = phylogenetic_signal, y = ratio)) +
      geom_point(aes(shape = method, color = method), size = 20, position = "jitter") +
      scale_color_manual(values=c("#56B4E9", "#CC79A7", "#009E73")) +
      scale_shape_manual(values=c(17, 15, 16)) +
      labs(title = plotTitle, y = "", x = "") +
      geom_hline(yintercept = 1, color = "gray", size = 3) +
      theme(axis.text.y = element_text(vjust = 0.6, size = 28, face = "bold"), axis.text.x = element_text(vjust = 0.6, size = 28, face = "bold")) + theme(legend.text = element_text(face = "bold")) + guides(shape = guide_legend(override.aes = list(size = 5)))
    # Print to tiff.
    print(plot)
    # Turn off dev.
    dev.off()
  }
}