# Functions for simulating missingness and imputing values in complete-case datasets. Returns lists of error rates for both numerical and categorical variables.

MeanModeSimputeMCAR <- function(data, vars, int = 100, missLevel = 0.1) {
  
  # Given a complete-case dataset, this function simulates missingness completely at random (MCAR) and imputes values using mean/mode replacement. Determines error rates for numerical (MSE) and categorical variables (PFC).
  # data = complete-case dataset containing trait data. Must also contain a column containing species name data = "species_name".
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # missLevel  = proportion of values to remove from dataset (e.g. 0.1). Performs simulation/interation for 10% itervals (e.g. 10%, 20%, etc.)
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Identify integer (count) traits, if any.   
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")
  # Ensure categorical traits are character type for mode imputation, including binary variables because we will be finding most frequent category (mode) value later on.
  data[, catTraits] <- lapply(data[, catTraits], as.character)
  # Determine original sample size for each trait.
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Imputation prep. ---
  # Convert missLevel into vector of missingness proportions.  
  missingness <- seq(0.1, missLevel, by = 0.1)
  # Create list to hold the dataframes with simulated missingness from each iteration and level of missingness.
  l_l_dfMiss <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration and level of missingness.
  l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create lists to hold the error rates for each iteration. There will be n = missLevel number of error rates for each iteration (pertaining to each level of missingness).
  # ** = l_l_Error (two levels because there is no parameter tuning for mean/mode).
  l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  
  # For every iteration...
  for(i in 1:int) {
    # Create lists to hold the error rates, dfMiss, and dfImp at each missingness level.
    # ** = l_* (one level for every iteration/no parameter tuning).
    l_dfMiss <- CreateNamedList(listLength = length(missingness), elementNames = missingness)
    l_dfImp <- CreateNamedList(listLength = length(missingness), elementNames = missingness)
    l_Error <- CreateNamedList(listLength = length(missingness), elementNames = missingness)
    # Log-transformation of numerical traits.
    dfLog <- LogTransform(data, contTraits)
    
    # For every level of missingness...
    for (s in 1:length(missingness)) { 
      # MCAR simulation. --- 
      # Make a copy of complete-case dataframe (transformed) for introducing NAs.
      dfMiss <- dfLog
      # Introduce NAs depending on the level of missingness. We also set the seed so that the NAs are added on top of the NAs already introduced at 0.1.
      set.seed(i)
      # For every level of missingness..
      for(l in 1:s) {
        # Remove 0.1 of the values.
        dfMiss <- prodNA(dfMiss[, vars], noNA = 0.1)
      }
      
      # Dataframe organization. ---
      # Add species_name information back.
      dfMiss$species_name <- dfLog$species_name
      # Bind missingness indicator columns and reorganize the dataframe.
      dfMiss <- BindAndOrganize(dfMiss, vars)
      # Make sure dfMiss and dfLog (the original data) are both ordered alphabetically.
      dfMiss <- dfMiss[order(dfMiss$species_name), ]
      dfLog <- dfLog[order(dfLog$species_name), ]
      
      # Imputation. ---
      # The ImputeMeanMode() function imputes values using mean/mode replacement for continuous and categorical variables respectively, and returns a list of error rates (MSE for continuous traits and PFC for categorical traits).
      impResult <- ImputeMeanMode(dfTrue = dfLog, dfMissing = dfMiss, cols = vars, cont = contTraits, cat = catTraits, inter = intTraits)

      # Result handling. ---
      # Extract dfImp.
      dfImp <- impResult[[1]]
      # Extract errorRates.
      errorRates <- impResult[[2]]
      
      # Back-transforming data. ---
      # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
      dfOrig <- data[data$species_name %in% dfLog$species_name, ]
      # Apply the Backtransform function to the continuous traits.
      dfImp <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMiss, cols = contTraits)
      # If there are any integer traits..
      if(length(intTraits) > 0){
        # Round to nearest whole number.
        dfImp[, intTraits] <- lapply(dfImp[, intTraits], function(x) as.integer(round(x)))
      }
      
      # Append results to sth element of lists. 
      l_dfMiss[[s]] <- dfMiss
      l_dfImp[[s]] <- dfImp
      l_Error[[s]] <- errorRates
      
    } ## s
    
    # Append to ith elements of lists.
    l_l_dfMiss[[i]] <- l_dfMiss
    l_l_dfImp[[i]] <- l_dfImp
    l_l_Error[[i]] <- l_Error

  } ## i
  
  # Create list to hold the results.
  l_results <- list(completeCaseData = data, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, sampleSizes = l_sampleSizes, missingness = missingness, reps = int, missingData = l_l_dfMiss, imputedData = l_l_dfImp, errorRates = l_l_Error)
  # Return l_results.
  return(l_results)
  
}

MeanModeSimputeMAR <- function(data, raw, vars, int = 100) {
  
  # Given a complete-case dataset, this function simulates missingness at random (MAR) and imputes values using mean/mode replacement. Determines error rates for numerical (mean squared error - MSE) and categorical variables (proportion falsely classified - PFC).
  # data = complete-case dataset containing trait data. Must also contain a column with species name information = "species_name"
  # raw = original dataset containing trait data with missing values. Used to build logistic regression models and simulate data that are MAR in the complete-case dataset. Must also contain a column with species name information = "species_name"
  # vars = names of columns for which to simulate missing data
  # int = number of iterations (missingness replicates). Default is 100
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Identify integer (count) traits, if any.
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")
  
  # Ensure categorical traits are factor type for logistic regression model building (data class required for glm function).
  data[, catTraits] <- lapply(data[, catTraits], as.factor)
  raw[, catTraits] <- lapply(raw[, catTraits], as.factor)
  # Determine original sample size for each trait.
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Logistic regression fitting. ---
  # Apply FitLogReg() function to raw data to identify significant predictors of missingness for each trait.
  l_models <- FitLogReg(data = raw, cols = vars)
  # Now, let's identify traits that CAN be simulated MAR versus those that cannot (i.e. they must be simulated MCAR). Traits that must be simulated MCAR either 1) have no missing values in the original data or 2) do not have any predictors of missingness in the dataset upon fitting of the logistic regression models. IDMissPattern() also refines the logistic regression models (i.e. drops insignificant terms) through use of the RefineModels() function.
  l_missPatt <- IDMissPattern(data = raw, vars = vars, models = l_models)
  # Extract MCAR traits.
  traitsMCAR <- l_missPatt[[1]]
  # Extract MAR traits.
  traitsMAR <- l_missPatt[[2]]
  # Extract final MAR models.
  l_MARfinalModels <- l_missPatt[[3]]
  
  # Imputation prep. ---
  # Create list to hold the dataframes with simulated missingness for each iteration.
  l_dfMiss <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes for each iteration.
  l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create lists to hold the error rates for each iteration.
  l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  
  # For every iteration...
  for(i in 1:int) {
    
    # MAR simulation. ---
    # Create copy of complete-case dataframe (untransformed) to introduce NAs into.
    dfMissOrig <- data
    # Create a list to hold the missingness proportion for each trait to be simulated MAR.
    l_missingness <- CreateNamedList(listLength = length(traitsMAR), elementNames = traitsMAR)
    
    # For every trait that can be simulated MAR..
    for(t in 1:length(traitsMAR)){
      # Take tth trait.
      trait <- traitsMAR[[t]]
      # Get the corresponding MAR model.
      index <- grep(pattern = trait, x = names(l_MARfinalModels))
      modelMAR <- l_MARfinalModels[[index]]
      # Simulate MAR data in complete-case dataset.
      res <- SimMAR(model = modelMAR, data = data)
      # Replace complete-case column in dfMiss with res.
      dfMissOrig[, trait] <- res
      # Get original sample size of trait.
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMiss after MAR simulation to determine actual number of NAs introduced.
      marN <- sum(is.na(dfMissOrig[[trait]])) - sum(is.na(data[[trait]]))
      # Divide by n to determine missingness percentage for trait and append to l_missingness.
      l_missingness[[t]] <- marN/n
    }
    
    # Log transformation of numerical traits prior to imputation.
    dfLog <- LogTransform(data, contTraits)
    # Introduce missing values from dfMissOrig into dfLog.
    dfMiss <- as.data.frame(mapply(function(x, y) ifelse(is.na(x), x, y), x = dfMissOrig[, c("species_name", contTraits)], y = dfLog[, c("species_name", contTraits)], SIMPLIFY = F))
    # Add categorical traits back.
    dfMiss <- merge(dfMiss, dfMissOrig[, c("species_name", catTraits)], by = "species_name")
    # Bind missingness indicator columns and reorganize the dataframe.
    dfMiss <- BindAndOrganize(dfMiss, vars)
    # Make sure dfMiss and dfLog (the original data) are both ordered alphabetically.
    dfMiss <- dfMiss[order(dfMiss$species_name), ]
    dfLog <- dfLog[order(dfLog$species_name), ]
    # Ensure categorical traits are factor type.
    dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.factor)
    dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.factor)
    # Finally, append traits that could not be simulated MAR (no NAs introduced based on fitted logistic regression models) to traitsMCAR.
    traitsMCAR <- AppendMCAR(dfMiss, cols = traitsMAR, varsMCAR = traitsMCAR)
    # Update traitsMAR based on this result.
    traitsMAR <- setdiff(traitsMAR, traitsMCAR)
    
    # Imputation. ---
    # Ensure categorical traits are character type for mode imputation, including binary variables because we will be finding most frequent category (mode) value later on.
    dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.character)
    dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.character)
    # The ImputeMeanMode() function imputes values using mean/mode replacement for continuous and categorical variables respectively, and returns a list of error rates (MSE for continuous traits and PFC for categorical traits).
    impResult <- ImputeMeanMode(dfTrue = dfLog, dfMissing = dfMiss, cols = vars, cont = contTraits, cat = catTraits, inter = intTraits)
    
    # Result handling. ---
    # Extract dfImp.
    dfImp <- impResult[[1]]
    # Extract errorRates.
    errorRates <- impResult[[2]]
    
    # Back-transforming data. ---
    # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfLog$species_name, ]
    # Apply the Backtransform function to the continuous traits.
    dfImp <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMiss, cols = contTraits)
    # If there are any integer traits..
    if(length(intTraits) > 0){
      # Round to nearest whole number.
      dfImp[, intTraits] <- lapply(dfImp[, intTraits], function(x) as.integer(round(x)))
    }
    
    # Append to ith elements of lists.
    l_dfMiss[[i]] <- dfMiss
    l_l_missingness[[i]] <- l_missingness
    l_dfImp[[i]] <- dfImp
    l_Error[[i]] <- errorRates
    
  }
  
  # Final MCAR check. ---
  # If there are any traits that could not be simulated MAR..
  if(length(traitsMCAR) > 0){
    # If the trait was identified in a previous iteration to be an MCAR variable, it is possible missing values were introduced in other iterations. These will likely be very low numbers of NAs and not enough for error rate analyses. So here we will subset l_l_Error to ensure it only contains values for varsMAR.
    l_Error <- lapply(l_Error, function(x) {
      # Identify error rates associated with traitsMAR.
      index <- which(names(x) %in% traitsMAR)
      # Remove traits that were identified as MCAR in later iterations.
      x <- x[index]
    })
    print("The following traits could not be simulated MAR:")
    print(unique(traitsMCAR))
  }
  
  # Average out the missingness for each trait. ---
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(traitsMAR, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  
  # Subset to sample sizes for traits that could be simulated MAR.
  l_MARn <- l_sampleSizes[names(l_sampleSizes) %in% traitsMAR]
  # Calculate the average number of NAs introduced for each trait.
  avgNAs <- mapply(function(x, y) x * y, y = l_MARn, x = avgMiss)
  
  # Create list to hold the results.
  l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, sampleSizes = l_sampleSizes, traitsNotSimulated = traitsMCAR, traitsMAR = traitsMAR, MARfinalModels = l_MARfinalModels, reps = int, missingData = l_dfMiss, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_dfImp, errorRates = l_Error)
  # Return l_results.
  return(l_results)
  
}

MeanModeSimputeMNAR <- function(data, vars, int = 100, quantiles = NULL, directions = NULL, absolutes = NULL, categories = NULL, ...) {
  
  # Given a complete-case dataset, this function simulates missingness not at random (MNAR) and imputes values using mean/mode replacement. Determines error rates for numerical (mean squared error - MSE) and categorical variables (proportion falsely classified - PFC).
  # data = complete-case dataset containing trait data. Must also contain a column with species name information = "species_name"
  # vars = names of columns for which to simulate missing data
  # int = number of iterations (missingness replicates). Default is 100
  # directions = for each numerical trait, named character vector of thresholds for which to introduce NAs
  # absolutes = for each numerical trait, named logical vector indicating whether to take the absolute values of the data when determining MNAR thresholds
  # categories = for each categorical trait, named character vector of categories in which to introduce NAs
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data = data, traitCols = vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Identify integer (count) traits, if any. **
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")
  # Determine original sample size for each trait.
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Imputation prep. ---
  # Create list to hold the dataframes with simulated missingness for each iteration.
  l_dfMissOrig <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes for each iteration.
  l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create lists to hold the error rates for each iteration.
  l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  
  # For every iteration...
  for(i in 1:int) {

    # MNAR simulation. ---
    # Create copy of complete-case dataframe (untransformed) to introduce NAs into.
    dfMissOrig <- data
    # Create a list to hold the missingness proportion for each trait to be simulated MNAR.
    l_missingness <- CreateNamedList(listLength = length(vars), elementNames = vars)
    
    # For every trait to be simulated MNAR..
    for(t in 1:length(vars)){
      # Take tth trait.
      trait <- vars[[t]]
      # If trait is numeric..
      if(trait %in% contTraits){
        # Take quantile that corresponds to trait.
        quant <- quantiles[[grep(trait, names(quantiles))]]
        # Take direction that corresponds to trait.
        dir <- directions[[grep(trait, names(directions))]]
        # Take absolute indicator that corresponds to trait.
        absInd <- absolutes[[grep(trait, names(absolutes))]]
        # Simulate MNAR data in trait.
        res <- SimContMNAR(var = dfMissOrig[[trait]], quantile = quant, direction = dir, absolute = absInd)
        # Replace complete-case column in dfMissOrig with res.
        dfMissOrig[, trait] <- res
      } else if(trait %in% catTraits){
        # Take category that corresponds to trait.
        cate <- categories[grep(trait, names(categories))]
        # Simulate MNAR data in trait. 
        res <- SimCatMNAR(var = dfMissOrig[[trait]], category = cate)
        # Replace complete-case column in dfMissOrig with res.
        dfMissOrig[, trait] <- res
      } 
      # Get original sample size of trait.
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMissOrig after MNAR simulation to determine actual number of NAs introduced.
      marN <- sum(is.na(dfMissOrig[[trait]])) - sum(is.na(data[[trait]]))
      # Divide by n to determine missingness percentage for trait and append to l_missingness.
      l_missingness[[t]] <- marN/n
    }
    
    # Dataframe organization. ---
    # Log transformation of numerical traits prior to imputation.
    dfLog <- LogTransform(data, contTraits)
    # Introduce missing values from dfMissOrig into dfLog.
    dfMiss <- as.data.frame(mapply(function(x, y) ifelse(is.na(x), x, y), x = dfMissOrig[, c("species_name", contTraits)], y = dfLog[, c("species_name", contTraits)], SIMPLIFY = F))
    # Add categorical traits back.
    dfMiss <- merge(dfMiss, dfMissOrig[, c("species_name", catTraits)], by = "species_name")
    # Bind missingness indicator columns and reorganize the dataframe.
    dfMiss <- BindAndOrganize(dfMiss, colnames(data)[-1])
    # Make sure dfMiss and dfLog (the original data) are both ordered alphabetically.
    dfMiss <- dfMiss[order(dfMiss$species_name), ]
    dfLog <- dfLog[order(dfLog$species_name), ]
    
    # Imputation. ---
    # Ensure categorical traits are character type for mode imputation, including binary variables because we will be finding most frequent category (mode) value later on.
    # If there is more than one categorical trait..
    if(length(catTraits) > 1){
      # Use lapply to convert to character.
      dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.character)
      dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.character)
    } else if(length(catTraits) == 1){
      # Convert trait to character.
      dfLog[[catTraits]] <- as.character(dfLog[[catTraits]])
      dfMiss[[catTraits]] <- as.character(dfMiss[[catTraits]])
    }
    # The ImputeMeanMode() function imputes values using mean/mode replacement for continuous and categorical variables respectively, and returns a list of error rates (MSE for continuous traits and PFC for categorical traits).
    impResult <- ImputeMeanMode(dfTrue = dfLog, dfMissing = dfMiss, cols = vars, cont = contTraits, cat = catTraits, inter = intTraits)
    
    # Result handling. ---
    # Extract dfImp.
    dfImp <- impResult[[1]]
    # Extract errorRates.
    errorRates <- impResult[[2]]
    
    # Back-transforming data. ---
    # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfLog$species_name, ]
    # Apply the Backtransform function to the continuous traits.
    dfImp <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMiss, cols = contTraits)
    # If there are any integer traits..
    if(length(intTraits) > 0){
      # Round to nearest whole number.
      dfImp[, intTraits] <- lapply(dfImp[, intTraits], function(x) as.integer(round(x)))
    }
    
    # Append to ith elements of lists.
    l_dfMissOrig[[i]] <- dfMissOrig
    l_l_missingness[[i]] <- l_missingness
    l_dfImp[[i]] <- dfImp
    l_Error[[i]] <- errorRates
    
  }
  
  # Average out the missingness for each trait. ---
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(vars, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  # Calculate the average number of NAs introduced for each trait.
  avgNAs <- mapply(function(x, y) x * y, y = l_sampleSizes, x = avgMiss)
  
  # Create list to hold the results.
  l_results <- list(completeCaseData = data, traits = vars, numeric = contTraits, categorical = catTraits, integer = intTraits, sampleSizes = l_sampleSizes, quantiles = quantiles, directions = directions, absolutes = absolutes, categories = categories, reps = int, missingData = l_dfMissOrig, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_dfImp, errorRates = l_Error)
  # Return l_results.
  return(l_results)
  
}

KNNSimputeMCAR <- function(data, vars, int = 100, missLevel = 0.1, k = 50, phyImp = F, tree = NULL) {
  
  # Given a complete-case dataset, this function simulates missingness completely at random (MCAR) and imputes values using the kNN() function in the "VIM" package for imputation. Determines error rates for numerical (MSE) and categorical variables (PFC).
  # Citations: Kowarik A, Templ M (2016). “Imputation with the R Package VIM.” Journal of Statistical Software, 74(7), 1–16. doi: 10.18637/jss.v074.i07.
  # https://cran.r-project.org/web/packages/VIM/VIM.pdf
  # data = complete-case dataset containing trait data. Must also contain a column containing species name data = "species_name".
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # missLevel* = proportion of values to remove from dataset (e.g. 0.1). Performs simulation/interation for 10% intervals (e.g. 10%, 20%, etc.)
  # k = maximum number of nearest neighbours to test
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Convert categorical traits to factors.
  data[, catTraits] <- lapply(data[, catTraits], as.factor)
  # Determine original sample size for each trait.
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data = data[, vars])
  
  # Imputation prep. ---
  # Convert missLevel into vector of missingness proportions.
  missingness <- seq(0.1, missLevel, by = 0.1)
  # Create list to hold the dataframes with simulated missingness from each iteration and level of missingness.
  l_l_dfMiss <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration and level of missingness.
  l_l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the error rates for each iteration. There will be n = missLevel number of error rates for each iteration (pertaining to each level of missingness).
  l_l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.
    l_l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # Log transformation of numerical traits prior to imputation.
    dfLog <- LogTransform(data, contTraits)
    # Create lists to hold the error rates, dfMiss, and dfImp at each missingness level.  
    l_dfMiss <- CreateNamedList(listLength = length(missingness), elementNames = missingness)
    l_l_dfImp <- CreateNamedList(listLength = length(missingness), elementNames = missingness)
    l_l_Error <- CreateNamedList(listLength = length(missingness), elementNames = missingness)
    # If phylogenetic imputation was chosen..  
    if(phyImp == T) {
      l_l_evs <- CreateNamedList(listLength = length(missingness), elementNames = missingness)
    }

    # For every level of missingness...
    for (s in 1:length(missingness)) {
  
      # MCAR simulation. ---  
      # Make a copy of complete-case dataframe (transformed) for introducing NAs.
      dfMiss <- dfLog
      # Introduce NAs depending on the level of missingness. We also set the seed so that the NAs are added on top of the NAs already introduced at 0.1.
      set.seed(i)
      # For every level of missingness..
      for(l in 1:s) {
        # Remove 0.1 of the values.
        dfMiss <- prodNA(dfMiss[, vars], noNA = 0.1)
      }
      
      # Dataframe organization. ---
      # Add species_name information back.
      dfMiss$species_name <- dfLog$species_name
      # Bind missingness indicator columns and reorganize the dataframe.
      dfMiss <- BindAndOrganize(dfMiss, vars)
      # Make sure dfMiss and dfLog (the original data) are both ordered alphabetically.
      dfMiss <- dfMiss[order(dfMiss$species_name), ]
      dfLog <- dfLog[order(dfLog$species_name), ]
      # Ensure categorical traits are factor type.
      dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.factor)
      dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.factor)
      
      # Phylogenetic eigenvector decomposition. ---
      if(phyImp == T) {
        # Append eigenvectors to dfMiss and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors. 
        l_evs <- AppendEigenvectors(data = dfMiss, vars = vars, tree = tree, predictors = l_predictors)
        # Extract list of dfMiss.
        l_dfMiss <- l_evs[[1]]
        # Extract updated list of predictors.
        l_EVPredictors <- l_evs[[2]]
      }
      
      # Identify outgroup(s) in trait datasets (this is the species that contains no trait data dfLog).
      outgroup <- dfLog$species_name[apply(dfLog[, traits], 1, function(x) all(is.na(x)))]
      # If found in dataframe..
      if(length(outgroup) > 0){
        # Remove outgroup from dataframes.
        dfLog <- dfLog[!dfLog$species_name %in% outgroup, ]
        dfMiss <- dfMiss[!dfMiss$species_name %in% outgroup, ]
        l_dfMiss <- lapply(l_dfMiss, function(x) x[!x$species_name %in% outgroup, ])
      }
 
      # Imputation. ---
      # The ImputeKNN() function entails a loop that uses different values of k (the number of nearest neighbours to use in the kNN algorithm) and returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.
      # If phylogenetic imputation was selected..
      if(phyImp == T) {
        # Impute data using eigenvectors.
        impResult <- ImputeKNN(dfTrue = dfLog, dfMissing = dfMiss, cols = vars, cont = contTraits, cat = catTraits, k = 50, predictors = l_EVPredictors, phyImp = T, l_dfMissing = l_dfMiss)
      } else if(phyImp == F){
        # Impute data only using trait data.
        impResult <- ImputeKNN(dfTrue = dfLog, dfMissing = dfMiss, cols = vars, cont = contTraits, cat = catTraits, k = 50, predictors = l_predictors)
      }

      # Result handling. ---
      # Extract l_dfImp.
      l_dfImp <- impResult[[1]]
      # Extract errorRates.
      l_Error <- impResult[[2]]
      
      # Back-transforming data. ---
      # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
      dfOrig <- data[data$species_name %in% dfLog$species_name, ]
      # For every imputed dataframe..
      for(d in 1:length(l_dfImp)){
        # Take the dth dfImp.
        dfImp <- l_dfImp[[d]]
        # For every continuous trait..
        for(t in 1:length(contTraits)){
          # Take the tth continuous trait.
          trait <- contTraits[[t]]
          # Apply the Backtransform function for the trait in question.
          dfImpBTrait <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMiss, cols = trait)
          # Replace data in dfImp with backtransformed data.
          dfImp[, trait] <- dfImpBTrait[, trait]
        }
        # Replace dfImp in l_dfImp with newly backtransformed dataset.
        l_dfImp[[d]] <- dfImp
      }
      
      # Append results to sth element of lists.
      l_dfMiss[[s]] <- dfMiss
      l_l_dfImp[[s]] <- l_dfImp
      l_l_Error[[s]] <- l_Error
      # If phylogenetic imputation was selected..
      if(phyImp == T) {
        l_l_evs[[s]] <- l_evs
      }

    } ## s
    
    # Append to ith elements of lists.
    l_l_dfMiss[[i]] <- l_dfMiss
    l_l_l_dfImp[[i]] <- l_l_dfImp
    l_l_l_Error[[i]] <- l_l_Error
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      l_l_l_evs[[i]] <- l_l_evs
    }
    
  } ## i
  
  # If phylogenetic imputation was selected..
  if(phyImp == T) { 
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, traits = l_traits, numeric = contTraits, categorical = catTraits, k = k, missingness = missingness, reps = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_l_evs, missingData = l_l_dfMiss, imputedData = l_l_l_dfImp, errorRates = l_l_l_Error)
  } else if(phyImp == F){
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, traits = l_traits, numeric = contTraits, categorical = catTraits, sampleSizes = l_sampleSizes, k = k, missingness = missingness, reps = int, predictors = l_predictors, missingData = l_l_dfMiss, imputedData = l_l_l_dfImp, errorRates = l_l_l_Error)
  }
  # Return l_results.
  return(l_results)
  
}

KNNSimputeMAR <- function(data, raw, vars, int = 100, k = 50, phyImp = F, tree = NULL) {
  
  # Given a complete-case dataset, this function simulates missingness at random (MAR) and imputes values using the kNN() function in the "VIM" package for imputation. Determines error rates for numerical (MSE) and categorical variables (PFC).
  # Citations: Kowarik A, Templ M (2016). “Imputation with the R Package VIM.” Journal of Statistical Software, 74(7), 1–16. doi: 10.18637/jss.v074.i07.
  # https://cran.r-project.org/web/packages/VIM/VIM.pdf
  
  # data = complete-case dataset containing trait data. Must also contain a column containing species name data = "species_name".
  # raw* = original dataset containing trait data with missing values. Used to build logistic regression models and simulate data that are MAR in the complete-case dataset. Must also contain a column with species name information = "species_name"
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # k = maximum number of nearest neighbours to test
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Ensure categorical traits are factor type for logistic regression model building (data class required for glm function).
  data[, catTraits] <- lapply(data[, catTraits], as.factor)
  raw[, catTraits] <- lapply(raw[, catTraits], as.factor)
  # Determine original sample size for each trait. &&
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data[, vars])
  
  # Logistic regression fitting. ---
  # Apply FitLogReg() function to raw data to identify significant predictors of missingness for each trait.
  l_models <- FitLogReg(data = raw, cols = vars)
  # Now, let's identify traits that CAN be simulated MAR versus those that cannot (i.e. they must be simulated MCAR). Traits that must be simulated MCAR either 1) have no missing values in the original data or 2) do not have any predictors of missingness in the dataset upon fitting of the logistic regression models. IDMissPattern() also refines the logistic regression models (i.e. drops insignificant terms) through use of the RefineModels() function.
  l_missPatt <- IDMissPattern(data = raw, vars = vars, models = l_models)
  # Extract MCAR traits.
  traitsMCAR <- l_missPatt[[1]]
  # Extract MAR traits.
  traitsMAR <- l_missPatt[[2]]
  # Extract final MAR models.
  l_MARfinalModels <- l_missPatt[[3]]
  
  # Imputation prep. ---
  # Create lists to hold the error rates for each iteration.
  l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the dataframes with simulated missingness from each iteration.
  l_dfMissOrig <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration.
  l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.
    l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # MAR simulation. ---
    # Create copy of complete-case dataframe (untransformed) to introduce NAs into.
    dfMissOrig <- data
    # Create a list to hold the missingness proportion for each trait to be simulated MAR.
    l_missingness <- CreateNamedList(listLength = length(traitsMAR), elementNames = traitsMAR)
    
    # For every trait that can be simulated MAR..
    for(t in 1:length(traitsMAR)){
      # Take tth trait.
      trait <- traitsMAR[[t]]
      # Get the corresponding MAR model.
      index <- grep(pattern = trait, x = names(l_MARfinalModels))
      modelMAR <- l_MARfinalModels[[index]]
      # Simulate MAR data in complete-case dataset.
      res <- SimMAR(model = modelMAR, data = data)
      # Replace complete-case column in dfMissOrig with res.
      dfMissOrig[, trait] <- res
      # Get original sample size of trait.
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMissOrig after MAR simulation to determine actual number of NAs introduced. &&
      marN <- sum(is.na(dfMissOrig[[trait]])) - sum(is.na(data[[trait]]))
      # Divide by n to determine missingness percentage for trait and append to l_missingness.
      l_missingness[[t]] <- marN/n
    }
    
    # Dataframe organization. ---
    # Log transformation of numerical traits prior to imputation.
    dfLog <- LogTransform(data = data, cols = contTraits)
    # Introduce missing values from dfMissOrig into dfLog.
    dfMiss <- as.data.frame(mapply(function(x, y) ifelse(is.na(x), x, y), x = dfMissOrig[, c("species_name", contTraits)], y = dfLog[, c("species_name", contTraits)], SIMPLIFY = F))
    # Add categorical traits back.
    dfMiss <- merge(dfMiss, dfMissOrig[, c("species_name", catTraits)], by = "species_name")
    # Bind missingness indicator columns and reorganize the dataframe.
    dfMiss <- BindAndOrganize(dfMiss, vars)
    # Make sure dfMiss and dfLog (the original data) are both ordered alphabetically.
    dfMiss <- dfMiss[order(dfMiss$species_name), ]
    dfLog <- dfLog[order(dfLog$species_name), ]
    # Ensure categorical traits are factor type.
    dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.factor)
    dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.factor)
    # Finally, append traits that could not be simulated MAR (no NAs introduced based on fitted logistic regression models) to traitsMCAR.
    traitsMCAR <- AppendMCAR(dfMiss, cols = traitsMAR, varsMCAR = traitsMCAR)
    # Update traitsMAR based on this result.  
    traitsMAR <- setdiff(traitsMAR, traitsMCAR)
    
    # Phylogenetic eigenvector decomposition. ---
    if(phyImp == T) {
      # Append eigenvectors to dfMiss and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors.
      l_evs <- AppendEigenvectors(data = dfMiss, vars = traitsMAR, tree = tree, predictors = l_predictors)
      # Extract list of dfMiss.
      l_dfMissEV <- l_evs[[1]]
      # Extract updated list of predictors.
      l_EVPredictors <- l_evs[[2]]
    }
    
    # Identify outgroup(s) in trait datasets (this is the species that contains no trait data dfLog). &&
    outgroup <- dfLog$species_name[apply(dfLog[, traits], 1, function(x) all(is.na(x)))]
    # If found in dataframe.. &&
    if(length(outgroup) > 0){
      # Remove outgroup from dataframes.
      dfLog <- dfLog[!dfLog$species_name %in% outgroup, ]
      dfMiss <- dfMiss[!dfMiss$species_name %in% outgroup, ]
      if(phyImp == T){
        # Also remove outgroup from l_dfMissEV.
        l_dfMissEV <- lapply(l_dfMissEV, function(x) x[!x$species_name %in% outgroup, ])
      }
    }
    
    # Imputation. ---   
    # The ImputeKNN() function entails a loop that uses different values of k (the number of nearest neighbours to use in the kNN algorithm) and returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      # Impute data using eigenvectors.
      impResult <- ImputeKNN(dfTrue = dfLog, dfMissing = dfMiss, cols = traitsMAR, cont = contTraits, cat = catTraits, k = 50, predictors = l_EVPredictors, phyImp = T, l_dfMissing = l_dfMissEV)
    } else if(phyImp == F){
      # Impute data only using trait data.
      impResult <- ImputeKNN(dfTrue = dfLog, dfMissing = dfMiss, cols = traitsMAR, cont = contTraits, cat = catTraits, k = 50, predictors = l_predictors)
    }
    
    # Result handling. ---  
    # Extract l_dfImp.  
    l_dfImp <- impResult[[1]]
    # Extract errorRates  .
    l_Error <- impResult[[2]]
    
    # Back-transforming data. ---
    # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfLog$species_name, ]
    
    # For every imputed dataframe..
    for(d in 1:length(l_dfImp)){
      # Take the dth dfImp.
      dfImp <- l_dfImp[[d]]
      # Apply BackTransform function to contTraits in dfImp.
      dfImp <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMiss, cols = traitsMAR[traitsMAR %in% contTraits])
      # Replace dfImp in l_dfImp with newly backtransformed dataset.
      l_dfImp[[d]] <- dfImp
    }
    
    # Append results to ith element of lists.
    l_dfMissOrig[[i]] <- dfMissOrig
    l_l_missingness[[i]] <- l_missingness
    l_l_dfImp[[i]] <- l_dfImp
    l_l_Error[[i]] <- l_Error
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      l_l_evs[[i]] <- l_evs
    }
    
  } ## i
  
  # Final MCAR check. ---
  # If there are any traits that could not be simulated MAR..
  if(length(traitsMCAR) > 0){
    # If the trait was identified in a previous iteration to be an MCAR variable, it is possible missing values were introduced in other iterations. These will likely be very low numbers of NAs and not enough for error rate analyses. So here we will subset l_l_Error to ensure it only contains values for varsMAR.
    l_l_Error <- lapply(l_l_Error, function(x) {
      # Identify error rates associated with traitsMAR.
      index <- which(names(x) %in% traitsMAR)
      # Remove traits that were identified as MCAR in later iterations.
      x <- x[index]
    })
    print("The following traits could not be simulated MAR:")
    print(unique(traitsMCAR))
  }
  
  # Average out the missingness for each trait. ---
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(traitsMAR, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  # Subset to sample sizes for traits that could be simulated MAR.
  l_MARn <- l_sampleSizes[names(l_sampleSizes) %in% traitsMAR]
  # Calculate the average number of NAs introduced for each trait.
  avgNAs <- mapply(function(x, y) x * y, y = l_MARn, x = avgMiss)
  
  # If phylogenetic imputation was selected..
  if(phyImp == T) { 
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, sampleSizes = l_sampleSizes, traitsNotSimulated = traitsMCAR, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, k = k, missingness = missingness, reps = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_evs, missingData = l_dfMissOrig, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  } else if(phyImp == F){
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, sampleSizes = l_sampleSizes, traitsNotSimulated = traitsMCAR, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, k = k, missingness = missingness, reps = int, predictors = l_predictors, missingData = l_dfMissOrig, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  }
  # Return results.  
  return(l_results)
  
}

KNNSimputeMNAR <- function(data, vars, int = 100, k = 50, quantiles = NULL, directions = NULL, absolutes = NULL, categories = NULL, phyImp = F, tree = NULL, ...) {
  
  # Given a complete-case dataset, this function simulates missingness not at random (MNAR) and imputes values using the kNN() function in the "VIM" package for imputation. Determines error rates for numerical (MSE) and categorical variables (PFC).
  # Citations: Kowarik A, Templ M (2016). “Imputation with the R Package VIM.” Journal of Statistical Software, 74(7), 1–16. doi: 10.18637/jss.v074.i07.
  # https://cran.r-project.org/web/packages/VIM/VIM.pdf
  
  # data = complete-case dataset containing trait data. Must also contain a column containing species name data = "species_name".
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # k = maximum number of nearest neighbours to test
  # quantiles = named numeric vector of quantiles at which to introduce NAs for each numerical trait  
  # directions = for each numerical trait, named character vector of directions ("greater" or "less") for which to introduce NAs with regards to the corresponding quantile  
  # absolutes = for each numerical trait, named logical vector indicating whether to take the absolute values of the data when determining MNAR thresholds  
  # categories = for each categorical trait, named character vector of categories in which to introduce NAs  
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---  
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Ensure categorical traits are factor type for model building (data class required for glm function).   
  data[, catTraits] <- lapply(data[, catTraits], as.factor)
  # Determine original sample size for each trait.
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---   
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data[, vars])
  
  # Imputation prep. ---  
  # Create lists to hold the error rates for each iteration.
  l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the dataframes with simulated missingness for each trait and each iteration.     
  l_l_dfMissTraitOrig <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration.
  l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..   
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.
    l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # MNAR simulation. ---  
    # Create a list of complete-case dataframe to introduce NAs into (one for each trait).
    l_dfMissTraitOrig <- lapply(1:length(vars), function(x) data)
    # Name l_dfMissTraitOrig according to traits.
    names(l_dfMissTraitOrig) <- vars
    # Create an empty list to hold the missingness proportion for each trait.
    l_missingness <- CreateNamedList(listLength = length(vars), elementNames = vars)
    
    # For every trait to be simulated MNAR..
    for(t in 1:length(vars)){
      # Take tth trait.
      trait <- vars[[t]]
      # Take the corresponding dfMissTrait for the trait.
      dfMissTraitOrig <- l_dfMissTraitOrig[[grep(trait, names(l_dfMissTraitOrig))]]
      # If trait is numeric..
      if(trait %in% contTraits){
        # Take quantile that corresponds to trait.
        quant <- quantiles[[grep(trait, names(quantiles))]]
        # Take direction that corresponds to trait.
        dir <- directions[[grep(trait, names(directions))]]
        # Take absolute indicator that corresponds to trait.
        absInd <- absolutes[[grep(trait, names(absolutes))]]
        # Simulate MNAR data in trait.
        res <- SimContMNAR(var = dfMissTraitOrig[[trait]], quantile = quant, direction = dir, absolute = absInd)
        # Replace complete-case column in dfMiss with res.
        dfMissTraitOrig[, trait] <- res
      } else if(trait %in% catTraits){
        # Take category that corresponds to trait.
        cate <- categories[[grep(trait, names(categories))]]
        # Simulate MNAR data in trait. 
        res <- SimCatMNAR(var = dfMissTraitOrig[[trait]], category = cate)
        # Replace complete-case column in dfMiss with res.
        dfMissTraitOrig[, trait] <- res
      } 
      # Get original sample size of trait.
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMiss after MNAR simulation to determine actual number of NAs introduced.
      marN <- sum(is.na(dfMissTraitOrig[[trait]])) - sum(is.na(data[[trait]]))
      # Divide by n to determine missingness percentage for trait and append to l_missingness.
      l_missingness[[t]] <- marN/n
      # Append dfMissTraitOrig with introduced NAs for the trait to l_dfMissTraitOrig.
      l_dfMissTraitOrig[[t]] <- dfMissTraitOrig
    }
    
    # Dataframe organization. ---
    # Log transformation of numerical traits prior to imputation.
    dfLog <- LogTransform(data = data, cols = contTraits)
    # For each numerical trait, introduce missing values from each dfMissOrig into copies of dfLog.
    l_dfMissTrait <- lapply(l_dfMissTraitOrig, function(x) {
      x <- mapply(function(x, y) ifelse(is.na(x), x, y), x = x[, c("species_name", contTraits)], y = dfLog[, c("species_name", contTraits)], SIMPLIFY = F)
    })
    # Convert all elements back to dataframe format. !!!
    l_dfMissTrait <- lapply(l_dfMissTrait, as.data.frame)
    # Merge categorical traits in each dataframe of l_dfMissTraitOrig back to each dataframe in l_dfMissTrait. !!!
    for(c in 1:length(l_dfMissTrait)){
      l_dfMissTrait[[c]] <- merge(l_dfMissTrait[[c]], l_dfMissTraitOrig[[c]][, c("species_name", catTraits)], by = "species_name")
    }
    # Bind missingness indicator columns and reorganize the dataframe.  
    l_dfMissTrait <- lapply(l_dfMissTrait, BindAndOrganize, vars = vars)
    # Make sure l_dfMissTrait and dfLog (the original data) are all ordered alphabetically by species_name.  
    l_dfMissTrait <- lapply(l_dfMissTrait, function(x) x[order(x$species_name), ])
    dfLog <- dfLog[order(dfLog$species_name), ]
    # Ensure categorical traits are factor type.
    dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.factor)
    l_dfMissTrait <- lapply(l_dfMissTrait, function(x){
      x[catTraits] <- lapply(x[catTraits], as.factor)
      return(x)})
    
    # Phylogenetic eigenvector decomposition. ---     
    if(phyImp == T) {
      # Append eigenvectors to dataframes in l_dfMissTrait and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors.
      l_evsTrait <- mapply(AppendEigenvectors, data = l_dfMissTrait, vars = vars, MoreArgs = list(tree = tree, predictors = l_predictors), SIMPLIFY = F)
      # Extract list of dfMissEV.
      l_dfMissEV <- sapply(l_evsTrait, function(x) x[[1]])
      names(l_dfMissEV) <- traits
      # Extract updated list of predictors.
      l_EVPredictors <- sapply(l_evsTrait, function(x) x[[2]])
    }
    # Identify outgroup(s) in trait datasets (this is the species that contains no trait data dfLog).   
    outgroup <- dfLog$species_name[apply(dfLog[, vars], 1, function(x) all(is.na(x)))]
    # If found in dataframe..
    if(length(outgroup) > 0){
      # Remove outgroup from dataframes.  
      dfLog <- dfLog[!dfLog$species_name %in% outgroup, ]
      l_dfMissTrait <- lapply(l_dfMissTrait, function(x) x[!x$species_name %in% outgroup, ])
      # If phylogenetic imputation was selected..
      if(phyImp == T){
        # Remove outgroup from dataframes in l_dfMissEV.  
        l_dfMissEV <- lapply(l_dfMissEV, function(x) x[!x$species_name %in% outgroup, ])
      }
    }
    
    # Imputation. ---
    
    # The ImputeKNN() function entails a loop that uses different values of k (the number of nearest neighbours to use in the kNN algorithm) and returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.
    
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      # Impute data using eigenvectors.
      impResult <- ImputeKNN(dfTrue = dfLog, dfMissing = l_dfMissEV[[1]], cols = vars, cont = contTraits, cat = catTraits, k = 50, predictors = l_EVPredictors, phyImp = T, l_dfMissing = l_dfMissEV)
      # Result handling. ---
      # Extract l_dfImp.
      l_dfImp <- impResult[[1]]
      # Extract errorRates.
      l_Error <- impResult[[2]]
    } else if(phyImp == F){
      # Impute data only using trait data.
      l_impResult <- mapply(ImputeKNN, dfMissing = l_dfMissTrait, cols = vars, MoreArgs = list(dfTrue = dfLog, cont = contTraits, cat = catTraits, k = 50, predictors = l_predictors), SIMPLIFY = F)
      # Since now we have a list of impResult, create new lists to hold the imputed dataframes and error rates.
      l_ImpTrait <- lapply(l_impResult, function(x) x[[1]])
      l_Error <- sapply(l_impResult, function(x) x[[2]])
      names(l_Error) <- vars
      # Creation of l_dfImp so that there is one dataframe per parameter.
      # Take species name and missingness information from dfMissing.
      dfImp <- as.data.frame(dfLog[, c("species_name")])
      names(dfImp) <- "species_name"
      # Make a list to hold the dataframes for each value of k.
      l_dfImp <- lapply(1:k, function(x) dfImp)
      names(l_dfImp) <- 1:k
      # For every trait...
      for (t in 1:length(vars)) { 
        # Take the tth trait.
        trait <- vars[[t]]
        # Extract the tth element in l_ImpTrait
        l_dfImputedTrait <- l_ImpTrait[[t]]
        # Get the name of the corresponding missingness info column.
        missCol <- grep(pattern = paste(trait, "_NA", sep = ""), x = colnames(l_dfImputedTrait[[1]]))
        missCol <- names(l_dfImputedTrait[[1]])[missCol]
        # For each value of k..
        for(k in 1:length(l_dfImputedTrait)){
          # Extract the kth dataframe in l_dfImputedTrait.
          dfImputedTrait <- l_dfImputedTrait[[k]]
          # Subset to only contain species_name, trait info, and missingness indicator column.
          dfImputedTrait <- dfImputedTrait[, c("species_name", trait, missCol)]
          # Extract the kth dataframe in l_dfImputedParam.
          dfImp <- l_dfImp[[k]]
          # Merge with dfImputedTrait.
          dfImp <- merge(dfImp, dfImputedTrait, by = "species_name")
          # Replace kth element of l_dfImputedParam.
          l_dfImp[[k]] <- dfImp
        }
      }
    }
    
    # Back-transforming data. ---
    # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfLog$species_name, ]
    # For every imputed dataframe..
    for(d in 1:length(l_dfImp)){
      # Take the dth dfImp.
      dfImp <- l_dfImp[[d]]
      # For every continuous trait..
      for(t in 1:length(contTraits)){
        # Take the tth continuous trait.
        trait <- contTraits[[t]]
        # Take the corresponding dfMissTrait.
        dfMissTrait <- l_dfMissTrait[[grep(trait, names(l_dfMissTrait))]]
        # Apply the Backtransform function for the trait in question.
        dfImpBTrait <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMissTrait, cols = trait)
        # Replace data in dfImp with backtransformed data.
        dfImp[, trait] <- dfImpBTrait[, trait]
      }
      # Replace dfImp in l_dfImp with newly backtransformed dataset.
      l_dfImp[[d]] <- dfImp
    }
    # Append results to ith element of lists.
    l_l_dfMissTraitOrig[[i]] <- l_dfMissTraitOrig
    l_l_missingness[[i]] <- l_missingness 
    l_l_dfImp[[i]] <- l_dfImp
    l_l_Error[[i]] <- l_Error
    # If phylogenetic imputation was selected..   
    if(phyImp == T) {
      l_l_evs[[i]] <- l_evsTrait
    }
  } ## i
  
  # Average out the missingness for each trait. ---  
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(vars, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  # Calculate the average number of NAs introduced for each trait.
  avgNAs <- mapply(function(x, y) x * y, y = l_sampleSizes, x = avgMiss)
  
  # If phylogenetic imputation was selected..     
  if(phyImp == T) { 
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, traits = vars, numeric = contTraits, categorical = catTraits, sampleSizes = l_sampleSizes, quantiles = quantiles, directions = directions, absolutes = absolutes, categories = categories, k = k, missingness = missingness, reps = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_evs, missingData = l_l_dfMissTraitOrig, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  } else if(phyImp == F){
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, traits = vars, numeric = contTraits, categorical = catTraits, sampleSizes = l_sampleSizes, quantiles = quantiles, directions = directions, absolutes = absolutes, categories = categories, k = k, missingness = missingness, reps = int, predictors = l_predictors, missingData = l_l_dfMissTraitOrig, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  }
  # Return results.
  return(l_results)
  
}

RFSimputeMCAR <- function(data, vars, int = 100, missLevel = 0.1, ntrees = c(100, 1000), phyImp = F, tree = NULL) {
  
  # Given a complete-case dataset, this function simulates missingness completely at random (MCAR) and imputes values using the missForest() function in the "missForest" package. Determines error rates for numerical (MSE) and categorical variables (PFC).
  # Citations: Daniel J. Stekhoven (2013). missForest: Nonparametric Missing Value Imputation using Random Forest. R package version 1.4.
  # Stekhoven D. J., & Buehlmann, P. (2012). MissForest - non-parametric missing value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.
  
  # data = complete-case dataset containing trait data. Must also contain a column containing species name called "species_name".
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # missLevel = proportion of values to remove from dataset (e.g. 0.1). Performs simulation/interation for 10% itervals (e.g. 10%, 20%, etc.)
  # ntrees = vector containing values of ntree to test
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Convert categorical traits to factors.
  data[, catTraits] <- lapply(data[, catTraits], as.factor)
  # Identify integer (count) traits, if any.   
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")
  # Determine original sample size for each trait. &&
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---   
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data[, vars])
  
  # Imputation prep. ---
  # Convert missLevel into vector of missingness proportions.  
  missingness <- seq(0.1, missLevel, by = 0.1)
  # Create list to hold the dataframes with simulated missingness from each iteration and level of missingness.    
  l_l_dfMiss <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration and level of missingness.    
  l_l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create lists to hold the error rates for each iteration. There will be n = missLevel number of error rates for each iteration (pertaining to each level of missingness).  
  l_l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..  
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.  
    l_l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # Log transformation of numerical traits prior to imputation.
    dfLog <- LogTransform(data, contTraits)
    # Create lists to hold the error rates, dfMiss, and dfImp at each missingness level.  
    l_dfMiss <- CreateNamedList(listLength = length(missingness), elementNames = missingness)  
    l_l_dfImp <- CreateNamedList(listLength = length(missingness), elementNames = missingness)  
    l_l_Error <- CreateNamedList(listLength = length(missingness), elementNames = missingness)
    # If phylogenetic imputation was chosen..  
    if(phyImp == T) {
      l_l_evs <- CreateNamedList(listLength = length(missingness), elementNames = missingness)    
    }
    
    # For every level of missingness...  
    for (s in 1:length(missingness)) {
      # MCAR simulation. ---  
      # Make a copy of complete-case dataframe (transformed) for introducing NAs.
      dfMiss <- dfLog
      # Introduce NAs depending on the level of miss. We also set the seed so that the NAs are added on top of the NAs already introduced at 0.1.
      set.seed(i)
      # For every level of missingness..
      for(l in 1:s) {
        # Remove 0.1 of the values.
        dfMiss <- prodNA(dfMiss[, vars], noNA = 0.1)
      }
      
      # Dataframe organization. ---
      # Add species_name information back.
      dfMiss$species_name <- dfLog$species_name
      # Bind missingness indicator columns and reorganize the dataframe.
      dfMiss <- BindAndOrganize(dfMiss, vars)
      # Make sure dfMiss and dfLog (the original data) are both ordered alphabetically.
      dfMiss <- dfMiss[order(dfMiss$species_name), ]
      dfLog <- dfLog[order(dfLog$species_name), ]
      # Ensure categorical traits are factor type.
      dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.factor)
      dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.factor)
      
      # Phylogenetic eigenvector decomposition. ---
      if(phyImp == T) {
        # Append eigenvectors to dfMiss and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors. 
        l_evs <- AppendEigenvectors(data = dfMiss, vars = vars, tree = tree, predictors = l_predictors)
        # Extract list of dfMiss.
        l_dfMissEV <- l_evs[[1]]
        # Extract updated list of predictors.
        l_EVPredictors <- l_evs[[2]]
      }
      
      # Identify outgroup(s) in trait datasets (this is the species that contains no trait data dfLog). &&
      outgroup <- dfLog$species_name[apply(dfLog[, traits], 1, function(x) all(is.na(x)))]
      # If found in dataframe.. &&
      if(length(outgroup) > 0){
        # Remove outgroup from dataframes.
        dfLog <- dfLog[!dfLog$species_name %in% outgroup, ]
        dfMiss <- dfMiss[!dfMiss$species_name %in% outgroup, ]
        if(phyImp == T){
          l_dfMissEV <- lapply(l_dfMissEV, function(x) x[!x$species_name %in% outgroup, ])
        }
      }
      
      # Imputation. ---
      # The ImputeRF() function entails a loop that uses different values of ntree (the number of trees grown in the forest) and returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.  
      # If phylogenetic imputation was selected..
      if(phyImp == T) {
        # Impute data using eigenvectors.
        impResult <- ImputeRF(dfTrue = dfLog, dfMissing = dfMiss, cols = vars, cont = contTraits, cat = catTraits, inter = intTraits, ntrees = c(100, 1000), predictors = l_EVPredictors, phyImp = T, l_dfMissing = l_dfMissEV)
      } else if(phyImp == F){
        # Impute data only using trait data.
        impResult <- ImputeRF(dfTrue = dfLog, dfMissing = dfMiss, cols = vars, cont = contTraits, cat = catTraits, inter = intTraits, ntrees = c(100, 1000), predictors = l_predictors)
      }
      
      # Result handling. ---  
      # Extract l_dfImp.  
      l_dfImp <- impResult[[1]]
      # Extract errorRates  .
      l_Error <- impResult[[2]]
      
      # Back-transforming data. ---
      # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
      dfOrig <- data[data$species_name %in% dfLog$species_name, ]
      # For every imputed dataframe..
      for(d in 1:length(l_dfImp)){
        # Take the dth dfImp.
        dfImp <- l_dfImp[[d]]
        # For every continuous trait..
        for(t in 1:length(contTraits)){
          # Take the tth continuous trait.
          trait <- contTraits[[t]]
          # Apply the Backtransform function for the trait in question.
          dfImpBTrait <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMiss, cols = trait)
          # Replace data in dfImp with backtransformed data.
          dfImp[, trait] <- dfImpBTrait[, trait]
        }
        # If there are any integer traits..
        if(length(intTraits) > 0){
          # Round to nearest whole number.
          dfImp[, intTraits] <- lapply(dfImp[, intTraits], function(x) as.integer(round(x)))
        }
        # Replace dfImp in l_dfImp with newly backtransformed dataset.
        l_dfImp[[d]] <- dfImp
      }
      
      # Append results to sth element of lists.
      l_dfMiss[[s]] <- dfMiss
      l_l_dfImp[[s]] <- l_dfImp
      l_l_Error[[s]] <- l_Error
      # If phylogenetic imputation was selected..
      if(phyImp == T) {
        l_l_evs[[s]] <- l_evs
      }
    } ## s
    # Append to ith elements of lists.
    l_l_dfMiss[[i]] <- l_dfMiss
    l_l_l_dfImp[[i]] <- l_l_dfImp
    l_l_l_Error[[i]] <- l_l_Error
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      l_l_l_evs[[i]] <- l_l_evs
    }
  } ## i
  
  # If phylogenetic imputation was selected..
  if(phyImp == T) { 
    # Create list to hold the results.
    l_results <- list(originalData = data, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, ntrees = ntrees, missingness = missingness, reps = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_l_evs, missingData = l_l_dfMiss, imputedData = l_l_l_dfImp, errorRates = l_l_l_Error)
  } else if(phyImp == F){
    # Create list to hold the results.
    l_results <- list(originalData = data, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, ntrees = ntrees, missingness = missingness, reps = int, predictors = l_predictors, missingData = l_l_dfMiss, imputedData = l_l_l_dfImp, errorRates = l_l_l_Error)
  }
  # Return l_results.
  return(l_results)
  
}

RFSimputeMAR <- function(data, raw, vars, int = 100, ntrees = c(100, 1000), phyImp = F, tree = NULL) {
  
  # Given a complete-case dataset, this function simulates missingness at random (MAR) and imputes values using the missForest() function in the "missForest" package. Determines error rates for numerical (MSE) and categorical variables (PFC).
  # Citations: Daniel J. Stekhoven (2013). missForest: Nonparametric Missing Value Imputation using Random Forest. R package version 1.4.
  # Stekhoven D. J., & Buehlmann, P. (2012). MissForest - non-parametric missing value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.
  
  # data = complete-case dataset containing trait data. Must also contain a column containing species name called "species_name".
  # raw* = original dataset containing trait data with missing values. Used to build logistic regression models and simulate data that are MAR in the complete-case dataset. Must also contain a column with species name information = "species_name"
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # ntrees = vector containing values of ntree to test
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Identify integer (count) traits, if any.
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")
  # Ensure categorical traits are factor type for logistic regression model building (data class required for glm function).
  data[, catTraits] <- lapply(data[, catTraits], as.factor)
  raw[, catTraits] <- lapply(raw[, catTraits], as.factor)
  # Determine original sample size for each trait. &&
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data[, vars])
  
  # Logistic regression fitting. ---
  # Apply FitLogReg() function to raw data to identify significant predictors of missingness for each trait.
  l_models <- FitLogReg(data = raw, cols = vars)
  # Now, let's identify traits that CAN be simulated MAR versus those that cannot (i.e. they must be simulated MCAR). Traits that must be simulated MCAR either 1) have no missing values in the original data or 2) do not have any predictors of missingness in the dataset upon fitting of the logistic regression models. IDMissPattern() also refines the logistic regression models (i.e. drops insignificant terms) through use of the RefineModels() function.
  l_missPatt <- IDMissPattern(data = raw, vars = vars, models = l_models)
  # Extract MCAR traits.
  traitsMCAR <- l_missPatt[[1]]
  # Extract MAR traits.
  traitsMAR <- l_missPatt[[2]]
  # Extract final MAR models.
  l_MARfinalModels <- l_missPatt[[3]]
  
  # Imputation prep. ---
  # Create lists to hold the error rates for each iteration.  
  l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the dataframes with simulated missingness from each iteration.    
  l_dfMissOrig <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.    
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration.    
  l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.    
    l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # MAR simulation. ---  
    # Create copy of complete-case dataframe (untransformed) to introduce NAs into.
    dfMissOrig <- data 
    # Create a list to hold the missingness proportion for each trait to be simulated MAR.  
    l_missingness <- CreateNamedList(listLength = length(traitsMAR), elementNames = traitsMAR)
    # For every trait that can be simulated MAR..
    for(t in 1:length(traitsMAR)){
      # Take tth trait.  
      trait <- traitsMAR[[t]]
      # Get the corresponding MAR model.
      index <- grep(pattern = trait, x = names(l_MARfinalModels))
      modelMAR <- l_MARfinalModels[[index]]
      # Simulate MAR data in complete-case dataset.
      res <- SimMAR(model = modelMAR, data = data)
      # Replace complete-case column in dfMissOrig with res.
      dfMissOrig[, trait] <- res
      # Get original sample size of trait. &&
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMissOrig after MAR simulation to determine actual number of NAs introduced. &&
      marN <- sum(is.na(dfMissOrig[[trait]])) - sum(is.na(data[[trait]]))
      # Divide by n to determine missingness percentage for trait and append to l_missingness. &&
      l_missingness[[t]] <- marN/n
    }
    
    # Dataframe organization. ---
    # Log transformation of numerical traits prior to imputation.  
    dfLog <- LogTransform(data, contTraits)
    # Introduce missing values from dfMissOrig into dfLog. !!!
    dfMiss <- as.data.frame(mapply(function(x, y) ifelse(is.na(x), x, y), x = dfMissOrig[, c("species_name", contTraits)], y = dfLog[, c("species_name", contTraits)], SIMPLIFY = F))
    # Add categorical traits back. !!!
    dfMiss <- merge(dfMiss, dfMissOrig[, c("species_name", catTraits)], by = "species_name")
    # Bind missingness indicator columns and reorganize the dataframe.
    dfMiss <- BindAndOrganize(dfMiss, vars)
    # Make sure dfMiss and dfLog (the original data) are both ordered alphabetically.
    dfMiss <- dfMiss[order(dfMiss$species_name), ]
    dfLog <- dfLog[order(dfLog$species_name), ]
    # Ensure categorical traits are factor type.
    dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.factor)
    dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.factor)
    # Finally, append traits that could not be simulated MAR (no NAs introduced based on fitted logistic regression models) to traitsMCAR.  
    traitsMCAR <- AppendMCAR(dfMiss, cols = traitsMAR, varsMCAR = traitsMCAR)
    # Update traitsMAR based on this result.  
    traitsMAR <- setdiff(traitsMAR, traitsMCAR)
    
    # Phylogenetic eigenvector decomposition. ---
    if(phyImp == T) {
      # Append eigenvectors to dfMiss and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors.
      l_evs <- AppendEigenvectors(data = dfMiss, vars = traitsMAR, tree = tree, predictors = l_predictors)
      # Extract list of dfMiss.
      l_dfMissEV <- l_evs[[1]]
      # Extract updated list of predictors.
      l_EVPredictors <- l_evs[[2]]
    }
    # Identify outgroup(s) in trait datasets (this is the species that contains no trait data dfLog). &&
    outgroup <- dfLog$species_name[apply(dfLog[, traits], 1, function(x) all(is.na(x)))]
    # If found in dataframe.. &&
    if(length(outgroup) > 0){
      # Remove outgroup from dataframes.
      dfLog <- dfLog[!dfLog$species_name %in% outgroup, ]
      dfMiss <- dfMiss[!dfMiss$species_name %in% outgroup, ]
      if(phyImp == T){
        # Also remove outgroup from l_dfMissEV.
        l_dfMissEV <- lapply(l_dfMissEV, function(x) x[!x$species_name %in% outgroup, ])
      }
    }
    
    # Imputation. ---   
    # The ImputeRF() function entails a loop that uses different values of ntree (the number of trees grown in the forest) and returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.
    # If phylogenetic imputation was selected..  
    if(phyImp == T) {
      # Impute data using eigenvectors.
      impResult <- ImputeRF(dfTrue = dfLog, dfMissing = dfMiss, cols = traitsMAR, cont = contTraits, cat = catTraits, inter = intTraits, ntrees = c(100, 1000), predictors = l_EVPredictors, phyImp = T, l_dfMissing = l_dfMissEV)
    } else if(phyImp == F){
      # Impute data only using trait data.
      impResult <- ImputeRF(dfTrue = dfLog, dfMissing = dfMiss, cols = traitsMAR, cont = contTraits, cat = catTraits, inter = intTraits, ntrees = c(100, 1000), predictors = l_predictors)
    }
    
    # Result handling. ---  
    # Extract l_dfImp.  
    l_dfImp <- impResult[[1]]
    # Extract errorRates  .
    l_Error <- impResult[[2]]
    
    # Back-transforming data. ---  
    # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfLog$species_name, ]
    # For every imputed dataframe..
    for(d in 1:length(l_dfImp)){
      # Take the dth dfImp.
      dfImp <- l_dfImp[[d]]
      # Apply BackTransform function to contTraits in dfImp.
      dfImp <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMiss, cols = traitsMAR[traitsMAR %in% contTraits])
      # If there are any integer traits..
      if(length(intTraits) > 0){
        # Round to nearest whole number.
        dfImp[, intTraits] <- lapply(dfImp[, intTraits], function(x) as.integer(round(x)))
      }
      # Replace dfImp in l_dfImp with newly backtransformed dataset.
      l_dfImp[[d]] <- dfImp
    }
    
    # Append results to ith element of lists.    
    l_dfMissOrig[[i]] <- dfMissOrig
    l_l_missingness[[i]] <- l_missingness      
    l_l_dfImp[[i]] <- l_dfImp
    l_l_Error[[i]] <- l_Error
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      l_l_evs[[i]] <- l_evs
    }
  } ## i
  
  # Final MCAR check. ---  
  # If there are any traits that could not be simulated MAR..
  if(length(traitsMCAR) > 0){
    # If the trait was identified in a previous iteration to be an MCAR variable, it is possible missing values were introduced in other iterations. These will likely be very low numbers of NAs and not enough for error rate analyses. So here we will subset l_l_Error to ensure it only contains values for varsMAR.
    l_l_Error <- lapply(l_l_Error, function(x) {
      # Identify error rates associated with traitsMAR.
      index <- which(names(x) %in% traitsMAR)
      # Remove traits that were identified as MCAR in later iterations.
      x <- x[index]
    })
    print("The following traits could not be simulated MAR:")
    print(unique(traitsMCAR))
  }
  
  # Average out the missingness for each trait. ---    
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(traitsMAR, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  # Subset to sample sizes for traits that could be simulated MAR. &&
  l_MARn <- l_sampleSizes[names(l_sampleSizes) %in% traitsMAR]
  # Calculate the average number of NAs introduced for each trait. &&
  avgNAs <- mapply(function(x, y) x * y, y = l_MARn, x = avgMiss)
  
  # If phylogenetic imputation was selected..    
  if(phyImp == T) { 
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, traitsNotSimulated = traitsMCAR, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, ntrees = ntrees, missingness = missingness, reps = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_evs, missingData = l_dfMissOrig, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  } else if(phyImp == F){
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, traitsNotSimulated = traitsMCAR, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, ntrees = ntrees, missingness = missingness, reps = int, predictors = l_predictors, missingData = l_dfMissOrig, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  }
  # Return l_results.  
  return(l_results)
  
}

RFSimputeMNAR <- function(data, vars, int = 100, ntrees = c(100, 1000), quantiles = NULL, directions = NULL, absolutes = NULL, categories = NULL, phyImp = F, tree = NULL) {
  
  # Given a complete-case dataset, this function simulates missingness not at random (MNAR) and imputes values using the missForest() function in the "missForest" package. Determines error rates for numerical (MSE) and categorical variables (PFC).
  # Citations: Daniel J. Stekhoven (2013). missForest: Nonparametric Missing Value Imputation using Random Forest. R package version 1.4.
  # Stekhoven D. J., & Buehlmann, P. (2012). MissForest - non-parametric missing value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.
  
  # data = complete-case dataset containing trait data. Must also contain a column containing species name called "species_name".
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # ntrees = vector containing values of ntree to test
  # quantiles = named numeric vector of quantiles at which to introduce NAs for each numerical trait  
  # directions = for each numerical trait, named character vector of directions ("greater" or "less") for which to introduce NAs with regards to the corresponding quantile  
  # absolutes = for each numerical trait, named logical vector indicating whether to take the absolute values of the data when determining MNAR thresholds  
  # categories = for each categorical trait, named character vector of categories in which to introduce NAs  
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---  
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Identify integer (count) traits, if any.   
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")
  # Ensure categorical traits are factor type for logistic regression model building (data class required for glm function).   
  data[, catTraits] <- lapply(data[, catTraits], as.factor)
  # Determine original sample size for each trait. &&
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data[, vars])
  
  # Imputation prep. ---
  # Create lists to hold the error rates for each iteration.
  l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the dataframes with simulated missingness for each trait and each iteration.     
  l_l_dfMissTraitOrig <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration.
  l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..   
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.
    l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # MNAR simulation. ---  
    # Create a list of complete-case dataframe to introduce NAs into (one for each trait).
    l_dfMissTraitOrig <- lapply(1:length(vars), function(x) data)
    # Name l_dfMissTraitOrig according to traits.
    names(l_dfMissTraitOrig) <- vars
    # Create an empty list to hold the missingness proportion for each trait.
    l_missingness <- CreateNamedList(listLength = length(vars), elementNames = vars)
    # For every trait to be simulated MNAR..
    for(t in 1:length(vars)){
      # Take tth trait.
      trait <- vars[[t]]
      # Take the corresponding dfMissTraitOrig for the trait.
      dfMissTraitOrig <- l_dfMissTraitOrig[[grep(trait, names(l_dfMissTraitOrig))]]
      # If trait is numeric..
      if(trait %in% contTraits){
        # Take quantile that corresponds to trait.
        quant <- quantiles[[grep(trait, names(quantiles))]]
        # Take direction that corresponds to trait.
        dir <- directions[[grep(trait, names(directions))]]
        # Take absolute indicator that corresponds to trait.
        absInd <- absolutes[[grep(trait, names(absolutes))]]
        # Simulate MNAR data in trait.
        res <- SimContMNAR(var = dfMissTraitOrig[[trait]], quantile = quant, direction = dir, absolute = absInd)
        # Replace complete-case column in dfMiss with res.
        dfMissTraitOrig[, trait] <- res
      } else if(trait %in% catTraits){
        # Take category that corresponds to trait.
        cate <- categories[[grep(trait, names(categories))]]
        # Simulate MNAR data in trait. 
        res <- SimCatMNAR(var = dfMissTraitOrig[[trait]], category = cate)
        # Replace complete-case column in dfMiss with res.
        dfMissTraitOrig[, trait] <- res
      } 
      # Get original sample size of trait.
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMiss after MNAR simulation to determine actual number of NAs introduced.
      marN <- sum(is.na(dfMissTraitOrig[[trait]])) - sum(is.na(data[[trait]]))
      # Divide by n to determine missingness percentage for trait and append to l_missingness.
      l_missingness[[t]] <- marN/n
      # Append dfMissTraitOrig with introduced NAs for the trait to l_dfMissTraitOrig.
      l_dfMissTraitOrig[[t]] <- dfMissTraitOrig
    }
    
    # Dataframe organization. ---  
    # Log transformation of numerical traits prior to imputation.  
    dfLog <- LogTransform(data = data, cols = contTraits)
    # For each numerical trait, introduce missing values from each dfMissOrig into copies of dfLog. !!!
    l_dfMissTrait <- lapply(l_dfMissTraitOrig, function(x) {
      x <- mapply(function(x, y) ifelse(is.na(x), x, y), x = x[, c("species_name", contTraits)], y = dfLog[, c("species_name", contTraits)], SIMPLIFY = F)
    })
    # Convert all elements back to dataframe format. !!!
    l_dfMissTrait <- lapply(l_dfMissTrait, as.data.frame)
    # To get original category data (and not numbers for each category as line above would do), merge categorical traits in each dataframe of l_dfMissTraitOrig back to each dataframe in l_dfMissTrait. !!!
    for(c in 1:length(l_dfMissTrait)){
      l_dfMissTrait[[c]] <- merge(l_dfMissTrait[[c]], l_dfMissTraitOrig[[c]][, c("species_name", catTraits)], by = "species_name")
    }
    # Bind missingness indicator columns and reorganize the dataframes.  
    l_dfMissTrait <- lapply(l_dfMissTrait, BindAndOrganize, vars = vars)
    # Make sure l_dfMissTrait and dfLog (the original data) are all ordered alphabetically by species_name.  
    l_dfMissTrait <- lapply(l_dfMissTrait, function(x) x[order(x$species_name), ])
    dfLog <- dfLog[order(dfLog$species_name), ]
    # Ensure categorical traits are factor type.
    dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.factor)
    l_dfMissTrait <- lapply(l_dfMissTrait, function(x){
      x[catTraits] <- lapply(x[catTraits], as.factor)
      return(x)})
    
    # Phylogenetic eigenvector decomposition. ---     
    if(phyImp == T) {
      # Append eigenvectors to dataframes in l_dfMissTrait and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors.
      l_evsTrait <- mapply(AppendEigenvectors, data = l_dfMissTrait, vars = vars, MoreArgs = list(tree = tree, predictors = l_predictors), SIMPLIFY = F)
      # Extract list of dfMissEV.
      l_dfMissEV <- sapply(l_evsTrait, function(x) x[[1]])
      names(l_dfMissEV) <- traits
      # Extract updated list of predictors.
      l_EVPredictors <- sapply(l_evsTrait, function(x) x[[2]])
    }
    # Identify outgroup(s) in trait datasets (this is the species that contains no trait data in dfLog).   
    outgroup <- dfLog$species_name[apply(dfLog[, vars], 1, function(x) all(is.na(x)))]
    # If found in dataframe..
    if(length(outgroup) > 0){
      # Remove outgroup from dataframes.  
      dfLog <- dfLog[!dfLog$species_name %in% outgroup, ]
      l_dfMissTrait <- lapply(l_dfMissTrait, function(x) x[!x$species_name %in% outgroup, ])
      # If phylogenetic imputation was selected..
      if(phyImp == T){
        # Remove outgroup from dataframes in l_dfMissEV.  
        l_dfMissEV <- lapply(l_dfMissEV, function(x) x[!x$species_name %in% outgroup, ])
      }
    }
    
    # Imputation. ---   
    # The ImputeRF() function entails a loop that uses different values of ntree (the number of trees grown in the forest) and returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.
    # If phylogenetic imputation was selected..  
    if(phyImp == T) {
      # Impute data using eigenvectors.
      impResult <- ImputeRF(dfTrue = dfLog, dfMissing = l_dfMissEV[[1]], cols = vars, cont = contTraits, cat = catTraits, inter = intTraits, ntrees = c(100, 1000), predictors = l_EVPredictors, phyImp = T, l_dfMissing = l_dfMissEV)
      # Result handling. ---  
      # Extract l_dfImp.
      l_dfImp <- impResult[[1]]
      # Extract errorRates.
      l_Error <- impResult[[2]]
    } else if(phyImp == F){
      # Impute data only using trait data.
      l_impResult <- mapply(ImputeRF, dfMissing = l_dfMissTrait, cols = vars, MoreArgs = list(dfTrue = dfLog, cont = contTraits, cat = catTraits, inter = intTraits, ntrees = c(100, 1000), predictors = l_predictors), SIMPLIFY = F)
      # Since now we have a list of impResult, create new lists to hold the imputed dataframes and error rates.  
      l_ImpTrait <- lapply(l_impResult, function(x) x[[1]])
      l_Error <- sapply(l_impResult, function(x) x[[2]])
      # Rename the list.
      names(l_Error) <- vars
      # Imputed dataset handling. ---  
      # Creation of l_dfImp so that there is one dataframe per parameter.  
      dfImp <- as.data.frame(dfLog[, c("species_name")])
      names(dfImp) <- "species_name"
      # Make a list to hold the dataframes for each value of ntree.     
      l_dfImp <- lapply(1:length(ntrees), function(x) dfImp)
      names(l_dfImp) <- ntrees
      # For every trait...
      for (t in 1:length(vars)) { 
        # Take the tth trait.
        trait <- vars[[t]]
        # Extract the tth element in l_ImpTrait
        l_dfImputedTrait <- l_ImpTrait[[t]]
        # Get the name of the corresponding missingness info column.
        missCol <- grep(pattern = paste(trait, "_NA", sep = ""), x = colnames(l_dfImputedTrait[[1]]))
        missCol <- names(l_dfImputedTrait[[1]])[missCol]
        # For each value of ntree..
        for(n in 1:length(l_dfImputedTrait)){
          # Extract the nth dataframe in l_dfImputedTrait.
          dfImputedTrait <- l_dfImputedTrait[[n]]
          # Subset to only contain species_name, trait info, and missingness indicator column.
          dfImputedTrait <- dfImputedTrait[, c("species_name", trait, missCol)]
          # Extract the nth dataframe in l_dfImputedParam.
          dfImp <- l_dfImp[[n]]
          # Merge with dfImputedTrait.
          dfImp <- merge(dfImp, dfImputedTrait, by = "species_name")
          # Replace nth element of l_dfImputedParam.
          l_dfImp[[n]] <- dfImp
        }
      }
    }
    
    # Back-transforming data. ---
    # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfLog$species_name, ]
    # For every imputed dataframe..
    for(d in 1:length(l_dfImp)){
      # Take the dth dfImp.
      dfImp <- l_dfImp[[d]]
      # For every continuous trait..
      for(t in 1:length(contTraits)){
        # Take the tth continuous trait.
        trait <- contTraits[[t]]
        # Take the corresponding dfMissTrait.
        dfMissTrait <- l_dfMissTrait[[grep(trait, names(l_dfMissTrait))]]
        # Apply the Backtransform function for the trait in question.
        dfImpBTrait <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMissTrait, cols = trait)
        # Replace data in dfImp with backtransformed data.
        dfImp[, trait] <- dfImpBTrait[, trait]
      }
      # If there are any integer traits..
      if(length(intTraits) > 0){
        # Round to nearest whole number.
        dfImp[, intTraits] <- lapply(dfImp[, intTraits], function(x) as.integer(round(x)))
      }
      # Replace dfImp in l_dfImp with newly backtransformed dataset.
      l_dfImp[[d]] <- dfImp
    }
    
    # Append results to ith element of lists.
    l_l_dfMissTraitOrig[[i]] <- l_dfMissTraitOrig
    l_l_missingness[[i]] <- l_missingness
    l_l_dfImp[[i]] <- l_dfImp
    l_l_Error[[i]] <- l_Error
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      l_l_evs[[i]] <- l_evsTrait
    }
  } ## i
  
  # Average out the missingness for each trait. ---  
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(vars, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  # Calculate the average number of NAs introduced for each trait.
  avgNAs <- mapply(function(x, y) x * y, y = l_sampleSizes, x = avgMiss)
  
  # If phylogenetic imputation was selected..
  if(phyImp == T) { 
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, traits = vars, numeric = contTraits, categorical = catTraits, integer = intTraits, sampleSizes = l_sampleSizes, quantiles = quantiles, directions = directions, absolutes = absolutes, categories = categories, ntrees = ntrees, missingness = missingness, int = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_evs, missingData = l_l_dfMissTraitOrig, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  } else if(phyImp == F){
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, traits = vars, numeric = contTraits, categorical = catTraits, integer = intTraits, sampleSizes = l_sampleSizes, quantiles = quantiles, directions = directions, absolutes = absolutes, categories = categories, ntrees = ntrees, missingness = missingness, reps = int, predictors = l_predictors, missingData = l_l_dfMissTraitOrig, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  }
  # Return l_results.
  return(l_results)
}

MICESimputeMCAR <- function(data, vars, int = 100, missLevel = 0.1, mSets = c(5, 10, 40), phyImp = F, tree = NULL) {
  
  # Given a complete-case dataset, this function simulates missingness completely at random (MCAR) and imputes values using the mice() function in the "mice" package. Determines error rates for numerical (MSE) and categorical variables (PFC).
  # Citations: van Buuren S, Groothuis-Oudshoorn K (2011). “mice: Multivariate Imputation by Chained Equations in R.” Journal of Statistical Software, 45(3), 1-67. https://www.jstatsoft.org/v45/i03/.
  # R package version 3.13.0. https://cran.r-project.org/web/packages/mice/mice.pdf
  # Vignettes consulted: Gerko Vink and Stef van Buuren. miceVignettes. https://www.gerkovink.com/miceVignettes/
  
  # data = complete-case dataset containing trait data. Must also contain a column containing species name called "species_name".
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # missLevel = proportion of values to remove from dataset (e.g. 0.1). Performs simulation/interation for 10% itervals (e.g. 10%, 20%, etc.)  
  # mSets = vector containing values of m to test (number of multiply imputed dataframes)
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Convert categorical traits to factors.
  data[, catTraits] <- lapply(data[, catTraits], as.factor)
  # Identify integer (count) traits, if any.   
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")
  # Determine original sample size for each trait. &&
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---   
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data[, vars])
  
  # Imputation prep. ---
  # Convert missLevel into vector of missingness proportions.  
  missingness <- seq(0.1, missLevel, by = 0.1)
  # Create list to hold the dataframes with simulated missingness from each iteration and level of missingness.    
  l_l_dfMiss <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration and level of missingness.    
  l_l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create lists to hold the error rates for each iteration. There will be n = missLevel number of error rates for each iteration (pertaining to each level of missingness).  
  l_l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..  
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.  
    l_l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
    l_l_dfMissEV <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # Log transformation of numerical traits prior to imputation.
    dfLog <- LogTransform(data, contTraits)
    # Create lists to hold the error rates, dfMiss, and dfImp at each missingness level.  
    l_dfMiss <- CreateNamedList(listLength = length(missingness), elementNames = missingness)    
    l_l_dfImp <- CreateNamedList(listLength = length(missingness), elementNames = missingness)    
    l_l_Error <- CreateNamedList(listLength = length(missingness), elementNames = missingness)  
    # If phylogenetic imputation was chosen..  
    if(phyImp == T) {
      l_l_evs <- CreateNamedList(listLength = length(missingness), elementNames = missingness)    
      l_dfMissEV <- CreateNamedList(listLength = length(missingness), elementNames = missingness)
    }
    
    # For every level of missingness...  
    for (s in 1:length(missingness)) {
      # MCAR simulation. ---  
      # Make a copy of complete-case dataframe (transformed) for introducing NAs.
      dfMiss <- dfLog
      # Introduce NAs depending on the level of miss. We also set the seed so that the NAs are added on top of the NAs already introduced at 0.1.
      set.seed(i)
      # For every level of missingness..
      for(l in 1:s) {
        # Remove 0.1 of the values.
        dfMiss <- prodNA(dfMiss[, vars], noNA = 0.1)
      }
      
      # Dataframe organization. ---
      # Add species_name information back.  
      dfMiss$species_name <- dfLog$species_name
      # Bind and organize missingness indicator columns.
      dfMiss <- BindAndOrganize(dfMiss, vars)
      # Make sure dfMiss and dfLog (the original data) are both ordered alphabetically.
      dfMiss <- dfMiss[order(dfMiss$species_name), ]
      dfLog <- dfLog[order(dfLog$species_name), ]
      # Ensure categorical traits are factor type.
      dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.factor)
      dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.factor)
      
      # Phylogenetic eigenvector decomposition. ---
      if(phyImp == T) {
        # Append eigenvectors to dfMiss and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors. 
        l_evs <- AppendEigenvectors(data = dfMiss, vars = vars, tree = tree, predictors = l_predictors)
        # Extract list of dfMiss.
        l_dfMissEV <- l_evs[[1]]
        # Extract updated list of predictors.
        l_EVPredictors <- l_evs[[2]]
        # Dataframe list collapse. ---
        # MICE requires a predictor matrix for imputation. So, we will first merge all of the dataframe in l_dfMissEV into one dataframe. This is so we can identify the eigenvector columns needed for imputation of each trait and update this in the predictor matrix. This should also same computation time since we are just imputing one dataframe and not a list of dataframes. 
        dfMissEV <- l_dfMissEV[[1]]
        # Identify the eigenvector columns and remove those dfMiss.
        evIndex <- grep(pattern = "V_", colnames(dfMissEV))
        dfMissEV <- dfMissEV[, -evIndex]
        # Also remove the missingness indicator columns.
        missIndex <- grep(pattern = "_NA", colnames(dfMissEV))
        dfMissEV <- dfMissEV[, -missIndex]
        # Now, let's merge dfMissEV with the corresponding eigenvector columns in each dataframe in l_dfMissEV.
        for(e in 1:length(l_dfMissEV)){
          # Take the eth l_dfMissEV.
          dfMisseth <- l_dfMissEV[[e]]
          # Identify eigenvector columns in eth dataframe.
          index <- grep("V_", colnames(dfMisseth))
          evCols <- colnames(dfMisseth)[index]
          # Merge these columns with dfMisseth (including species_name).
          dfMissEV <- merge(dfMissEV, dfMisseth[, c("species_name", evCols)], by = "species_name")
        }
      }
      # Identify outgroup(s) in trait datasets (this is the species that contains no trait data dfLog). &&
      outgroup <- dfLog$species_name[apply(dfLog[, traits], 1, function(x) all(is.na(x)))]
      # If found in dataframe.. &&
      if(length(outgroup) > 0){
        # Remove outgroup from dataframes.
        dfLog <- dfLog[!dfLog$species_name %in% outgroup, ]
        dfMiss <- dfMiss[!dfMiss$species_name %in% outgroup, ]
        # If phylogenetic imputation was selected..
        if(phyImp == T) {
          l_dfMissEV <- lapply(l_dfMissEV, function(x) x[!x$species_name %in% outgroup, ])
          dfMissEV <- dfMissEV[!dfMissEV$species_name %in% outgroup, ]
        }
      }
      
      # Predictor matrix creation. ---   
      # If phylogenetic imputation wasn't chosen..
      if(phyImp == F){
        # Create predictor matrix for use in MICE imputation (based on results of our trait predictor screening).
        matPredictors <- CreatePredictorMatrix(cols = vars, predictors = l_predictors, dfMissing = dfMiss)
        # If phylogenetic imputation was chosen..
      } else if(phyImp == T){
        # Create predictor matrix for use in MICE imputation with appended eigenvectors. We can indicate the column names of dfMissEV (excluding species_name) as names of variables to consider in imputation process.
        matPredictors <- CreatePredictorMatrix(dfMissing = dfMissEV, cols = colnames(dfMissEV)[-1], predictors = l_EVPredictors)
      }
      
      # Imputation. ---   
      # The ImputeMICE() function entails a loop that uses different values of m (the number of multiply imputed datasets) and returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.  
      # If phylogenetic imputation wasn't selected..
      if(phyImp == F) {
        # Impute data only using trait data.
        impResult <- ImputeMICE(dfTrue = dfLog, dfMissing = dfMiss, cols = vars, cont = contTraits, cat = catTraits, inter = intTraits, mSets = c(5, 10, 40), matPredictors = matPredictors)
        # If phylogenetic imputation was selected..
      } else if(phyImp == T){
        # Bind back the missing columns.   
        dfMissEV <- bind_shadow(dfMissEV, only_miss = T)
        # Impute data using eigenvectors.
        impResult <- ImputeMICE(dfTrue = dfLog, dfMissing = dfMissEV, cols = vars, cont = contTraits, cat = catTraits, inter = intTraits, mSets = c(5, 10, 40), matPredictors = matPredictors)
      }
      
      # Result handling. ---
      # Extract l_dfImp.
      l_dfImp <- impResult[[1]]
      # Extract errorRates  .
      l_Error <- impResult[[2]]
      
      # Back-transforming data. ---  
      # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
      dfOrig <- data[data$species_name %in% dfLog$species_name, ]
      # For every imputed dataframe..
      for(d in 1:length(l_dfImp)){
        # Take the dth dfImp.
        dfImp <- l_dfImp[[d]]
        # For every continuous trait..
        for(t in 1:length(contTraits)){
          # Take the tth continuous trait.
          trait <- contTraits[[t]]
          # Apply the Backtransform function for the trait in question.
          dfImpBTrait <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMiss, cols = trait)
          # Replace data in dfImp with backtransformed data.
          dfImp[, trait] <- dfImpBTrait[, trait]
        }
        # If there are any integer traits..
        if(length(intTraits) > 0){
          # Round to nearest whole number.
          dfImp[, intTraits] <- lapply(dfImp[, intTraits], function(x) as.integer(round(x)))
        }
        # Replace dfImp in l_dfImp with newly backtransformed dataset.
        l_dfImp[[d]] <- dfImp
      }
  
      # Append results to sth element of lists.
      l_dfMiss[[s]] <- dfMiss
      l_l_dfImp[[s]] <- l_dfImp
      l_l_Error[[s]] <- l_Error
      # If phylogenetic imputation was selected..
      if(phyImp == T) {
        l_l_evs[[s]] <- l_evs
        l_dfMissEV[[s]] <- dfMissEV
      }
    } ## s
    
    # Append to ith elements of lists.    
    l_l_dfMiss[[i]] <- l_dfMiss
    l_l_l_dfImp[[i]] <- l_l_dfImp
    l_l_l_Error[[i]] <- l_l_Error
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      l_l_l_evs[[i]] <- l_l_evs
      l_l_dfMissEV[[i]] <- l_dfMissEV
    }
  } ## i
  
  # If phylogenetic imputation was selected..
  if(phyImp == T) { 
    # Create list to hold the results.
    l_results <- list(originalData = data, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, m = mSets, missingness = missingness, int = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_l_evs, missingData = l_l_dfMiss, eigenTraitData = l_l_dfMissEV, imputedData = l_l_l_dfImp, errorRates = l_l_l_Error)
  } else if(phyImp == F){
    # Create list to hold the results.
    l_results <- list(originalData = data, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, m = mSets, predictorMatrix = matPredictors, missingness = missingness, reps = int, predictors = l_predictors, missingData = l_l_dfMiss, imputedData = l_l_l_dfImp, errorRates = l_l_l_Error)
  }
  # Return l_results.
  return(l_results)
  
}

MICESimputeMAR <- function(data, raw, vars, int = 100, mSets = c(5, 10, 40), phyImp = F, tree = NULL) {
  
  # Given a complete-case dataset, this function simulates missingness at random (MAR) and imputes values using the mice() function in the "mice" package. Determines error rates for numerical (MSE) and categorical variables (PFC).
  # Citations: van Buuren S, Groothuis-Oudshoorn K (2011). “mice: Multivariate Imputation by Chained Equations in R.” Journal of Statistical Software, 45(3), 1-67. https://www.jstatsoft.org/v45/i03/.
  # R package version 3.13.0. https://cran.r-project.org/web/packages/mice/mice.pdf
  # Vignettes consulted: Gerko Vink and Stef van Buuren. miceVignettes. https://www.gerkovink.com/miceVignettes/
  
  # data = complete-case dataset containing trait data. Must also contain a column containing species name called "species_name".
  # raw = original dataset containing trait data with missing values. Used to build logistic regression models and simulate data that are MAR in the complete-case dataset. Must also contain a column with species name information = "species_name"  
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # mSets = vector containing values of m to test (number of multiply imputed dataframes)
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Identify integer (count) traits, if any.   
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")
  # Ensure categorical traits are factor type for logistic regression model building (data class required for glm function).  
  data[, catTraits] <- lapply(data[, catTraits], as.factor)
  raw[, catTraits] <- lapply(raw[, catTraits], as.factor)
  # Determine original sample size for each trait. &&
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data[, vars])
  
  # Logistic regression fitting. ---  
  # Apply FitLogReg() function to raw data to identify significant predictors of missingness for each trait.
  l_models <- FitLogReg(data = raw, cols = vars)
  # Now, let's identify traits that CAN be simulated MAR versus those that cannot (i.e. they must be simulated MCAR). Traits that must be simulated MCAR either 1) have no missing values in the original data or 2) do not have any predictors of missingness in the dataset upon fitting of the logistic regression models. IDMissPattern() also refines the logistic regression models (i.e. drops insignificant terms) through use of the RefineModels() function.
  l_missPatt <- IDMissPattern(data = raw, vars = vars, models = l_models)
  # Extract MCAR traits.
  traitsMCAR <- l_missPatt[[1]]
  # Extract MAR traits.
  traitsMAR <- l_missPatt[[2]]
  # Extract final MAR models.
  l_MARfinalModels <- l_missPatt[[3]]
  
  # Imputation prep. ---
  # Create lists to hold the error rates for each iteration.  
  l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the dataframes with simulated missingness from each iteration.    
  l_dfMissOrig <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.    
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration.    
  l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.    
    l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
    l_dfMissEV <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # MAR simulation. ---  
    # Create copy of complete-case dataframe (untransformed) to introduce NAs into.
    dfMissOrig <- data
    # Create a list to hold the missingness proportion for each trait to be simulated MAR.  
    l_missingness <- CreateNamedList(listLength = length(traitsMAR), elementNames = traitsMAR)
    # For every trait that can be simulated MAR..
    for(t in 1:length(traitsMAR)){
      # Take tth trait.  
      trait <- traitsMAR[[t]]
      # Get the corresponding MAR model.
      index <- grep(pattern = trait, x = names(l_MARfinalModels))
      modelMAR <- l_MARfinalModels[[index]]
      # Simulate MAR data in complete-case dataset.
      res <- SimMAR(model = modelMAR, data = data)
      # Replace complete-case column in dfMissOrig with res.
      dfMissOrig[, trait] <- res
      # Get original sample size of trait.
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMissOrig after MAR simulation to determine actual number of NAs introduced.
      marN <- sum(is.na(dfMissOrig[[trait]])) - sum(is.na(data[[trait]]))
      # Divide by n to determine missingness percentage for trait and append to l_missingness.
      l_missingness[[t]] <- marN/n
      
    }
    
    # Dataframe organization. ---
    # Log transformation of numerical traits prior to imputation.  
    dfLog <- LogTransform(data, contTraits)
    # Introduce missing values from dfMissOrig into dfLog. !!!
    dfMiss <- as.data.frame(mapply(function(x, y) ifelse(is.na(x), x, y), x = dfMissOrig[, c("species_name", contTraits)], y = dfLog[, c("species_name", contTraits)], SIMPLIFY = F))
    # Add categorical traits back. !!!
    dfMiss <- merge(dfMiss, dfMissOrig[, c("species_name", catTraits)], by = "species_name")
    # Bind missingness indicator columns and reorganize the dataframe.
    dfMiss <- BindAndOrganize(dfMiss, vars)
    # Make sure dfMiss and dfLog (the original data) are both ordered alphabetically.
    dfMiss <- dfMiss[order(dfMiss$species_name), ]
    dfLog <- dfLog[order(dfLog$species_name), ]
    # Ensure categorical traits are factor type.
    dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.factor)
    dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.factor)
    # Finally, append traits that could not be simulated MAR (no NAs introduced based on fitted logistic regression models) to traitsMCAR.  
    traitsMCAR <- AppendMCAR(dfMiss, cols = traitsMAR, varsMCAR = traitsMCAR)
    # Update traitsMAR based on this result.  
    traitsMAR <- setdiff(traitsMAR, traitsMCAR)
    
    # Phylogenetic eigenvector decomposition. ---
    if(phyImp == T) {
      # Append eigenvectors to dfMiss and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors.
      l_evs <- AppendEigenvectors(data = dfMiss, vars = vars, tree = tree, predictors = l_predictors)
      # Extract list of dfMiss.
      l_dfMissEV <- l_evs[[1]]
      # Extract updated list of predictors.
      l_EVPredictors <- l_evs[[2]]
      # Dataframe list collapse. ---   
      # MICE requires a predictor matrix for imputation. So, we will first merge all of the dataframe in l_dfMissEV into one dataframe. This is so we can identify the eigenvector columns needed for imputation of each trait and update this in the predictor matrix. This should also same computation time since we are just imputing one dataframe and not a list of dataframes. 
      dfMissEV <- l_dfMissEV[[1]]
      # Identify the eigenvector columns and remove those dfMiss.
      evIndex <- grep(pattern = "V_", colnames(dfMissEV))
      dfMissEV <- dfMissEV[, -evIndex]
      # Also remove the missingness indicator columns.
      missIndex <- grep(pattern = "_NA", colnames(dfMissEV))
      dfMissEV <- dfMissEV[, -missIndex]
      # Now, let's merge dfMissEV with the corresponding eigenvector columns in each dataframe in l_dfMissEV.
      for(e in 1:length(l_dfMissEV)){
        # Take the eth l_dfMissEV.
        dfMisseth <- l_dfMissEV[[e]]
        # Identify eigenvector columns in eth dataframe.
        index <- grep("V_", colnames(dfMisseth))
        evCols <- colnames(dfMisseth)[index]
        # Merge these columns with dfMisseth (including species_name).
        dfMissEV <- merge(dfMissEV, dfMisseth[, c("species_name", evCols)], by = "species_name")
      }
    }
    # Identify outgroup(s) in trait datasets (this is the species that contains no trait data dfLog). &&
    outgroup <- dfLog$species_name[apply(dfLog[, vars], 1, function(x) all(is.na(x)))]
    # If found in dataframe.. &&
    if(length(outgroup) > 0){
      # Remove outgroup from dataframes.
      dfLog <- dfLog[!dfLog$species_name %in% outgroup, ]
      dfMiss <- dfMiss[!dfMiss$species_name %in% outgroup, ]
      # If phylogenetic imputation was selected.. &&   
      if(phyImp == T) {
        l_dfMissEV <- lapply(l_dfMissEV, function(x) x[!x$species_name %in% outgroup, ])
        dfMissEV <- dfMissEV[!dfMissEV$species_name %in% outgroup, ]
      }
    }
    
    # Predictor matrix creation. ---   
    # If phylogenetic imputation wasn't chosen..
    if(phyImp == F){
      # Create predictor matrix for use in MICE imputation (based on results of our trait predictor screening).
      matPredictors <- CreatePredictorMatrix(cols = vars, predictors = l_predictors, dfMissing = dfMiss)
      # If phylogenetic imputation was chosen..
    } else if(phyImp == T){
      # Create predictor matrix for use in MICE imputation with appended eigenvectors. We can indicate the column names of dfMissEV (excluding species_name) as names of variables to consider in imputation process.
      matPredictors <- CreatePredictorMatrix(dfMissing = dfMissEV, cols = colnames(dfMissEV)[-1], predictors = l_EVPredictors)
      # We can set all of the rows for the eigenvectors to 0 as they themselves will not be imputed.
      eigens <- grep(pattern = "V_", rownames(matPredictors))
      matPredictors[eigens, ] <- 0
    }
    
    # Imputation. ---   
    # The ImputeMICE() function entails a loop that uses different values of m (the number of multiply imputed datasets) and returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.
    # If phylogenetic imputation wasn't selected..  
    if(phyImp == F) {
      # Impute data only using trait data.
      impResult <- ImputeMICE(dfTrue = dfLog, dfMissing = dfMiss, cols = traitsMAR, cont = contTraits, cat = catTraits, inter = intTraits, mSets = c(5, 10, 40), matPredictors = matPredictors)
      # If phylogenetic imputation was selected..
    } else if(phyImp == T){
      # Bind back the missing columns.   
      dfMissEV <- bind_shadow(dfMissEV, only_miss = T)
      # Impute data using eigenvectors.
      impResult <- ImputeMICE(dfTrue = dfLog, dfMissing = dfMissEV, cols = traitsMAR, cont = contTraits, cat = catTraits, inter = intTraits, mSets = c(5, 10, 40), matPredictors = matPredictors)
    }
    # Result handling. ---  
    # Extract l_dfImp.  
    l_dfImp <- impResult[[1]]
    # Extract errorRates  .
    l_Error <- impResult[[2]]
    
    # Back-transforming data. ---  
    # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfLog$species_name, ]
    # For every imputed dataframe..
    for(d in 1:length(l_dfImp)){
      # Take the dth dfImp.
      dfImp <- l_dfImp[[d]]
      # Apply BackTransform function to contTraits in dfImp.
      dfImp <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMiss, cols = traitsMAR[traitsMAR %in% contTraits])
      # If there are any integer traits..
      if(length(intTraits) > 0){
        # Round to nearest whole number.
        dfImp[, intTraits] <- lapply(dfImp[, intTraits], function(x) as.integer(round(x)))
      }
      # Replace dfImp in l_dfImp with newly backtransformed dataset.
      l_dfImp[[d]] <- dfImp
    }
    
    # Append results to ith element of lists.    
    l_dfMissOrig[[i]] <- dfMiss
    l_l_missingness[[i]] <- l_missingness      
    l_l_dfImp[[i]] <- l_dfImp
    l_l_Error[[i]] <- l_Error
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      l_l_evs[[i]] <- l_evs
      l_dfMissEV[[i]] <- dfMissEV
    }
  } ## i
  
  # Final MCAR check. ---  
  # If there are any traits that could not be simulated MAR..
  if(length(traitsMCAR) > 0){
    # If the trait was identified in a previous iteration to be an MCAR variable, it is possible missing values were introduced in other iterations. These will likely be very low numbers of NAs and not enough for error rate analyses. So here we will subset l_l_Error to ensure it only contains values for varsMAR.
    l_l_Error <- lapply(l_l_Error, function(x) {
      lapply(x, function(e){
        # Identify error rates associated with traitsMAR.
        index <- which(names(e) %in% traitsMAR)
        # Remove traits that were identified as MCAR in later iterations.
        e <- e[index]
      })
    })
    print("The following traits could not be simulated MAR:")
    print(unique(traitsMCAR))
  }
  
  # Average out the missingness for each trait. ---    
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(traitsMAR, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  # Subset to sample sizes for traits that could be simulated MAR.
  l_MARn <- l_sampleSizes[names(l_sampleSizes) %in% traitsMAR]
  # Calculate the average number of NAs introduced for each trait.
  avgNAs <- mapply(function(x, y) x * y, y = l_MARn, x = avgMiss)
  
  # If phylogenetic imputation was selected..
  if(phyImp == T) { 
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, traitsNotSimulated = traitsMCAR, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, m = mSets, predictorMatrix = matPredictors, missingness = missingness, reps = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_evs, missingData = l_dfMissOrig, eigenTraitData = l_dfMissEV, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  } else if(phyImp == F){
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, traitsNotSimulated = traitsMCAR, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, m = mSets, predictorMatrix = matPredictors, missingness = missingness, reps = int, predictors = l_predictors, missingData = l_dfMissOrig, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  }
  # Return l_results.  
  return(l_results)
  
}

MICESimputeMNAR <- function(data, vars, int = 100, mSets = c(5, 10, 40), quantiles = NULL, directions = NULL, absolutes = NULL, categories = NULL, phyImp = F, tree = NULL) {
  
  # Given a complete-case dataset, this function simulates missingness not at random (MNAR) and imputes values using the mice() function in the "mice" package. Determines error rates for numerical (MSE) and categorical variables (PFC).
  # Citations: van Buuren S, Groothuis-Oudshoorn K (2011). “mice: Multivariate Imputation by Chained Equations in R.” Journal of Statistical Software, 45(3), 1-67. https://www.jstatsoft.org/v45/i03/.
  # R package version 3.13.0. https://cran.r-project.org/web/packages/mice/mice.pdf
  # Vignettes consulted: Gerko Vink and Stef van Buuren. miceVignettes. https://www.gerkovink.com/miceVignettes/
  
  # data = complete-case dataset containing trait data. Must also contain a column containing species name called "species_name".
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # mSets = vector containing values of m to test (number of multiply imputed dataframes)
  # quantiles = named numeric vector of quantiles at which to introduce NAs for each numerical trait  
  # directions = for each numerical trait, named character vector of directions ("greater" or "less") for which to introduce NAs with regards to the corresponding quantile  
  # absolutes = for each numerical trait, named logical vector indicating whether to take the absolute values of the data when determining MNAR thresholds  
  # categories = for each categorical trait, named character vector of categories in which to introduce NAs  
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---  
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Identify integer (count) traits, if any.   
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")
  # Ensure categorical traits are factor type for logistic regression model building (data class required for glm function).   
  data[, catTraits] <- lapply(data[, catTraits], as.factor)
  # Determine original sample size for each trait. &&
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---   
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data[, vars])
  
  # Imputation prep. ---  
  # Create lists to hold the error rates for each iteration.  
  l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the dataframes with simulated missingness for each trait and each iteration.  
  l_l_dfMissTraitOrig <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.  
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration.    
  l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..   
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.
    l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
    l_dfMissEV <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # MNAR simulation. ---  
    # Create a list of complete-case dataframe to introduce NAs into (one for each trait).
    l_dfMissTraitOrig <- lapply(1:length(vars), function(x) data)
    # Name l_dfMissTraitOrig according to traits.
    names(l_dfMissTraitOrig) <- vars
    # Create an empty list to hold the missingness proportion for each trait.
    l_missingness <- CreateNamedList(listLength = length(vars), elementNames = vars)
    # For every trait to be simulated MNAR..
    for(t in 1:length(vars)){
      # Take tth trait.
      trait <- vars[[t]]
      # Take the corresponding dfMissTraitOrig for the trait.
      dfMissTraitOrig <- l_dfMissTraitOrig[[grep(trait, names(l_dfMissTraitOrig))]]
      # If trait is numeric..
      if(trait %in% contTraits){
        # Take quantile that corresponds to trait.
        quant <- quantiles[[grep(trait, names(quantiles))]]
        # Take direction that corresponds to trait.
        dir <- directions[[grep(trait, names(directions))]]
        # Take absolute indicator that corresponds to trait.
        absInd <- absolutes[[grep(trait, names(absolutes))]]
        # Simulate MNAR data in trait.
        res <- SimContMNAR(var = dfMissTraitOrig[[trait]], quantile = quant, direction = dir, absolute = absInd)
        # Replace complete-case column in dfMiss with res.
        dfMissTraitOrig[, trait] <- res
      } else if(trait %in% catTraits){
        # Take category that corresponds to trait.
        cate <- categories[[grep(trait, names(categories))]]
        # Simulate MNAR data in trait. 
        res <- SimCatMNAR(var = dfMissTraitOrig[[trait]], category = cate)
        # Replace complete-case column in dfMiss with res.
        dfMissTraitOrig[, trait] <- res
      } 
      # Get original sample size of trait.
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMiss after MNAR simulation to determine actual number of NAs introduced.
      marN <- sum(is.na(dfMissTraitOrig[[trait]])) - sum(is.na(data[[trait]]))
      # Divide by n to determine missingness percentage for trait and append to l_missingness.
      l_missingness[[t]] <- marN/n
      # Append dfMissTraitOrig with introduced NAs for the trait to l_dfMissTraitOrig.
      l_dfMissTraitOrig[[t]] <- dfMissTraitOrig
    }
    
    # Dataframe organization. ---  
    # Log transformation of numerical traits prior to imputation.
    dfLog <- LogTransform(data = data, cols = contTraits)
    # For each numerical trait, introduce missing values from each dfMissOrig into copies of dfLog.  
    l_dfMissTrait <- lapply(l_dfMissTraitOrig, function(x) {
      x <- mapply(function(x, y) ifelse(is.na(x), x, y), x = x[, c("species_name", contTraits)], y = dfLog[, c("species_name", contTraits)], SIMPLIFY = F)
    })
    # Convert all elements back to dataframe format.  
    l_dfMissTrait <- lapply(l_dfMissTrait, as.data.frame)
    # To get original category data (and not numbers for each category as line above would do), merge categorical traits in each dataframe of l_dfMissTraitOrig back to each dataframe in l_dfMissTrait.  
    for(c in 1:length(l_dfMissTrait)){
      l_dfMissTrait[[c]] <- merge(l_dfMissTrait[[c]], l_dfMissTraitOrig[[c]][, c("species_name", catTraits)], by = "species_name")
    }
    # Bind missingness indicator columns and reorganize the dataframe.  
    l_dfMissTrait <- lapply(l_dfMissTrait, BindAndOrganize, vars = vars)
    # Make sure l_dfMissTrait and dfLog (the original data) are all ordered alphabetically by species_name.  
    l_dfMissTrait <- lapply(l_dfMissTrait, function(x) x[order(x$species_name), ])
    dfLog <- dfLog[order(dfLog$species_name), ]
    # Ensure categorical traits are factor type.
    dfLog[, catTraits] <- lapply(dfLog[, catTraits], as.factor)
    l_dfMissTrait <- lapply(l_dfMissTrait, function(x){
      x[catTraits] <- lapply(x[catTraits], as.factor)
      return(x)})
    
    # Phylogenetic eigenvector decomposition. ---     
    if(phyImp == T) {
      # Append eigenvectors to dataframes in l_dfMissTrait and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors.
      l_evsTrait <- mapply(AppendEigenvectors, data = l_dfMissTrait, vars = vars, MoreArgs = list(tree = tree, predictors = l_predictors), SIMPLIFY = F)
      # Extract list of dfMissEV.  
      l_dfMissEV <- sapply(l_evsTrait, function(x) x[[1]])
      names(l_dfMissEV) <- traits
      # Extract updated list of predictors.
      l_EVPredictors <- sapply(l_evsTrait, function(x) x[[2]])
      names(l_EVPredictors) <- traits
    }
    # Outgroup removal. ---   
    # Identify outgroup(s) in trait datasets (this is the species that contains no trait data dfLog).   
    outgroup <- dfLog$species_name[apply(dfLog[, vars], 1, function(x) all(is.na(x)))]
    # If found in dataframe..  
    if(length(outgroup) > 0){
      # Remove outgroup from dataframes.
      dfLog <- dfLog[!dfLog$species_name %in% outgroup, ]
      l_dfMissTrait <- lapply(l_dfMissTrait, function(x) x[!x$species_name %in% outgroup, ])
      # If phylogenetic imputation was selected..   
      if(phyImp == T) {
        # Remove outgroup from dataframes in l_dfMissEV and dfMissEV.
        l_dfMissEV <- lapply(l_dfMissEV, function(x) x[!x$species_name %in% outgroup, ])
      }
    }
    
    # Predictor matrix creation. ---     
    # If phylogenetic imputation wasn't chosen..
    if(phyImp == F){
      # Create predictor matrix for use in MICE imputation (based on results of our trait predictor screening).   l_dfMissTrait[[1]] in MNAR
      matPredictors <- CreatePredictorMatrix(cols = vars, predictors = l_predictors, dfMissing = l_dfMissTrait[[1]])
      # If phylogenetic imputation was chosen..
    } else if(phyImp == T){
      # Create a predictor matrix for each dataframe in l_dfMissEV for use in MICE imputation with appended eigenvectors. We can indicate the column names of each dfMissEV (excluding species_name) as names of variables to consider in imputation process.  
      # Create list of column names for each dataframe in l_dfMissEV.  
      l_colsEV <- lapply(l_dfMissEV, function(x) {
        # Extract column names.
        colNames <- colnames(x)
        # Remove species_name.
        colNames <- colNames[!colNames %in% "species_name"]
        # Remove the missingness indicator columns.
        missIndex <- grep(pattern = "_NA", colNames)
        colNames <- colNames[-missIndex] 
      })
      # Modify list of l_EVPredictors so that only the trait in question has associated eigenvectors in the imputation.
      # Create list to hold updated predictors.
      l_l_EVPredictors <- lapply(1:length(vars), function(x) l_EVPredictors)
      names(l_l_EVPredictors) <- vars
      # For every trait..
      for(t in 1:length(vars)){
        # Take tth trait. 
        trait <- vars[[t]]
        # Take tth l_l_EVPredictors.
        l_EVTrait <- l_l_EVPredictors[[t]]
        # For every set of predictors..
        for(p in 1:length(l_EVTrait)){
          # If these are not the predictors for the tth trait..
          if(!names(l_EVTrait[p]) %in% trait){
            # Identify eigenvector predictors in l_EVPredictors.
            index <- grep(pattern = "V_", l_EVTrait[[p]])
            # Remove eigenvectors from predictor list as these are not in the corresponding dfMissEV. 
            updatedPredictors <- l_EVTrait[[p]][-index]
            # Replace l_EVTrait[[p]] with updated predictors.
            l_EVTrait[[p]] <- updatedPredictors
          }
        }
        # Append updated l_EVTrait to l_l_EVPredictors.
        l_l_EVPredictors[[t]] <- l_EVTrait
      }
      # Apply CreatePredictorMatrix using l_dfMissEV, l_colsEV, and l_l_EVPredictors to create a list of predictor matrices for use in MICE imputation.
      l_matPredictors <- mapply(CreatePredictorMatrix, dfMissing = l_dfMissEV, cols = l_colsEV, predictors = l_l_EVPredictors)
      names(l_matPredictors) <- vars
    }
    
    # Imputation. ---   
    # The ImputeMICE() function entails a loop that uses different values of m (the number of multiply imputed datasets) and returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.
    # If phylogenetic imputation was selected..  
    if(phyImp == T) {
      # Impute data using eigenvectors.
      l_impResult <- mapply(ImputeMICE, dfMissing = l_dfMissEV, matPredictors = l_matPredictors, cols = vars, MoreArgs = list(dfTrue = dfLog, cont = contTraits, cat = catTraits, inter = intTraits, mSets = c(5, 10, 40)), SIMPLIFY = F)
      # Since now we have a list of impResult, create new lists to hold the imputed dataframes and error rates.
      l_ImpTrait <- lapply(l_impResult, function(x) x[[1]])
      l_Error <- lapply(l_impResult, function(x) x[[2]])
      # Imputed dataset handling. ---  
      # Creation of l_dfImp so that there is one dataframe per parameter.  
      # Take species name and missingness information from dfMissing.  
      dfImp <- as.data.frame(dfLog[, c("species_name")])
      names(dfImp) <- "species_name"
      # Make a list to hold the dataframes for each value of mSets.     
      l_dfImp <- lapply(1:length(mSets), function(x) dfImp)
      names(l_dfImp) <- mSets
      # For every trait...
      for (t in 1:length(vars)) { 
        # Take the tth trait.
        trait <- vars[[t]]
        # Extract the tth element in l_ImpTrait
        l_dfImputedTrait <- l_ImpTrait[[t]]
        # Get the name of the corresponding missingness info column.
        missCol <- grep(pattern = paste(trait, "_NA", sep = ""), x = colnames(l_dfImputedTrait[[1]]))
        missCol <- names(l_dfImputedTrait[[1]])[missCol]
        # For each value of mSet..
        for(m in 1:length(l_dfImputedTrait)){
          # Extract the mth dataframe in l_dfImputedTrait.
          dfImputedTrait <- l_dfImputedTrait[[m]]
          # Subset to only contain species_name, trait info, and missingness indicator column.
          dfImputedTrait <- dfImputedTrait[, c("species_name", trait, missCol)]
          # Extract the mth dataframe in l_dfImputedParam.
          dfImp <- l_dfImp[[m]]
          # Merge with dfImputedTrait.
          dfImp <- merge(dfImp, dfImputedTrait, by = "species_name")
          # Replace mth element of l_dfImputedParam.
          l_dfImp[[m]] <- dfImp
        }
      }
      
    } else if(phyImp == F){
      # Impute data only using trait data.
      l_impResult <- mapply(ImputeMICE, dfMissing = l_dfMissTrait, cols = vars, MoreArgs = list(dfTrue = dfLog, cont = contTraits, cat = catTraits, inter = intTraits, mSets = c(5, 10, 40), matPredictors = matPredictors), SIMPLIFY = F)
      # Since now we have a list of impResult, create new lists to hold the imputed dataframes and error rates.  
      l_ImpTrait <- lapply(l_impResult, function(x) x[[1]])
      l_Error <- lapply(l_impResult, function(x) x[[2]])
      names(l_Error) <- vars
      # Imputed dataset handling. ---  
      # Creation of l_dfImp so that there is one dataframe per parameter.  
      # Take species name and missingness information from dfMissing.  
      dfImp <- as.data.frame(dfLog[, c("species_name")])
      names(dfImp) <- "species_name"
      # Make a list to hold the dataframes for each value of mSets.     
      l_dfImp <- lapply(1:length(mSets), function(x) dfImp)
      names(l_dfImp) <- mSets
      # For every trait...
      for (t in 1:length(vars)) { 
        # Take the tth trait.
        trait <- vars[[t]]
        # Extract the tth element in l_ImpTrait
        l_dfImputedTrait <- l_ImpTrait[[t]]
        # Get the name of the corresponding missingness info column.
        missCol <- grep(pattern = paste(trait, "_NA", sep = ""), x = colnames(l_dfImputedTrait[[1]]))
        missCol <- names(l_dfImputedTrait[[1]])[missCol]
        # For each value of ntree..
        for(n in 1:length(l_dfImputedTrait)){
          # Extract the kth dataframe in l_dfImputedTrait.
          dfImputedTrait <- l_dfImputedTrait[[n]]
          # Subset to only contain species_name, trait info, and missingness indicator column.
          dfImputedTrait <- dfImputedTrait[, c("species_name", trait, missCol)]
          # Extract the nth dataframe in l_dfImputedParam.
          dfImp <- l_dfImp[[n]]
          # Merge with dfImputedTrait.
          dfImp <- merge(dfImp, dfImputedTrait, by = "species_name")
          # Replace nth element of l_dfImputedParam.
          l_dfImp[[n]] <- dfImp
        }
      }
    }
    
    # Back-transforming data. ---  
    # First, match data species (original species in complete-case) to species now in dfLog in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfLog$species_name, ]
    # For every imputed dataframe..
    for(d in 1:length(l_dfImp)){
      # Take the dth dfImp.
      dfImp <- l_dfImp[[d]]
      # For every continuous trait..
      for(t in 1:length(contTraits)){
        # Take the tth continuous trait.
        trait <- contTraits[[t]]
        # Take the corresponding dfMissTrait.
        dfMissTrait <- l_dfMissTrait[[grep(trait, names(l_dfMissTrait))]]
        # Apply the Backtransform function for the trait in question.
        dfImpBTrait <- BackTransform(origData = dfOrig, tfData = dfImp, missData = dfMissTrait, cols = trait)
        # Replace data in dfImp with backtransformed data.
        dfImp[, trait] <- dfImpBTrait[, trait]
      }
      # If there are any integer traits..
      if(length(intTraits) > 0){
        # Round to nearest whole number.
        dfImp[, intTraits] <- lapply(dfImp[, intTraits], function(x) as.integer(round(x)))
      }
      # Replace dfImp in l_dfImp with newly backtransformed dataset.
      l_dfImp[[d]] <- dfImp
    }
    
    # Append results to ith element of lists. ---
    l_l_dfMissTraitOrig[[i]] <- l_dfMissTraitOrig
    l_l_missingness[[i]] <- l_missingness
    l_l_dfImp[[i]] <- l_dfImp 
    l_l_Error[[i]] <- l_Error
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      l_l_evs[[i]] <- l_evsTrait
    }
  } ## i
  
  # Average out the missingness for each trait. ---
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(vars, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  # Calculate the average number of NAs introduced for each trait.
  avgNAs <- mapply(function(x, y) x * y, y = l_sampleSizes, x = avgMiss)
  
  # Result object creation. ---
  # If phylogenetic imputation was selected..
  if(phyImp == T) { 
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, traits = vars, numeric = contTraits, categorical = catTraits, integer = intTraits, sampleSizes = l_sampleSizes, quantiles = quantiles, directions = directions, absolutes = absolutes, categories = categories, m = mSets, predictorMatrix = l_matPredictors, missingness = missingness, reps = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_evs, missingData = l_l_dfMissTraitOrig, eigenTraitData = l_dfMissEV, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  } else if(phyImp == F){
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, traits = traits, numeric = contTraits, categorical = catTraits, integer = intTraits, sampleSizes = l_sampleSizes, quantiles = quantiles, directions = directions, absolutes = absolutes, categories = categories, m = mSets, predictorMatrix = matPredictors, missingness = missingness, reps = int, predictors = l_predictors, missingData = l_l_dfMissTraitOrig, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error)
  }
  # Return l_results.
  return(l_results)
  
}
