# Functions for imputation scripts.

AppendMCAR <- function(data, cols, varsMCAR){
  
  # Function for identifying additional variables that could not be simulated MAR (variables with no missing values introduced upon application of the logistic regression models).
  # data = dataframe with missingness indicator columns
  # cols = columns to check for associated missingness indicator columns
  # varsMCAR = previously identified variables that can only be simulated MCAR
  
  # Get the names of the missingness indicator columns as we will needs these later on.
  missingCols <- colnames(data)[grep(pattern = "_NA", x = colnames(data))]
  
  # For each column..
  for(t in 1:length(cols)) {
    # Take the name of the tth column.
    variable <- cols[[t]]
    # Get the name of the corresponding missingness info column.
    missCol <- grep(pattern = variable, x = missingCols)
    missCol <- missingCols[missCol]
    # If no missing values were introduced (i.e. no associated missingness indicator with the variable)..
    if(length(missCol) == 0){
      # Append variable name to varsMCAR.
      varsMCAR <- c(variable, varsMCAR)
    }
  }
  # Return updated list of MCAR variables.
  return(varsMCAR)
}

AverageErrors <- function(results, data, vars, missLevel = 0.1, method, paramTrack = F, treeName = NULL) {
  
  # Function for averaging out error rates of the results from the imputation functions.
  
  # results = list object containing error rates for either numerical or categorical traits
  # df = original dataframe containing trait information
  # vars = vector of trait names
  # method = one of "MeanMode", "KNN", "RF", or "MICE"
  # paramTrack = Whether to track the parameter tuning (e.g. optimal K for KNN, optimal ntree for RF, optimal m for MICE)
  
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  
  # Convert missLevel into vector of missingness proportions.
  missingness <- seq(0.1, missLevel, by = 0.1)
  # Determine number of missingness replicates.
  reps <- length(results)
  # Make a list to hold the average error rates and SEs for the vars at each missingness level.
  l_AVGFinal <- CreateNamedList(listLength = length(vars), elementNames = vars)
  l_SEFinal <- CreateNamedList(listLength = length(vars), elementNames = vars)
  
  # If parameter tracking is on...
  if(paramTrack == T) {
    # Create empty list to hold the optimal parameter for the vars at each missingness level.
    l_ParamFinal <- CreateNamedList(listLength = length(vars), elementNames = vars)
  }
  
  # For each trait...
  for(t in 1:length(vars)) {
    # Take the tth trait.
    trait <- vars[[t]] 
    # Create empty vector to hold the error rates and SEs for each missingness level.
    finalRates <- vector(mode = "numeric", length = length(missingness))
    finalSE <- vector(mode = "numeric", length = length(missingness))
    # If parameter tracking is on...
    if(paramTrack == T) {
      # Create empty list to hold the optimal parameter for each missingness level.
      finalParam <- vector(mode = "numeric", length = length(missingness))
    }
    
    # For every level of missingness...
    for (s in 1:length(missingness)) {
      if(method == "MeanMode"|method == "KNN"|method == "RF"){
        # Create empty vector to hold the error rates for each replicate.
        repRates <- vector(mode = "list", length = reps)
        # For every replicate...
        for (i in 1:reps) { 
          # Putting the results into dataframe format.
          dfRes <- as.data.frame(results[[i]])
          # Subset by column for the trait of interest.
          traitCol <- grep(pattern = trait, x = colnames(dfRes))
          dfResTrait <- dfRes[traitCol]
          # Extract the column that corresponds to the missingness level.
          missingCol <- dfResTrait[, s]
          names(missingCol) <- rownames(dfResTrait) # ***
          # Append to repRates.
          repRates[[i]] <- missingCol
        }
        # Convert repRates to dataframe format.
        dfRep <- as.data.frame(repRates)
        # Name the columns.
        colnames(dfRep) <- 1:reps
        # Calculate the average error rate across the replicates.
        dfRep$average <- rowMeans(dfRep)
        # Calculate the standard error.
        dfRep$SE <- apply(dfRep[, 1:reps], MARGIN = 1, FUN = std.error)
        # If parameter tracking is on...
        if(paramTrack == T) {
          if(method == "KNN"){
            # Identify value of parameter that minimizes the imputation error rate.
            bestParam <- which.min(dfRep$average) 
          } else if(method == "RF"){
            # Identify value of parameter that minimizes the imputation error rate.
            bestParam <- which.min(dfRep$average) 
            # Get name of parameter value.
            bestParam <- as.numeric(rownames(dfRep[bestParam, ])) ## ***
          }
        }
        
        # If the method is MeanMode...
        if(method == "MeanMode") {
          # Subset out the average error rate (we only have value here for Mean/Mode).
          avgError <- as.numeric(dfRep["average"])
          # Subset out the SE.
          SE <- as.numeric(dfRep["SE"])
        } else {
          # Identify which parameter value results in the lowest error rate.
          bestContK <- which.min(dfRep$average)
          # Subset out the lowest average error rate.
          avgError <- dfRep[bestContK, "average"]
          # Subset out the corresponding SE.
          SE <- dfRep[bestContK, "SE"]
        }
      } else if(method == "MICE"){
        # Create a list to hold the error rates for each replicate.
        l_dfRep <- vector(mode = "list", length = reps)
        # For every replicate...
        for(i in 1:reps) {
          # Take the results for the ith replicate.
          l_results <- results[[i]]
          # Subset list for missingness level.
          result <- l_results[[s]]
          dfRep <- as.data.frame(result)
          # Subset for trait in question.
          traitRow <- grep(trait, rownames(dfRep))
          dfRep <- dfRep[traitRow, ]
          # Transpose dataframe.
          dfRep <- t(dfRep)
          # Append to l_dfRep.
          l_dfRep[[i]] <- dfRep
        }
        # Combine l_dfResTrait.
        dfRepAll <- as.data.frame(do.call("cbind", l_dfRep))
        # Name columns according to int.
        colnames(dfRepAll) <- 1:reps
        # Calculate the average MSE for each m across the iterations.
        dfRepAll$average <- rowMeans(dfRepAll)
        # Calculate the standard error.
        dfRepAll$SE <- apply(dfRepAll[, 1:reps], MARGIN = 1, FUN = std.error)
        # Which value of m minimizes the error for the trait?
        bestParam <- which.min(dfRepAll$average)
        # Subset out the lowest error.
        avgError <- dfRepAll[bestParam, "average"]
        SE <- dfRepAll[bestParam, "SE"]
        # Convert bestParam to name of parameter value instead of row name.
        bestParam <- rownames(dfRepAll[bestParam, ])
        # Remove X from name and convert to numeric class. ***
        bestParam <- as.numeric(gsub("X", "", bestParam))
      }
      # Append value to finalRates.
      finalRates[[s]] <- avgError
      # Append value to finalSE.
      finalSE[[s]] <- SE
      # If parameter tracking is on...
      if(paramTrack == T) {
        # Append value to finalParam.
        finalParam[[s]] <- bestParam
      }
    }
    # Append to the list that holds rates for the trait at each missingness level.
    l_AVGFinal[[t]] <- finalRates
    l_SEFinal[[t]] <- finalSE
    # If parameter tracking is on...
    if(paramTrack == T) {
      # Append value to finalParam.
      l_ParamFinal[[t]] <- finalParam
    }
  }
  
  # Create a vector of missingness values.
  missingness_level <- missingness
  # Create a dataframe that contains information about missingness levels.
  dfError <- as.data.frame(missingness_level)
  # Make a list to hold the error dataframes.
  l_dfError <- lapply(1:length(vars), function(x) dfError)
  names(l_dfError) <- vars
  # If parameter tracking is on...
  if(paramTrack == T) {
    # Create a dataframe that contains information about levels of K at different missingness levels.
    dfParam <- as.data.frame(missingness_level)
    # Make a list to hold the K dataframes.
    l_dfParam <- lapply(1:length(vars), function(x) dfParam)
    names(l_dfParam) <- vars
  }
  
  # For each trait...
  for(t in 1:length(vars)) {
    # Take the tth trait.
    trait <- vars[[t]]
    # Unlist the error rates.
    rates <- unlist(l_AVGFinal[[t]])
    # Unlist the standard errors.
    se <- unlist(l_SEFinal[[t]])
    # Take the tth dfError.
    dfTrait <- l_dfError[[t]]
    # If trait is numerical...
    if(trait %in% contTraits) {
      # Append the average MSE and standard errors to the corresponding dataframe.
      dfTrait[, paste(method, "MSE", sep = "_")] <- rates
      dfTrait[, paste(method, "SE", sep = "_")] <- se
      # If trait is categorical...
    } else if(trait %in% catTraits){
      # Append the average PFC and standard errors to the corresponding dataframe.
      dfTrait[, paste(method, "PFC", sep = "_")] <- rates
      dfTrait[, paste(method, "SE", sep = "_")] <- se
    }
    # Append to l_dfError.
    l_dfError[[t]] <- dfTrait
    if(is.null(treeName)){
      # Write dataframe to file according to name of taxa and imputation method.
      WriteErrorRates(dfTrait, fileName = paste(trait, method, "MCAR_ErrorRates.csv", sep = "_"))
    } else if(!is.null(treeName)){
      # Write dataframe to file according to name of taxa and imputation method.
      WriteErrorRates(dfTrait, fileName = paste(trait, method, treeName, "MCAR_ErrorRates.csv", sep = "_"))
    }
    # If parameter tracking is on...
    if(paramTrack == T) {
      # Unlist the parameters.
      param <- unlist(l_ParamFinal[[t]])
      # Take the tth dataframe.
      dfParam <- l_dfParam[[t]]
      # Append the param values to the corresponding dataframe.
      dfParam$param <- param
      # Append to l_dfParam.
      l_dfParam[[t]] <- dfParam
      if(is.null(treeName)){
        # Write dataframe to file according to name of taxa and imputation method.
        WriteErrorRates(dfParam, fileName = paste(trait, method, "MCAR_Parameters.csv", sep = "_"))
      } else if (!is.null(treeName)){
        # Write dataframe to file according to name of taxa and imputation method.
        WriteErrorRates(dfParam, fileName = paste(trait, method, treeName, "MCAR_Parameters.csv", sep = "_"))
      }
    }
  }
}

AverageErrorsBias <- function(results, data, cont, cat, method, type = "MAR", paramTrack = F, treeName = NULL) {
  
  # Function for averaging out error rates of the results from the SimputeMAR and SimputeMNAR functions.
  
  # results = list object containing error rates for either numerical or categorical traits
  # data = original dataframe containing trait information
  # cont = names of continuous variables
  # cat = names of categorical variables
  # method = imputation method. One of "MeanMode", "KNN", "RF", or "MICE"
  # type = type of missingness pattern. One of "MAR" or "MNAR"
  # paramTrack = Whether to track the parameter tuning (e.g. optimal K for KNN, optimal ntree for RF, optimal m for MICE)
  
  # If the method is MeanMode, KNN, or RF...
  if(method == "MICE" & type == "MAR"){
    # Extract variable names from the first element of results[[1]][[1]]. This is because MICE result is structured by parameter value and not trait.
    vars <- names(results[[1]][[1]])
    # Otherwise..
  } else {
    # Extract variable names using the first element of results[[1]].
    vars <- names(results[[1]])
  }
  
  # Determine number of missingness replicates.
  reps <- length(results)
  # Make a list to hold the average error rates and SEs for the vars.
  AVGFinal <- CreateNamedList(listLength = length(vars), elementNames = vars)
  SEFinal <- CreateNamedList(listLength = length(vars), elementNames = vars)
  
  # If parameter tracking is on...
  if(paramTrack == T) {
    # Create empty list to hold the optimal parameter for the vars at each missingness level.
    finalParam <- CreateNamedList(listLength = length(vars), elementNames = vars)
  }
  
  # For each trait...
  for(t in 1:length(vars)) {
    # Take the tth trait.
    trait <- vars[[t]] 
    # If the method is MeanMode, KNN, or RF...
    if(method == "MeanMode"|method == "KNN"|method == "RF"){
      # Create empty vector to hold the error rates for each replicate.
      repRates <- vector(mode = "list", length = reps)
      # For every replicate...
      for (i in 1:reps) { 
        # Putting the results into dataframe format.
        dfRes <- as.data.frame(results[[i]])
        # Subset by column for the trait of interest.
        traitCol <- grep(pattern = trait, x = colnames(dfRes))
        missingCol <- dfRes[, traitCol]
        names(missingCol) <- rownames(dfRes) ### ***
        # Append to repRates.
        repRates[[i]] <- missingCol
      }
      
      # Convert repRates to dataframe format.
      dfRep <- as.data.frame(repRates)
      # Name the columns.
      colnames(dfRep) <- 1:reps
      # Calculate the average error rate across the replicates.
      dfRep$average <- rowMeans(dfRep)
      # Calculate the standard error.
      dfRep$SE <- apply(dfRep[, 1:reps], MARGIN = 1, FUN = std.error)
      # If parameter tracking is on...
      if(paramTrack == T) {
        if(method == "KNN"){
          # Identify value of parameter that minimizes the imputation error rate.
          bestParam <- which.min(dfRep$average) 
        } else if(method == "RF"){
          # Identify value of parameter that minimizes the imputation error rate.
          bestParam <- which.min(dfRep$average) 
          # Get name of parameter value.
          bestParam <- as.numeric(rownames(dfRep[bestParam, ])) ## ***
        }
      }
      # If the method is MeanMode...
      if(method == "MeanMode") {
        # Subset out the average error rate (we only have one value here for Mean/Mode).
        avgError <- as.numeric(dfRep["average"])
        # Subset out the SE.
        SE <- as.numeric(dfRep["SE"])
      } else {
        # Identify which parameter value results in the lowest error rate.
        bestP <- which.min(dfRep$average)
        # Subset out the lowest average error rate.
        avgError <- dfRep[bestP, "average"]
        # Subset out the corresponding SE.
        SE <- dfRep[bestP, "SE"]
      }
    } else if(method == "MICE" & type == "MAR"){
      # Create a list to hold the error rates for each replicate.
      l_dfRep <- vector(mode = "list", length = reps)
      
      # For every replicate...
      for(i in 1:reps) {
        # Take the results for the ith replicate.
        result <- results[[i]]
        # Convert to dataframe format.
        dfRep <- as.data.frame(result)
        # Subset for trait in question.
        traitRow <- grep(trait, rownames(dfRep))
        dfRep <- dfRep[traitRow, ]
        # Transpose dataframe.
        dfRep <- t(dfRep)
        # Append to l_dfRep.
        l_dfRep[[i]] <- dfRep
      }
      
      # Combine l_dfRep.
      dfRepAll <- as.data.frame(do.call("cbind", l_dfRep))
      # Name columns according to int.
      colnames(dfRepAll) <- 1:reps
      # Calculate the average MSE for each m across the iterations.
      dfRepAll$average <- rowMeans(dfRepAll)
      # Calculate the standard error.
      dfRepAll$SE <- apply(dfRepAll[, 1:reps], MARGIN = 1, FUN = std.error)
      # Which value of m minimizes the error for the trait?
      bestParam <- which.min(dfRepAll$average)
      # Subset out the lowest error.
      avgError <- dfRepAll[bestParam, "average"]
      SE <- dfRepAll[bestParam, "SE"]
      # Convert bestParam to name of parameter value instead of row name.
      bestParam <- rownames(dfRepAll[bestParam, ])
      # Remove X from name and convert to numeric class. ***
      bestParam <- as.numeric(gsub("X", "", bestParam))
      
    } else if(method == "MICE" & type == "MNAR"){
      # Create a list to hold the error rates for each replicate.
      l_dfRep <- vector(mode = "list", length = reps)
      
      # For every replicate...
      for(i in 1:reps) {
        # Take the results for the ith replicate.
        result <- results[[i]]
        # Convert to dataframe format.
        dfRep <- as.data.frame(result)
        # Transpose dataframe.
        dfRep <- t(dfRep)
        # Subset for trait in question.
        traitRow <- grep(trait, rownames(dfRep))
        dfRep <- dfRep[traitRow, ]
        # Append to l_dfRep.
        l_dfRep[[i]] <- dfRep
      }
      
      # Combine l_dfRep.
      dfRepAll <- as.data.frame(do.call("cbind", l_dfRep))
      # Name columns according to int.
      colnames(dfRepAll) <- 1:reps
      # Calculate the average error for each m across the iterations.
      dfRepAll$average <- rowMeans(dfRepAll)
      # Calculate the standard error.
      dfRepAll$SE <- apply(dfRepAll[, 1:reps], MARGIN = 1, FUN = std.error)
      # Which value of m minimizes the error for the trait?
      bestParam <- which.min(dfRepAll$average)
      # Subset out the lowest error.
      avgError <- dfRepAll[bestParam, "average"]
      SE <- dfRepAll[bestParam, "SE"]
      # Convert bestParam to name of parameter value instead of row name.
      bestParam <- rownames(dfRepAll[bestParam, ])
      # Remove trait from name and convert to numeric class. ***
      bestParam <- as.numeric(gsub(paste(trait, ".", sep = ""), "", bestParam))
    }
    
    # Append to the list that holds rates for the trait.
    AVGFinal[[t]] <- avgError
    SEFinal[[t]] <- SE
    # If parameter tracking is on...
    if(paramTrack == T) {
      # Append value to finalParam.
      finalParam[[t]] <- bestParam
    }
  }
  
  # Create a dataframe that contains information about missingness levels.
  dfError <- as.data.frame(type)
  colnames(dfError) <- "missingness_level"
  # Make a list to hold the error dataframes.
  l_dfError <- lapply(1:length(vars), function(x) dfError)
  names(l_dfError) <- vars
  # If parameter tracking is on...
  if(paramTrack == T) {
    # Create a dataframe that contains information about levels of parameter at different missingness levels.
    dfParam <- as.data.frame(type)
    colnames(dfParam) <- "missingness_level"
    # Make a list to hold the K dataframes.
    l_dfParam <- lapply(1:length(vars), function(x) dfParam)
    names(l_dfParam) <- vars
  }
  
  # For each trait...
  for(t in 1:length(vars)) {
    # Take the tth trait.
    trait <- vars[[t]]
    # Extract the error rate.
    rate <- AVGFinal[[t]]
    # Extract the SE.
    se <- SEFinal[[t]]
    # Take the tth dfError.
    dfTrait <- l_dfError[[t]]
    # If trait is numerical...
    if(trait %in% cont) {
      # Append the average MSE and standard errors to the corresponding dataframe.
      dfTrait[, paste(method, "MSE", sep = "_")] <- rate
      dfTrait[, paste(method, "SE", sep = "_")] <- se
      # If trait is categorical...
    } else if(trait %in% cat){
      # Append the average PFC and standard errors to the corresponding dataframe.
      dfTrait[, paste(method, "PFC", sep = "_")] <- rate
      dfTrait[, paste(method, "SE", sep = "_")] <- se
    }
    # Append to l_dfError.
    l_dfError[[t]] <- dfTrait
    if(is.null(treeName)){
      # Write dataframe to file according to name of taxa, imputation method, and missingness pattern.
      WriteErrorRates(dfTrait, fileName = paste(trait, method, type, "ErrorRates.csv", sep = "_"))
    } else if(!is.null(treeName)){
      # Write dataframe to file according to name of taxa and imputation method.
      WriteErrorRates(dfTrait, fileName = paste(trait, method, treeName, type, "ErrorRates.csv", sep = "_"))
    }
    # If parameter tracking is on...
    if(paramTrack == T) {
      # Unlist the parameters.
      param <- unlist(finalParam[[t]])
      # Take the tth dataframe.
      dfParam <- l_dfParam[[t]]
      # Append the param values to the corresponding dataframe.
      dfParam$param <- param
      # Append to l_dfParam.
      l_dfParam[[t]] <- dfParam
      if(is.null(treeName)){
        # Write dataframe to file according to name of taxa, imputation method, and missingness pattern.
        WriteErrorRates(dfParam, fileName = paste(trait, method, type, "Parameters.csv", sep = "_"))
      } else if (!is.null(treeName)){
        # Write dataframe to file according to name of taxa, imputation method, and missingness pattern.
        WriteErrorRates(dfParam, fileName = paste(trait, method, treeName, type, "Parameters.csv", sep = "_"))
      }
    }
  }
}

BackTransform <- function(origData, tfData, missData, cols) {
  
  # Function for applying back-transformations to a dataframe's continuous variables. This function takes into account whether a variable contains negative values originally.
  
  # origData = dataframe containing untransformed data
  # tfData = dataframe containing transformed continuous data
  # missData = dataframe containing missingness information
  # cols = column names in tfData to backtransform
  
  # Ensure they are all in dataframe format.
  origData <- as.data.frame(origData)
  tfData <- as.data.frame(tfData)
  missData <- as.data.frame(missData)
  
  # For continuous traits, check which ones contain negative values in original data. 
  negTest <- lapply(origData[, cols, drop = F], function(x) sum(na.omit(x) < 0) > 0)
  negCols <- cols[which(negTest == T)]
  # Apply exp function to all other columns.
  expCols <- setdiff(cols, negCols)
  # If there's more than one expCols
  if(length(expCols) > 1){
    # Lapply exp function.
    tfData[, expCols] <- lapply(tfData[, expCols], exp)
  } else {
    # Apply exp function to one column.
    tfData[, expCols] <- exp(tfData[, expCols])
  }
  # If a column with negative data was identified..
  if(length(negCols) > 0){
    # For every column in negCols..
    for(i in 1:length(negCols)){
      # Take the ith negCol.
      negCol <- negCols[[i]]
      # Apply ExpNegative function using imputed and original data.
      btCol <- ExpNegative(variable = tfData[, negCol], original = origData[, negCol], missing = missData[, negCol])
      # Replace column in tfData.
      tfData[, negCol] <- btCol
    }
  }
  # Return dataframe with backtransformed values.
  return(tfData)
}

BindAndOrganize <- function(data, vars){
  
  # Function that binds missingness indicator columns (bind_shadow function from the "naniar" package) and reorganizes the columns in the dataframe.
  # data = dataframe with missing values and column called species_name
  # vars = names of columns with trait variables
  
  # We need to keep track of which observations were imputed by adding a shadow matrix (bind_shadow() function).
  dfMiss <- bind_shadow(data, only_miss = T)
  # Get the names of the missing columns as we will needs these later on.
  missingCols <- colnames(dfMiss)[grep(pattern = "_NA", x = colnames(dfMiss))]
  # Reorganize the dataframe.
  dfMiss <- dfMiss[, c("species_name", vars, missingCols)]
  # Return reorganized dataframe.
  return(dfMiss)
}

CombineMIDataframes <- function(midata, method, m, contVars, catVars) {
  
  # Function for combining imputed dataframes from an object of class mids.
  
  # midata = object containing multiply imputed dataframes
  # method = "MICE", "Amelia", or "mi"
  # m = number of multiply imputed datasets
  # contVars = column names for continuous variables
  # catVars = column names for categorical variables
  
  # If MICE was used for imputation...
  if(method == "MICE") {
    # Create a list to hold the complete datasets.
    l_completeDatasets <- vector(mode = "list", length = m)
    # For every dataset...
    for(i in 1:m) {
      # Get the complete dataset.
      dfComplete <- complete(midata, i)
      # Add it to l_completeDatasets.
      l_completeDatasets[[i]] <- dfComplete
    }
    # Bind the datasets together by row.
    dfCombined <- rbindlist(l_completeDatasets)
    # Set to data.table.
    setDT(dfCombined)
    # Create a new dataframe to hold our averaged info.
    dfAveraged <- as.data.frame(unique(dfCombined$species_name))
    colnames(dfAveraged)[1] <- "species_name"
    # For continuous variables, get the mean value for each ID (should be m number of IDs).
    dfAvgCont <- dfCombined[, lapply(.SD, mean), by = species_name, .SDcols = contVars]
    # For categorical variables, get the mode value for each ID (should be m number of IDs).
    dfAvgCat <- dfCombined[, lapply(.SD, function(x) names(which.max(table(x)))), by = species_name, .SDcols = catVars]
    # Merge with dfAveraged.
    dfAveraged <- merge(dfAveraged, dfAvgCont, by = "species_name")
    dfAveraged <- merge(dfAveraged, dfAvgCat, by = "species_name")
    # Ensure it is returned as dataframe format.
    dfAveraged <- as.data.frame(dfAveraged)
    return(dfAveraged)
  } else if (method == "Amelia") {
    l_completeDatasets <- list()
    for(i in 1:m) {
      dfComplete <- midata$imputations[[i]]
      l_completeDatasets[[i]] <- dfComplete
    }
    dfCombined <- rbindlist(l_completeDatasets)
    dfCombined <- dfCombined[, lapply(.SD, mean), by = species_name]
    dfCombined <- as.data.frame(dfCombined)
  }
  
}

CreatePredictorMatrix <- function(dfMissing, cols, predictors){
  
  # Function for creating a predictor matrix for use in MICE imputation.
  # dfMissing = dataframe with missing values
  # cols = names of variables to be imputed
  # predictors = named list of predictors to be used to impute each variable
  
  # We need to create a predictor matrix for the traits. To do this, we must run mice (setting maxit to 0) to gain access to the predictor matrix.
  imputedMICE <- mice(dfMissing[, c("species_name", cols)], maxit = 0, print = F)
  # Access the predictor matrix.
  predMatrix <- imputedMICE$pred
  # Set all elements to 0.
  predMatrix[,] <- 0
  # Now we need to modify the predictor matrix so that only the traits previously identified as predictors will be used in the imputation process.
  # First, identify trait names from the list of predictors.
  traits <- names(predictors)
  # For every trait...
  for(t in 1:length(traits)){
    # Get the trait.
    trait <- traits[[t]]
    # Get the predictors for the trait.
    predInd <- grep(pattern = trait, x = names(predictors))
    preds <- predictors[[predInd]]
    # The predictors to 1 in the predictor matrix so they are used to impute the target trait.
    predMatrix[trait, preds] <- 1
  }
  # Return predictor matrix.
  return(predMatrix)
  
}

ExpNegative <- function(variable, original, missing) {
  
  # Function for back-transforming data with negative values.
  # variable = numeric vector that were transformed using the NormalizeNegative function
  # original = numeric vector of original data
  # missing = missingness indicator vector
  
  # Apply exp function to variable.
  untransformed <- exp(variable)
  # Add constant back that was originally used to log-transform negative values.
  finaltransformed <- untransformed - abs(min(original, na.rm = T))
  # Identify position of minimum value in original data.
  index <- which.min(original)
  # If the minimum value wasn't imputed..
  if(!is.na(missing[index])){
    # Subtract 1 from minimum value so it matches original. # TODO: Can remove this line once we add one in NormalizeNegative function.
    finaltransformed[index] <- finaltransformed[index] - 1 
  }
  # Return finaltransformed variable.
  return(finaltransformed)
  
}

FitLogReg <- function(data, cols) {
  
  # Function for performing logistic regression analyses on a dataframe to identify significant predictors of missingness for a particular column.
  # data = dataframe 
  # cols = columns with missing values (NAs) to consider for regression analyses
  
  # Convert data to dataframe format.
  data <- as.data.frame(data)
  # Replace blanks with NAs.
  data[data == ""] <- NA
  # For each trait with missing values (only_miss = T), create an indicator column for missingness using the bind_shadow() function.
  dfShadow <- bind_shadow(data, only_miss = T)
  # Get the names of the missingness indicator columns.
  missCols <- colnames(dfShadow)[grep(pattern = "_NA", x = colnames(dfShadow))]
  # For missCols, convert !NA to 1 (present) and NA to 0 (absent).
  dfShadow[, missCols] <- lapply(dfShadow[, missCols], function(x) ifelse(x == "!NA", 1, 0))
  # Create a list to hold the logistic regression models for each missingness indicator column in the dataframe.
  l_models <- vector(mode = "list", length = length(missCols))
  # Name l_models according to the names of the missingness indicator columns.
  names(l_models) <- missCols
  
  # For every missingness indicator column...
  for(i in 1:length(missCols)) {
    # Take the name of the ith missingness indicator column.
    response <- missCols[i]
    # Since we are simulating MAR and don't need the original response data as a covariate, remove the corresponding trait from the covariates.
    rmTrait <- gsub("_NA", "", response) ## Remove _NA from the response name to get name of trait to remove.
    candidates <- cols[!cols %in% rmTrait] ## Remove from cols vector.
    # Create empty vector to hold predictors for the response variable in question.
    predictors <- vector(mode = "character")
    
    # For every candidate predictor...
    for (j in 1:length(candidates)) {
      # Take the jth candidate.
      candidate <- candidates[j]
      # Create a subset of data that contains only response and candidate columns.
      dfSubset <- dfShadow[, c(response, candidate)]
      # Convert to dataframe format.
      dfSubset <- as.data.frame(dfSubset)
      # Removing missing values in candidate column.
      dfSubset <- na.omit(dfSubset)
      # If the candidate is numeric or there are 2 categories in the candidate (binary)...
      if(is.numeric(dfSubset[, candidate]) | length(unique(dfSubset[, candidate])) == 2) {
        # Fit glm model with binomial family specified. 
        fit <- glm(dfSubset[, response] ~ dfSubset[, candidate], family = "binomial", maxit = 100)
        # Extract the p-values from the model.
        pVal <- coef(summary(fit))[, 4]
        # Remove the intercept.
        pVal <- pVal[-1]
        # If there are more than 2 categories in the candidate (multi-categorical)..
      } else if(length(unique(dfSubset[, candidate])) > 2) {
        # Fit the alternative model.
        fitALT <- glm(dfSubset[, response] ~ dfSubset[, candidate], family = "binomial", maxit = 100)
        # Run an anova to compare the models.
        result <- anova(fitALT, test = "LRT")
        # Extract the p-value from the model.
        pVal <- result$`Pr(>Chi)`[2]
      }
      if(pVal < 0.05) {
        # Add the name of the significant predictor to the list of predictors.
        predictors[[j]] <- candidate
      }
    }
    # Remove NA values from predictors.
    predictors <- na.omit(predictors)
    # Convert dfShadow to dataframe.
    dfShadow <- as.data.frame(dfShadow)
    # If there are predictors of missingness for the trait..
    if(length(predictors) > 0){
      # Create regression formula.
      regForm <- as.formula(paste(response, "~", paste(predictors, collapse = "+")))
      # Fit a logistic regression model.
      fit <- glm(regForm, data = dfShadow, family = "binomial", na.action = na.omit)
    } else if(length(predictors) == 0){
      fit <- "No significant predictors of missingness!"
    }
    # Append to the list of predictors.
    l_models[[i]] <- fit
  }
  # Return list of predictors for each response variable.
  return(l_models)
}

IDMissPattern <- function(data, vars, models){
  
  # Function for identifying variables with no missing values in the raw data (these variables cannot be simulated missing at random (MAR) and so must be simulated missing completely at random (MCAR)).
  # data = dataframe with missing values
  # vars = column names to consider
  
  # First, identify variables with no missing values in the raw data (sum of NAs = 0).
  index <- sapply(data[, vars], function(x) sum(is.na(x)) == 0)
  # Subset out these variables.
  varsMCAR <- vars[index]
  
  # Now, out of the variables for which logistic regression models were fitted, identify the variables with no significant predictors of missingness (character vector result after application of the FitLogReg() function).
  index <- sapply(models, function(x) is.character(x) == T)
  # Separate from the variables that can be simulated MAR.
  varsMCAR_2 <- names(models)[index]
  # Remove "_NA" from names of varsMCAR_2 so we can get the original variable names.
  varsMCAR_2 <- gsub("_NA", "", varsMCAR_2)
  # Append to varsMCAR variable.
  varsMCAR <- c(varsMCAR, varsMCAR_2)
  # Subset out MAR models from models list.
  l_MARmodels <- models[!index]
  # Build the final models (e.g. by dropping insignificant terms in the model using the RefineModel() function).
  l_MARfinalModels <- lapply(l_MARmodels, RefineModel, data = data)
  # Remove "_NA" from names of l_MARfinalmodels so we can get the original variable names.
  names(l_MARfinalModels) <- gsub("_NA", "", names(l_MARfinalModels))
  # Separate from the variables that can be simulated MCAR.
  varsMAR <- names(l_MARfinalModels)
  # Return list of variables that can be simulated MAR vs. those that cannot (MCAR). Also return refined MAR models.
  return(list(MCAR = varsMCAR, MAR = varsMAR, refinedModels = l_MARfinalModels))
  
}

ImputeKNN <- function(dfTrue, dfMissing, cols, cont, cat, k, predictors, phyImp = F, l_dfMissing = NULL){
  
  # Function for imputing missing values using the kNN() function in the "VIM" package for imputation. Returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.
  # Citations: Kowarik A, Templ M (2016). “Imputation with the R Package VIM.” Journal of Statistical Software, 74(7), 1–16. doi: 10.18637/jss.v074.i07.
  # https://cran.r-project.org/web/packages/VIM/VIM.pdf
  
  # dfTrue = dataframe containing true (known) observations
  # dfMissing = dataframe containing missing values
  # cols = names of traits (column names)
  # cont = names of numerical traits
  # cat = names of categorical traits
  # k = maximum value of k (the number of nearest neighbours to use in the kNN algorithm) to consider
  # predictors = list of predictors for each trait
  # phyImp = whether data are to be imputed using phylogenetic information
  # l_dfMissing = if phyImp == T, a named list of dataframes with missing values corresponding to each trait
  
  # Get the names of the missing columns as we will need these later on.
  missingCols <- colnames(dfMissing)[grep(pattern = "_NA", x = colnames(dfMissing))]
  # Create a list of the values of k we are considering.
  neighbours <- seq(1, k, by = 1)
  
  # Create an empty list to hold the error rate results for each trait and value of k. ^
  l_ErrorParams <- CreateNamedList(listLength = length(cols), elementNames = cols)
  # Create an empty list to hold the imputed data for each trait and value of k. ^
  l_l_dfImputedTrait <- CreateNamedList(listLength = length(cols), elementNames = cols)
  
  # For each trait..
  for(t in 1:length(cols)) {
    # Take the name of the tth trait.
    trait <- cols[[t]]
    # Get the name of the corresponding missingness info column.
    missCol <- grep(pattern = trait, x = missingCols)
    missCol <- missingCols[missCol]
    # Take the corresponding predictors for that trait.
    predIndex <- grep(pattern = trait, x = names(predictors))
    preds <- predictors[[predIndex]]
    # Identify species with missing values in original dataset. &&
    origMiss <- dfTrue$species_name[is.na(dfTrue[, trait] == T)]
    if(phyImp == T){
      # Get the name of the corresponding dfMissing.
      index <- grep(pattern = trait, x = names(l_dfMissing))
      # Replace original dfMissing with dfMissing with appended eigenvectors for the trait.
      dfMissing <- l_dfMissing[[index]]
    }
    # Create a vector to hold the error rate results for each value of k.
    paramRes <- vector(mode = "numeric", length = length(neighbours))
    # Create a list to hold the imputed trait data for each value of k. ^
    l_dfImputedTrait <- CreateNamedList(listLength = length(neighbours), elementNames = neighbours)
    
    # Impute the trait using different values of k.
    for(k in 1:length(neighbours)) {
      # Take the number of neighbours.
      NN <- neighbours[[k]]
      # Impute the fold using NN. ^
      dfImputedTrait <- kNN(dfMissing, variable = trait, dist_var = preds, k = NN)
      # If the trait is numerical...
      if(trait %in% cont) {
        # Subset to only contain species_name, trait of interest, and information about whether the trait has been imputed. ^
        dfSubset <- dfImputedTrait[, c("species_name", trait, missCol)]
        # Subset to only contain rows that were imputed.
        dfSubset <- dfSubset[which(dfSubset[, missCol] == "NA"), ]
        # Rename column.
        colnames(dfSubset)[2] <- "Imputed_Value"
        # Merge with the true data.
        dfError <- merge(dfSubset, dfTrue, by = "species_name")
        # If origMiss contains at least one species..
        if(length(origMiss) > 0){ 
          # Remove from dfError. &
          dfError <- dfError[!dfError$species_name %in% origMiss, ]
        }
        # Calculate the mean squared error (MSE).
        errorRate <- mse(dfError[, trait], dfError[, "Imputed_Value"])
        # If the trait is categorical...
      } else if(trait %in% cat){
        # If origMiss contains at least one species.. &&
        if(length(origMiss) > 0){ 
          # Remove origMiss from dfTrue, dfImputed, and dfMissing. &&
          dfTruePFC <- dfTrue[!dfTrue$species_name %in% origMiss, ]
          dfImputedPFC <- dfImputedTrait[!dfImputedTrait$species_name %in% origMiss, ] # ** name is different
          dfMissingPFC <- dfMissing[!dfMissing$species_name %in% origMiss, ]
          # Calculate the PFC (percent falsely classified) with modified dataframes.
          errorRate <- pfc(x = dfTruePFC[[trait]], y = dfImputedPFC[[trait]], m = is.na(dfMissingPFC[, trait]))
          # If origMiss is empty..
        } else if(length(origMiss) == 0){
          # Calculate the PFC (percent falsely classified) using original dataframes.
          errorRate <- pfc(x = dfTrue[[trait]], y = dfImputedTrait[[trait]], m = is.na(dfMissing[, trait]))
        }
      }
      # Name the errorRate according to the value of k.
      names(errorRate) <- NN
      # Append results to the kth element in the lists.
      paramRes[[k]] <- errorRate
      l_dfImputedTrait[[k]] <- dfImputedTrait
    }
    # Append to the tth elements in the lists.
    l_ErrorParams[[t]] <- paramRes
    l_l_dfImputedTrait[[t]] <- l_dfImputedTrait
  }
  
  
  # Imputed dataframe handling. --- ^
  # Because we imputed data on a trait-by-trait basis, let's merge dataframes together for each parameter value so results are a bit easier to view. ^
  # Take species name and missingness information from dfMissing.
  dfImputedParam <- dfMissing[, c("species_name", missingCols)]
  # Make a list to hold the dataframes for each value of k.
  l_dfImputedParam <- lapply(1:length(neighbours), function(x) dfImputedParam)
  names(l_dfImputedParam) <- neighbours
  
  # For every trait...
  for (t in 1:length(l_l_dfImputedTrait)) { 
    # Take the tth trait.
    trait <- names(l_l_dfImputedTrait)[[t]]
    # Extract the tth element.
    l_dfImputedTrait <- l_l_dfImputedTrait[[t]]
    # For each value of k..
    for(k in 1:length(l_dfImputedTrait)){
      # Extract the kth dataframe in l_dfImputedTrait.
      dfImputedTrait <- l_dfImputedTrait[[k]]
      # Subset to only contain species_name and trait info.
      dfImputedTrait <- dfImputedTrait[, c("species_name", trait)]
      # Extract the kth dataframe in l_dfImputedParam.
      dfImputedParam <- l_dfImputedParam[[k]]
      # Merge with dfImputedTrait.
      dfImputedParam <- merge(dfImputedParam, dfImputedTrait, by = "species_name")
      # Replace kth element of l_dfImputedParam.
      l_dfImputedParam[[k]] <- dfImputedParam
    }
  }
  # Rearrange columns in l_dfImputedParam.
  l_dfImputedParam <- lapply(l_dfImputedParam, function(x) x[, c("species_name", cols, missingCols)])
  # Return list of error rates for each trait and parameter value and list of imputed dataframes for each parameter value. ^
  return(list(l_dfImputedParam = l_dfImputedParam, l_ErrorParams = l_ErrorParams))
  
}

ImputeMeanMode <- function(dfTrue, dfMissing, cols, cont, cat, inter = NULL){
  
  # Function for imputing missing values using the mean/mode replacement. Returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.
  
  # dfTrue = dataframe containing true (known) observations
  # dfMissing = dataframe containing missing values and missingness indicator columns
  # cols = names of traits (column names) to impute
  # cont = names of numerical traits
  # cat = names of categorical traits
  # inter = names of integer traits
  
  # Get the names of the missing columns as we will need these later on.
  missingCols <- colnames(dfMissing)[grep(pattern = "_NA", x = colnames(dfMissing))]
  
  # Imputation. ---
  # Make a copy of dfMissing called dfImputed which we will use to impute values.
  dfImputed <- dfMissing
  
  # For each trait..
  for(t in 1:length(cols)) {
    # Take the name of the tth trait.
    trait <- cols[[t]]
    # Identify missing values in trait.
    index <- which(is.na(dfImputed[trait]))
    # If trait is numeric..
    if(trait %in% cont){
      # Replace all NA values with the mean of the known observations for the continuous variables.
      dfImputed[index, trait] <- mean(dfImputed[[trait]], na.rm = T)
    } else if(trait %in% cat){
      # Replace all NA values with the mode of the known observations for the categorical variables.
      dfImputed[index, trait] <- names(which.max(table(dfImputed[[trait]])))
    }
  }
  
  # Error rate calculation. ---
  
  # Create an empty list to hold the error rate results for each trait.
  l_Error <- CreateNamedList(listLength = length(cols), elementNames = cols)
  
  # For each trait..
  for(t in 1:length(cols)) {
    # Take the name of the tth trait.
    trait <- cols[[t]]
    # Identify species with missing values in original dataset. &&
    origMiss <- dfTrue$species_name[is.na(dfTrue[, trait] == T)]
    # If the trait is numerical...
    if(trait %in% cont){
      # Get the name of the corresponding missingness info column.
      missCol <- grep(pattern = trait, x = missingCols)
      missCol <- missingCols[missCol]
      # Subset to only contain species_name, trait of interest, and info about whether the trait has been imputed.
      dfSubset <- dfImputed[, c("species_name", trait, missCol)]
      # Subset to only contain rows that were imputed.
      dfSubset <- dfSubset[which(dfSubset[, missCol] == "NA"), ]
      # Rename column.
      colnames(dfSubset)[2] <- "Imputed_Value"
      # Merge with the known data.
      dfError <- merge(dfSubset, dfTrue, by = "species_name")
      # If origMiss contains at least one species.. &&
      if(length(origMiss) > 0){ 
        # Remove from dfError. &
        dfError <- dfError[!dfError$species_name %in% origMiss, ]
      }
      # Calculate the mean squared error.
      errorRate <- mse(dfError[, "Imputed_Value"], dfError[, trait])
      # If the trait is categorical...
    } else if (trait %in% cat){
      # If origMiss contains at least one species.. &&
      if(length(origMiss) > 0){ 
        # Remove origMiss from dfTrue, dfImputed, and dfMissing. &&
        dfTruePFC <- dfTrue[!dfTrue$species_name %in% origMiss, ]
        dfImputedPFC <- dfImputed[!dfImputed$species_name %in% origMiss, ]
        dfMissingPFC <- dfMissing[!dfMissing$species_name %in% origMiss, ]
        # Calculate the PFC (percent falsely classified) with modified dataframes.
        errorRate <- pfc(x = dfTruePFC[[trait]], y = dfImputedPFC[[trait]], m = is.na(dfMissingPFC[, trait]))
        # If origMiss is empty..
      } else if(length(origMiss) == 0){
        # Calculate the PFC (percent falsely classified) using original dataframes.
        errorRate <- pfc(x = dfTrue[[trait]], y = dfImputed[[trait]], m = is.na(dfMissing[, trait]))
      }
    }
    
    # Name the errorRate according to the trait.
    names(errorRate) <- trait
    # Append to the list. 
    l_Error[[t]] <- errorRate
  }
  # Return list of error rates and dfImputed.
  return(list(dfImputed, l_Error))
}

ImputeMICE <- function(dfTrue, dfMissing, cols, cont, cat, inter = NULL, mSets, matPredictors){
  
  # Function for imputing missing values using the mice() function in the "MICE" package. Returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.
  # Citations: van Buuren S, Groothuis-Oudshoorn K (2011). “mice: Multivariate Imputation by Chained Equations in R.” Journal of Statistical Software, 45(3), 1-67. https://www.jstatsoft.org/v45/i03/.
  # R package version 3.13.0. https://cran.r-project.org/web/packages/mice/mice.pdf
  
  # dfTrue = dataframe containing true (known) observations
  # dfMissing = dataframe containing missing values and missingness indicator columns
  # cols = names of traits (column names)
  # cont = names of numerical traits
  # cat = names of categorical traits
  # inter = names of integer traits
  # mSets = vector containing values of m to test (number of multiply imputed dataframes)
  # matPredictors = predictor matrix indicated which traits to use in the imputation process (argument for mice)
  # phyImp = whether data are to be imputed using phylogenetic information
  
  # Ensure dfMissing is a dataframe.
  dfMissing <- as.data.frame(dfMissing)
  # Get the names of the missing columns as we will needs these later on.
  missingCols <- colnames(dfMissing)[grep(pattern = "_NA", x = colnames(dfMissing))]
  
  # Create an empty list to hold the imputed dataset for each value of m. ^
  l_dfImputedParam <- CreateNamedList(listLength = length(mSets), elementNames = mSets)
  # Create an empty list to hold the error rate results for each value of m and each trait.
  l_ErrorParams <- CreateNamedList(listLength = length(mSets), elementNames = mSets)
  
  # For each value of m...
  for(m in 1:length(mSets)) {
    # Take the number of mSets.
    M <- mSets[[m]]
    # Impute the datasets using M and matPredictors. 
    imputedMICE <- mice(dfMissing[, c(colnames(matPredictors))], predictorMatrix = matPredictors, m = M, maxit = 10, print = FALSE)
    # Combine the imputed data into one dataframe.
    dfImputed <- CombineMIDataframes(midata = imputedMICE, method = "MICE", m = M, contVars = cont, catVars = cat)
    # Make sure the variables are of the correct type.
    dfImputed[cont] <- lapply(dfImputed[cont], as.numeric)
    dfImputed[cat] <- lapply(dfImputed[cat], as.factor)
    # Merge dfImputed with the columns containing the missing info.
    dfImputed <- merge(dfImputed, dfMissing[, c("species_name", missingCols)], by = "species_name")
    # Merge with the true data.
    dfError <- merge(dfImputed, dfTrue, by = "species_name")
    # Create a vector to hold the results for each trait.
    param <- vector(mode = "numeric", length = length(cols))
    # Name according to traits.
    names(param) <- cols
    
    # For every trait...
    for(t in 1:length(cols)) {
      # Take the name of the tth trait.
      trait <- cols[[t]]
      # Identify index of imputed column.
      impCol <- grep(paste(trait, ".x", sep = ""), x = colnames(dfError))
      # Identify index of true column.
      trueCol <- grep(paste(trait, ".y", sep = ""), x = colnames(dfError))
      # Get the name of the corresponding missingness info column.
      missCol <- grep(pattern = trait, x = missingCols)
      missCol <- missingCols[missCol]
      # Identify species with missing values in original dataset. &&
      origMiss <- dfTrue$species_name[is.na(dfTrue[, trait])]
      # If the trait is numerical...
      if(trait %in% cont) {
        # Subset to only contain rows that were imputed.
        dfSubset <- dfError[which(dfError[, missCol] == "NA"), ]
        # If origMiss contains at least one species.. &&
        if(length(origMiss) > 0){ 
          # Remove from dfSubset. && name change for MICE
          dfSubset <- dfSubset[!dfSubset$species_name %in% origMiss, ]
        }
        # If trait is in inter (i.e. is count data)..
        if(trait %in% inter) {
          # Back-transformed the imputed value.
          dfSubset[, impCol] <- exp(dfSubset[, impCol])
          # Round to the nearest whole number, as random forest treated these data as continuous.
          dfSubset[, impCol] <- as.integer(round(dfSubset[, impCol]))
          # Log-transform once again so we can compare to the original log-transformed data.
          dfSubset[, impCol] <- log(dfSubset[, impCol])
        }
        # Calculate the mean squared error.
        errorRate <- mse(dfSubset[, trueCol], dfSubset[, impCol])
      } else if(trait %in% cat){
        # If origMiss contains at least one species..
        if(length(origMiss) > 0){ 
          # Remove origMiss from dfTrue, dfImputed, and dfMissing.
          dfTruePFC <- dfTrue[!dfTrue$species_name %in% origMiss, ]
          dfImputedPFC <- dfImputed[!dfImputed$species_name %in% origMiss, ]
          dfMissingPFC <- dfMissing[!dfMissing$species_name %in% origMiss, ]
          # Calculate the PFC (percent falsely classified) with modified dataframes.
          errorRate <- pfc(x = dfTruePFC[[trait]], y = dfImputedPFC[[trait]], m = is.na(dfMissingPFC[, trait]))
          # If origMiss is empty.. &
        } else if(length(origMiss) == 0){
          # Calculate the PFC (percent falsely classified) using original dataframes. (dfImputed in MICE)
          errorRate <- pfc(x = dfTrue[[trait]], y = dfImputed[[trait]], m = is.na(dfMissing[, trait]))
        }
      }
      # Name the errorRate according to the value of ntree.
      names(errorRate) <- trait
      # Append to the list. 
      param[[t]] <- errorRate
    }
    # Append to the mth elements in the lists. ^
    l_dfImputedParam[[m]] <- dfImputed
    l_ErrorParams[[m]] <- param
  }
  # Return list of error rates for each trait and parameter value and list of imputed dataframes for each parameter value. ^
  return(list(l_dfImputedParam = l_dfImputedParam, l_ErrorParams = l_ErrorParams))
  
  
}

ImputeRF <- function(dfTrue, dfMissing, cols, cont, cat, inter = NULL, ntrees, predictors, phyImp = F, l_dfMissing = NULL){
  # Function for imputing missing values using the missForest() function in the "missForest" package for imputation. Returns a list of error rates (MSE for continuous traits and PFC for categorical traits) for each parameter value tested.
  # Citations: Daniel J. Stekhoven (2013). missForest: Nonparametric Missing Value Imputation using Random Forest. R package version 1.4.
  # Stekhoven D. J., & Buehlmann, P. (2012). MissForest - non-parametric missing value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.

  # dfTrue = dataframe containing true (known) observations
  # dfMissing = dataframe containing missing values
  # cols = names of traits (column names)
  # cont = names of numerical traits
  # cat = names of categorical traits
  # inter = names of integer traits
  # ntrees = vector of ntree (number of trees to grow in the forest) values to consider
  # predictors = list of predictors for each trait
  # phyImp = whether data are to be imputed using phylogenetic information
  # l_dfMissing = if phyImp == T, a named list of dataframes with missing values corresponding to each trait
  
  # Get the names of the missing columns as we will need these later on.
  missingCols <- colnames(dfMissing)[grep(pattern = "_NA", x = colnames(dfMissing))]
  # Create an empty list to hold the error rate results for each trait and value of k. ^
  l_ErrorParams <- CreateNamedList(listLength = length(cols), elementNames = cols)
  # Create an empty list to hold the imputed data for each trait and value of k. ^
  l_l_dfImputedTrait <- CreateNamedList(listLength = length(cols), elementNames = cols)
  
  # For each trait..
  for(t in 1:length(cols)) {
    
    # Take the name of the tth trait.
    trait <- cols[[t]]
    # Get the name of the corresponding missingness info column.
    missCol <- grep(pattern = trait, x = missingCols)
    missCol <- missingCols[missCol]
    # Take the corresponding predictors for that trait.
    predIndex <- grep(pattern = trait, x = names(predictors))
    preds <- predictors[[predIndex]]
    # Identify species with missing values in original dataset.
    origMiss <- dfTrue$species_name[is.na(dfTrue[, trait])]
    if(phyImp == T){
      # Get the name of the corresponding dfMissing.
      index <- grep(pattern = trait, x = names(l_dfMissing))
      # Replace original dfMissing with dfMissing with appended eigenvectors for the trait.
      dfMissing <- l_dfMissing[[index]]
    }
    # Create a vector to hold the results for each value of ntree.
    paramRes <- vector(mode = "numeric", length = length(ntrees))
    # Name according to ntrees.
    names(paramRes) <- ntrees
    # Create a list to hold the imputed trait data for each value of ntree. ^
    l_dfImputedTrait <- CreateNamedList(listLength = length(ntrees), elementNames = ntrees)
    
    # Impute the trait using different values of ntree.
    for(n in 1:length(ntrees)) {
      # Take the number of ntree.
      N <- ntrees[[n]]
      # Impute the dataset (which contains the trait in question and its predictors) using missForest.
      imputedRF <- missForest(as.data.frame(dfMissing[, c(trait, preds)]), ntree = N)
      # Access the imputed dataframe.
      dfImputedTrait <- imputedRF$ximp
      # Add species_name info back.
      dfImputedTrait$species_name <- dfMissing$species_name
      # If the trait is numerical...
      if(trait %in% cont) {
        # Subset dfImputedTrait to only contain species_name, trait of interest, and information about whether the trait has been imputed.
        dfSubset <- dfImputedTrait[, c("species_name", trait)]
        # Rename column.
        colnames(dfSubset)[2] <- "Imputed_Value"
        # Merge with dfMissing to obtain missingness indicator column.
        dfSubset <- merge(dfSubset, dfMissing[c("species_name", missCol)], by = "species_name")
        # Subset to only contain rows that were imputed.
        dfSubset <- dfSubset[which(dfSubset[, missCol] == "NA"), ]
        # Merge with the true data.
        dfError <- merge(dfSubset, dfTrue, by = "species_name")
        # If origMiss contains at least one species..
        if(length(origMiss) > 0){ 
          # Remove from dfError.
          dfError <- dfError[!dfError$species_name %in% origMiss, ]
        }
        # If trait is in inter (i.e. is count data)..
        if(trait %in% inter) {
          # Back-transform the imputed value.
          dfError$Imputed_Value <- exp(dfError$Imputed_Value)
          # Round to the nearest whole number, as random forest treated these data as continuous.
          dfError$Imputed_Value <- as.integer(round(dfError$Imputed_Value))
          # Log-transform once again so we can compare to the original log-transformed data.
          dfError$Imputed_Value <- log(dfError$Imputed_Value)
        }
        # Calculate the mean squared error (MSE).
        errorRate <- mse(dfError[, trait], dfError[, "Imputed_Value"])
        # If the trait is categorical...
      } else if(trait %in% cat){
        # If origMiss contains at least one species..
        if(length(origMiss) > 0){ 
          # Remove origMiss from dfTrue, dfImputed, and dfMissing.
          dfTruePFC <- dfTrue[!dfTrue$species_name %in% origMiss, ]
          dfImputedPFC <- dfImputedTrait[!dfImputedTrait$species_name %in% origMiss, ] # ** name is different
          dfMissingPFC <- dfMissing[!dfMissing$species_name %in% origMiss, ]
          # Calculate the PFC (percent falsely classified) with modified dataframes.
          errorRate <- pfc(x = dfTruePFC[[trait]], y = dfImputedPFC[[trait]], m = is.na(dfMissingPFC[, trait]))
          # If origMiss is empty..
        } else if(length(origMiss) == 0){
          # Calculate the PFC (percent falsely classified) using original dataframes.
          errorRate <- pfc(x = dfTrue[[trait]], y = dfImputedTrait[[trait]], m = is.na(dfMissing[, trait]))
        }
      }
      # Name the errorRate according to the value of ntree.
      names(errorRate) <- N
      # Append results to the nth element in the lists. ^
      paramRes[[n]] <- errorRate
      l_dfImputedTrait[[n]] <- dfImputedTrait
    }
    # Append to the tth elements in the lists. ^
    l_ErrorParams[[t]] <- paramRes
    l_l_dfImputedTrait[[t]] <- l_dfImputedTrait
  }
  
  # Imputed dataframe handling. ---
  # Because we imputed data on a trait-by-trait basis, let's merge dataframes together for each parameter value so results are a bit easier to view.
  # Take species name and missingness information from dfMissing.
  dfImputedParam <- dfMissing[, c("species_name", missingCols)]
  # Make a list to hold the dataframes for each value of ntree.
  l_dfImputedParam <- lapply(1:length(ntrees), function(x) dfImputedParam)
  names(l_dfImputedParam) <- ntrees
  
  # For every trait...
  for (t in 1:length(l_l_dfImputedTrait)) { 
    # Take the tth trait.
    trait <- names(l_l_dfImputedTrait)[[t]]
    # Extract the tth element.
    l_dfImputedTrait <- l_l_dfImputedTrait[[t]]
    # For each value of ntree..
    for(n in 1:length(l_dfImputedTrait)){
      # Extract the nth dataframe in l_dfImputedTrait.
      dfImputedTrait <- l_dfImputedTrait[[n]]
      # Subset to only contain species_name and trait info.
      dfImputedTrait <- dfImputedTrait[, c("species_name", trait)]
      # Extract the nth dataframe in l_dfImputedParam.
      dfImputedParam <- l_dfImputedParam[[n]]
      # Merge with dfImputedTrait.
      dfImputedParam <- merge(dfImputedParam, dfImputedTrait, by = "species_name")
      # Replace nth element of l_dfImputedParam.
      l_dfImputedParam[[n]] <- dfImputedParam
    }
  }
  # Rearrange columns in l_dfImputedParam.
  l_dfImputedParam <- lapply(l_dfImputedParam, function(x) x[, c("species_name", cols, missingCols)])
  # Return list of error rates for each trait and parameter value and list of imputed dataframes for each parameter value.
  return(list(l_dfImputedParam = l_dfImputedParam, l_ErrorParams = l_ErrorParams))
}

LogTransform <- function(data, cols) {
  
  # Function for applying log-transformation to a dataframe's continuous variables. This function takes into account whether a variable contains negative values. If so, it will apply the NormalizeNegative function instead of log to apply a constant during the transformation.
  
  # data = dataframe containing continuous data
  # cols = column names containing continuous data
  
  # For continuous traits, check which ones contain negative values. 
  negTest <- lapply(data[, cols, drop = F], function(x) sum(na.omit(x) < 0) > 0)
  negTraits <- cols[which(negTest == T)]
  # Apply the NormalizeNegative function for log transformation for these traits
  data[, negTraits] <-  lapply(data[, negTraits, drop = F], NormalizeNegative)
  # Apply log-transformation to all other continuous traits.
  logTraits <- setdiff(cols, negTraits)
  data[, logTraits] <- lapply(data[, logTraits, drop = F], log)
  # Replace -Inf values with 0.
  data[data == -Inf] <- 0
  # Return dataframe with transformed values.
  return(data)
  
}

NormalizeNegative <- function(variable) {
  
  # Function for normalizing data with negative values by adding a constant (the absolute minimum value in the dataset) and log transforming the data.
  # variable = numeric vector
  
  # Add constant to each value.
  transformed <- variable + abs(min(variable, na.rm = T))
  # Add one to deal with 0 values. # TODO: Add 1 for zero values.
  transformed <- log(transformed) 
  # Return transformed variable.
  return(transformed)
  
}

RefineModel <- function(model, data){
  
  # Function for dropping insignificant terms from a glm object.
  # model = fitted glm object
  # data = dataframe with missing values used to fit glm
  
  # For each trait, create an indicator column for missingness using the bind_shadow() function. Setting only_miss = T so only those variables with missing values have corresponding missing indicator columns created.
  dfShadow <- bind_shadow(data, only_miss = T)
  # Get the names of the missingness indicator columns.
  missCols <- colnames(dfShadow)[grep(pattern = "_NA", x = colnames(dfShadow))]
  # For missCols, convert !NA to 1 (present) and NA to 0 (absent).
  dfShadow[, missCols] <- lapply(dfShadow[, missCols], function(x) ifelse(x == "!NA", 1, 0))
  
  # Extract the response variable from the model.
  response <- all.vars(model$formula[[2]])
  # Extract the covariates from the model.
  covariates <- all.vars(model$formula[[3]])
  # Extract the p-values from the model.
  pVals <- coef(summary(model))[, 4]
  # Remove the intercept.
  pVals <- pVals[-1]
  
  # Take the complete-case dataset for response and covariates. We need this for comparison purposes (e.g. fitting an anova) when we fit an alternative glm without a particular term.
  dfShadow <- na.omit(dfShadow[, c(response, covariates)])
  
  # Identify factor variables.
  catTraits <- colnames(dfShadow)[sapply(dfShadow, is.factor)]
  # Determine number of levels in the factor variables.
  factorLevels <- sapply(dfShadow[, catTraits], function(x) length(levels(x)))
  # Identify the names of the binary variables.
  binCat <- names(which(factorLevels == 2))
  # Identify the names of the multi-categorical variables.
  multiCat <- names(which(factorLevels > 2))
  
  # If there are factor variables..
  if(length(catTraits) > 0){
    # For every factor variable..
    for(m in 1:length(factorLevels)){
      # Get name of variable.
      FV <- names(factorLevels)[[m]]
      # If the variable is binary..
      if(FV %in% binCat) {
        # Identify p-value associated with FV.
        index <- grep(FV, names(pVals))
        # Since the name of the category compared to the reference level is appended to the name of the p-value, we are replacing the modified name with the name of the original variable (because we are eventually modifying the glm formula using the original variable name).
        names(pVals)[index] <- FV
        # If the variable is multi-categorical..
      } else if(FV %in% multiCat) {
        # Remove m from covariates.
        modCov <- covariates[!covariates %in% FV]
        # Create a new formula using modCov.
        modForm <- as.formula(paste(response, "~", paste(modCov, collapse = "+")))
        # Fit a logistic regression model without m.
        modFit <- glm(modForm, data = dfShadow, family = "binomial", na.action = na.omit)
        # Compare to original model.
        result <- anova(modFit, model, test = "LRT")
        # Extract the p-value from the result. This is necessary to obtain a p-value for the overall factor variable because it is multi-categorical.
        modP <- result$`Pr(>Chi)`[2]
        # Name according to FV.
        names(modP) <- FV
        # Identify p-values associated with different levels of multicategorical variable.
        index <- grep(FV, names(pVals))
        # Remove from original pVals.
        pVals <- pVals[-index]
        # Append new p-value for multi-categorical variable to pVals.
        pVals <- c(pVals, modP)
      }
    }
  }

  # If there are any non-significant terms in the model..
  if(any(pVals > 0.05)){
    # Set ALLSIG == F. This is an indicator variable to check that all terms are significant in the glm model. This is set to T when all terms are significant (p < 0.05) and the loop will end.
    ALLSIG <- F
    # Order p-values largest to smallest.
    pVals <- sort(pVals, decreasing = T)
    # Create a vector to hold the original number of covariates.
    covariates <- names(pVals)
    # Create a vector to track the updated number of covariates to include in the model.
    tempCovariates <- covariates
    while(ALLSIG == F){
      # For every covariate..
      for(c in 1:length(covariates)){
        # Take name of cth covariate.
        covariate <- covariates[[c]]
        # Match covariate to corresponding p-value.
        p <- pVals[grep(pattern = covariate, covariates)]
        # If the term is insignificant..
        if(p > 0.05){
          # Update the tempCovariates to exclude the term.
          newCands <- tempCovariates[!tempCovariates %in% covariate]
          # Create new formula without the term.
          newForm <- as.formula(paste(response, "~", paste(newCands, collapse = "+")))
          # Fit a new logistic regression model.
          newFit <- glm(newForm, data = dfShadow, family = "binomial", na.action = na.omit)
          # Extract the p-values from the new model.
          newP <- coef(summary(newFit))[, 4]
          # Remove the intercept.
          newP <- newP[-1]
          # If all the terms in the new model are significant...
          if(all(newP < 0.05)){
            # Set ALLSIG == TRUE.
            ALLSIG <- TRUE
            # Return newly fitted model.
            return(newFit)
          } else {
            # Update tempCovariates vector to include the most recent set of covariates.
            tempCovariates <- newCands
          }
        } ## if
      } ## for
    } ## WHILE
  } else {
    print("No terms to remove!")
    # Return the original model.
    return(model)
  }
  
}

SelectPredictors <- function(data) {
  
  # Function for performing regression analyses between dataframe columns, according to data class.
  # data = dataframe
  # cols = column numbers to consider for regression analyses

  # Create a list to hold the significant predictors for each response variable (column) in the dataframe.
  l_predictors <- vector(mode = "list", length = ncol(data))
  # Name l_predictors according to the names of the columns.
  names(l_predictors) <- colnames(data)
  
  # For every column in the dataframe...
  for(i in 1:ncol(data)) {
    # Take the name of the ith column.
    resName <- colnames(data)[i]
    # Take the ith column.
    resClass <- data[, i]
    # Remove the response variable from data.
    dfPredictors <- data[, -i]
    # Create empty vector to hold predictors for the response variable in question.
    predictors <- vector(mode = "character")
    
    # If it is a continuous response variable...
    if(is.double(resClass)) {
      # For every column in dfPredictors...
      for (j in 1:ncol(dfPredictors)) {
        # Get the name of the predictor.
        predictorName <- colnames(dfPredictors[j])
        # Subset to get complete observations for response and pred.
        dfSubset <- na.omit(data[, c(resName, predictorName)])
        # Take the response column.
        response <- dfSubset[, resName]
        # Take the predictor column.
        pred <- dfSubset[, predictorName]
        # If the predictor is numeric (0 levels) OR there are 2 categories in the predictor (binary)...
        if(length(levels(pred)) == 0 | length(levels(pred)) == 2) {
          # Fit glm model with gaussion family specified.
          fit <- glm(response ~ pred, family = "gaussian")
          # Extract the p-values from the model.
          pVal <- coef(summary(fit))[, 4]
          # Remove the intercept.
          pVal <- pVal[-1]
          # If there are more than 2 categories in the predictor...
        } else if(length(levels(pred)) > 2) {
          # Fit the null model.
          fitNULL <- glm(response ~ 1, family = "gaussian")
          # Fit the alternative model.
          fitALT <- glm(response ~ pred, family = "gaussian")
          # Run an anova to compare the models.
          result <- anova(fitNULL, fitALT, test = "F")
          # Extract the p-value from the model.
          pVal <- result$`Pr(>F)`[2]
        }
        if(pVal < 0.15) {
          # Add the name of the significant predictor to the list of predictors.
          predictors[[j]] <- predictorName
        } 
        # Remove empty elements from the predictors vector.
        predictors <- na.omit(predictors)
        # If predictors is empty or there is only 1, add all column names as predictors.
        if(length(predictors) == 0) {
          predictors <- colnames(dfPredictors)
        }
        # Append to the list of predictors.
        l_predictors[[i]] <- predictors
      }
      
      # If it is an integer class response variable (count data)...
    } else if(class(resClass) == "integer") {
      # For every column in dfPredictors...
      for (j in 1:ncol(dfPredictors)) {
        # Take the jth column.
        pred <- dfPredictors[, j]
        # Get the name of the predictor.
        predictorName <- colnames(dfPredictors[j])
        # Subset to get complete observations for response and pred.
        dfSubset <- na.omit(data[, c(resName, predictorName)])
        # Take the response column.
        response <- dfSubset[, resName]
        # Take the predictor column.
        pred <- dfSubset[, predictorName]
        # If the predictor is numeric (0 levels) OR there are 2 categories in the predictor (binary)...
        if(length(levels(pred)) == 0 | length(levels(pred)) == 2) {
          # Fit glm model with poisson family specified.
          fit <- glm(response ~ pred, family = "poisson")
          # Extract the p-values from the model.
          pVal <- coef(summary(fit))[, 4]
          # Remove the intercept.
          pVal <- pVal[-1]
          # If there are more than 2 categories in the predictor...
        } else if(length(levels(pred)) > 2) {
          # Fit the null model.
          fitNULL <- glm(response ~ 1, family = "poisson")
          # Fit the alternative model.
          fitALT <- glm(response ~ pred, family = "poisson")
          # Run an anova to compare the models.
          result <- anova(fitNULL, fitALT, test = "Chisq")
          # Extract the p-value from the model.
          pVal <- result$`Pr(>Chi)`[2]
        }
        if(pVal < 0.15) {
          # Add the name of the significant predictor to the list of predictors.
          predictors[[j]] <- predictorName
        } 
        # Remove empty elements from the predictors vector.
        predictors <- na.omit(predictors)
        # If predictors is empty, add all column names as predictors.
        if(length(predictors) == 0) {
          predictors <- colnames(dfPredictors)
        }
        # Append to the list of predictors.
        l_predictors[[i]] <- predictors
      }
      
    } # If it is a categorical response variable...
    else if(class(resClass) == "factor") {
      # If it is a binary response variable...
      if(length(levels(resClass)) == 2) {
        # For every column in dfPredictors...
        for (j in 1:ncol(dfPredictors)) {
          # Get the name of the predictor.
          predictorName <- colnames(dfPredictors[j])
          # Subset to get complete observations for response and pred.
          dfSubset <- na.omit(data[, c(resName, predictorName)])
          # Take the response column.
          response <- dfSubset[, resName]
          # Take the predictor column.
          pred <- dfSubset[, predictorName]
          # If the predictor is numeric (0 levels) OR there are 2 categories in the predictor (binary)...
          if(length(levels(pred)) == 0 | length(levels(pred)) == 2) {
            # Fit glm model with binomial family specified. Every other variable in the dataframe will be specified as the predictor variable in turn.
            fit <- glm(response ~ pred, family = "binomial")
            # # Extract the p-values from the model.
            pVal <- coef(summary(fit))[, 4]
            # # Remove the intercept.
            pVal <- pVal[-1]
            # If there are more than 2 categories in the predictor...
          } else if(length(levels(pred)) > 2) {
            # Fit the alternative model.
            fitALT <- glm(response ~ pred, family = "binomial")
            # Run an anova to compare the models.
            result <- anova(fitALT, test = "LRT")
            # Extract the p-value from the model.
            pVal <- result$`Pr(>Chi)`[2]
          }
          if(pVal < 0.15) {
            # Add the name of the significant predictor to the list of predictors.
            predictors[[j]] <- predictorName
          }
        }
        
        # Remove NA values from predictors.
        predictors <- na.omit(predictors)
        # If predictors is empty, add all column names as predictors.
        if(length(predictors) == 0) {
          predictors <- colnames(dfPredictors)
        }
        # Append to the list of predictors.
        l_predictors[[i]] <- predictors
        
        # If it is a multicategorical response variable...
      } else if(length(levels(resClass)) > 2) {
        # For every column in dfPredictors...
        for (j in 1:ncol(dfPredictors)) {
          # Get the name of the predictor.
          predictorName <- colnames(dfPredictors[j])
          # Subset to get complete observations for response and pred.
          dfSubset <- na.omit(data[, c(resName, predictorName)])
          # Take the response column.
          response <- dfSubset[, resName]
          # Take the predictor column.
          pred <- dfSubset[, predictorName]
          # Fit multinom model. Every other variable in the dataframe will be specified as the predictor variable in turn.
          # Fit the null model.
          fitNULL <- multinom(response ~ 1)
          # Fit the alternative model.
          fitALT <- multinom(response ~ pred)
          # Run a likelihood ratio test.
          result <- lrtest(fitNULL, fitALT) 
          # Extract the p-value.
          pVal <- result$`Pr(>Chisq)`[2]
          if(pVal < 0.15) {
            # Add the name of the significant predictor to the list of predictors.
            predictors[[j]] <- predictorName
          }
        }
        # Remove NA values from predictors.
        predictors <- na.omit(predictors)
        # If predictors is empty, add all column names as predictors.
        if(length(predictors) == 0) {
          predictors <- colnames(dfPredictors)
        }
        # Append to the list of predictors.
        l_predictors[[i]] <- predictors
      }
    }
  }
  # Return list of predictors for each response variable.
  return(l_predictors)
  
}

SimMAR <- function(model, data){
  
  # Function for simulating missing at random (MAR) data in a complete-case dataset. Combines predict.glm() and rbinom() functions and returns a vector with NAs introduced MAR.
  # model = fitted glm object
  # data = dataframe containing covariate data in model$formula
  
  # Create dataframe for introducing NAs.
  dfMissing <- data
  # Extract response variable and remove "_NA" from name.
  response <- gsub(x = all.vars(model$formula[[2]]), pattern = "_NA", replacement = "")
  # Extract covariates from model.
  covariates <- all.vars(model$formula[[3]])
  # If there's only one covariate..
  if(length(covariates) == 1) {
    # Subset dataframe to only include covariates. Also ensure it is dataframe format.
    dfCov <- as.data.frame(data[, covariates])
    # Set column name to name of covariate (it does not do this automatically when only subsetting a single column).
    colnames(dfCov) <- covariates
    # If there's more than one covariate..
  } else if(length(covariates) > 1) {
    # Subset dataframe to only include covariates. Also ensure it is dataframe format.
    dfCov <- as.data.frame(data[, covariates])
  }
  # Using predict.glm function, use the final logistic regression model to estimate probability that a value is present for a certain variable. newdata argument to set to the subsetted complete-case dataframe (dfCov) as we are predicting probability of missingness in the response using this data. type argument is set to "response" as we are predicting probability of a value being present in the response variable of the model.
  indNA <- predict.glm(model, newdata = dfCov, type = "response")
  # Using rbinom function to generate values following a binomial distribution. n is set to 1 observation, size is set to 1 trial. Using ifelse function to skip over existing NA values.
  preds <- sapply(indNA, function(x) ifelse(test = is.na(x), x, rbinom(x, n = 1, size = 1)))
  # Check that missingness proportion is at least 8% of trait data so that we that sample size for error rate calculations are reliable.
  zerosProp <- sum(na.omit(preds) == 0)/length(na.omit(preds))
  # If missingness proportion is less than 0.08 (0s in preds vector)..
  if(zerosProp < 0.08){
    # While the proportion is less than 0.08..
    while (zerosProp < 0.08) {
      # Reduce intercept in model by half.
      model$coefficients[1] <- (model$coefficients[1])/2
      # Estimate new probabilities. 
      indNA <- predict.glm(model, newdata = dfCov, type = "response")
      # Generate new values following a binomial distribution
      preds <- sapply(indNA, function(x) ifelse(test = is.na(x), x, rbinom(x, n = 1, size = 1)))
      # Recalculate zerosProp.
      zerosProp <- sum(na.omit(preds) == 0)/length(na.omit(preds))
    }
  }
  # If the response is a factor (categorical) type..
  if(class(dfMissing[, response]) == "factor"){
    # Convert to character type so we can retain the original names of the categories and it's not converted to integer type.
    dfMissing[, response] <- as.character(dfMissing[, response])
  }
  # Replace 0s with NAs and 1s with original value in response.
  resMiss <- ifelse(preds == 0, yes = NA, no = dfMissing[, response])
  # Return the response variable with introduced NAs.
  return(resMiss)
  
}

SimCatMNAR <- function(var, category, proportion = 0.15){
  
  # Function for simulating data that are missing not at random (MNAR) in categorical variables.
  # var = categorical (factor) vector in which to introduce NAs MNAR
  # category = category of var in which to introduce missing values
  # proportion = proportion of category to remove
  
  # Identify observations in categories
  index <- which(var == category)
  # Determine how many values are 15% of var.
  numNAs <- round(proportion*(length(na.omit(var))))
  # Randomly sample numNAs from index.
  rmThese <- sample(index, numNAs)
  # Remove from var.
  var[rmThese] <- NA
  # Return var.
  return(var)
  
}

SimContMNAR <- function(var, quantile = 0.9, direction = "greater", absolute = F){
  
  # Function for simulating data that are missing not at random (MNAR) in numeric variables.
  # var = numeric vector in which to introduce NAs MNAR
  # quantile = threshold at which to introduce NAs
  # direction = either "greater" or "less". Remove values that exceed this threshold.
  # absolute = if TRUE, take absolute values when estimating thresholds
  
  if(absolute == F){
    # Determine quantile.
    extreme <- quantile(var, probs = quantile, na.rm = T)
  } else if(absolute == T){
    # Determine quantile using absolute values.
    extreme <- quantile(abs(var), probs = quantile, na.rm = T)
  }
  # Identify values that exceed this threshold.
  # If greater than..
  if(direction == "greater"){
    if(absolute == F){
      # If the extreme is equal to the maximum (i.e. in case of skewed count data..)
      if(extreme == max(na.omit(var))){
        # Identify observations equal to extreme.
        indExt <- which(var == extreme)
        # Determine how many values are 10% of var.
        numNAs <- round(0.10*(length(na.omit(var))))
        # Randomly sample numNAs from indExt.
        index <- sample(indExt, numNAs)
      } else {
        # Take values that exceed the value selected for extreme.
        index <- which(var > extreme)
      }
    } else if(absolute == T){
      # Take absolute values that exceed the value selected for extreme.
      index <- which(abs(var) > extreme)
    }
    # If less than..
  } else if(direction == "less"){
    if(absolute == F){
      # If the extreme is equal to the min (i.e. in case of skewed count data..)
      if(extreme == min(na.omit(var))){
        # Identify observations equal to extreme.
        indExt <- which(var == extreme)
        # Determine how many values are 10% of var.
        numNAs <- round(0.10*(length(na.omit(var))))
        # Randomly sample numNAs from indExt.
        index <- sample(indExt, numNAs)
      } else {
        # Identify observations that exceed or match the value selected for extreme.
        indExt <- which(var <= extreme)
        # Determine how many values are 10% of var.
        numNAs <- round(0.10*(length(na.omit(var))))
        # Randomly sample numNAs from indExt.
        index <- sample(indExt, numNAs)
      }
    } else if(absolute == T){
      # Take absolute values that are less than the value selected for extreme.
      index <- which(abs(var) < extreme)
    }
  }
  # Remove indices from var.
  var[index] <- NA
  # Return var.
  return(var)
  
}

WriteErrorRates <- function(data, fileName) {
  
  # Function for writing dataframes with error rate info to file.
  # df = dataframe error rate info.
  # fileName = name to give to csv file
  
  # Write the dataframe to file.
  write.csv(data, fileName)
  print(paste("Wrote records to ", fileName, sep = ""))
  
}
