# Data handling functions for preparing trait and sequence data for imputation runs.

BreakIntoTypes <- function(data, traitCols) {
  
  # Function for identifying trait types of columns in a dataframe.
  
  # data = dataframe containing trait and taxonomic data
  # traitsCols = columns containing trait data
  
  # Identify categorical traits.
  catTraits <- GetTraitNames(data = data[, traitCols], class = "character")
  # Identify binary traits and append them to catTraits as they will also be treated as factors.
  binTraits <- GetTraitNames(data = data[, traitCols], class = "binary")
  catTraits <- c(binTraits, catTraits)
  # Convert categorical traits to factors to differentiate them from numeric traits.
  # If there is more than 1 categorical trait..
  if(length(catTraits) > 1){
    # Use lapply.
    data[, catTraits] <- lapply(data[, catTraits], as.factor)
    # if there's only one categorical trait..
  } else if(length(catTraits) == 1){
    data[[catTraits]] <- as.factor(data[[catTraits]])
  }
  # Identify continuous traits.
  contTraits <- GetTraitNames(data = data[, traitCols], class = "numeric")
  
  # Combine contTraits and catTraits into a list.
  l_traits <- list(contTraits, catTraits)
  # Name according to trait type.
  names(l_traits) <- c("Numerical", "Categorical")
  
  # Return l_traits.
  return(l_traits)
  
}

CreateNamedList <- function(listLength, elementNames) {
  
  # Function for creating a named list.
  # listLength = length of list
  # elementNames = vector element containing names for the list
  
  # Create empty list with length = listLength.
  l_list <- vector(mode = "list", length = listLength)
  # Name the list with provided vector.
  names(l_list) <- elementNames
  
  # Return named list.
  return(l_list)
  
}

CreateTraitSubset <- function(data, traits, taxonomy, speciesColumn, critNumber) {
  
  # This function identifies the traits that maximize the number of complete cases and returns a complete-case dataframe.
  
  # data = dataframe containing trait information
  # traits = column names of trait data to consider (should be names of traits)
  # taxonomy = column names that contain taxonomy information
  # speciesColumn = character vector that contains name of species column in dataframe
  # critNumber = the minimum sample size required for the complete-case subset
  
  # Use the IdentifyTraits to identify how many complete cases we get for each trait subset.
  traitNum <- IdentifyTraits(data, variables = traits)
  # How many of the traits create subsets greater than critNumber?
  goodTraits <- names(which(traitNum > critNumber))
  
  # If goodTraits contains more than 3 traits...
  if(length(goodTraits) >= 3) {
    
    # Add the taxonomy info back.
    goodTraits <- c(taxonomy, goodTraits)
    # Get the complete-case dataset using the TakeCompleteCases function.
    dfComplete <- TakeCompleteCases(data, goodTraits)
    # Getting it ready for merging with sequence data. Renaming Species column and adding underscores to species names.
    speciesIndex <- which(colnames(dfComplete) == speciesColumn)
    colnames(dfComplete)[speciesIndex] <- "species_name"
    dfComplete$species_name <- gsub(" ", "_", dfComplete$species_name)
    
    # Return the complete-case subset.
    return(dfComplete)
    
  } else if(length(goodTraits) < 3) {
    
    print(paste("Not enough trait data at this sample size:", critNumber, sep = " "))
    
  }
  
}

CritNumberLoop <- function(critNumbers =  seq(100, 5000, by = 100), data, traitCols, taxCols) {
  
  # Function to create complete-case datasets for trait data, considering different criterion for the required sample size.
  
  # critNumbers = numeric vector containing sample size thresholds to consider, minimum must be 100
  # data = dataframe containing trait data
  # traitCols = column names containing trait data
  # taxCols = column names containing taxonomy data
  
  # For each value of critNumber...
  for(i in 1:length(critNumbers)) {
    
    # Take the ith critNumber.
    num <- critNumbers[[i]]
    # Create complete-case subsets for each taxon using the CreateTraitSubset() function.
    dfCompleteCase <- CreateTraitSubset(data, traits = traitCols, taxonomy = taxCols, speciesColumn = "species_name", critNumber = num)
    # If sample size criterion has been met (dfCompleteCase is a dataframe)..
    if(class(dfCompleteCase) == "data.frame"){
      # Write the dataframe to file using the WriteCompleteCases() function.
      WriteCompleteCases(dfCompleteCase)
    }
    
  }
  
}

GetCategoricalInfo <- function(vector) {
  
  # Function for providing descriptive information for a categorical (character) vector.
  # vector = character vector
  
  # Create empty list to hold results.
  result <- vector(mode = "list", length = 5)
  # Name the list.
  names(result) <- c("n", "NumberOfCategories", "Frequencies", "Proportions")
  
  # Get the sample size (n) of the data.
  result[[1]] <- length(na.omit(vector))
  # Get the number of categories (excluding NAs).
  result[[2]] <- length(unique(na.omit(vector)))
  # Get the categorical frequencies including NAs.
  result[[3]] <- table(vector, useNA = "ifany")
  # Get the categorical proportions including NAs.
  result[[4]] <- table(vector, useNA = "ifany")/length(vector)
  
  # Return the list of results.
  return(result)
  
}

GetNumericalInfo <- function(vector) {
  
  # Function for providing descriptive information for a numerical vector.
  # vector = numerical vector
  
  # Create empty list to hold results.
  result <- vector(mode = "list", length = 7)
  # Name the list.
  names(result) <- c("n", "Mean", "Median", "Range", "NumberNAs", "ProportionNAs", "Plot")
  
  # Get the sample size (n) of the data.
  result[[1]] <- length(na.omit(vector))
  # Get the mean of the data.
  result[[2]] <- mean(vector, na.rm = T)
  # Get the median of the data.
  result[[3]] <- median(vector, na.rm = T)
  # Get the range of the data.
  result[[4]] <- range(vector, na.rm = T)
  # Get the numer of NAs in the data.
  result[[5]] <- sum(is.na(vector))
  # Get the proportion of NAs in the data.
  result[[6]] <- sum(is.na(vector))/length(vector)
  # Plot the data.
  result[[7]] <- hist(vector, col = "skyblue")
  
  # Return the list of results.
  return(result)
  
}

GetTraitNames <- function(data, class = "numeric") { 
  
  ## TODO: Fix warning message if single trait provided.
  
  # Function for extracting trait names from a dataframe that are of a certain class.
  
  # data = dataframe containing trait data
  # class = type of class for which to extract names of trait columns (either "numeric", "binary", or "character").
  
  if(class == "numeric") {
    # Apply is.numeric() function across columns.
    cons <- sapply(data, is.numeric)
    # Get the names of the columns that are TRUE.
    traits <- names(which(cons == TRUE))
  } else if(class == "binary") {
    # Apply is.binary() function from the "Information" package across columns.
    bin <- sapply(data, Information::is.binary)
    # Get the names of the columns that are TRUE.
    traits <- names(which(bin == TRUE))
  } else if(class == "integer") {
    # Apply is.character() function across columns.
    ints <- sapply(data, is.integer)
    # Get the names of the columns that are TRUE.
    traits <- names(which(ints == TRUE))
  } else if(class == "character") {
    # Apply is.character() function across columns.
    chars <- sapply(data, is.character)
    # Get the names of the columns that are TRUE.
    traits <- names(which(chars == TRUE))
  }
  
  # Return names of the traits.
  return(traits)
  
}

IdentifyTraits <- function(data, variables) {
  # Function for identifying the traits that maximize the number of complete-cases.
  
  # data = dataframe containing species and trait data
  # variables = vector of trait names (should match column names in data)
  
  # Convert data to data.table format.
  dfMultivariable <- as.data.table(data[, variables])
  # Order the columns by the amount of missing data (NA values).
  suppressWarnings(dfTraitsNA <- sort(dfMultivariable[, lapply(.SD, function(x) sum(is.na(x)))]))
  # Reorder the dataset so that the columns with the least amount of NA values are now first.
  setcolorder(dfMultivariable, names(dfTraitsNA))
  
  # Now I want to loop through the traits, removing one column (trait) at a time and count the number of complete cases. This will provide us some information as to which traits would provide an adequate sample size for downstream analysis.
  # Take the number of columns in dfMultivariable.
  len <- ncol(dfMultivariable)
  # Create a numeric vector to hold the results of the loop.
  all.cc <- NULL
  # Start the loop:
  for (i in 1:len) {
    # Works best if you set dfMultivariable back to a dataframe.
    x <- as.data.frame(dfMultivariable)
    # x is the placeholder dataframe in the loop.
    x <- x[, 1:len]
    # Determine which rows are "complete" using the "len" subset of traits.
    x <- complete.cases(x)
    # Complete rows of data will be "TRUE".
    x <- which(x == "TRUE")
    # Find the number of complete cases.
    x <- length(x)
    # Add it to the all.cc variable that's holding all of the results of the loop.
    all.cc[i] <- x
    # Minus 1 from tempLen so we can check the next subset of traits (we started at the last column because the columns were ordered by number of NA values).
    len <- len - 1
  }
  # Now, decide where to cut the datatable. (i.e. pick an adequate subset of traits that maximize sample size).
  names(all.cc) <- rev(colnames(dfMultivariable))
  # Look at the results.
  all.cc
  
  # Return the vector that contains information about how to maximize complete cases.
  return(all.cc)
  
}

OutlierCheck <- function(data, col){
  
  # Function for identifying outliers based on the interquartile range (IQR) of the data. Returns a subsetted dataframe of the outliers so user can check the values for potential error.
  # data = dataframe. Make sure dataframe also contains a column called "species_name" with species name information.
  # col = column name for which to identify outliers
  
  # Code for calculating thresholds adapted from: 
  # Author: https://stackoverflow.com/users/1312519/by0.
  # https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r.
  
  # Determine the 25% quantile.
  lowerQuantile <- quantile(data[, col], na.rm = T)[2]
  # Determine the 75% quantile.
  upperQuantile <- quantile(data[, col], na.rm = T)[4]
  # Calculate the IQR.
  iqr <- upperQuantile - lowerQuantile
  # Calculate upper threshold ((3 x the IQR) + upperQuantile).
  upperThreshold <- (iqr * 3) + upperQuantile
  # Calculate lower threshold (lowerQuantile - (3 x the IQR)).
  lowerThreshold <- lowerQuantile - (iqr * 3)
  
  # Identify outliers based on whether they exceed the upper or lower thresholds.
  outliers <- which(data[, col] > upperThreshold | data[, col] < lowerThreshold)
  
  # Subset the outliers. Ensure species name information is also kept.
  dfOutliers <- data[outliers, c("species_name", col)]
  
  if(nrow(dfOutliers) == 0){
    print(paste("No outliers detected for:", col, sep = " "))
  } else {
    print(paste("Outliers detected for:", col, sep = " "))
    # Return the dataframe of outliers and their corresponding extreme values.
    return(dfOutliers)
  }
  
}

PasteColNames <- function(df, colNames, string, ...) {
  
  # Function for pasting a string onto a subset of column names in a dataframe.
  # df = dataframe
  # colNames = column names to paste string onto
  # string = string to paste onto column names
  # ... = arguments for paste
  
  # Get trait col indices.
  index <- which(colnames(df) %in% colNames)
  # Paste col names so we can identify complete data.
  colnames(df)[index] <- paste(colnames(df[, colNames]), string, ...)
  
  # Return dataframe with new column names.
  return(df)
  
}

ScreenCategories <- function(variable, threshold = 0.10) {
  
  # Function for determining whether a variable contains categories with observations below a certain threshold.
  
  # variable = character/factor type vector
  # varName = names of variable
  # threshold = numeric threshold (e.g. 0.10)
  
  # Count the number of observations per category.
  counts <- plyr::count(variable)
  # Calculate the proportions for each category.
  proportions <- counts$freq/sum(counts$freq)
  # Name the proportions according to category.
  names(proportions) <- counts$x
  # Create an empty list to hold categories that meet the required threshold.
  goodCat <- vector(mode = "numeric")
  
  print("Results:")
  
  # # For every category...
  for(i in 1:length(proportions)) {
    # Take the ith category.
    category <- proportions[i]
    # If the category comprises less than 10% of the dataset...
    if(category < threshold) {
      print(paste(category, "is less than", threshold, "of dataset. Remove", names(category), "!"))
    } else if(category > threshold) {
      print(paste(category, "is greater than", threshold, "of dataset. Keep", names(category), "!"))
      # Append to goodCat.
      goodCat[[i]] <- category
      # Name according to category.
      names(goodCat)[[i]] <- names(category)
    }
  }
  
  # Remove NAs from goodCat.
  goodCat <- na.omit(goodCat)
  
  # Check if there is more than one level in the variable. Otherwise the whole trait should be removed because it is invariant.
  if(length(goodCat) == 1) {
    print("This variable only contains one category. Removal from dataset is recommended.")
  } else {
    print(paste("This variable contains", length(goodCat), "categories. It can remain in the dataset."))
  }
  
  # Return the proportions for each category.
  return(goodCat)
  
}

TakeCompleteCases <- function(data, variables) {
  # Function for subsetting dataset based on column names providing and removing NAs.
  # data = dataframe containing species and trait information
  # variables = column names to subset
  
  # Subset the column names in variables.
  dfSubset <- subset(data, select = variables)
  # Take only the complete-cases.
  dfSubset <- na.omit(dfSubset)
  
  return(dfSubset)
  
}

WriteCompleteCases <- function(data) {
  
  # Function for writing complete case dataset to file.
  # data = dataframe containing complete case trait information. Should be at the order level.
  
  # Create file name.
  fileName <- paste(nrow(data), ".csv", sep = "")
  
  # Write the dataframe to file.
  write.csv(data, fileName)
  print(paste("Wrote records to ", fileName, sep = ""))
  
}
