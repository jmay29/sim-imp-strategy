### OPTIONAL: Creation of "nearly" complete-case dataset. ----

# Use the following code if you want to allow for some missingness in your complete-case dataset. I.e. a "nearly" complete case dataset. If using phylogenetic imputation, you may wish to do this after you have already matched your species names in the trait dataset to those in the tree(s) to ensure you have enough data (e.g. n > 100). It is recommended here to not exceed 10% missingness in any individual trait so that simulation/imputations have enough observed data to provide accurate results.

# Load libraries and functions.
library(data.table)
source("R/Functions/DataHandling_Functions.R")

# Read in trait dataset you wish to impute.
fileName <- file.choose()
dfTraits <- read.csv(fileName)
# Ensure blanks are replaced with NAs.
dfTraits[dfTraits == ""] <- NA
# Create vector of column names that contain taxonomic/misc information. For example:
taxCols <- "species_name"
# Extract trait names for this taxon (everything other than taxCols). I.e. the traits you want to consider for imputation.
traits <- setdiff(colnames(dfTraits), taxCols)

# Make copy of dfTraits.
dfTraitsCount <- dfTraits
# Create a new column called trait count.
dfTraitsCount$trait_count <- apply(dfTraitsCount[, traits], MARGIN = 1, function(x) sum(!is.na(x)))

# Create empty list to hold subsetted dataframes (recommend including at least 3 traits for imputing).
l_dataframes <- CreateNamedList(listLength = length(3:max(dfTraitsCount$trait_count)), elementNames = 3:max(dfTraitsCount$trait_count))

# Loop over different trait missingness thresholds.
for(t in max(dfTraitsCount$trait_count):3){
  
  print(paste("Results for", t, "trait minimum:", sep = " "))
  # Subset to only include rows with data for at least t traits.
  dfSubset <- dfTraitsCount[dfTraitsCount$trait_count >= t,]
  # Subset to only keep traits of interest and taxonomic information.
  dfSubset <- dfSubset[, c(taxCols, traits)]
  
  print("Missingness proportion in each trait column:")
  print(apply(dfSubset[, traits], MARGIN = 2, function(x) sum(is.na(x))/nrow(dfSubset)))
  
  print("Sample size:")
  print(nrow(dfSubset))
  
  # Append to l_dataframes.
  index <- grep(t, names(l_dataframes))
  l_dataframes[[index]] <- dfSubset
  
}

# Extract the dataframe you want to use for your "nearly-complete" dataset. For example, if we wanted to keep 5 traits:
myNC <- which(names(l_dataframes) %in% 5)
dfNC <- l_dataframes[[myNC]]
# Write dataframe to file.
fwrite(dfNC, "NearlyCC.csv")
