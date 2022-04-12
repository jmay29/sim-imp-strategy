# Functions related to phylogenetic analyses.

AppendEigenvectors <- function(data, vars, tree, predictors){
  
  # Function for appending eigenvectors to dataframe containing trait data with missing values. Returns a list of dataframes with corresponding eigenvectors for each trait.
  
  # data = dataframe with missing values
  # vars = names of traits (columns in data)
  # tree = phylo object
  # predictors = list of predictors for each trait
  
  # Create a list to hold dataframes of predictors and eigenvectors for each trait.
  l_dfMiss <- CreateNamedList(listLength = length(vars), elementNames = vars)
  # Create a list to hold names of updated predictors.
  l_newPreds <- CreateNamedList(listLength = length(vars), elementNames = vars)
  
  # For every trait..
  for(t in 1:length(vars)){
    # Take the tth trait.
    trait <- vars[[t]]
    # Make a complete-case subset of trait data.
    dfTrait <- na.omit(data[, c("species_name", trait)])
    # Decompose phylogenetic distrance matrix into orthogonal vectors (phylogenetic eigenvectors).
    PEMtrait <- DecomposeTree(tree = tree, data = dfTrait, trait = trait, method = "PEM", evselection = "cutoff", cutoff = 65)
    # Update the PEM using derived parameters and extract the eigenvectors so we also have them for species with missing data (cite: Johnson et al. 2021).
    dfTraitEVs <- UpdatePEM(tree, alpha = PEMtrait$a, rate = PEMtrait$psi)
    # Subset to those eigenvectors selected according to the cutoff threshold.
    dfTraitEVs <- dfTraitEVs[, names(dfTraitEVs) %in% names(PEMtrait[[1]])]
    # Identify column numbers of the eigenvectors.
    index <- grep("V_", colnames(dfTraitEVs))
    # Append name of trait to eigenvector names.
    colnames(dfTraitEVs)[index] <- paste(colnames(dfTraitEVs)[index], trait, sep = "_")
    # Assign eigenvector names to a variable (without species_name).
    eigenvectors <- colnames(dfTraitEVs[, -ncol(dfTraitEVs)])
    # Identify the predictors for the trait in question.
    index <- grep(trait, names(predictors))
    # Update the names of the predictors.
    newPreds <- c(predictors[[index]], eigenvectors)
    # Append the eigenvectors to data.
    l_dfMiss[[t]] <- merge(data, dfTraitEVs, by = "species_name")
    # Append the updated predictors.
    l_newPreds[[t]] <- newPreds
  }
  # Return list of dataframes and modified predictors.
  return(list(l_dfMiss, l_newPreds))
  
}

DecomposeTree <- function(tree, method = "PVR", evselection = "cutoff", cutoff = 65, data, trait = NULL) {
  
  # Function for decomposing phylogenetic distance matrix into orthogonal vectors using the PVRdecomp (package "PVR") or PEM.fitSimple (package "MPSEM") functions.
  
  # Code adapted from:
  # PVR documentation: Santos T. 2018. https://cran.r-project.org/web/packages/PVR/PVR.pdf
  # MPSEM tutorial: G Geunard. 2019. A phylogenetic modelling tutorial using Phylogenetic Eigenvector Maps (PEM) as implemented in R package MPSEM (0.3-6).
  
  # tree = object of class "phylo". Branch lengths should be in expected # of subs/site/unit of time (ultrametric tree).
  # method = one of "PVR" 'phylogenetic eigenvector regression' or "PEM" phylogenetic eigenvector mapping
  # evselection = one of "cutoff" or "stepwise" methods for selecting eigenvector number. If cutoff is selected, a threshold variability method is used; if "stepwise" is selected, stepwise regressiod is used (function lmforwardsequentialAICc from the MPSEM package)
  # cutoff = number indicating the variability cutoff for the eigenvalues for eigenvector selection (default is 65)
  # data = dataframe containing trait data
  # trait = when PEM is specified, indicates target trait (response variable)
  
  # Make sure the trait data and tree tips match.
  tree <- drop.tip(phy = tree, tip = tree$tip.label[!tree$tip.label %in% data$species_name])
  # Make sure the order matches between the phylogeny and dataframe.
  data <- data[match(tree$tip.label, data$species_name), ]
  # Ensure data is in data.frame format.
  data <- as.data.frame(data)
  
  if(method == "PVR") {
    # Decomposing phylogenetic distance matrix into orthogonal vectors.
    matrixPhylo <- PVRdecomp(tree)
    # Extract the eigenvalues so we can determine how many eigenvectors we should be using.
    dfEigenvalues <- as.data.frame(matrixPhylo@Eigen$values) ## Place eigenvalues into dataframe.
    # How much variability does the eigenvalue explain?
    dfEigenvalues$variability <- dfEigenvalues$`matrixPhylo@Eigen$values`/sum(dfEigenvalues$`matrixPhylo@Eigen$values`)
    # Get the percentage.
    dfEigenvalues$percentage <- dfEigenvalues$variability*100
    # Get the cumulative percentage for each eigenvalues.
    dfEigenvalues$cumsum <- cumsum(dfEigenvalues$percentage)
    # Plot the cumulative percentages.
    plot(dfEigenvalues$cumsum)
    # Add a horizontal line at cutoff.
    abline(h = cutoff, col = "red")
    # How many eigenvalues explain the variability?
    evs <- which(dfEigenvalues$cumsum <= cutoff)
    EVN <- length(evs)
    # Now we need to extract the eigenvectors for each species.
    # Extract the eigenvalues so we can determine how many eigenvectors we should be using.
    dfEigenvectors <- as.data.frame(matrixPhylo@Eigen$vectors) ## Place eigenvalues into dataframe.
    # Add the species names to this dataframe so we can merge with the trait data later on.
    dfEigenvectors$species_name <- matrixPhylo@phylo$tip.label
    # Subset according to the number of eigenvalues we determined previously.
    dfEigenvectors <- dfEigenvectors[, c(1:EVN, ncol(dfEigenvectors))]
    # Return dataframe of eigenvectors.
    return(dfEigenvectors)
    
  } else if(method == "PEM") {
    # Convert phylogenetic tree to graph-class object needed to build phylogenetic eigenvector map (PEM).
    graphPhylo <- Phylo2DirectedGraph(tree)
    # Calculate phylogenetic eigenvector maps (PEMs). Specifying PEM.fitSimple allows to estimate values of a (steepness parameter) and psi (relative evolution rate along edges).
    PEM <- PEM.fitSimple(
      y = data[, trait],
      x = NULL,
      w = graphPhylo)
    # If the cutoff method for eigenvector selection was selected..
    if(evselection == "cutoff"){
      # Extract the singular values so we can determine how many eigenvectors we should be using. 
      dfEigenvalues <- as.data.frame(PEM$d) ## Place singular values into dataframe.
      # Square root to convert to obtains eigenvalues.
      dfEigenvalues$`PEM$d` <- sqrt(dfEigenvalues$`PEM$d`)
      # How much variability does the eigenvalue explain?
      dfEigenvalues$variability <- dfEigenvalues$`PEM$d`/sum(dfEigenvalues$`PEM$d`)
      # Get the percentage.
      dfEigenvalues$percentage <- dfEigenvalues$variability*100
      # Get the cumulative percentage for each eigenvalues.
      dfEigenvalues$cumsum <- cumsum(dfEigenvalues$percentage)
      # Plot the cumulative percentages.
      plot(dfEigenvalues$cumsum)
      # Add a horizontal line at cutoff.
      abline(h = cutoff, col = "red")
      # How many eigenvalues explain the variability?
      evs <- which(dfEigenvalues$cumsum <= cutoff)
      EVN <- length(evs)
      # Now we need to extract the eigenvectors for each species.
      dfEigenvectors <- as.data.frame(PEM)
      dfEigenvectors$species_name <- rownames(dfEigenvectors)
      # Subset according to the number of eigenvalues we determined previously using PVR.
      dfEigenvectors <- dfEigenvectors[, c(1:EVN, ncol(dfEigenvectors))]
      rownames(dfEigenvectors) <- NULL
      
      # If the stepwise method was selected..
    } else if(evselection == "stepwise"){
      # Apply lmforwardsequentialAICc function.
      lm1 <- lmforwardsequentialAICc(y = data[, trait], object = PEM)
      # Extract eigenvectors.
      EVN <- names(lm1$coefficients)[-1]
      # Now we need to extract the eigenvectors for each species.
      dfEigenvectors <- as.data.frame(PEM)
      dfEigenvectors$species_name <- rownames(dfEigenvectors)
      # Rearrange dataframe.
      dfEigenvectors <- dfEigenvectors[, c(EVN, "species_name")]
      # Remove row names.
      rownames(dfEigenvectors) <- NULL
    }
    # Combine dfEigenvectors into a list with parameters estimated for PEM graph (a = steepness parameter; psi = relative rate of evolution along branch lengths)
    result <- list(eigenvectors = dfEigenvectors, a = unique(PEM$a), psi = unique(PEM$psi))
    # Return result.
    return(result)
  }
  
}

DropAndMatch <- function(tree, data){
  
  # Function for matching tree tip labels to a dataframe. Also ensures the order is the same.
  # tree = phylo object with tips as species names
  # data = dataframe with column entitled "species_name"
  
  # Make sure the trait data and tree tips match.
  tree <- drop.tip(phy = tree, tip = tree$tip.label[!tree$tip.label %in% data$species_name])
  # Make sure the order matches between the phylogeny and dataframe.
  data <- data[match(tree$tip.label, data$species_name), ]
  
  # Return dataframe and tree in list format.
  return(list(phylo = tree, dataframe = data))
  
}

TestPhyloSig <- function(df, colName, phylo, sig = "lambda") {
  
  # Function that estimates phylogenetic signal of the trait given a phylogenetic tree.
  # df = dataframe containing trait information. Must also contain column called "species_name" with species name info
  # colName = column name of trait of interest
  # phylo = tree of class phylo
  # sig = one of "K" (Blomberg's K) or "lambda" (Pagel's lambda)
  
  # Subset the trait from the dataframe.
  trait <- df[, colName]
  # If the trait is numeric...
  if(is.numeric(trait)) {
    # Ensure df is data.frame format.
    df <- as.data.frame(df)
    # Ensure that the order of the data matches the tree.
    df <- df[match(phylo$tip.label, df$species_name), ]
    # Name it according to the tips of the provided tree.
    names(trait) <- phylo$tip.label
    # Estimate Blomberg's K (Blomberg, 2003) or Pagel's lambda using the phylosig function from the "phytools" package. A test of significance is also performed by setting test = TRUE.
    sig <- phylosig(phylo, trait, test = TRUE, method = sig)
    # If the trait is factor (categorical)...
  } else if(is.factor(trait)) {
    # Ensure df is data.frame format.
    df <- as.data.frame(df)
    # First prune the tree so that only tips in the dataframe are present. This is necessary for the phylo.d function to work properly.
    prunedPhylo <- drop.tip(phy = phylo, tip = phylo$tip.label[!phylo$tip.label %in% df$species_name])
    # Make sure the dataframe matches the order of the tree's tip labels.
    df <- df[match(prunedPhylo$tip.label, df$species_name), ]
    # Estimate the D metric (measure of phylogenetic signal for binary traits) by Fritz and Purvis (2010) using the phylo.d function from the "caper" package.
    # Temporarily rename column as "colName" so that phylo.d function works (as it takes colName as literal name of column).
    names(df)[names(df) == colName] <- "colName"
    sig <- phylo.d(df, prunedPhylo, names.col = species_name, binvar = colName) 
  }
  # Return the estimate of phylogenetic signal.
  return(sig)
  
}

UpdatePEM <- function(tree, alpha, rate){
  
  # Function for building phylogenetic eigenvector map (PEM) using empirically derived parameters.
  # tree = object of class "phylo". Branch lengths should be in expected # of subs/site/unit of time (ultrametric tree).
  # alpha = alpha (steepness) parameter
  # rate = relative evolution rate of the trait
  
  # Code adapted from:
  # MPSEM tutorial: G Geunard. 2019. A phylogenetic modelling tutorial using Phylogenetic Eigenvector Maps (PEM) as implemented in R package MPSEM (0.3-6).
  
  # Convert phylogenetic tree to graph-class object needed to build phylogenetic eigenvector map (PEM).
  graphPhylo <- Phylo2DirectedGraph(tree)
  # Create new PEM using derived parameter values.
  updatedPEM <- PEM.build(graphPhylo, a = alpha, psi = rate)
  # Extract the eigenvectors into a dataframe.
  dfEigenvectors <- as.data.frame(updatedPEM)
  # Create species_name column.
  dfEigenvectors$species_name <- rownames(dfEigenvectors)
  # Remove rownames.
  rownames(dfEigenvectors) <- NULL
  # Return dataframe of updated eigenvectors.
  return(dfEigenvectors)
  
}
