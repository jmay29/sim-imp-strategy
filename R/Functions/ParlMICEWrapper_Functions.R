# These are functions to call parlmice from within another function (needed due to a glitch). Credit to users JoshuaSimon and vwrobel on https://github.com/amices/mice/issues/189 (January 2022)

# Original code from MICE package: 
# mice() function in the "mice" package for imputation.
# Citations: van Buuren S, Groothuis-Oudshoorn K. MICE: Multivariate Imputation by Chained Equations in R. J Stat Softw. 2011;45(3):1–67. R package v. 3.13.0.
# https://cran.r-project.org/web/packages/mice/mice.pdf.

check.dataform <- function(data) {
  if (!(is.matrix(data) || is.data.frame(data)))
    stop("Data should be a matrix or data frame", call. = FALSE)
  if (ncol(data) < 2)
    stop("Data should contain at least two columns", call. = FALSE)
  data <- as.data.frame(data)
  mat <- sapply(data, is.matrix)
  if (any(mat)) stop("Cannot handle columns with class matrix: ", 
                     colnames(data)[mat])
  
  dup <- duplicated(colnames(data))
  if (any(dup)) stop("Duplicate names found: ", 
                     paste(colnames(data)[dup], collapse = ", "))
  
  return(data)
}

parlmiceMOD <- function (data, m = 5, seed = NA, cluster.seed = NA, n.core = NULL, 
                         n.imp.core = NULL, maxit = 10, predictorMatrix, cl.type = "PSOCK", ...) 
  
  # JM EDIT: I added maxit and predictorMatrix as arguments.
  
{
  data <- check.dataform(data)
  m <- check.m(m)
  if (sum(is.na(data)) == 0) {
    stop("Data has no missing values")
  }
  if (!is.null(n.core)) {
    if (n.core > parallel::detectCores()) {
      stop("Number of cores specified is greater than the number of logical cores in your CPU")
    }
  }
  if (!is.null(n.core) & is.null(n.imp.core)) {
    n.imp.core <- m
    warning(paste("Number of imputations per core not specified: n.imp.core = m =", 
                  m, "has been used"))
  }
  if (is.null(n.core) & !is.null(n.imp.core)) {
    n.core <- parallel::detectCores() - 1
    warning(paste("Number of cores not specified. Based on your machine a value of n.core =", 
                  parallel::detectCores() - 1, "is chosen"))
  }
  if (is.null(n.core) & is.null(n.imp.core)) {
    specs <- match.cluster(n.core = parallel::detectCores() - 
                             1, m = m)
    n.core <- specs$cores
    n.imp.core <- specs$imps
  }
  if (!is.na(seed)) {
    if (n.core > 1) {
      warning("Be careful; the specified seed is equal for all imputations. Please consider specifying cluster.seed instead.")
    }
  }
  args <- match.call(mice, expand.dots = TRUE)
  args[[1]] <- NULL
  args$m <- n.imp.core
  cl <- parallel::makeCluster(n.core, type = cl.type)
  ## JM EDIT: I added maxit and predictorMatrix to varlist.
  parallel::clusterExport(cl, varlist = c("data", "m", 
                                          "seed", "cluster.seed", "n.core", "n.imp.core", 
                                          "cl.type", "maxit", "predictorMatrix", ls(parent.frame())), envir = environment())
  parallel::clusterExport(cl, varlist = "do.call")
  parallel::clusterEvalQ(cl, library(mice))
  if (!is.na(cluster.seed)) {
    parallel::clusterSetRNGStream(cl, cluster.seed)
  }
  imps <- parallel::parLapply(cl = cl, X = 1:n.core, function(x) do.call(mice, 
                                                                         as.list(args), envir = environment()))
  parallel::stopCluster(cl)
  imp <- imps[[1]]
  if (length(imps) > 1) {
    for (i in 2:length(imps)) {
      imp <- ibind(imp, imps[[i]])
    }
  }
  for (i in 1:length(imp$imp)) {
    colnames(imp$imp[[i]]) <- 1:imp$m
  }
  imp
}


ParlMiceWrapper <- function(data, cluster.seed, n.core, n.imp.core, maxit, predictorMatrix, ...) {
  
  # This is a function to call parlmice from within another function. Credit to user JoshuaSimon on https://github.com/amices/mice/issues/189 (January 2022).

  imputedMICE <- parlmiceMOD(data = data, cluster.seed = cluster.seed,
                          n.core = n.core, n.imp.core = n.imp.core, 
                          maxit = maxit, predictorMatrix = predictorMatrix)
  return(imputedMICE)
}

check.data <- function(data, method) {
  check.dataform(data)
}


check.dataform <- function(data) {
  if (!(is.matrix(data) || is.data.frame(data))) {
    stop("Data should be a matrix or data frame", call. = FALSE)
  }
  if (ncol(data) < 2) {
    stop("Data should contain at least two columns", call. = FALSE)
  }
  data <- as.data.frame(data)
  mat <- sapply(data, is.matrix)
  df <- sapply(data, is.data.frame)
  if (any(mat)) {
    stop(
      "Cannot handle columns with class matrix: ",
      colnames(data)[mat]
    )
  }
  if (any(df)) {
    stop(
      "Cannot handle columns with class data.frame: ",
      colnames(data)[df]
    )
  }
  
  dup <- duplicated(colnames(data))
  if (any(dup)) {
    stop(
      "Duplicate names found: ",
      paste(colnames(data)[dup], collapse = ", ")
    )
  }
  data
}


check.m <- function(m) {
  m <- m[1L]
  if (!is.numeric(m)) {
    stop("Argument m not numeric", call. = FALSE)
  }
  m <- floor(m)
  if (m < 1L) {
    stop("Number of imputations (m) lower than 1.", call. = FALSE)
  }
  m
}


check.cluster <- function(data, predictorMatrix) {
  # stop if the cluster variable is a factor
  isclassvar <- apply(predictorMatrix == -2, 2, any)
  for (j in colnames(predictorMatrix)) {
    if (isclassvar[j] && lapply(data, is.factor)[[j]]) {
      stop("Convert cluster variable ", j, " to integer by as.integer()")
    }
  }
  TRUE
}

check.ignore <- function(ignore, data) {
  if (is.null(ignore)) {
    return(rep(FALSE, nrow(data)))
  }
  if (!is.logical(ignore)) {
    stop("Argument ignore not a logical.")
  }
  if (length(ignore) != nrow(data)) {
    stop(
      "length(ignore) (", length(ignore),
      ") does not match nrow(data) (", nrow(data), ")."
    )
  }
  if (sum(!ignore) < 10L) {
    warning(
      "Fewer than 10 rows for fitting the imputation model. Are you sure?",
      call. = FALSE
    )
  }
  ignore
}

check.newdata <- function(newdata, data) {
  if (is.null(newdata)) {
    stop("No newdata found.")
  }
  if (!is.data.frame(newdata)) {
    stop("newdata not a data.frame.")
  }
  newdata
}
