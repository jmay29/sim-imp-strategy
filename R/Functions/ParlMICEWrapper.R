ParlMiceWrapper <- function(data, cluster.seed, n.core, n.imp.core, maxit, predictorMatrix, ...) {
  
  # This is a function to call parlmice from within another function. Credit to user JoshuaSimon on https://github.com/amices/mice/issues/189 (January 2022).
  # predictorMatrix <- predictorMatrix
  # maxit <- maxit
  imputedMICE <- parlmiceMOD(data = data, cluster.seed = cluster.seed,
                          n.core = n.core, n.imp.core = n.imp.core, 
                          maxit = maxit, predictorMatrix = predictorMatrix)
  return(imputedMICE)
}
