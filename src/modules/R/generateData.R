#
# R interface to the SAS mixed model simulator
# 
# Author: Sarah Kreidler
# Created: 8/14/2013
#
#

library(MASS)
library(plyr)

#
# generate a single mixed data set 
#
#
generateMixedDataSetGivenSigma <- function(n, mu, Sigma) {
  tmpData = mvrnorm(n, mu, Sigma)
  
}

# generateMixedGivenSigma
# 
# Generate data sets with the specified covariance matrix.
# Assumes no missing data, although fixed patterns of imbalance
# may be specified in the covariance matrix
# 
# Arguments:
#   replicates - the number of replicates to produce for simulation
#   blockSize - number of replicates to include in each data set written to disk
#   X - the design matrix for fixed effects
#   Beta - matrix containing choices for regression coefficients for fixed effects
#   SigmaS - the stacked covariance (block diagonal with one block
#             per independent sampling unit)
#
generateMixedGivenSigma <- function(replicates=10000, blockSize=1000, 
                                    datasetPrefix="genMixedData", mu, Sigma, X) {
  
  ### determine the number of sets and the size of each set ###
  
  # when the blocksize does not divide evenly into the
  # total replicates, the size of the last set is the remainder
  lastSetSize = replicates %% blockSize;
  numSets = floor(replicates / blockSize)
  
  tmp = sapply(1:numSets, function(setNumber, X) {
    # generate the data in multivariate format
    data = data.frame(mvrnorm(n = blockSize, rep(0,3), diag(3)))
    startIter = ((setNumber-1)*blockSize) + 1
    endIter = startIter + blockSize - 1;
    data$setID = startIter:endIter
    
    # create an empty aggregation block
    dataBlock = data.frame(setID=c('A'),y=(0),X)
    dataBlock <- dataBlock[which(is.na(dataBlock$y))]
                      
    # split by row, transpose

                      
    ddply(data, 1, function(row, X, dataBlock){
      tmp = data.frame(setID=c(row$setID), y=row, X)
      dataBlock = rbind(dataBlock, tmp)
    }, X, dataBlock)
    
  }, X)
  
  for(i in 1:numSets) {
    startIter = ((i-1)*blockSize) + 1
    endIter = startIter + blockSize - 1;
    
    data = mvrnorm(n = blockSize, rep(0,3), diag(3))
    d_ply
  }
  
  if (lastSetSize > 0) {
    # build last data set
    startIter = ((i-1)*blockSize) + 1
    endIter = replicates;
    
    
  }
  
  
  d_ply(mvrnorm(n=10,rep(0,3),diag(3)), 1, function(row){
    
  })
  
}


dd<-data.frame(matrix(rnorm(216),72,3),c(rep("A",24),rep("B",24),rep("C",24)),c(rep("J",36),rep("K",36)))
colnames(dd) <- c("v1", "v2", "v3", "dim1", "dim2")

#ddply is the plyr function
ddply(dd, 1, function(df)df$dim1)


