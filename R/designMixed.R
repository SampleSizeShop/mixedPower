# 
#  Package glmmPower calculates power for the general linear 
#  multivariate model
#  Copyright (C) 2013 Sarah Kreidler.
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
# 
# Provides the design.mixed class which defines matrices
# for the general linear mixed model with fixed predictors.
#
# For notation and theoretical details, see
# 
# 1. Muller KE, Stewart PW. Linear model theory: univariate, multivariate, and mixed models. 
# Hoboken, New Jersey: John Wiley and Sons; 2006.
#
#

#
# design.mixed
#
# Class describing the general linear model with fixed predictors and
# one or more Gaussian covariates
#
setClass (
  "design.mixed",
  representation ( name = "character",
                   description = "character",
                   XEssence = "matrix",
                   perGroupN = "numeric",
                   Beta = "matrix",
                   SigmaISU = "matrix",
                   observed = "numeric"
  ),
  prototype ( name ="",
              description ="",
              XEssence = diag(4),
              perGroupN = 10,
              Beta = matrix(c(1,0,0,0),nrow=4),
              SigmaISU = matrix(c(1,0.2,0.2,1)),
              observed = rep(1,40)
  ),
  validity = function(object) {
    #
    # Note, the class definition will already enforce that the
    # matrices are non-null 
    #
    
    # make sure that Sigma ISU is square, symmetric and positive definite
    if (!isSymmetric(object@SigmaISU)) {
      stop("SigmaISU matrix is not symmetric")
    } else {
      sigmaISUEigenValues = eigen(object@SigmaISU)
      if (sum(sigmaISUEigenValues$values > 0) < length(sigmaISUEigenValues$values)) {
        stop("SigmaISU matrix is not positive definite")
      }
    }
    
    # check matrix conformance
    if (ncol(object@XEssence) != nrow(object@Beta)) {
      stop("The number of columns in the essence X matrix does not match the 
           number of columns in the Beta matrix")
    } else if (ncol(object@XEssence) %% nrow(object@SigmaISU) != 0) {
      stop("The columns of XEssence must be divisble by the rows of SigmaISU")    
    } else if (length(object@observed) != nrow(object@XEssence) * object@perGroupN) {
      stop("The length of the observed array must equal the total number of observations (rows XEssence * perGroupN)")    
    }
    
    return(TRUE)
  }
)

test = new("design.mixed", name="Repeated measures with missing data",
          description="Mixed model with 2 treatment groups, 
           4 planned observations, missing last observation in 5 ISUs in each group", 
           XEssence=diag(8), 
           perGroupN=10,
           Beta=matrix(c(1,1,1,1,0,0,0,0), nrow=8),
           SigmaISU=matrix(c(1, 0.2, 0.1, 0.05, 0.2, 1, 0.2, 0.1, 0.1, 0.2, 1, 0.2, 0.05, 0.1, 0.2, 1),
                           nrow=4,byrow=TRUE),
           observed=c(rep(c(1,1,1,1),5), rep(c(1,1,1,0),5),rep(c(1,1,1,1),5), rep(c(1,1,1,0),5))
)

##
# simulateData: simulate data sets for the specified study design
#
# Args:
#  design (required) - the glmmFG study design object
#  replicates (optional) - the total number of data sets
#  blockSize (optional) - the number of data sets to include in each file
#  outputDir (required) - directory in which to write the data sets
#  filePrefix (optional) - prefix added to filenames
#  xNames (optional) - column names for the predictors
#  yNames (optional) - column names for the outcomes
#
# Outputs:
#  writes simulated data to disk in CSV format.  
#  Filnames follow the pattern <filePrefix><start>To<end>.csv 
#  (if filePrefix unspecified, then simulatedData<start>To<end>.csv)
#  Files have the format
# 
#  dataSetID | Y (outcomes) | X (predictors)
#
#
setGeneric("simulateData", function(design, replicates=10000, blockSize=1000,
                                    outputDir=".", filePrefix="simulatedData",
                                    xNames=NA, yNames=NA, ...) standardGeneric("simulateData"))
setMethod("simulateData", "design.glmmFG", 
          function(design, replicates=1000, blockSize=100,
                   outputDir=".", filePrefix="simulatedData",
                   xNames=NA, yNames=NA, realizations=1000) {
            if (is.na(design@perGroupN)) {
              stop("Per group sample size not specified in study design")
            }
            if (!is.na(xNames) && length(xNames) != (ncol(design@XEssence) + ncol(design@SigmaG))) {
              stop("The number of values in xNames does not match the number of columns in XEssence")
            }
            if (!is.na(yNames) && length(yNames) != ncol(design@Beta)) {
              stop("The number of values in yNames does not match the number of columns in Beta")
            }
            
            ### calculate the mean, XB for fixed predictors ###
            # first, calculate the full X matrix
            X = matrix(rep(1,design@perGroupN), nrow=design@perGroupN) %x% design@XEssence
            # calculate XB
            XB = X %*% design@Beta
            # total sample size
            totalN = nrow(X) 
            
            ### determine the number of sets and the size of each set ###
            
            # when the blocksize does not divide evenly into the
            # total replicates, the size of the last set is the remainder
            lastSetSize = replicates %% blockSize;
            numSets = floor(replicates / blockSize)
            if (lastSetSize > 0) {
              numSets = numSets + 1
            }
            
            ### generate column names for the data sets ###
            yPredef = sapply(1:ncol(XB), function(x) {paste("Y.",x,sep="",collapse="")})
            xPredef = c(sapply(1:ncol(X), function(x) {paste("XF.",x,sep="",collapse="")}),
                        sapply(1:ncol(design@SigmaG), function(x) {paste("XG.",x,sep="",collapse="")}))
            if (!is.na(yNames) && !is.na(xNames)) {
              dataSetNames = c("realizationID", "setID", yNames, xNames)
            } else if (!is.na(yNames) && is.na(xNames)) {
              dataSetNames = c("realizationID", "setID", yNames, xPredef)
            } else if (is.na(yNames) && !is.na(xNames)) {
              dataSetNames = c("realizationID", "setID", yPredef, xNames)
            } else {
              dataSetNames = c("realizationID", "setID", yPredef, xPredef)
            }
            
            #
            # Calculate the covariance of errors per Glueck and Muller
            #  
            SigmaError = design@SigmaY - design@SigmaYG %*% solve(design@SigmaG) %*% t(design@SigmaYG)
            
            ### generate the realizations ###
            sapply(1:realizations, function(realizationID) {
              
              XG = mvrnorm(n = totalN, 
                           mu=rep(0, nrow(design@SigmaG)), 
                           design@SigmaG)
              ### buil the data set blocks ###
              sapply(1:numSets, function(setNumber) {
                # set start and end iteration numbers for this block
                startIter = ((setNumber-1)*blockSize) + 1
                if (setNumber == numSets && lastSetSize > 0 && lastSetSize != blockSize) {
                  numDataSets = lastSetSize
                } else {
                  numDataSets = blockSize
                }
                endIter = startIter + numDataSets - 1;
                cat("Generating data sets ", startIter, " to ", endIter, "\n")
                
                # generate numDataSets number of data sets
                dataSet = do.call("rbind", 
                                  lapply(1:numDataSets, function(blockNumber) {
                                    errorMatrix = 
                                      mvrnorm(n = totalN, 
                                              mu=rep(0, nrow(SigmaError)), 
                                              SigmaError)
                                    yData = XB + errorMatrix
                                    dataBlock = data.frame(realizationID=realizationID,
                                                           setID=((setNumber-1)*blockSize+blockNumber),
                                                           Y=yData, X=X, XG=XG)
                                    return(dataBlock)
                                    
                                  }))
                # write to the output directory
                names(dataSet) = dataSetNames
                filename = paste(outputDir,"/", filePrefix,"Realization", realizationID, "Iter", 
                                 startIter,"to",endIter,".csv",sep="")
                write.csv(dataSet, row.names=FALSE, file=filename)
                return(filename)
              }) # end sapply
            })
          }
)
