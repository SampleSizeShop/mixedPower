# 
#  Package mixedPower calculates power for the linear mixed model
#  Copyright (C) 2013 University of Colorado Denver.
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
# Define classes required to calculate mixed model power
#
#
# defineClasses.R
#
# Defines the following classes and related methods
#  missingDataPattern - class describing the set of responses observed for a group
#                       of independent sampling units.
#  design.mixed - class describing a mixed model study design
#  glh.mixed - class describing a mixed model hypothesis related to fixed effects
#
#

#'
#' missingDataPattern
#'
#' Class describing the pattern of observations for a group of
#' independent sampling units within a mixed model design. 
#'
#' @slot group The group (i.e. treatment assignment) for the sampling unit.
#' @slot observations The list of indices for which 
#' @slot size The number of independent sampling units with this pattern
#'    of observations
#' An optional \code{character} string specifying the name of the study design.
#' @slot designMatrix the design matrix for independent sampling units with
#'    the described pattern of observations.     
#' @seealso design.mixed
#' 
#' @name missingDataPattern 
#' @rdname missingDataPattern
#' @note For theoretical details, please see
#' 
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
#' A Power Approximation for Longitudinal Studies Using the 
#' Kenward and Roger Wald Test in the Linear Mixed Model, In review.
#'
setClass(
  "missingDataPattern",
  representation (
    group = "numeric",
    observations = "numeric",
    size = "numeric",
    designMatrix = "matrix"
  ),
  prototype (
    group = 1,
    observations = c(1),
    size = 10,
    designMatrix = diag(1)
  )  
) 

#'
#' design.mixed
#'
#' Class describing a mixed model study design for use in power analysis
#' 
#' @slot name An optional \code{character} string specifying the name of the study design.
#' @slot description An optional \code{character} string specifying the brief description of the study design.
#' @slot xPatternList list of missing data patterns for all planned participants  
#' @slot beta The \code{matrix} of regression coefficients.
#' @slot Sigma The residual covariance for a complete data case \code{matrix}.
#' 
#' @name design.mixed 
#' @rdname design.mixed
#' @note For theoretical details, please see
#' 
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
#' A Power Approximation for Longitudinal Studies Using the 
#' Kenward and Roger Wald Test in the Linear Mixed Model, In review.
#'
setClass (
  "design.mixed",
  representation ( name = "character",
                   description = "character",
                   xPatternList = "list",
                   beta = "matrix",
                   Sigma = "matrix"
  ),
  prototype ( name ="",
              description ="",
              xPatternList = c(
                new("missingDataPattern", group=1, observations=c(1,2,3), size=20,
                    designMatrix=matrix(cbind(diag(3), matrix(rep(0,9), nrow=3)))),
                new("missingDataPattern", group=1, observations=c(1,2), size=20,
                    designMatrix=matrix(cbind(diag(2), matrix(rep(0,6), nrow=2)))),
                new("missingDataPattern", group=2, observations=c(1,2,3), size=20,
                    designMatrix=matrix(cbind(matrix(rep(0,9), nrow=3), diag(3)))),
                new("missingDataPattern", group=2, observations=c(1,2), size=15,
                    designMatrix=matrix(cbind(matrix(rep(0,6), nrow=2), diag(2)))) 
              ),
              beta = matrix(c(1,1,1,0,0,0),ncol=1),
              Sigma = matrix(c(1,0.3,0.3,0.3,1,0.3,0.3,0.3,1), nrow=3)
  ),
  validity = function(object) {
    # TODO - check max observations. make sure no one has more listed
    # also make sure all observations are in range
    
    return(TRUE)
  }
)

#'
#' glh.mixed
#'
#' Class describing the general linear hypothesis for fixed effects in the mixed model
#' 
#' @slot alpha The Type I error rate
#' @slot fixedContrast The \code{matrix} of fixed effects contrasts.
#' @slot thetaNull The \code{matrix} of between participant contrasts.
#' @slot test A \code{character} string indicating the statistical test. At present,
#' only the Kenward and Roger Wald test is supported with value "Wald, KR ddf".
#' 
#' @name glh.mixed 
#' @rdname glh.mixed
#' @note For theoretical details, please see
#' 
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
#' A Power Approximation for Longitudinal Studies Using the 
#' Kenward and Roger Wald Test in the Linear Mixed Model, In review.
#'
setClass (
  "glh.mixed",
  representation ( alpha = "numeric",
                   fixedContrast = "matrix",
                   thetaNull = "matrix",
                   test = "character"
  ),
  prototype ( alpha = 0.05,
              fixedContrast = matrix(c(1/3,1/3,1/3,-1/3,-1/3,-1/3), nrow=1),
              thetaNull = matrix(c(0)),
              test = "Wald, KR ddf"
  ),
  validity = function(object) {
    # make sure thetaNull conforms with the between and within contrasts
    if (nrow(object@fixedContrast) != nrow(object@thetaNull)) {
      stop("The number of rows in the between contrast must match the number of rows of thetaNull")
    }
    if (ncol(object@thetaNull) > 1) {
      stop("theta null must be a vector")
    }
    
    return(TRUE)
  }
)


#
#
#
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
setMethod("simulateData", "design.mixed", 
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

