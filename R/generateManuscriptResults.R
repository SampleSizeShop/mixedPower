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

#' convertToMultivariateDesign
#' 
#' Create a multivariate design (design.glmmF object) equivalent to the specified
#' mixed model design if all sampling units had complete data.
#' Please note that this is NOT a general purpose function,
#' but makes assumptions specific to the validation study in the
#' manuscript
#' 
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
#' A Power Approximation for Longitudinal Studies Using the 
#' Kenward and Roger Wald Test in the Linear Mixed Model, In review.
#' 
convertToMultivariateDesign <- function(mixedDesign) {
  
  # determine the number of unique treatment groups
  numGroups = length(unique(sapply(mixedDesign@xPatternList, function(x) {
    return(x@group)
  })))
  
  # get the number of sampling units per group
  perGroupN = sum(sapply(mixedDesign@xPatternList, function(pattern) {
    if (pattern@group == 1) {
      return(pattern@size)
    } else {
      return(0)
    }
  }))
  
  Beta = matrix(rep(0,numGroups*5), nrow=numGroups)
  Beta[1,1] = mixedDesign@beta[1,1]
  return(new("design.glmmF", name=mixedDesign@name,
             description=mixedDesign@description,
             XEssence=diag(numGroups),
             perGroupN=perGroupN,
             Beta=Beta,
             SigmaError=mixedDesign@Sigma))
}

#' convertToMultivariateGLH
#' 
#' Create a multivariate linear hypothesis (glh object) equivalent to the specified
#' mixed model hypothesis if all sampling units had complete data.
#' Please note that this is NOT a general purpose function,
#' but makes assumptions specific to the validation study in the
#' manuscript
#' 
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
#' A Power Approximation for Longitudinal Studies Using the 
#' Kenward and Roger Wald Test in the Linear Mixed Model, In review.
#' 
convertToMultivariateGLH <- function(mixedGLH) {
  # assumes a time x treatment interaction with 5 repeated measures
  numGroups = ncol(mixedGLH@fixedContrast) / 5
  if (numGroups == 2) {
    betweenContrast = matrix(c(1,-1),nrow=1)
  } else {
    # 4 groups
    betweenContrast = cbind(matrix(rep(1,3),nrow=3), -1*diag(3))
  }
  withinContrast = t(cbind(matrix(rep(1,4),nrow=4), -1*diag(4)))
  
  thetaNull.rows = nrow(betweenContrast)
  thetaNull.columns = ncol(withinContrast)

  return(new("glh", alpha=mixedGLH@alpha,
             betweenContrast=betweenContrast, withinContrast=withinContrast,
             thetaNull=matrix(rep(0,thetaNull.rows*thetaNull.columns), nrow=thetaNull.rows)))
}

#' heterogeneousCS
#' 
#' generate a heterogeneous compound symmetric covariance matrix.
#' The size of the matrix is determined by the number of
#' variance values in the sigmaList. Note that
#' when all sigma values are equal, a compound symmetric 
#' covariance is produced.  
#'  
#' @param sigmaList list of variances
#' @param rho correlation parameter
#' @return covariance \code{matrix} with heterogeneous compound symmetric
#' structure
#' 
heterogeneousCSMatrix <- function(sigmaList,rho) {
  #sigmaList = c(1,rep(0.1,4))
  size = length(sigmaList)
  Sigma = matrix(rep(0,size*size),nrow=size)
  for(r in 1:nrow(Sigma)) {
    for(c in 1:ncol(Sigma)) {
      if (r==c) {
        Sigma[r,c] = sigmaList[r]*sigmaList[c]
      } else {
        Sigma[r,c] = sigmaList[r]*sigmaList[c] * rho
      }
    }
  }
  return(Sigma)
}


#' learMatrix
#' 
#' generate a linear exponent auto-regressive correlation matrix.
#' Assumes equal spacing of measurements.
#'   
#' @param size row/column dimension of the correlation matrix
#' @param rho max correlation
#' @param delta rate of correlation decay
#' @return correlation \code{matrix} with LEAR structure
#'
#' @note Implements the methods of 
#' Simpson, S. L., Edwards, L. J., Muller, K. E., Sen, P. K., & Styner, M. A. (2010). 
#' A linear exponent AR(1) family of correlation structures. Statistics in Medicine, 
#' 29(17), 1825–1838. doi:10.1002/sim.3928
#' 
learMatrix = function(size, rho, delta) {
  dmin=1
  dmax=size-1
  lear = diag(size)
  for(r in 1:(size-1)) {
    for(c in (r+1):size) {
      value = rho^(dmin + delta*((c-r-dmin)/(dmax-dmin)))
      lear[r,c] = value
      lear[c,r] = value
    }
  }
  return(lear)
}

#' ar1Matrix
#' 
#' generate an auto-regressive correlation matrix.
#'   
#' @param size row/column dimension of the correlation matrix
#' @param rho correlation between measurements 1 unit apart
#' @return correlation \code{matrix} with auto-regressive structure
#' 
ar1Matrix = function(size, rho) {
  dmin=1
  dmax=size-1
  ar1 = diag(size)
  for(r in 1:(size-1)) {
    for(c in (r+1):size) {
      value = rho^(c-r)
      ar1[r,c] = value
      ar1[c,r] = value
    }
  }
  return(ar1)
}

#' 
#' getBetaScaleByPower
#' 
#' Identify the beta scale such that mixedPower returns the 
#' requested target power.
#'   
#' @param design a \code{design.mixed} object
#' @param glh a \code{glh.mixed} object
#' @param targePower desired power value
#' @param lower lower bound for beta scale search
#' @param upper upper bound for beta scale search
#'
getBetaScaleByPower <- function(design, glh, targetPower=0.90, 
                                lower=0.0001, upper=1000) {
  betaScale = uniroot(function(x) {
    design@beta = design@beta * x
    return(mixedPower(design, glh) - targetPower)
  }, c(lower, upper))
  return(betaScale$root)
}

#' 
#' getGlh
#' 
#' Create a time by treatment interaction hypothesis object for 
#' a design with the specified number of groups and repeated
#' measures
#'   
#' @param numGroups number of between participant groups
#' @param maxObs number of repeated measures for a complete data case
#'
getGlh = function(numGroups, maxObs) {
  betweenContrast = cbind(matrix(rep(1,numGroups-1)), -1*diag(numGroups-1))
  withinContrast = cbind(matrix(rep(1,maxObs-1)), -1*diag(maxObs-1))
  fixedContrast = betweenContrast %x% withinContrast
  return (new("glh.mixed",
              alpha = 0.05,
              fixedContrast = fixedContrast,
              thetaNull = matrix(rep(0,nrow(fixedContrast)), nrow=nrow(fixedContrast)),
              test = "Wald, KR ddf"))
}

#' 
#' generateLongitudinalDesign
#' 
#' Create a logitudinal design.mixed object and a corresponding
#' glh.mixed object which tests the time by treatment interaction
#'   
#' @param params input parameters describing the designs including
#' number of groups, per group N, percent of ISUs with missing data,
#' missing pattern type, covariance structure, and target power.
#'
generateLongitudinalDesign = function(params) {
  
  name = "longitudinal"
  description = paste(c(params$numGroups, 
                        " group longitudinal design, test of time by treatment interaction. ",
                        "Per group N = ", params$perGroupN, 
                        ", percent missing = ", params$missingPercent,
                        ", missing pattern type = ", params$missingType,
                        ", covariance = ", params$covariance,
                        ", target power = ", params$targetPower), collapse="")
  # number of observations for a complete data case
  maxObservations = params$maxObservations
  
  # build the pattern list
  patternList = list()
  pattern = 1
  for(grp in 1:params$numGroups) {
    # between effects design matrix
    designBetween = matrix(diag(params$numGroups)[grp,], nrow=1)
    
    # calculate the number of ISUs with missing data
    numComplete = floor(params$perGroupN * (1-params$missingPercent))
    numIncomplete = params$perGroupN - numComplete
    
    # build missing data patterns
    if (numIncomplete == 0) {
      # complete case
      patternList[[pattern]] = 
        new("missingDataPattern", group=grp, observations=1:maxObservations, 
            designMatrix = (designBetween %x% diag(maxObservations)),
            size=numComplete)
    } else if (params$missingType == 'monotone') {
      ## monotone dropout pattern - once missing, never come back
      # complete case
      patternList[[pattern]] = 
        new("missingDataPattern", group=grp, observations=1:maxObservations, 
            designMatrix = (designBetween %x% diag(maxObservations)),
            size=numComplete)
      pattern = pattern + 1
      # missing last 2 observations
      missingPattern2 = 1:(maxObservations-2)
      patternList[[pattern]] = 
        new("missingDataPattern", group=grp, observations=missingPattern2, 
            designMatrix = (designBetween %x% 
                              matrix(diag(maxObservations)[missingPattern2,], 
                                     nrow=length(missingPattern2))),
            size=numIncomplete)
      
    } else {
      ## Non-monotone dropout pattern 
      ## - delete observations 2 and 4
      # complete case
      patternList[[pattern]] = 
        new("missingDataPattern", group=grp, observations=1:maxObservations, 
            designMatrix = (designBetween %x% diag(maxObservations)),
            size=numComplete)
      pattern = pattern + 1
      
      # missing 2nd (and 4th) observation(s)
      incompleteObs = c(1,3,5)
      patternList[[pattern]] = 
        new("missingDataPattern", group=grp, observations=incompleteObs, 
            designMatrix = (designBetween %x% 
                              matrix(diag(maxObservations)[incompleteObs,], 
                                     nrow=length(incompleteObs))),
            size=numIncomplete)
    }
    pattern = pattern + 1
  }
  
  # build sigma
  rho = 0.04
  sigmaSq = 1
  if (params$covariance == 'CS') {
    Sigma = heterogeneousCSMatrix(rep(1,maxObservations),rho)
  } else if (params$covariance == 'CSH') {
    Sigma = heterogeneousCSMatrix(c(1,0.5,0.3,0.1,0.1),rho)
  } else {
    Sigma = sigmaSq * 
      ar1Matrix(maxObservations,rho)
  }
  
  # build the design
  beta = matrix(c(1,rep(0,(maxObservations-1)),
                  rep(0,maxObservations*(params$numGroups-1))))
  design = new("design.mixed", name = name, description = description,
               xPatternList = patternList,
               beta = beta,
               Sigma = Sigma
  )
  # get the appropriate hypothesis
  glh = getGlh(params$numGroups, params$maxObservations)
  # update the beta scale
  betaScale = getBetaScaleByPower(design, glh, targetPower=params$targetPower)  
  design@beta = beta * betaScale
  return(design)
  
}


#' 
#' generateDesignsForManuscript
#' 
#' Generate the longitudinal study designs used in the validation
#' experiment in the manuscript
#' 
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
#' A Power Approximation for Longitudinal Studies Using the 
#' Kenward and Roger Wald Test in the Linear Mixed Model, In review.
#'   
#' @param output.data.dir directory to which design/hypothesis information
#' and design input parameter files are written
#'
generateDesignsForManuscript = function(output.data.dir=".") {
  #
  # For each longitudinal randomized design, we
  # vary the following parameters
  #
  # number of treatment groups
  numGroupsList = c(2, 4)
  # total ISUs per treatment group
  perGroupNList = c(50, 100)
  # missing data pattern (either monotone or non-monotone)
  missingTypeList = c("monotone", "non-monotone")
  # percent missing
  missingPercentList = c(0, 0.2, 0.4)
  # covariance
  covarianceList = c("CS", "CSH", "AR(1)")
  # in all cases, we select the scale factor 
  # for beta to achieve the following power
  targetPowerList = c(0.2, 0.5, 0.8)
  # maximum number of observations 
  maxObservationsList = c(5)
  
  # generate parameters
  paramList = list(targetPower=targetPowerList,
                   covariance=covarianceList,
                   missingType=missingTypeList,
                   missingPercent=missingPercentList,
                   perGroupN=perGroupNList,
                   numGroups=numGroupsList,
                   maxObservations=maxObservationsList)
  paramComboList = data.frame(expand.grid(paramList))
  
  #
  # Calculate the appropriate betaScale values
  # and build the list of designs
  #
  longitudinalDesignList = list()
  betaScaleList = vector()
  for(i in 1:length(paramComboList$numGroups)) {
    cat("Case ", i, "\n")
    params = paramComboList[i,]
    longitudinalDesignList[[i]] = list(generateLongitudinalDesign(params), 
                                       getGlh(params$numGroups, params$maxObservations))                                     
    betaScaleList[i] = longitudinalDesignList[[i]][[1]]@beta[1,1]
  }
  paramComboList$betaScale = betaScaleList
  
  # write the parameter data to a csv file and an RData file
  write.csv(paramComboList, 
            file=paste(c(output.data.dir,"longitudinalParams.csv"),collapse="/"),
            row.names=FALSE)
  save(paramComboList, 
       file=paste(c(output.data.dir,"longitudinalParams.RData"),collapse="/"))
  
  # write the designs to an Rdata file
  save(longitudinalDesignList, 
       file=paste(c(output.data.dir,"longitudinalDesigns.RData"),collapse="/"))
  
  return(longitudinalDesignList)
}


#' 
#' calculateEmpiricalPowerWithSAS
#' 
#' Calls a SAS program included in the mixedPower package to
#' generate the empirical power results for the validation experiment
#' in the manuscript
#' 
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
#' A Power Approximation for Longitudinal Studies Using the 
#' Kenward and Roger Wald Test in the Linear Mixed Model, In review.
#'   
#' @param output.data.dir directory to which empirical power data files are written
#'
calculateEmpiricalPowerWithSAS <- function(output.data.dir) {
  # generate the common.sas file
  
  # call the sas code to run empirical power
  result <- system(paste(c("sas.exe -i ", 
                           paste(c(path.package("mixedPower"), 
                                   "inst/sas/calculateEmpiricalPower.sas"), 
                                 collapse="/")), collapse=""), 
                   intern = TRUE, show.output.on.console = TRUE)
  # remove common.sas file
  
  
  # load the empirical power results
  empiricalPowerData = read.csv(
    paste(c(output.data.dir, "longitudinalEmpiricalPower.csv"), collapse="/"),
    stringsAsFactors=FALSE)

  # save as a Rdata file
  save(empiricalPowerData, file=paste(c(output.data.dir, "longitudinalEmpiricalPower.RData"), collapse="/"))
  
  return(empiricalPowerData)
}

#' 
#' calculateExemplaryPowerWithSAS
#' 
#' Calls a SAS program included in the mixedPower package to
#' generate the exemplary power results for the validation experiment
#' in the manuscript
#' 
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
#' A Power Approximation for Longitudinal Studies Using the 
#' Kenward and Roger Wald Test in the Linear Mixed Model, In review.
#'   
#' @param output.data.dir directory to which exemplary power data files are written
#'
calculateExemplaryPowerWithSAS <- function(output.data.dir) {
  # generate the common.sas file
  
  # call the sas code to run empirical power
  result <- system(paste(c("sas.exe -i ", 
                           paste(c(path.package("mixedPower"), 
                                   "inst/sas/calculateStroupPower.sas"), 
                                 collapse="/")), collapse=""), 
                   intern = TRUE, show.output.on.console = TRUE)
  # remove common.sas file
  
  
  # load the empirical power results
  exemplaryPowerData = read.csv(
    paste(c(output.data.dir, "longitudinalExemplaryPower.csv"), collapse="/"),
    stringsAsFactors=FALSE)
  
  # save as a Rdata file
  save(exemplaryPowerData, file=paste(c(output.data.dir, "longitudinalExemplaryPower.RData"), collapse="/"))
  
  return(exemplaryPowerData)
}

#' calculateApproximatePowerForDesignList
#' 
#' Calculate approximate power values for a list of design/glh pairs.
#'   
#' @param designList list of pairs of design.glmmFG and glh objects
#' @param output.data.dir directory to which data sets are written
#' @return data frame containing approximate power results
#' @note This function takes about 15-30 minutes to run
#' 
calculateApproximatePowerForDesignList <- function(designList, output.data.dir=".") {

  # add power using the method described by 
  #
  # Helms, R. W. (1992). Intentionally incomplete longitudinal designs: 
  # I. Methodology and comparison of some full span designs. 
  # Statistics in Medicine, 11(14-15), 1889–1913.
  #
  helmsPowerList = sapply(designList, function(designAndGlh) {
    return(mixedPower.helms(designAndGlh[[1]], designAndGlh[[2]]))
  })
  
  # add power assuming complete, balanced data. This method is described by 
  #
  # Muller, K. E., Lavange, L. M., Ramey, S. L., & Ramey, C. T. (1992). 
  # Power Calculations for General Linear Multivariate Models Including 
  # Repeated Measures Applications. Journal of the American Statistical 
  # Association, 87(420), 1209–1226.
  #
  multivariatePowerList=sapply(designList, function(x) {
    multivariateDesign = convertToMultivariateDesign(x[[1]])
    multivariateHypothesis = convertToMultivariateGLH(x[[2]])
    return(glmmPower.fixed(multivariateDesign, multivariateHypothesis)) 
  })
  
  # add power using the method described by
  #
  # Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
  # A Power Approximation for Longitudinal Studies Using the 
  # Kenward and Roger Wald Test in the Linear Mixed Model, In review.
  #
  kreidlerPowerList = sapply(designList, function(designAndGlh) {
    return(mixedPower(designAndGlh[[1]], designAndGlh[[2]]))
  })
  
  approxPowerData = data.frame(helmsPower=helmsPowerList,
                               multivariatePower=multivariatePowerList,
                               kreidlerPower=kreidlerPowerList)
  ## write the approximate power data to disk     
  save(approxPowerData,
       file=paste(c(output.data.dir, "longitudinalApproximatePower.RData"), collapse="/"))
  
  return(approxPowerData)
}

#' summarizeResults
#' 
#' Generate box plots which summarize the deviations from empirical
#' power for the approximate methods.
#'   
#' @param output.data.dir directory to which data sets are written
#' @return output.figures.dir directory to which figures are written
#'
summarizeResults = function(output.data.dir=".", output.figures.dir=".") {
  
  # load the data - creates a variable 'approximateAndEmpiricalPowerData'
  load(paste(c(output.data.dir, "approximateAndEmpiricalPower.RData"), collapse="/"))
  
  # calculate the deviations
  approximateAndEmpiricalPowerData$id = 1:length(approximateAndEmpiricalPowerData$targetPower)
  approximateAndEmpiricalPowerData$diff.kreidler = approximateAndEmpiricalPowerData$kreidlerPower - 
    approximateAndEmpiricalPowerData$empiricalPower 
  approximateAndEmpiricalPowerData$diff.helms = approximateAndEmpiricalPowerData$helmsPower - 
    approximateAndEmpiricalPowerData$empiricalPower 
  approximateAndEmpiricalPowerData$diff.multivariate = approximateAndEmpiricalPowerData$multivariatePower - 
    approximateAndEmpiricalPowerData$empiricalPower
  approximateAndEmpiricalPowerData$diff.stroup = approximateAndEmpiricalPowerData$exemplaryPower - 
    approximateAndEmpiricalPowerData$empiricalPower
  
  # get max absolute deviations
  maxDeviationData = data.frame(method=c("Kreidler", "Helms", "Multivariate", "Stroup"),
                                maxDeviation=c(max(abs(approximateAndEmpiricalPowerData$diff.kreidler)),
                                               max(abs(approximateAndEmpiricalPowerData$diff.helms)),
                                               max(abs(approximateAndEmpiricalPowerData$diff.multivariate)),
                                               max(abs(approximateAndEmpiricalPowerData$diff.stroup))))
  # save the max deviation info to a file
  save(maxDeviationData,
       file=paste(c(output.data.dir, "maxDeviationsByMethod.RData"), collapse="/"))
  
  # convert to long with factor identifying power method
  powerDataLong = reshape(approximateAndEmpiricalPowerData, 
                          varying=c(
                            "diff.kreidler",
                            "diff.helms",
                            "diff.multivariate",
                            "diff.stroup"
                          ),
                          timevar="method",
                          times = c("Kreidler", "Helms", "Multivariate", "Stroup"),
                          idvar="id", direction="long")
  # make 'method' into a factor so we can sort the boxplots
  powerDataLong$method = factor(powerDataLong$method, 
                                levels=c("Kreidler", "Helms", "Multivariate", "Stroup"),
                                labels=c("Kreidler", "Helms", "Multivariate", "Stroup"))
  
    
  # Plot deviation from empirical across all designs
  pdf(file=paste(c(output.figures.dir, "PowerBoxPlot_Overall.pdf"), collapse="/"), family="Times")
  par(lab=c(3,3,7))
  boxplot(diff ~ method, data=powerDataLong, las=1, ylim=c(-0.2,0.2),
          ylab="Deviation from Empirical Power")
  abline(h=0,lty=3)
  dev.off()
  
  # number of groups
  pdf(file=paste(c(output.figures.dir, "PowerBoxPlot_NumGroups.pdf"), collapse="/"), family="Times")
  par(mfrow=c(2,1), oma=c(5,1,1,1), mar=c(1,4,0,0), lab=c(3,3,7))
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$numGroups==2,],
          xaxt='n', ylim=c(-0.2, 0.2), las=1,
          ylab="2 treatments")
  abline(h=0,lty=3)
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$numGroups==4,],
          ylim=c(-0.2, 0.2), las=1,
          ylab="4 treatments")
  abline(h=0,lty=3)
  dev.off()
  
  # plot by small and large sample size
  pdf(file=paste(c(output.figures.dir, "PowerBoxPlot_PerGroupN.pdf"), collapse="/"), family="Times")
  par(mfrow=c(2,1), oma=c(5,1,1,1), mar=c(1,4,0,0), lab=c(3,3,7))
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$perGroupN==50,],
          xaxt='n', ylim=c(-0.2, 0.2), las=1,
          ylab="Per Group N = 50")
  abline(h=0,lty=3)
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$perGroupN==100,],
          ylim=c(-0.2, 0.2), las=1,
          ylab="Per Group N = 100")
  abline(h=0,lty=3)
  dev.off()
  
  # plot by amount of missing data
  pdf(file=paste(c(output.figures.dir, "PowerBoxPlot_MissingPercent.pdf"), collapse="/"), family="Times")
  par(mfrow=c(3,1), oma=c(5,1,1,1), mar=c(1,4,0,0), lab=c(3,3,7))
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$missingPercent==0,],
          xaxt='n', ylab="Complete data", las=1,
          ylim=c(-0.2, 0.2))
  abline(h=0,lty=3)
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$missingPercent==0.2,],
          xaxt='n', ylab="20% Missing", las=1,
          ylim=c(-0.2, 0.2))
  abline(h=0,lty=3)
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$missingPercent==0.4,],
          xaxt='n', ylab="40% Missing", las=1,
          ylim=c(-0.2, 0.2))
  abline(h=0,lty=3)
  dev.off()
  
  
  # plot by covariance
  pdf(file=paste(c(output.figures.dir, "PowerBoxPlot_MissingPercent.pdf"), collapse="/"), family="Times")
  par(mfrow=c(3,1), oma=c(5,1,1,1), mar=c(1,4,0,0), lab=c(3,3,7))
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$covariance=="CS",],
          xaxt='n', ylab="Compound Symmetry", las=1,
          ylim=c(-0.2, 0.2))
  abline(h=0,lty=3)
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$covariance=="CSH",],
          xaxt='n', ylab="Heterogeneous CS", las=1,
          ylim=c(-0.2, 0.2))
  abline(h=0,lty=3)
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$covariance=="AR(1)",],
          xaxt='n', ylab="Auto-regressive", las=1,
          ylim=c(-0.2, 0.2))
  abline(h=0,lty=3)
  dev.off()
  
  
  par(mfrow=c(3,1), oma=c(5,1,1,1), mar=c(1,4,0,0), lab=c(3,3,7))
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$covariance=="CS" & powerDataLong$missingPercent==0.4
                                            & powerDataLong$numGroups==4,],
          xaxt='n', ylab="Compound Symmetry", las=1,
          ylim=c(-0.1, 0.1))
  abline(h=0,lty=3)
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$covariance=="CSH" & powerDataLong$missingPercent==0.4
                                            & powerDataLong$numGroups==4,],
          xaxt='n', ylab="Heterogeneous CS", las=1,
          ylim=c(-0.1, 0.1))
  abline(h=0,lty=3)
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$covariance=="AR(1)" & powerDataLong$missingPercent==0.4
                                            & powerDataLong$numGroups==4,],
          xaxt='n', ylab="Auto-regressive", las=1,
          ylim=c(-0.1, 0.1))
  abline(h=0,lty=3)
  
}

#' runSimulationStudy
#' 
#' This function reproduces the simulation study results for the manuscript:\cr
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
#' A Power Approximation for Longitudinal Studies Using the 
#' Kenward and Roger Wald Test in the Linear Mixed Model, In review.
#'
#' @param study.seed the random number seed (defaults to 7634)
#' @param study.data.dir the directory into which data files are written (defaults to
#' current working directory)
#' @param study.figures.dir the directory into which pdf figures are written (defaults
#' to the current working directory)
#' @param study.runEmpirical if true, empirical power values will be recalculated. If false,
#' existing empirical power values will be loaded from the R package.
#' 
#' @note 
#' The empirical power calculations may take several hours to run
#'
runSimulationStudy <- function(study.seed=5896, study.data.dir=".", study.figures.dir=".",
                               study.runEmpirical=TRUE, study.runExemplary=TRUE) {
  # set the random seed
  set.seed(study.seed)
  
  # generate the designs for the validation study
  cat("### Generating designs and hypotheses\n")
  designList = generateDesignsForManuscript()
  
  # calculate empirical power
  if (study.runEmpirical) {
    cat("### Running empirical power calculations\n")
    # calculate empirical power for each design
    # !! Requires several hours to run !!
    empiricalPowerData = calculateEmpiricalPowerWithSAS(study.data.dir)
    
  } else {
    cat("### Loading existing empirical power calculations\n")
    # load the existing empirical data 
    data("empiricalPower", package="mixedPower")
  }
  
  # calculate approximate power using the exemplary data method
  # described by Stroup et al.
  if (study.runExemplary) {
    cat("### Running exemplary data approximate power calculations\n")
    # calculate empirical power for each design
    exemplaryPowerData = calculateExemplaryDataPowerWithSAS(study.data.dir)
    
  } else {
    cat("### Loading existing exemplary power data\n")
    # load the existing empirical data 
    data("exemplaryPower", package="mixedPower")
  }
  
  # calculate approximate power - runs in about 10-30 minutes
  cat("### Running approximate power calculations\n")
  approxPowerData = calculateApproximatePowerForDesignList(designList, study.data.dir)
  
  # combine the data into a single data set and write to disk
  cat("### Combining empirical and approximate values into common data set\n")
  load(paste(c(output.data.dir,"longitudinalParams.RData"),collapse="/"))
  approximateAndEmpiricalPowerData = data.frame(paramComboList, approxPowerData, 
                                                exemplaryPower=exemplaryPowerData$exemplaryPower,
                                                empiricalPower=empiricalPowerData$empiricalPower)
  # save to disk
  save(approximateAndEmpiricalPowerData,
       file=paste(c(study.data.dir, "approximateAndEmpiricalPower.RData"), collapse="/"))
  
  # Produce summary figures and calculate max deviations for each method
  cat("### Summarizing results and generating figures\n")
  summarizeResults(study.data.dir, study.figures.dir)
}



