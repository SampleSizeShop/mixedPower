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


calculateEmpiricalPowerForDesignList <- function(output.data.dir) {
  # call the sas code to run empirical power
  result <- system(paste(c("sas.exe -i ", 
                           paste(c(path.package("mixedPower"), 
                                   "inst/sas/calculateEmpiricalPower.sas"), 
                                 collapse="/")), collapse=""), 
                   intern = TRUE, show.output.on.console = TRUE)
  
  # load the empirical power results
  empiricalPowerData = read.csv(
    paste(c(path.package("mixedPower"), "inst/sas/empiricalPower.csv"), collapse="/"),
    stringsAsFactors=FALSE)

  # save as a Rdata file
  save(empiricalPowerData, file=paste(c(output.data.dir, "empiricalPower.RData"), collapse="/"))
}

#' calculateApproximatePowerForDesignList
#' 
#' Calculate approximate power values for a list of design/glh pairs.
#'   
#' @param designList list of pairs of design.glmmFG and glh objects
#' @param output.data.dir directory to which data sets are written
#' @return design.glmmFG object with a single covariate
#' @note This function takes about 15-30 minutes to run
#' 
calculateApproximatePowerForDesignList <- function(designList, output.data.dir=".") {

  # add power using the method described by 
  #
  # Helms, R. W. (1992). Intentionally incomplete longitudinal designs: 
  # I. Methodology and comparison of some full span designs. 
  # Statistics in Medicine, 11(14-15), 1889–1913.
  #
  approxPowerList = sapply(designList, function(designAndGlh) {
    return(mixedPower.helms(designAndGlh[[1]], designAndGlh[[2]]))
  })
  
  # add power using the method described by (pp. )
  #
  # Littell, P. D., Russell C., Milliken, P. D., George A., 
  # Stroup, P. D., Walter W., Wolfinger, P. D., Russell D., & Schabenberger, 
  # P. D., Oliver. (2006). SAS for Mixed Models, Second Edition (2nd ed.). SAS Institute.
  #
  approxPowerList = sapply(longitudinalDesignList, function(designAndGlh) {
    return(mixedPower.stroup(designAndGlh[[1]], designAndGlh[[2]]))
  })
  
  # add power assuming complete, balanced data. This method is described by 
  #
  # Muller, K. E., Lavange, L. M., Ramey, S. L., & Ramey, C. T. (1992). 
  # Power Calculations for General Linear Multivariate Models Including 
  # Repeated Measures Applications. Journal of the American Statistical 
  # Association, 87(420), 1209–1226.
  #
  approxPowerData$power.fixedOnly=sapply(designList, function(x) {
    multivariateDesign = x[[3]]
    multivariateHypothesis = x[[4]]
    return(glmmPower.fixed(newDesign, newHypothesis)) 
  })
  
  # add power using the method described by
  #
  # Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
  # A Power Approximation for Longitudinal Studies Using the 
  # Kenward and Roger Wald Test in the Linear Mixed Model, In review.
  #
  approxPowerList = sapply(longitudinalDesignList, function(designAndGlh) {
    return(mixedPower(designAndGlh[[1]], designAndGlh[[2]]))
  })
  
  ## write the approximate power data to disk     
  save(approxPowerData,
       file=paste(c(output.data.dir, "approximatePower.RData"), collapse="/"))
  
  return(approxPowerData)
}

#' generateDesignsForManuscript 
#' 
#' Generate the list of test cases (design.glmmFG and glh objects) for
#' the manuscript:\cr
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
#' A Power Approximation for Longitudinal Studies Using the 
#' Kenward and Roger Wald Test in the Linear Mixed Model, In review.
#' 
#' @return list of pairs of design.glmmFG and glh objects
#' @keywords internal
#' 
generateDesignsForManuscript <- function() {
  
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
                               study.runEmpirical=TRUE, study.sasEmpirical=FALSE) {
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
    empiricalPowerData = calculateEmpiricalPowerForDesignList(designList, study.data.dir)
    
  } else {
    cat("### Loading existing empirical power calculations\n")
    # load the existing empirical data 
    data("empiricalPower", package="rPowerlib")
  }
  
  # calculate approximate power - runs in about 10-30 minutes
  cat("### Running approximate power calculations\n")
  approxPowerData = calculateApproximatePowerForDesignList(designList, study.data.dir)
  
  # combine the data into a single data set and write to disk
  cat("### Combining empirical and approximate values into common data set\n")
  approximateAndEmpiricalPowerData = data.frame(approxPowerData, 
                                                empiricalPower=empiricalPowerData$empiricalPower,
                                                empiricalTime=empiricalPowerData$time)
  # save to disk
  save(approximateAndEmpiricalPowerData,
       file=paste(c(study.data.dir, "approximateAndEmpiricalPower.RData"), collapse="/"))
  
  # Produce summary figures and calculate max deviations for each method
  cat("### Summarizing results and generating figures\n")
  summarizeResults(study.data.dir, study.figures.dir)
}



