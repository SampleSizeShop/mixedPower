# 
#  Package mixedPower calculates power for the linear mixed model
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
# Example power calculations for longitudinal designs
# 
#
source("../R/mixedPower.R")

#
# Convenience routine for generate paths
# to files in the data directory
#
dataFile = function(filename) {
  return(paste(c("../data/", filename), collapse=""))
}

#
# Calculate a lear correlation matrix of the specified
# size, correlation, and decay rate
# Assume equal spacing
#
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

#
# Identify the beta scale such that mixedPower returns the 
# requested target power
#
getBetaScaleByPower <- function(design, glh, targetPower=0.90, lower=0.0001, upper=1000) {
  betaScale = uniroot(function(x) {
    design@beta = design@beta * x
    return(mixedPower(design, glh) - targetPower)
  }, c(lower, upper))
  return(betaScale$root)
}

#
# Get time by treatment interaction hypothesis object
#
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



#
# Generate longitudinal design
#
generateLongitudinalDesign = function(params) {
  
  name = "longitudinal"
  description = paste(c(params$numGroups, " group longitudinal design, test of time by treatment interaction. ",
                        "Per group N = ", params$perGroupN, 
                        ", max observations = ", params$maxObservations,
                        ", monotone missing? = ", (params$monotone==1),
                        ", target power = ", params$targetPower), collapse="")

  
  # build the pattern list
  patternList = list()
  pattern = 1
  for(grp in 1:params$numGroups) {
    # between effects design matrix
    designBetween = matrix(diag(params$numGroups)[grp,], nrow=1)
    
    if (params$monotone) {
      ## monotone dropout pattern - once missing, never come back
      # complete case
      patternList[[pattern]] = 
        new("missingDataPattern", group=grp, observations=1:params$maxObservations, 
            designMatrix = (designBetween %x% diag(params$maxObservations)),
            size=floor(params$perGroupN*0.5))
      # missing last observation
      missingPattern1 = 1:(params$maxObservations-1)
      patternList[[pattern+1]] = 
        new("missingDataPattern", group=grp, observations=missingPattern1, 
            designMatrix = (designBetween %x% 
                              matrix(diag(params$maxObservations)[missingPattern1,], 
                                     nrow=length(missingPattern1))),
            size=floor(params$perGroupN*0.3))
      # missing last 2 observations
      missingPattern2 = 1:(params$maxObservations-2)
      patternList[[pattern+2]] = 
        new("missingDataPattern", group=grp, observations=missingPattern2, 
            designMatrix = (designBetween %x% 
                              matrix(diag(params$maxObservations)[missingPattern2,], 
                                     nrow=length(missingPattern2))),
            size=floor(params$perGroupN*0.2))

    } else {
      ## Non-monotone dropout pattern 
      ## - delete 2nd observations in 30%, delete 3rd in 20%
      ## - if max observations > 3, also delete 4th
      # complete case
      patternList[[pattern]] = 
        new("missingDataPattern", group=grp, observations=1:params$maxObservations, 
            designMatrix = (designBetween %x% diag(params$maxObservations)),
            size=floor(params$perGroupN*0.5))
      # missing 2nd (and 4th) observation(s)
      if (params$maxObservations == 3) {
        missingPattern1 = c(1,3)
      } else {
        missingPattern1 = c(1,3,5:params$maxObservations)
      }
      patternList[[pattern+1]] = 
        new("missingDataPattern", group=grp, observations=missingPattern1, 
            designMatrix = (designBetween %x% 
                              matrix(diag(params$maxObservations)[missingPattern1,], 
                                     nrow=length(missingPattern1))),
            size=floor(params$perGroupN*0.3))
      # missing 3rd (and 4th) observation(s)
      if (params$maxObservations == 3) {
        missingPattern2 = c(1,2)
      } else {
        missingPattern2 = c(1,2,5:params$maxObservations)
      }
      patternList[[pattern+2]] = 
        new("missingDataPattern", group=grp, observations=missingPattern2, 
            designMatrix = (designBetween %x% 
                              matrix(diag(params$maxObservations)[missingPattern2,], 
                                     nrow=length(missingPattern2))),
            size=floor(params$perGroupN*0.2))
      
    }
    pattern = pattern + 3
  }
  
  # build the design
  design = new("design.mixed", name = name, description = description,
               xPatternList = patternList,
               beta = matrix(c(1,rep(0,(params$maxObservations-1)),rep(0,params$maxObservations*(params$numGroups-1)))),
               Sigma = learMatrix(params$maxObservations,params$rho,params$delta)
  )
  # get the appropriate hypothesis
  glh = getGlh(params$numGroups, params$maxObservations)
  # update the beta scale
  betaScale = getBetaScaleByPower(design, glh, targetPower=params$targetPower)  
  design@beta = matrix(c(1,rep(0,params$numGroups-1))) * betaScale  
  
  return(design)
  
}



generateDesigns.longitudinal = function() {
  #
  # For each longitudinal randomized design, we
  # vary the following parameters
  #
  # number of treatment groups
  numGroupsList = c(2, 4)
  # total ISUs per treatment group
  perGroupNList = c(30, 60)
  # max observations for each participants
  maxObservationsList = c(5, 10)
  # missing data pattern (either monotone or non-monotone)
  monotoneList = c(1, 0)
  # in all cases, we select the scale factor 
  # for beta to achieve the following power
  targetPowerList = c(0.2, 0.5, 0.8)
  # Lear parameters
  rhoList = c(0.4)
  deltaList = c(0.5)
  # sigma squared
  sigmaSqList = c(2)
  
  # generate parameters
  paramList = list(monotone=monotoneList, 
                   maxObservations=maxObservationsList,
                   perGroupN=perGroupNList, 
                   targetPower=targetPowerList,
                   rho=rhoList,
                   delta=deltaList,
                   sigmaSq=sigmaSqList,
                   numGroups=numGroupsList)
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
  
  # write the parameter data to a csv file
  write.csv(paramComboList, file=dataFile("longitudinalParams.csv"),
            row.names=FALSE, eol="\r\n")
  # write the designs to an Rdata file
  save(longitudinalDesignList, file=dataFile("longitudinalDesigns.RData"))
  
}

calculatePower.longitudinal = function(runEmpirical=FALSE) {
  if (!file.exists(dataFile("longitudinalParams.csv")) ||
        !file.exists(dataFile("longitudinalDesigns.RData"))) {
    generateDesigns.longitudinal()
  }
  
  if (runEmpirical) {
    # exec SAS file to run empirical power for longitudinal designs
    # requires SAS installation
  }
  
  empiricalFile = dataFile("longitudinalEmpirical.csv");
  if (!file.exists(empiricalFile)) {
    stop(paste(c("Missing empirical power file: ", empiricalFile), collapse=""))
  }
  # load the empirical values
  powerResults = read.csv(empiricalFile)
  
  # load the designs and calculate power
  load(dataFile("longitudinalDesigns.RData"))
  for(i in 1:length(powerResults$targetPower)) {
    
  }
  
  # combine with the empirical set and save to disk
  powerResults$approxPower = approxPowerList
  write.csv(dataFile("longitudinalResults.csv"))
  
}

#
# Build summary tables and power curves
#
summarizeResults.longitudinal = function() {
  powerResults = read.csv(dataFile("longitudinalResults.csv"))
  
  
}


generateDesigns.longitudinal()

