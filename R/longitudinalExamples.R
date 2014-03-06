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
# Example power calculations for cluster randomized trials
# with 4 treatments
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
# Identify the beta scale such that mixedPower returns the 
# requested target power
#
getBetaScaleByPower <- function(design, glh, targetPower=0.90, lower=0, upper=1000) {
  betaScale = uniroot(function(x) {
    design@beta = design@beta * x
    return(mixedPower(design, glh) - targetPower)
  }, c(lower, upper))
  return(betaScale$root)
}

getGlhByNumGroups = function(numGroups) {
  return (new("glh.mixed",
              alpha = 0.05,
              fixedContrast = cbind(matrix(rep(1,numGroups-1)), -1*diag(numGroups-1)),
              thetaNull = matrix(rep(0,numGroups-1)),
              test = "Wald, KR ddf"))
}

#
# Generate longitudinal design
#
generateLongitudinalDesign = function(params) {
  
  name = "longitudinal"
  description = paste(c(params$numGroups, " group longitudinal design, test of time by treatment interaction. ",
                        "Per group N = ", params$perGroupN, 
                        ", cluster size = ", params$clusterSize,
                        ", missing % (in half of ISUs) = ", params$missingPercent,
                        ", target power = ", params$targetPower), collapse="")
  incompleteSize = floor(params$clusterSize * (1-params$missingPercent))
  
  # build the pattern list
  patternList = list()
  pattern = 1
  for(grp in 1:params$numGroups) {
    patternList[[pattern]] = 
      new("missingDataPattern", group=grp, observations=1:params$clusterSize, 
          designMatrix = matrix(rep(1,params$clusterSize)) %x% matrix(diag(params$numGroups)[grp,],nrow=1),
          size=params$perGroupN/2)
    patternList[[pattern+1]] = 
      new("missingDataPattern", group=grp, observations=1:incompleteSize, 
          designMatrix = matrix(rep(1,incompleteSize)) %x% matrix(diag(params$numGroups)[grp,],nrow=1),
          size=params$perGroupN/2)
    pattern = pattern + 2
  }
  
  # build the design
  cluster = matrix(rep(1,params$clusterSize))
  design = new("design.mixed", name = name, description = description,
               xPatternList = patternList,
               beta = matrix(c(1,rep(0,params$numGroups-1))),
               Sigma = params$sigmaSq * (params$icc * (cluster %*% t(cluster)) + 
                                           diag(params$clusterSize) * (1 - params$icc))
  )
  # get the appropriate hypothesis
  glh = getGlhByNumGroups(params$numGroups)
  # update the beta scale
  betaScale = getBetaScaleByPower(design, glh, targetPower=params$targetPower)  
  design@beta = matrix(c(1,rep(0,params$numGroups-1))) * betaScale  
  
  return(design)
  
}



generateDesigns.longitudinal = function() {
  #
  # For each cluster randomized design, we
  # vary the following parameters
  #
  # number of treatment groups
  numGroupsList = c(2, 4)
  # total clusters per treatment group
  perGroupNList = c(10, 40)
  # total participants per cluster
  numObservations = c(3, 5, 10)
  # percent missing (in half of the clusters)
  monotone = c(1, 0)
  # in all cases, we select the scale factor 
  # for beta to achieve the following power
  targetPowerList = c(0.2, 0.5, 0.8)
  # Lear parameters
  rhoList = c(0.4)
  deltaList = c(0.5)
  # sigma squared
  sigmaSqList = c(2)
  
  # generate parameters
  paramList = list(missingPercent=missingPercentList, 
                   clusterSize=clusterSizeList,
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
  clusterDesignList = list()
  betaScaleList = vector()
  for(i in 1:length(paramComboList$numGroups)) {
    params = paramComboList[i,]
    clusterDesignList[[i]] = list(generateClusterRandomizedDesign(params), 
                                  getGlhByNumGroups(params$numGroups))                                     
    betaScaleList[i] = clusterDesignList[[i]][[1]]@beta[1,1]
  }
  paramComboList$betaScale = betaScaleList
  
  # write the parameter data to a csv file
  write.csv(paramComboList, file=dataFile("clusterRandomizedParams.csv"),
            row.names=FALSE, eol="\r\n")
  # write the designs to an Rdata file
  save(clusterDesignList, file=dataFile("clusterRandomizedDesigns.RData"))
  
}

calculatePower.clusterRandomized = function(runEmpirical=FALSE) {
  if (!file.exists(dataFile("clusterRandomizedParams.csv")) ||
        !file.exists(dataFile("clusterRandomizedDesigns.RData"))) {
    generateDesigns.clusterRandomized()
  }
  
  if (runEmpirical) {
    # exec SAS file to run empirical power for 4 group cluster designs
    # requires SAS installation
  }
  
  empiricalFile = dataFile("clusterRandomizedEmpirical.csv");
  if (!file.exists(empiricalFile)) {
    stop(paste(c("Missing empirical power file: ", empiricalFile), collapse=""))
  }
  # load the empirical values
  powerResults = read.csv(empiricalFile)
  
  # load the designs and calculate power
  load(dataFile("clusterRandomizedDesigns.RData"))
  for(i in 1:length(powerResults$targetPower)) {
    
  }
  
  # combine with the empirical set and save to disk
  powerResults$approxPower = approxPowerList
  write.csv(dataFile("clusterRandomizedResults.csv"))
  
}

#
# Build summary tables and power curves
#
summarizeResults.clusterRandomized = function() {
  powerResults = read.csv(dataFile("cluster4GroupPowerResults.csv"))
  
  
}


generateDesigns.clusterRandomized()

