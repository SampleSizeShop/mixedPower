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
# Generate a single 4 group cluster design
#
generateClusterRandomizedDesign = function(params) {

  name = "clusterRandomized"
  description = paste(c(params$numGroups, " group cluster design, test of main effect of treatment. ",
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
  betaScale = getBetaScaleByPower(design, glh, targetPower=params$targetPower, 
                                  lower=0.001, upper=100)  
  design@beta = matrix(c(1,rep(0,params$numGroups-1))) * betaScale  
  
  return(design)

}



generateDesigns.clusterRandomized = function() {
  #
  # For each cluster randomized design, we
  # vary the following parameters
  #
  # number of treatment groups
  numGroupsList = c(2, 4)
  # total clusters per treatment group
  perGroupNList = c(10, 40)
  # total participants per cluster
  clusterSizeList = c(5, 50)
  # percent missing (in half of the clusters)
  missingPercentList = c(0, 0.20, 0.40)
  # in all cases, we select the scale factor 
  # for beta to achieve the following power
  targetPowerList = c(0.2, 0.5, 0.8)
  # intracluster correlation
  iccList = c(0.04)
  # sigma squared
  sigmaSqList = c(2)
  
  # generate parameters
  paramList = list(missingPercent=missingPercentList, 
                   clusterSize=clusterSizeList,
                   perGroupN=perGroupNList, 
                   targetPower=targetPowerList,
                   sigmaSq=sigmaSqList,
                   icc=iccList,
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
    cat("Generating design ", i, "\n")
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
  approxPowerList = sapply(clusterDesignList, function(designAndGlh) {
    return(mixedPower(designAndGlh[[1]], designAndGlh[[2]]))
  })

  # combine with the empirical set and save to disk
  powerResults$approxPower = approxPowerList
  write.csv(powerResults, dataFile("clusterRandomizedResults.csv"))
  
}

#
# Build summary tables and power curves
#
summarizeResults.clusterRandomized = function() {
  powerResults = read.csv(dataFile("clusterRandomizedResults.csv"))
  
  # remove designs which violate assumption of Nd > q + pd + 1
  powerResults = powerResults[powerResults$clusterSize != 50 | 
                                powerResults$perGroupN !=10,]
  powerResults = powerResults[powerResults$clusterSize != 5 | 
                                powerResults$numGroups != 2,]
  
  powerResults$deviation = powerResults$approxPower - powerResults$empiricalPower
  mean(powerResults$deviation)
  range(powerResults$deviation)
  
  pdf(file="../inst/figures/ClusterPowerBoxPlots.pdf", height=5)
  par(mfrow=c(1,3), oma=c(5,5,5,5), mar=c(5,2,1,1))
  boxplot(powerResults$deviation ~ powerResults$numGroups, ylim=c(-0.1,0.1),
          xlab="Number of Treatment Groups")
  abline(h=0, lty=3)
  boxplot(powerResults$deviation ~ powerResults$clusterSize, ylim=c(-0.1,0.1),
          xlab="Cluster Size")
  abline(h=0, lty=3)
  boxplot(powerResults$deviation ~ powerResults$missingPercent, ylim=c(-0.1,0.1),
          xlab="Percent Missing (half of clusters)")
  abline(h=0, lty=3)  
  dev.off()
  
}

#
# Applied example for the manuscript
#
# A proposed cluster-randomized trial in which
# worksites are randomized to 2 smoking cessation programs
#
hennrikusExample = function() {
  
  completeSize = 30
  clusterComplete = matrix(rep(1,completeSize))
  incompleteSize = 20
  clusterIncomplete = matrix(rep(1,incompleteSize))
  
  group1X = matrix(c(1,0), nrow=1)
  group2X = matrix(c(0,1), nrow=1)
  
  design = new("design.mixed", name = "Hennrikus Examle", 
               description = "Cluster randomized trial of smoking cessation intervention",
               xPatternList = c(
                 new("missingDataPattern", group=1, observations=1:completeSize, size=25,
                     designMatrix=clusterComplete %x% group1X),
                 new("missingDataPattern", group=1, observations=1:incompleteSize, size=15,
                     designMatrix=clusterIncomplete %x% group1X),
                 new("missingDataPattern", group=1, observations=1:completeSize, size=25,
                     designMatrix=clusterComplete %x% group2X),
                 new("missingDataPattern", group=1, observations=1:incompleteSize, size=15,
                     designMatrix=clusterIncomplete %x% group2X)
                 ),
               beta = matrix(c(25,0)),
               Sigma = (125^2) * (0.04 * (matrix(rep(1,30)) %*% t(matrix(rep(1,30)))) + 
                                           diag(30) * (1 - 0.04))
  )
  # get the appropriate hypothesis
  glh = new("glh.mixed",
            alpha = 0.05,
            fixedContrast = matrix(c(1,-1), nrow=1),
            thetaNull = matrix(0),
            test = "Wald, KR ddf")
  
  return(mixedPower(design, glh))
  
}

#generateDesigns.clusterRandomized()