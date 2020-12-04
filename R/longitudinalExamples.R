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

heterogeneousCS <- function(sigmaList,rho) {
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
# Generate an auto-regressive covariance matrix
#
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

#
# Identify the beta scale such that mixedPower returns the 
# requested target power
#
getBetaScaleByPower <- function(design, glh, targetPower=0.90, 
                                lower=0.0001, upper=1000) {
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
    Sigma = heterogeneousCS(rep(1,maxObservations),rho)
  } else if (params$covariance == 'CSH') {
    Sigma = heterogeneousCS(c(1,0.5,0.3,0.1,0.1),rho)
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

generateDesigns.longitudinal = function(data.dir=getwd()) {
  #
  # For each longitudinal randomized design, we
  # vary the following parameters
  #
  # number of treatment groups
  numGroupsList = c(2, 4)
  # total ISUs per treatment group
  perGroupNList = c(50)
  # max observations for each participants
  maxObservationsList = c(5)
  # missing data pattern (either monotone or non-monotone)
  monotoneList = c(1, 0)
  # percent missing
  missingPercentList = c(0, 0.2, 0.4)
  # in all cases, we select the scale factor 
  # for beta to achieve the following power
  targetPowerList = c(0.2, 0.5, 0.8)
  # Lear parameters
  rhoList = c(0.4)
  deltaList = c(1)
  # sigma squared
  sigmaSqList = c(1)
  # covariance 
  covarianceList = c("AR(1)")
   
  # generate parameters
  paramList = list(monotone=monotoneList, 
                   maxObservations=maxObservationsList,
                   missingPercent=missingPercentList,
                   perGroupN=perGroupNList, 
                   targetPower=targetPowerList,
                   rho=rhoList,
                   delta=deltaList,
                   sigmaSq=sigmaSqList,
                   numGroups=numGroupsList,
                   covariance=covarianceList)
  paramComboList = data.frame(expand.grid(paramList))
  paramComboList$missingType = ifelse(paramComboList$monotone == 1, "monotone", "non-monotone")
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
  write.csv(paramComboList, 
            file=file.path(data.dir,"longitudinalParams.csv"),
            row.names=FALSE, eol="\r\n")
  # write the designs to an Rdata file
  save(longitudinalDesignList, 
       file=file.path(data.dir,"longitudinalDesigns.RData"))
  
}

# 
# generateDesigns.longitudinal = function(output.data.dir=getwd()) {
#   #
#   # For each longitudinal randomized design, we
#   # vary the following parameters
#   #
#   # number of treatment groups
#   numGroupsList = c(2, 4)
#   # total ISUs per treatment group
#   perGroupNList = c(50)
#   # missing data pattern (either monotone or non-monotone)
#   missingTypeList = c("monotone", "non-monotone")
#   # percent missing
#   missingPercentList = c(0, 0.2, 0.4)
#   # covariance
#   covarianceList = c("AR(1)")
#   # in all cases, we select the scale factor 
#   # for beta to achieve the following power
#   targetPowerList = c(0.2, 0.5, 0.8)
#   # maximum number of observations 
#   maxObservationsList = c(5)
#   
#   # generate parameters
#   paramList = list(targetPower=targetPowerList,
#                    covariance=covarianceList,
#                    missingType=missingTypeList,
#                    missingPercent=missingPercentList,
#                    perGroupN=perGroupNList,
#                    numGroups=numGroupsList,
#                    maxObservations=maxObservationsList)
#   paramComboList = data.frame(expand.grid(paramList))
#   
#   #
#   # Calculate the appropriate betaScale values
#   # and build the list of designs
#   #
#   longitudinalDesignList = list()
#   betaScaleList = vector()
#   for(i in 1:length(paramComboList$numGroups)) {
#     cat("Case ", i, "\n")
#     params = paramComboList[i,]
#     longitudinalDesignList[[i]] = list(generateLongitudinalDesign(params), 
#                                   getGlh(params$numGroups, params$maxObservations))                                     
#     betaScaleList[i] = longitudinalDesignList[[i]][[1]]@beta[1,1]
#   }
#   paramComboList$betaScale = betaScaleList
#   
#   # write the parameter data to a csv file
#   write.csv(paramComboList, 
#             file=file.path(data.dir,"longitudinalParams.csv"),
#             row.names=FALSE, eol="\r\n")
#   # write the designs to an Rdata file
#   save(longitudinalDesignList, 
#        file=file.path(data.dir,"longitudinalDesigns.RData"))
#   
# }

calculatePower.longitudinal = function(data.dir=getwd(), figures.dir=getwd(),
                                       runEmpirical=FALSE) {
  if (!file.exists(file.path(data.dir,"longitudinalParams.csv")) ||
        !file.exists(file.path(data.dir,"longitudinalDesigns.RData"))) {
    generateDesigns.longitudinal()
  }
  
  if (runEmpirical) {
    # exec SAS file to run empirical power for longitudinal designs
    # requires SAS installation
    empiricalFile = file.path(data.dir,"longitudinalEmpirical.csv");
    
  } else {
    empiricalFile = file.path(path.package("mixedPower"),"data","longitudinalEmpiricalPower.csv")

  }
  
  if (!file.exists(empiricalFile)) {
    stop(paste(c("Missing empirical power file: ", empiricalFile), collapse=""))
  }
  
  # load the empirical values
  powerResults = read.csv(empiricalFile)
  
  # load the designs and calculate power
  load(file.path(data.dir,"longitudinalDesigns.RData"))
  approxPowerList = sapply(longitudinalDesignList, function(designAndGlh) {
    return(mixedPower(designAndGlh[[1]], designAndGlh[[2]]))
  })
  
  # combine with the empirical set and save to disk
  powerResults$approxPower = approxPowerList
  write.csv(powerResults, file=file.path(data.dir,"longitudinalResults.csv"))
  
}

#
# Build summary tables and power curves
#
summarizeResults.longitudinal = function(data.dir=getwd(), figures.dir=getwd()) {
  powerResults = read.csv(file.path(data.dir,"longitudinalResults.csv"))
  
  powerResults$deviation = powerResults$approxPower - powerResults$empiricalPower
  mean(powerResults$deviation)
  range(powerResults$deviation)
  fivenum(powerResults$deviation)
  
  pdf(file=file.path(figures.dir, "LongitudinalPowerBoxPlots.pdf"), height=5)
  par(mfrow=c(1,3), oma=c(5,5,5,5), mar=c(5,2,1,1))
  boxplot(powerResults$deviation ~ powerResults$numGroups, ylim=c(-0.1,0.1),
          xlab="Total Treatment Groups")
  abline(h=0, lty=3)
  powerResults$pattern = factor(powerResults$monotone, 
                              levels=c(1,0),
                              labels=c("Monotone", "Non"))
  boxplot(powerResults$deviation ~ powerResults$pattern, ylim=c(-0.1,0.1),
          xlab="Missing Data Pattern")
  abline(h=0, lty=3)
  boxplot(powerResults$deviation ~ powerResults$missingPercent, ylim=c(-0.1,0.1),
          xlab="Total Incomplete Sampling Units")
  abline(h=0, lty=3)
  
  dev.off()
  
  tiff(file=file.path(figures.dir, "LongitudinalPowerBoxPlots.tiff"), units="in", 
                      width=5, height=5, res=300)
  par(mfrow=c(1,3), oma=c(5,5,5,5), mar=c(5,2,1,1))
  boxplot(powerResults$deviation ~ powerResults$numGroups, ylim=c(-0.1,0.1),
          xlab="Total Treatment Groups")
  abline(h=0, lty=3)
  powerResults$pattern = factor(powerResults$monotone, 
                                levels=c(1,0),
                                labels=c("Monotone", "Non"))
  boxplot(powerResults$deviation ~ powerResults$pattern, ylim=c(-0.1,0.1),
          xlab="Missing Data Pattern")
  abline(h=0, lty=3)
  boxplot(powerResults$deviation ~ powerResults$missingPercent, ylim=c(-0.1,0.1),
          xlab="Total Incomplete Sampling Units")
  abline(h=0, lty=3)
  
  dev.off()
}

# 
# generateDesigns.longitudinal()
# summarizeResults.longitudinal()
