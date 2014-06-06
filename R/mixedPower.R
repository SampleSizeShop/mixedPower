#
# Paper: Kenward and Roger Power
# Author: Sarah Kreidler
# Date: 4/19/2013
#
# calculatePowerKenwardRoger.sxs
#
# SAS/IML Module to calculate power for the Wald test
# with denominator degrees of freedom as described by Kenward and Roger
# 
#
#


#
#
#
mixedPower.stroup <- function(design, glh) {
  # create an exemplary data set
  
  # write the exemplary data to disk
  
  # write the SAS code to fit the exemplary data and calculate power
  
  
  # run the SAS code
  
  # read the resulting power result
  
  # cleanup the temporary files
  
  # return the power value
    
}

#
# Helms power approximation for mixed models
#
mixedPower.helms = function(design, glh) {
  
  # get the max number of planned observations 
  maxObs = nrow(design@Sigma)
  
  designMatrixList = list()
  SigmaList = list()
  isu = 1
  for(i in 1:length(design@xPatternList)) {
    pattern = design@xPatternList[[i]]
    # get sigma matrix for this pattern
    deletionMatrix = matrix(diag(maxObs)[pattern@observations,], 
                            nrow=length(pattern@observations))
    SigmaD = (deletionMatrix %*% design@Sigma %*% t(deletionMatrix))
    
    # add to the list
    for(j in 1:pattern@size) {
      designMatrixList[[isu]] = pattern@designMatrix
      SigmaList[[isu]] = SigmaD
      isu = isu + 1
    }
  }
  
  # calculate the approximate covariance of thetaDiff = C*beta-thetaNull
  sumXtSigmaInvX = Reduce("+", lapply(1:length(designMatrixList), function(isu) {
    return(t(designMatrixList[[isu]]) %*% solve(SigmaList[[isu]]) %*% 
             designMatrixList[[isu]])
  }))
  thetaDiffHat.Sigma = (glh@fixedContrast %*%
                          solve(sumXtSigmaInvX) %*%
                          t(glh@fixedContrast))
  # calculate the mean
  thetaDiffHat.mu = glh@fixedContrast %*% design@beta - glh@thetaNull
  omega = t(thetaDiffHat.mu) %*% solve(thetaDiffHat.Sigma) %*% thetaDiffHat.mu
  
  # numerator degrees of freedom
  a = nrow(glh@fixedContrast)
  # denominator degrees of freedom
  # get the total number of independent sampling units
  N = sum(sapply(design@xPatternList, function(pattern) {
    return(pattern@size)
  }))
  # get the rank of the design matrix (assumes X is full rank)
  rankX = ncol(design@xPatternList[[1]]@designMatrix)
  ddf = N - rankX
  
  # get the critical value
  Fcrit = qf(1-glh@alpha,a,ddf)
  # calculate power;
  power = 1 - pf(Fcrit, a, ddf, omega);
}


#
# Returns lambda, omega, and ddf for the distribution
# of the Kenward-Roger statistic
#
getKRParams = function(a, moments) {
  # calculate rho
  rho = moments$var.alt / (2*((moments$mu.null)^2)) 
  # calculate noncentrality parameter
  omega = a*((moments$mu.alt/moments$mu.null)-1)
  # denominator degrees of freedom
  ddf = 4 + (2*(a + 2*omega) + (a + omega)^2)/(rho*(a^2) - a - 2*omega)
  # scale factor for KR adjustment
  lambda = ddf / ((ddf - 2)*moments$mu.null)
  
  return(data.frame(lambda=lambda, ddf=ddf, omega=omega))
}

#
# Get the approximate moments of the Wald
# statistic given the mixed model design
#
#
getWaldMoments = function(design, glh, homoscedastic=TRUE) {
  
  # determine the number of unique treatment groups
  numGroups = length(unique(sapply(design@xPatternList, function(x) {
    return(x@group)
  })))
  
  # get the max number of planned observations 
  maxObs = nrow(design@Sigma)
  
  # generate the inverse Wishart list 
  patternDfList = vector()
  for(pattern in design@xPatternList) {
    patternString = paste(sort(pattern@observations), collapse="_")
    if (is.na(patternDfList[patternString])) {
      patternDfList[patternString] = pattern@size
    } else {
      patternDfList[patternString] = patternDfList[patternString] + pattern@size
    }
  }
  
  # generate the corresponding design matrix list.
  # also generate the list of Sigma matrices for each ISU
  invWishartList = list()
  designMatrixList = list()
  SigmaList = list()
  isu = 1
  for(i in 1:length(design@xPatternList)) {
    pattern = design@xPatternList[[i]]
    # get sigma matrix for this pattern
    deletionMatrix = matrix(diag(maxObs)[pattern@observations,], 
                            nrow=length(pattern@observations))
    SigmaD = (deletionMatrix %*% design@Sigma %*% t(deletionMatrix))
              
    # add to the list
    for(j in 1:pattern@size) {
      designMatrixList[[isu]] = pattern@designMatrix
      df = as.numeric(patternDfList[paste(sort(pattern@observations), collapse="_")]) - numGroups
      invWishartList[[isu]] = 
        invert(new("wishart", 
                   df=df,
                   covariance=(SigmaD/df)))
      SigmaList[[isu]] = SigmaD
      isu = isu + 1
    }
  }
  class(invWishartList) = "inverseWishart"
  
  # calculate the approximating inverse Wishart for X' * SigmaInv * X
  distXtSigmaInvX = approximateInverseWishart(invWishartList, 
                                              designMatrixList,
                                              method="trace")
  # now invert to get a Wishart
  distXtSigmaInvXInv = invert(distXtSigmaInvX)
  # scale by C matrices and form corresponding wishart
  distCXtSigmaInvXInvCt = new("wishart", df=distXtSigmaInvXInv@df,
                              covariance=(glh@fixedContrast %*% 
                                            distXtSigmaInvXInv@covariance %*% 
                                            t(glh@fixedContrast))) 
  
  # calculate the approximate covariance of thetaDiff = C*beta-thetaNull
  sumXtSigmaInvX = Reduce("+", lapply(1:length(designMatrixList), function(isu) {
    return(t(designMatrixList[[isu]]) %*% solve(SigmaList[[isu]]) %*% 
             designMatrixList[[isu]])
  }))
  thetaDiffHat.Sigma = (glh@fixedContrast %*%
                          solve(sumXtSigmaInvX) %*%
                          t(glh@fixedContrast))
  # calculate the mean
  thetaDiffHat.mu = glh@fixedContrast %*% design@beta - glh@thetaNull

  # form the appropriate F distribution
  a = nrow(glh@fixedContrast)  
  if (a == 1) {
    # special case when contrast has a single row
    ndf = 1
    ddf = distCXtSigmaInvXInvCt@df
    omega = t(thetaDiffHat.mu) %*% solve(thetaDiffHat.Sigma) %*% thetaDiffHat.mu
    Fscale = ((1/(ddf)) * 
                (as.numeric(thetaDiffHat.Sigma) / 
                   as.numeric(distCXtSigmaInvXInvCt@covariance)))
    
    
  } else {
    # a > 1 case
    precision = solve(distCXtSigmaInvXInvCt@covariance)
    covarRatio = precision %*% thetaDiffHat.Sigma
    tmpOmega = t(thetaDiffHat.mu) %*% precision %*% thetaDiffHat.mu
    deltaStar = (
      ((tmpOmega * sum(diag(covarRatio))) + 2 * tmpOmega^2) / 
        (sum(diag(covarRatio %*% covarRatio)) + 
           2 * (t(thetaDiffHat.mu) %*% precision %*% thetaDiffHat.Sigma %*% 
                  precision %*% thetaDiffHat.mu))
    )
    nStar = deltaStar * sum(diag(covarRatio)) / tmpOmega
    lambdaStar = tmpOmega / deltaStar
    
    ndf = nStar
    ddf = distCXtSigmaInvXInvCt@df + a - 1
    omega = deltaStar
    
    Fscale = sum(diag(covarRatio))/((ddf) * a)
    
  }
  
  # Calculate the moments of the F under the null 
  mu.null = Fscale * (ddf / (ddf - 2))
  var.null = Fscale^2 * (2 * ddf^2 * (ndf + ddf - 2)) / (ndf * (ddf - 2)^2 * (ddf - 4))
  # Calculate the moments of the F under the alternative
  mu.alt = Fscale * (ddf * (ndf + omega)) / (ndf * (ddf - 2))
  var.alt = Fscale^2 * (2 * (ddf / ndf)^2 * 
                          ((ndf + omega)^2 + (ndf + 2 * omega) * (ddf - 2)) / 
                          ((ddf - 2)^2 * (ddf - 4)))
  
  return(data.frame(mu.null=mu.null, var.null=var.null, mu.alt=mu.alt, var.alt=var.alt))
}

#
# Calculate power for the KR test of fixed effects
# in the mixed model
#
mixedPower = function(design, glh, homoscedastic=TRUE) {
  
  # get the approximate moments of the Wald statistic
  moments = getWaldMoments(design, glh, homoscedastic)
  
  # calculate the parameters of the KR statistic
  params = getKRParams(nrow(glh@fixedContrast), moments)
  
  a = nrow(glh@fixedContrast)
  # get the critical F
  Fcrit = qf(1-glh@alpha,a,params$ddf)
#   cat("mu=(", moments$mu.null, ", ", moments$mu.alt, ")",
#       Fcrit, " ", params$lambda, "\n")
#   cat("params(l=", params$lambda, ", ddf=", params$ddf, 
#       ", omega=", params$omega, ")\n")
  # calculate power;
  power = 1 - pf(Fcrit, a, params$ddf, params$omega);
  
  return(power)
  
}


csMatrix = function(size, rho) {
  ones = matrix(rep(1,size))
  return(ones %*% t(ones) * rho + diag(size) * (1 - rho))
}


