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
library(magic)
library(invWishartSum)

#
# Class which describes a specific data pattern
# (either complete or with missing observations)
# for a mixed model design
#
setClass(
  "missingDataPattern",
  representation (
    group = "numeric",
    observations = "numeric",
    size = "numeric"
    ),
  prototype (
    group = 1,
    observations = c(1),
    size = 10
    )  
  ) 
  
#
# Mixed model design for use in power analysis
#
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
                new("missingDataPattern", group=1, observations=c(1,2,3), size=20),
                new("missingDataPattern", group=1, observations=c(1,2), size=20),
                new("missingDataPattern", group=2, observations=c(1,2,3), size=20),
                new("missingDataPattern", group=2, observations=c(1,2), size=15) 
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

#
# glh
#
# Class describing the general linear hypothesis
# for fixed effects in the mixed model
#
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
# Returns lambda, omega, and ddf for the distribution
# of the Kenward-Roger statistic
#
getKRParams = function(a, moments) {
  # calculate rho
  rho = moments$var.alt / (2*moments$mu.null^2) 
  # calculate noncentrality parameter
  omega = a*((moments$mu.alt/moments$mu.null)-1)
  # denominator degrees of freedom
  ddf = 4 + (2*(a + 2*omega) + (a + omega)^2)/(rho*a^2 - a - 2*omega)
  # scale factor for KR adjustment
  lambda = ddf / ((ddf - 2)*moments$mu.null)
  
  return(data.frame(lambda=lambda, ddf=ddf, omega=omega))
}

#
# Get the approximate moments of the Wald
# statistic given the mixed model design
#
#
getWaldMoments = function(design, glh, homoscedastic) {
  
  # determine the number of unique treatment groups
  numGroups = length(unique(sapply(design@xPatternList, function(x) {
    return(x@group)
  })))
  
  # get the max number of planned observations 
  maxObs = nrow(design@Sigma)
  
  # build the weight matrix for C*beta-thetaNull.  We weight by the
  # number of deletion classes containing the ith element of beta
  weightMatrix = diag(unlist(lapply(1:numGroups, function(group) {
    groupPatterns = Filter(function(x) { x@group==group }, design@xPatternList)
    return(sapply(1:maxObs, function(obs, groupPatterns) {
      return(1/sum(sapply(groupPatterns, 
                          function(pattern, obs) {
        return(as.numeric(obs %in% pattern@observations))
      }, obs)))
    }, groupPatterns))
  })))

  
  # calculate the approximate distribution of the stacked Sigma matrix
  # by estimating the distribution of the sum of quadratic 
  # forms in inverse Wishart matrices
  # 
  # first, build the list of Wisharts and corresponding design matrices
  invWishartList = list()
  designMatrixList = list()
  betaCovarSum = matrix(rep(0, (numGroups*maxObs)^2), nrow=(numGroups*maxObs))
  j = 1
  for(i in 1:length(design@xPatternList)) {
    pattern = design@xPatternList[[i]]
    # get the between/within portions of the design
    designBetween = matrix(diag(numGroups)[pattern@group,], nrow=1)
    designWithin = diag(maxObs)[pattern@observations,]
    
    # calculate Sigma for the deletion class
    SigmaD = designWithin %*% design@Sigma %*% t(designWithin)
    
    ##
    ## Setup the inputs to calculate the approximate Wishart for
    ## the middle term of the Wald statistic
    ##
    # build the design matrix - kronecker
    X = designBetween %x% designWithin
    # build the corresponding inverse Wisharts
    invWishart = invert(new("wishart", df=pattern@size, covariance=SigmaD ))
    # add to the list
    for(h in 1:pattern@size) {
      designMatrixList[[j]] = X
      invWishartList[[j]] = invWishart
      j = j+1
    }
    
    ##
    ## add in the next term for the beta covariance
    ##
    betaCovarSum = (betaCovarSum + 
                   ((1/pattern@size) * 
                      (t(designWithin) %*% SigmaD %*% designWithin) %x%
                      (t(designBetween) %*% designBetween)
                    ))
    
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
  thetaDiffHat.Sigma = (glh@fixedContrast %*% weightMatrix %*% 
                      betaCovarSum %*% weightMatrix %*% t(glh@fixedContrast))
  # calculate the mean
  thetaDiffHat.mu = glh@fixedContrast %*% design@beta - glh@thetaNull
                    
  qp = numGroups * maxObs
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
    tmpOmega = t(thetaDiffHat.mu) %*% solve(thetaDiffHat.Sigma) %*% thetaDiffHat.mu
    deltaStar = (
      2*((tmpOmega * sum(diag(covarRatio))) + 2 * tmpOmega^2) / 
        (sum(diag(covarRatio))^2 + 
           2 * (t(thetaDiffHat.mu) %*% precision %*% thetaDiffHat.Sigma %*% 
                  precision %*% thetaDiffHat.mu))
      )
    nStar = deltaStar * sum(diag(covarRatio)) / tmpOmega
    lambdaStar = tmpOmega / deltaStar
    
    ndf = nStar
    ddf = distCXtSigmaInvXInvCt@df
    omega = deltaStar

    Fscale = (lambdaStar * nStar)/((ddf) * a)
    
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
  
  # scale it down by lambda
  #Fcrit = params$lambda * Fcrit
  
  # calculate power;
  power = 1 - pf(Fcrit, a, params$ddf, params$omega);
  
  return(power)
  
}



scaleDesign = function(betaScale) {
  return (
    new("design.mixed",
        name ="3 treatments, 3 repeated measures, unbalanced",
        description ="",
        xPatternList = c(
          new("missingDataPattern", group=1, observations=c(1,2,3), size=20),
          new("missingDataPattern", group=1, observations=c(1,2), size=20),
          new("missingDataPattern", group=2, observations=c(1,2,3), size=20),
          new("missingDataPattern", group=2, observations=c(1,2), size=15),
          new("missingDataPattern", group=3, observations=c(1,2,3), size=30),
          new("missingDataPattern", group=3, observations=c(1,2), size=15) 
        ),
        beta = betaScale*matrix(c(1,1,1,0,0,0,0,0,0),ncol=1),
        Sigma = matrix(c(1,0.3,0.3,0.3,1,0.3,0.3,0.3,1), nrow=3)
    )
    )
}

## 
mixed1 = new("design.mixed",
             name ="3 treatments, 3 repeated measures, unbalanced",
             description ="",
             xPatternList = c(
               new("missingDataPattern", group=1, observations=c(1,2,3), size=20),
               new("missingDataPattern", group=1, observations=c(1,2), size=20),
               new("missingDataPattern", group=2, observations=c(1,2,3), size=20),
               new("missingDataPattern", group=2, observations=c(1,2), size=15),
               new("missingDataPattern", group=3, observations=c(1,2,3), size=30),
               new("missingDataPattern", group=3, observations=c(1,2), size=15) 
             ),
             beta = 1*matrix(c(1,1,1,0,0,0,0,0,0),ncol=1),
             Sigma = matrix(c(1,0.3,0.3,0.3,1,0.3,0.3,0.3,1), nrow=3)
             )

mixed.balanced = new("design.mixed",
             name ="3 treatments, 3 repeated measures, balanced",
             description ="",
             xPatternList = c(
               new("missingDataPattern", group=1, observations=c(1,2,3), size=400),
               new("missingDataPattern", group=2, observations=c(1,2,3), size=400),
               new("missingDataPattern", group=3, observations=c(1,2,3), size=400)
             ),
             beta = matrix(c(1,1,1,0,0,0,0,0,0),ncol=1),
             Sigma = matrix(c(1,0.3,0.3,0.3,1,0.3,0.3,0.3,1), nrow=3)
)

csMatrix = function(size, rho) {
  ones = matrix(rep(1,size))
  return(ones %*% t(ones) * rho + diag(size) * (1 - rho))
}


glh1= new("glh.mixed",
          alpha = 0.05,
          fixedContrast = matrix(c(1/3,1/3,1/3,-1/3,-1/3,-1/3,0,0,0,
                                   1/3,1/3,1/3,0,0,0,-1/3,-1/3,-1/3), nrow=2, byrow=TRUE),
          thetaNull = matrix(c(0,0), nrow=2),
          test = "Wald, KR ddf"
          )

glh2= new("glh.mixed",
          alpha = 0.05,
          fixedContrast = matrix(c(1/3,1/3,1/3,-1/3,-1/3,-1/3,0,0,0), nrow=1, byrow=TRUE),
          thetaNull = matrix(c(0), nrow=1),
          test = "Wald, KR ddf"
)

