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
                new("missingDataPattern", group=1, observations=c(1,2,3), size=10),
                new("missingDataPattern", group=1, observations=c(1,2), size=10),
                new("missingDataPattern", group=2, observations=c(1,2,3), size=10),
                new("missingDataPattern", group=2, observations=c(1,2), size=5)   
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
    if (nrow(object@betweenContrast) != nrow(object@thetaNull)) {
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
getKRParams = function(a, waldMeanNull, waldMeanAlt, waldVarAlt) {
  # calculate rho
  rho = waldVarAlt / (2*waldMeanNull^2) 
  # calculate noncentrality parameter
  omega = a*((waldMeanAlt/waldMeanNull)-1)
  # denominator degrees of freedom
  ddf = 4 + (2*(a + 2*omega) + (a + omega)^2)/(rho*a^2 - a - 2*omega)
  # scale factor for KR adjustment
  lambda = ddf / ((ddf - 2)*waldMeanNull)
  
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
    designMatrixList[[i]] = designBetween %x% designWithin
    
    # build the corresponding Wisharts
    wishart = new("wishart", df=pattern@size, covariance=SigmaD )
    invWishartList[[i]] = invert(wishart)
    
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
  thetaDiffCovar = (glh@fixedContrast %*% weightMatrix %*% 
                      betaCovarSum %*% weightMatrix %*% t(glh@fixedContrast))
  # calculate the mean
  thetaDiffMu = glh@fixedContrast %*% design@beta - glh@thetaNull
                    

  thetaDiff
  

  # Form the parameters of the approximate F
  #
  ndf = nrow(glh@fixedContrast)
  ddf = 16  
  omega=5

  # Calculate the moments of the F under the null 
  mu.null = ddf / (ddf - 2)
  var.null = (2 * ddf^2 * (ndf + ddf - 2)) / (ndf * (ddf - 2)^2 * (ddf - 4))
  # Calculate the moments of the F under the alternative
  mu.alt = (ddf * (ndf + omega)) / (ndf * (ddf - 2))
  var.alt = (2 * (ddf / ndf)^2 * 
               ((ndf + omega)^2 + (ndf + 2 * omega) * (ddf - 2)) / 
               ((ddf - 2)^2 * (ddf - 4)))
    
    
  
  return(data.frame(mu.null=mu.null, var.null=var.null, muAlt=muAlt, varAlt=varAlt))
}

#
# Calculate power for the KR test of fixed effects
# in the mixed model
#
mixedPower = function(design, glh, homoscedastic=TRUE) {
  
  # get the approximate moments of the Wald statistic
  moments = getWaldMoments(design, glh, homoscedastic)
  
  # calculate the parameters of the KR statistic
  params = getKRParams(a, moments$muNull, moments$muAlt, moments$varAlt)
  
  # get Wald statistic
  wald
  
  # get the critical F
  Fcrit = qf(1-alpha,a,params$ddf)
  
  # scale it down by lambda
  Fcrit = params$lambda * Fcrit
  
  # calculate power;
  power = 1 - pf(Fcrit, a, params$ddf, params$omega);
  
  return(power)
  
}

test = data.frame(treatment=c(1,1,1,2,2,2))
test$N = c(5,5,12,13,4,5) 
test$observations = list(c(1), c(2), c(1,2,3), c(1,2), c(1,2,3), c(1))

, 
                  observations=list(c(1), c(2), c(1,2,3), c(1,2), c(1,2,3), c(1)),
                  N=c(5,5,12,13,4,5))


#
# getObservedF
#
# Get the Wald statistic that would be obtained if
# we observed Beta and SigmaS exactly as specified
#
# Arguments:
#  X - design matrix for expected value
#  Beta - choices for parameters related to the mean
#  C - linear contrast for fixed effects
#  SigmaS - covariance structure for the complete model
#  ThetaNull - matrix of null hypothesis
#
# Returns:
#  observed F value
#
getObservedF(X, Beta, C, SigmaS, ThetaNull);
# get the row dimension of C;
a = NROW(C);

# calculate the middle matrix in the Wald statistic;
M = C*INV(X`*INV(SigmaS)*X)*C`;
          
          # calculate the obsered theta and difference from the null hypothesis;
          thetaObs = C*Beta;
          thetaObsDiff = thetaObs - ThetaNull;
          
          Fobserved = thetaObsDiff`*INV(M)*thetaObsDiff / a;
          
          return(Fobserved);
          finish;
          
          /#
            #
            # Calculate the denominator degrees of freedom for
            # the approximate noncentral F used for power
            #
            # Arguments
            #  Fobserved - f statistic that would be obtained if we
            #              observed the Beta and SigmaS as specified
            #  a - numerator degrees of freedom, also the rows of C
            #  nuW - the denominator degrees of freedom for the true
            #        distribution of the Wald test (power will vary
            #        depending on the choice of this value)
          #
          # Returns
          #  ddf - the denominator degrees of freedom for the
          #        approximate noncentral F
          #
          #/
          start getDenominatorDegreesOfFreedom(Fobserved, a, nuW);
          
          # calculate rho, the ratio of the variance of the true noncentral Wald distribution;
          # to 2 times the square of the expectation of the true central Wald distribution;
          rhoNumerator = a*(1+Fobserved)*(1+Fobserved) + (1+(2*Fobserved))*(nuW - 2);
          rho = rhoNumerator / (a*(nuW-4));
          
          # calculate the numerator for ddf ;
          numerator = 2*(1+2*Fobserved)+ a*(1+Fobserved)*(1+Fobserved);
          denominator = rho*a-1-2*Fobserved;
          
          # calculate the degrees of freedom for the approximate F;
          ddf = 4 + (numerator / denominator);
          return (ddf);
          
          finish;
          
          /#
            # calculatePowerKenwardRoger
            #
            # Calculate power for the linear mixed model Wald test
            # with denominator degrees of freedom as described by 
            # Kenward and Roger (1997)
            #
            # Implements the noncentral F approximation developed
            # by Kreidler et al.
            #
            # Arguments:
            #  X - design matrix for expected value
          #  Beta - choices for parameters related to the mean
          #  C - linear contrast for fixed effects
          #  SigmaS - covariance structure for the complete model
          #  ThetaNull - matrix of null hypothesis
          #  alpha - desired Type I error rate
          #
          # Returns:
          #  power
          #/
          start calculatePowerKenwardRoger(X, Beta, C, SigmaS, ThetaNull, alpha, numISU);
          
          # calculate the total sample size;
          N = numISU;
          
          # calculate the rank of X;
          rankX=round(trace(ginv(X)*X));
          
          # choose value for true ddf;
          ddfTrue = N - rankX; 
          
          # get the F we would obtain were we to observe Beta and SigmaS;
          Fobserved = getObservedF(X, Beta, C, SigmaS, ThetaNull);
          
          # calculate numerator degrees of freedom;
          ndf = NROW(C);
          
          # calculate the denominator degrees of freedom;
          ddf = getDenominatorDegreesOfFreedom(Fobserved, ndf, ddfTrue);
          
          # get the critical F;
          Fcrit = FINV(1-alpha,ndf,ddf);
          
          # calculate lambda;
          lambda = (ddf*(ddfTrue-2))/(ddfTrue*(ddf-2));
          # scale F crit;
          Fcrit = Fcrit * lambda;
          
          # calculate the noncentrality;
          omegaP = ndf * Fobserved;
          
          # calculate power;
          power = 1 - CDF("F", FCRIT, ndf, ddf, omegaP);
          
          # build the results, including F distribtuion info and power;
          results = Fcrit || Fobserved || ndf || ddf || omegaP || power;
          
          # calculate power;
          return (results);
          
          finish;
          
          
          
          
          
          
          
          
          
          ###### Test case ######
          
          Xessence = diag(2)
          
          X = Xessence %x% matrix(rep(1,35));
          
          ISU = matrix(
            c(rep(1,5), rep(2,10), rep(3,10), rep(4,5), rep(5,5),
              rep(6,5), rep(7,10), rep(8,10), rep(9,5), rep(10,5)))
          
          X = cbind(ISU, X)
          
          names(X) = c("clusterId", "A", "B")
          
          Beta = matrix(c(1,0), nrow=1)
          
          SigmaSmall = 2*matrix(rep(1,25),nrow=5)*0.3+diag(5)*0.7;
          SigmaBig = 2*matrix(rep(1,100),nrow=10)*0.3+diag(10)*0.7;
          
          SigmaS = adiag(SigmaSmall,SigmaBig,SigmaBig,SigmaSmall,SigmaSmall,
                         SigmaSmall,SigmaBig,SigmaBig,SigmaSmall,SigmaSmall)
          
          C = matrix(c(1 -1), nrow=1);
          thetaNull = matrix(c(0));
          alpha = 0.05;
          