#
#
#
source("../R/mixedPower.R")
source("../R/empiricalPower.R")

#
# Get the approximate moments of the Wald
# statistic given the mixed model design
#
#
getWaldMoments2 = function(design, glh, invWishartList, designMatrixList, X, SigmaS,
                           fiddle=1) {
  
  # determine the number of unique treatment groups
  numGroups = length(unique(sapply(design@xPatternList, function(x) {
    return(x@group)
  })))
  
  # get the max number of planned observations 
  maxObs = nrow(design@Sigma)
  
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
  thetaDiffHat.Sigma = (glh@fixedContrast %*%
                          solve(t(as.matrix(X)) %*% solve(SigmaS) 
                                %*% as.matrix(X)) %*%
                          t(glh@fixedContrast))*fiddle
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
    tmpOmega = t(thetaDiffHat.mu) %*% precision %*% thetaDiffHat.mu
    deltaStar = (
      ((tmpOmega * sum(diag(covarRatio))) + 2 * tmpOmega^2) / 
        (sum(diag(covarRatio))^2 + 
           2 * (t(thetaDiffHat.mu) %*% precision %*% thetaDiffHat.Sigma %*% 
                  precision %*% thetaDiffHat.mu))
    )
    nStar = deltaStar * sum(diag(covarRatio)) / tmpOmega
    lambdaStar = tmpOmega / deltaStar
    
    ndf = nStar
    ddf = distCXtSigmaInvXInvCt@df
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
mixedPower2 = function(design, glh, invWishartList, designMatrixList, X, SigmaS,
                       fiddle=1) {
  
  # get the approximate moments of the Wald statistic
  moments = getWaldMoments2(design, glh, invWishartList, designMatrixList, X, SigmaS,
                            fiddle)
  
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



generateClusterDesign = function(numGroups, perGroupN, clusterN.full, clusterN.deleted,
                                  groupDiff=1, icc=0.04, sigmaSq=1,
                                  name="", description="") {
  
  print(numGroups)
  
  patternList = list()
  pattern = 1
  for(grp in 1:numGroups) {
    patternList[pattern] = 
      new("missingDataPattern", group=grp, observations=1:clusterN.full, size=perGroupN)
    patternList[pattern+1] = 
      new("missingDataPattern", group=grp, observations=1:clusterN.deleted, size=perGroupN)
    pattern = pattern + 2
  }

  cluster = matrix(rep(1,clusterN.full))
  return(new("design.mixed", name = name, description = description,
                 xPatternList = patternList,
                 beta = matrix(c(groupDiff,rep(0,numGroups-1))),
                 Sigma = sigmaSq * (icc * (cluster %*% t(cluster)) + 
                                      diag(clusterN.full) * (1 - icc))
             ))

}

generateClusterHypothesis = function(design, alpha=0.05, thetaNull=NULL) {
  numGroups = length(unique(sapply(design@xPatternList, function(x) { return(x@group) })))
  
  if (is.null(thetaNull)) {
    thetaNull = matrix(rep(0,numGroups-1))
  }

  return(new("glh.mixed",
            alpha = alpha,
            fixedContrast = cbind(matrix(rep(1,numGroups-1)), -1*diag(numGroups-1)),
            thetaNull = thetaNull,
            test = "Wald, KR ddf"
  ))
}

generateEmpiricalPowerCode = function(design, glh, filename, dataPrefix) {
  numGroups = length(unique(sapply(design@xPatternList, function(x) { return(x@group) })))
  
  mixedCall = paste(c(
    "proc mixed data=&datasetName;\n",
    "model y = trt1 trt2 trt / noint solution ddfm=KR;",
    "random int / subject=clusterID;",
    "by setID;",
    "contrast \"trt\" trt1 1 trt2 -1, ;",
    "run;"), collapse="\n")

  genSASCode(filename, design, glh, mixedCall,
           dataPrefix)
}


cluster1 = generateClusterDesign(2, 10, 5, 3, groupDiff=1, icc=0.04, sigmaSq=1,
                                            name="cluster, 2grp", description="")
glh1 = generateClusterHypothesis(cluster1, alpha=0.05) 
generateEmpiricalPowerCode(cluster1, glh1, "clusterRandomized1.sas", "clusterEx1")



#
# Terribly programmed example
#
Xessence = diag(2)
X = data.frame(Xessence %x% matrix(rep(1,40)))
names(X) = c("A", "B")
X$id = c(rep(1,10), rep(2,10), rep(3,10), rep(4,5), rep(5,5),
         rep(6,10), rep(7,10), rep(8,10), rep(9,5), rep(10, 5))

SigmaSmall = 2*(matrix(rep(1,25), nrow=5)*0.3+diag(5)*0.7)
SigmaBig = 2*(matrix(rep(1,100), nrow=10)*0.3+diag(10)*0.7)
SigmaS = as.matrix(bdiag(SigmaBig, SigmaBig, SigmaBig, SigmaSmall, SigmaSmall,
                         SigmaBig, SigmaBig, SigmaBig, SigmaSmall, SigmaSmall))
w5 = new("wishart", df=6, covariance=SigmaSmall)
w10 = new("wishart", df=5, covariance=SigmaBig)
iw5 = invert(w5)
iw10 = invert(w10)

invWishartList = c(iw10, iw10, iw10, iw5, iw5,
                   iw10, iw10, iw10, iw5, iw5)

scaleMatrixList = list(
  matrix(c(rep(1,10),rep(0,10)), nrow=10),
  matrix(c(rep(1,10),rep(0,10)), nrow=10),
  matrix(c(rep(1,10),rep(0,10)), nrow=10),
  matrix(c(rep(1,5),rep(0,5)), nrow=5),
  matrix(c(rep(1,5),rep(0,5)), nrow=5),
  matrix(c(rep(0,10),rep(1,10)), nrow=10),
  matrix(c(rep(0,10),rep(1,10)), nrow=10),
  matrix(c(rep(0,10),rep(1,10)), nrow=10),
  matrix(c(rep(0,5),rep(1,5)), nrow=5),
  matrix(c(rep(0,5),rep(1,5)), nrow=5) 
  )


mDist = approximateInverseWishart(invWishartList, 
                                            scaleMatrixList,
                                            method="trace")

C = matrix(c(1,-1),nrow=1)
beta=(matrix(c(1,0), nrow=2))
# now invert to get a Wishart
distXtSigmaInvXInv = invert(mDist)
# scale by C matrices and form corresponding wishart
distCXtSigmaInvXInvCt = new("wishart", df=distXtSigmaInvXInv@df,
                            covariance=(C %*% 
                                          distXtSigmaInvXInv@covariance %*% 
                                          t(C))) 


# calculate the approximate covariance of thetaDiff = C*beta-thetaNull
thetaDiffHat.Sigma = (glh@fixedContrast %*% weightMatrix %*% 
                        betaCovarSum %*% weightMatrix %*% t(glh@fixedContrast))

thetaDiffHat.Sigma = (C %*% 
                        solve(t(as.matrix(X[,1:2])) %*% solve(SigmaS) 
                              %*% as.matrix(X[,1:2])) %*%
                        t(C))

thetaDiffHat.Sigma = distCXtSigmaInvXInvCt@covariance

# calculate the mean
thetaDiffHat.mu = C %*% beta

# special case when contrast has a single row
ndf = 1
ddf = distCXtSigmaInvXInvCt@df
omega = t(thetaDiffHat.mu) %*% solve(thetaDiffHat.Sigma) %*% thetaDiffHat.mu
Fscale = ((1/(ddf)) * 
            (as.numeric(thetaDiffHat.Sigma) / 
               as.numeric(distCXtSigmaInvXInvCt@covariance)))
#Fscale = 1

# Calculate the moments of the F under the null 
mu.null = Fscale * (ddf / (ddf - 2))
var.null = Fscale^2 * (2 * ddf^2 * (ndf + ddf - 2)) / (ndf * (ddf - 2)^2 * (ddf - 4))
# Calculate the moments of the F under the alternative
mu.alt = Fscale * (ddf * (ndf + omega)) / (ndf * (ddf - 2))
var.alt = Fscale^2 * (2 * (ddf / ndf)^2 * 
                        ((ndf + omega)^2 + (ndf + 2 * omega) * (ddf - 2)) / 
                        ((ddf - 2)^2 * (ddf - 4)))


params = getKRParams(nrow(C), 
                     data.frame(mu.null=mu.null, var.null=var.null, 
                                mu.alt=mu.alt, var.alt=var.alt))



# get the critical F
Fcrit = qf(1-0.05,1,params$ddf)

# calculate power;
power = 1 - pf(Fcrit, 1, params$ddf, params$omega);



print X;
Beta = {1,0};

SigmaSmall = 2#(J(5,5,1)*0.3+I(5)*0.7);
SigmaBig = 2#(J(10,10,1)*0.3+I(10)*0.7);
SigmaS = block(SigmaSmall,SigmaBig,SigmaBig,SigmaSmall,SigmaSmall,
               SigmaSmall,SigmaBig,SigmaBig,SigmaSmall,SigmaSmall);
print SigmaS;
C = {1 -1};
thetaNull = {0};
alpha = {0.05};



#
# a = 1, balanced
#
Xessence = diag(2)
X = data.frame(Xessence %x% matrix(rep(1,50)))
names(X) = c("A", "B")
X$id = c(sapply(1:20, function(x) { return(rep(x,5)) }))

SigmaSmall = 2*(matrix(rep(1,25), nrow=5)*0.3+diag(5)*0.7)
SigmaBig = 2*(matrix(rep(1,100), nrow=10)*0.3+diag(10)*0.7)
SigmaS = diag(20) %x% SigmaSmall

iw5 = invert(new("wishart", df=20, covariance=SigmaSmall))

invWishartList = list()
scaleMatrixList = list()
for(i in 1:20) {
  invWishartList[[i]] = iw5
  if (i <= 10) {
    scaleMatrixList[[i]] = matrix(c(rep(1,5),rep(0,5)), nrow=5)
    
  } else {
    scaleMatrixList[[i]] = matrix(c(rep(0,5),rep(1,5)), nrow=5)
    
  }
}
class(invWishartList) = "inverseWishart"

design.a1bal = new("design.mixed",
                   name ="a=1, balanced",
                   description ="",
                   xPatternList = c(
                     new("missingDataPattern", group=1, observations=c(1:5), size=10),
                     new("missingDataPattern", group=2, observations=c(1:5), size=10)
                   ),
                   beta = matrix(c(1,0), nrow=2),
                   Sigma = matrix(c(1,0.3,0.3,0.3,1,0.3,0.3,0.3,1), nrow=3)
)
glh.a1bal = new("glh.mixed",
            alpha = 0.05,
            fixedContrast = matrix(c(1,-1),nrow=1),
            thetaNull = matrix(c(0), nrow=1),
            test = "Wald, KR ddf")

mixedPower2(design.a1bal, glh.a1bal, invWishartList, 
            scaleMatrixList, X[,1:2], SigmaS, fiddle=1.07)


#
# a > 1, balanced
#
Xessence = diag(3)
X = data.frame(Xessence %x% matrix(rep(1,50)))
names(X) = c("A", "B", "C")
X$id = c(sapply(1:30, function(x) { return(rep(x,5)) }))

SigmaSmall = 2*csMatrix(5, 0.3)
SigmaS = diag(30) %x% SigmaSmall

iw5 = invert(new("wishart", df=30, covariance=SigmaSmall))

invWishartList = list()
scaleMatrixList = list()
for(i in 1:30) {
  invWishartList[[i]] = iw5
  if (i <= 10) {
    scaleMatrixList[[i]] = matrix(c(rep(1,5),rep(0,10)), nrow=5)
    
  } else if (i <= 20) {
    scaleMatrixList[[i]] = matrix(c(rep(0,5),rep(1,5), rep(0,5)), nrow=5)
    
  } else {
    scaleMatrixList[[i]] = matrix(c(rep(0,10),rep(1,5)), nrow=5)
    
  }
}
class(invWishartList) = "inverseWishart"

design.a2bal = new("design.mixed",
                   name ="a>1, balanced",
                   description ="",
                   xPatternList = c(
                     new("missingDataPattern", group=1, observations=c(1:5), size=10),
                     new("missingDataPattern", group=2, observations=c(1:5), size=10),
                     new("missingDataPattern", group=3, observations=c(1:5), size=10)
                   ),
                   beta = matrix(c(1,0,0), nrow=3),
                   Sigma = 2*csMatrix(5,0.3)
)
glh.a2bal = new("glh.mixed",
                alpha = 0.05,
                fixedContrast = matrix(c(1,-1, 0,1,0,-1),nrow=2, byrow=TRUE),
                thetaNull = matrix(c(0,0), nrow=2),
                test = "Wald, KR ddf")

mixedPower2(design.a2bal, glh.a2bal, invWishartList, 
            scaleMatrixList, X[,1:3], SigmaS)


# determine the number of unique treatment groups
numGroups = length(unique(sapply(design.a2bal@xPatternList, function(x) {
  return(x@group)
})))

# get the max number of planned observations 
maxObs = nrow(design.a2bal@Sigma)

# calculate the approximating inverse Wishart for X' * SigmaInv * X
distXtSigmaInvX = approximateInverseWishart(invWishartList, 
                                            scaleMatrixList,
                                            method="trace")
# now invert to get a Wishart
distXtSigmaInvXInv = invert(distXtSigmaInvX)
# scale by C matrices and form corresponding wishart
distCXtSigmaInvXInvCt = new("wishart", df=distXtSigmaInvXInv@df,
                            covariance=(glh.a2bal@fixedContrast %*% 
                                          distXtSigmaInvXInv@covariance %*% 
                                          t(glh.a2bal@fixedContrast))) 

# calculate the approximate covariance of thetaDiff = C*beta-thetaNull
thetaDiffHat.Sigma = (glh.a2bal@fixedContrast %*%
                        solve(t(as.matrix(X[,1:3])) %*% solve(SigmaS) 
                              %*% as.matrix(X[,1:3])) %*%
                        t(glh.a2bal@fixedContrast))
# calculate the mean
thetaDiffHat.mu = glh.a2bal@fixedContrast %*% design.a2bal@beta - glh.a2bal@thetaNull

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
  a = nrow(glh.a2bal@fixedContrast)
  precision = solve(distCXtSigmaInvXInvCt@covariance)
  covarRatio = precision %*% thetaDiffHat.Sigma
  tmpOmega = t(thetaDiffHat.mu) %*% precision %*% thetaDiffHat.mu
  deltaStar = (
    ((tmpOmega * sum(diag(covarRatio))) + 2 * tmpOmega^2) / 
      (sum(diag(covarRatio))^2 + 
         2 * (t(thetaDiffHat.mu) %*% precision %*% thetaDiffHat.Sigma %*% 
                precision %*% thetaDiffHat.mu))
  )
  nStar = deltaStar * sum(diag(covarRatio)) / tmpOmega
  lambdaStar = tmpOmega / deltaStar
  
  ndf = nStar
  ddf = distCXtSigmaInvXInvCt@df
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



#
# a > 1, unbalanced
#
design = mixed1
X = getXStacked(design)
SigmaS = getSigmaStacked(design)
SigmaSmall = diag(3)[1:2,] %*% design@Sigma %*% t(diag(3)[1:2,])
mixed1@beta = mixed1@beta *0.5
# determine the number of unique treatment groups
numGroups = length(unique(sapply(design@xPatternList, function(x) {
  return(x@group)
})))

# get the max number of planned observations 
maxObs = nrow(design@Sigma)

iwFull = invert(new("wishart", df=70, covariance=design@Sigma))
iwSmall = invert(new("wishart", df=50, covariance=SigmaSmall))

invWishartList = list()
scaleMatrixList = list()
isu = 1
for(pattern in design@xPatternList) {
  for(i in 1:pattern@size) {
    designBetween = matrix(diag(numGroups)[pattern@group,], nrow=1)
    designWithin = diag(maxObs)[pattern@observations,]
    
    scaleMatrixList[[isu]] = designBetween %x% designWithin

    if (length(pattern@observations) == 3) {
      invWishartList[[isu]] = iwFull
    } else {
      invWishartList[[isu]] = iwSmall
    }
    isu = isu + 1
  }
}
class(invWishartList) = "inverseWishart"

glh.test = glh1
glh.test@fixedContrast = matrix(glh1@fixedContrast[1,], nrow=1)
glh.test@thetaNull = matrix(glh1@thetaNull[1,], nrow=1)
mixedPower2(mixed1, glh.test, invWishartList, 
            scaleMatrixList, X[,2:ncol(X)], SigmaS)









