

#
# Asymptotic behavior for Wald with two-sample T
#
X = diag(2)
X

C = matrix(c(1,-1),nrow=1)
C

B = matrix(c(1,0),nrow=2)
B

thetaObs = C%*%B
thetaObs

thetaNull = matrix(c(0),nrow=1)
thetaNull

thetaDiff = C%*%B-thetaNull

sigmaSq = matrix(c(2))
sigmaSqInv = 0.5

calcWald <- function(N) {
  Xfull = X%x%matrix(rep(1,N))
  sigInvFull = sigmaSqInv*diag(nrow(Xfull))
  Xpx = t(Xfull)%*%sigInvFull%*%Xfull
  F = t(thetaDiff)%*%solve(C%*%solve(Xpx)%*%t(C))%*%thetaDiff*(1/nrow(C))
  return(F)                       
                        
}

N = 2:400
Y = lapply(N,calcWald)

plot(N,unlist(Y))
lines(N,N)

#
# Anova
#
X = diag(4)
X

C = matrix(c(1,-1,0,0,1,0,-1,0,1,0,0,-1),nrow=3,byrow=TRUE)
C

B = matrix(c(1,0,0,0),nrow=4)
B

thetaObs = C%*%B
thetaObs

thetaNull = matrix(c(0,0,0),nrow=3)
thetaNull

thetaDiff = C%*%B-thetaNull

sigmaSq = matrix(c(2))
sigmaSqInv = 0.5

N = 2:400
Y = lapply(N,calcWald)

plot(N,unlist(Y))
lines(N,N)


#
# Where does rho go?
#
r = 4
a = nrow(C)
N = 10:400
Y = lapply(N,calcWald)
Fobs = unlist(Y)
rho = (a*(1+Fobs)^2 + (1+2*Fobs)*(N-r-2)) / (a*(N-r-4))
plot(N, Fobs)
plot(N,rho)

#
# Where does nu_p go?
#
nu_p = 4 + (2*(1+2*Fobs)+a*(1+Fobs)^2)/(rho*a-1-2*Fobs)

plot(N,nu_p,"l")
lines(N,N,col="green")