#
# Calculate power using the Stroup exemplary data method
#
stroupPower <- function(X, N, beta, sigma, C, alpha) {
  ndf = nrow(C)
  CB = C %*% b
  Fvalue = t(CB) %*% solve(C %*% solve(t(X) %*% solve(sigma) %*% X) %*% t(C)) %*% CB / ndf
  
  ddf = N - ncol(X)
  noncentrality = ndf*Fvalue
  FCrit=qf(1 - alpha, ndf, ddf, 0)
  Power=1-pf(FCrit,ndf,ddf,noncentrality)
}

# build design
alpha = 0.05
b = matrix(c(20,25,25),nrow=3)
sigma = 9 * diag(15)
ones = matrix(rep(1,5),nrow=5)
X = diag(3) %x% ones
# power holder
powerList = data.frame(contrast=character(), power=numeric())
names(powerList) = c("contrast", "power")
# omnibus test
C = matrix(c(1,-1,0,
             1,0,-1), nrow=2, byrow=TRUE)
contrastName = "Omnibus"
power = stroupPower(X, nrow(X), b, sigma, C, alpha)
powerList = rbind(powerList, data.frame(contrast=contrastName, power=power))

# ctrl vs. exp
C = matrix(c(2,-1,-1), nrow=1, byrow=TRUE)
contrastName = "ctrl vs. trt"
power = stroupPower(X, nrow(X), b, sigma, C, alpha)
powerList = rbind(powerList, data.frame(contrast=contrastName, power = power))

# pairwise
C = matrix(c(1,-1,0), nrow=1, byrow=TRUE)
contrastName = "ctrl vs. trt 1"
power = stroupPower(X, nrow(X), b, sigma, C, alpha)
powerList = rbind(powerList, data.frame(contrast=contrastName, power = power))

C = matrix(c(1,0,-1), nrow=1, byrow=TRUE)
contrastName = "ctrl vs. trt 2"
power = stroupPower(X, nrow(X), b, sigma, C, alpha)
powerList = rbind(powerList, data.frame(contrast=contrastName, power = power))

C = matrix(c(0,1,-1), nrow=1, byrow=TRUE)
contrastName = "trt 1 vs. trt 2"
power = stroupPower(X, nrow(X), b, sigma, C, alpha)
powerList = rbind(powerList, data.frame(contrast=contrastName, power = power))

