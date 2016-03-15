rm(list = ls())

library("SpatialGEVBMA")

l.dot <- function (tau, alpha, lambda, D, a, b) 
  {
      res <- c(NA,NA)
      E.l <- exp(-1/lambda * D)
      diag(E.l) <- diag(E.l) + 1e-05
      E.inv <- solve(E.l)
      F.l <- 1/lambda^2 * D * E.l
      M.l <- E.inv %*% (-F.l) %*% E.inv
      res[1] <- -0.5 * sum(diag(E.inv %*% F.l)) - 0.5 * alpha * t(tau) %*%  M.l %*% tau - b + (a - 1)/lambda
      G.l <- -2/lambda^3 * (D * E.l) + 1/lambda^2 * (D * F.l)
      L.l <- M.l %*% F.l + E.inv %*% G.l
      N.l <- M.l %*% (-F.l) %*% E.inv +
          E.inv %*% (-G.l) %*% E.inv +
              E.inv %*% (-F.l) %*% M.l
      res[2] <- -0.5 * sum(diag(L.l)) - 0.5 * alpha * t(tau) %*% N.l %*% tau - (a - 1)*lambda^(-2)
      return(res)
  }

data(norway)
attach(norway)

set.seed(1)
G <- gev.init(Y.list, X, S, NULL, FALSE, NULL)
for(i in 1:5)G <- gev.update(G)

##------ Unload --------
tau <- G$tau.mu
alpha <- G$alpha.mu
lambda <- G$lambda.mu
D <- G$D
a <- G$prior$mu$lambda.a
b <- G$prior$mu$lambda.b

E <- exp(-D/lambda)
diag(E) <- diag(E) + 1e-5
E.inv <- solve(E)
F <- 1/lambda^2 * D * E
M <- E.inv %*% (-F) %*% E.inv
G <- -2/lambda^3 * (D * E) + 1/lambda^2 * (D * F)
L <- M %*% F + E.inv %*% G
N <- M %*% (-F) %*% E.inv +
    E.inv %*% (-G) %*% E.inv +
        E.inv %*% (-F) %*% M

for(i in 1:1e1)
    {
        A.1 <- l.dot(tau, alpha, lambda, D, a, b)
        A <- ldot(tau, alpha, lambda, D, a, b)
  }

l.prime(tau, alpha, lambda, D, a, b)
l.double.prime(tau, alpha, lambda, D, a, b)



##----------------------

### Timing:

ptm=proc.time()
for (i in 1:(2*10^3))
{
  A.11=l.prime(tau, alpha, lambda, D, a, b)
  A.12=l.double.prime(tau, alpha, lambda, D, a, b)
  }
proc.time()-ptm
#user  system elapsed 
#10.896   0.012  11.206

ptm=proc.time()
for (i in 1:(2*10^3))
{
  A.1=l.dot(tau, alpha, lambda, D, a, b)
}
proc.time()-ptm
#user  system elapsed 
#8.44    0.00    8.72

ptm=proc.time()
for (i in 1:(2*10^3))
{
  A <- ldot(tau, alpha, lambda, D, a, b)
}
proc.time()-ptm
#user  system elapsed 
#6.356   0.012   6.762 

