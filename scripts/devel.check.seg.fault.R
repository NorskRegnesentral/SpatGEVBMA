rm(list = ls())

make.D <- function (x.1, x.2) 
{
    n.1 <- dim(x.1)[1]
    n.2 <- dim(x.2)[1]
    d <- dim(x.1)[2]
    D <- matrix(0, n.1, n.2)
    for (j in 1:n.2) {
        D[, j] <- sqrt(colSums((t(x.1) - x.2[j, ])^2))
    }
    return(D)
}

S <- matrix(rnorm(69*2),69,2)
alpha <- .59
tau <- rnorm(69)
lambda <- 1
D <- make.D(S,S)
C <- exp(-1/lambda * D)/alpha
##C <- matrix(0, 69,69)
diag(C) <- diag(C) + 1e-5
solve(C)

library(SpatGEVBMA)

solve(C)

library(SpatGEVBMA)

data(norway)
attach(norway)
a <- spatial.gev.bma(Y.list,X,as.matrix(S),2)
