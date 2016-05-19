A <- list()

A[[1]] <- matrix(1:6,3,2)
A[[2]] <- matrix(1:6 + 6,3,2)
A[[3]] <- matrix(1:6 + 12,3,2)
A[[4]] <- matrix(1:6 + 18,3,2)

B <- array(data = unlist(A), dim = c(3,2,4))
