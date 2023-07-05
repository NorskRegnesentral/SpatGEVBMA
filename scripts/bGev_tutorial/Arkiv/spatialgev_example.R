#install.packages("TMB")
install.packages("SpatialGEV")

library(TMB)
library(INLA)

library(SpatialGEV)
library(evd)


a <- simulatedData2$a
logb <- simulatedData2$logb
logs <- simulatedData2$logs
locs <- simulatedData2$locs
n_loc <- nrow(locs)
y <- Map(evd::rgev, n=sample(50:70, n_loc, replace=TRUE),
         loc=a, scale=exp(logb), shape=exp(logs))

filled.contour(unique(locs$x), unique(locs$y), matrix(a, ncol=sqrt(n_loc)), 
               color.palette = terrain.colors, xlab="Longitude", ylab="Latitude", 
               main="Spatial variation of a",
               cex.lab=1,cex.axis=1)
filled.contour(unique(locs$x), unique(locs$y), matrix(exp(logb), ncol=sqrt(n_loc)), 
               color.palette = terrain.colors, xlab="Longitude", ylab="Latitude", 
               main="Spatial variation of b",
               cex.lab=1,cex.axis=1)

filled.contour(unique(locs$x), unique(locs$y), matrix(exp(logs), ncol=sqrt(n_loc)), 
               color.palette = terrain.colors, xlab="Longitude", ylab="Latitude", 
               main="Spatial variation of s",
               cex.lab=1,cex.axis=1)


barplot(sapply(y, length), 
        xlab = "Location", ylab = "Number of observations at each location",
        main = "Summary of number of observations per location")


set.seed(123)
n_sam <-8
sam_inds <- sample(1:n_loc, n_sam, replace=FALSE)
par(mfrow=c(2, n_sam/2))
for (i in sam_inds){
  obs_i <- y[[i]]
  hist(obs_i, breaks=8,
       xlab="Observation value", main=paste("Observations at location", i))
}


fit <- spatialGEV_fit(y = y, locs = locs, random = "abs",
                      init_param = list(a = rep(60, n_loc),
                                        log_b = rep(2,n_loc),
                                        s = rep(-3,n_loc),
                                        beta_a = 60, beta_b = 2, beta_s = -2,
                                        log_sigma_a = 1.5, log_ell_a = 5,
                                        log_sigma_b = 1.5, log_ell_b = 5,
                                        log_sigma_s = -1, log_ell_s = 5),
                      reparam_s = "positive", kernel="exp", silent = T)

print(fit)

#posterior sampling:
sam <- spatialGEV_sample(model = fit, n_draw = 2000, observation = T)
print(sam)
pos_summary <- summary(sam)
pos_summary$param_summary[c(1:5, 
                            (n_loc+1):(n_loc+5), 
                            (2*n_loc+1):(2*n_loc+5),
                            (3*n_loc):(nrow(pos_summary$param_summary))),]



par(mfrow=c(1,3))
plot(a, pos_summary$param_summary[1:n_loc,"mean"], 
     xlab="True a",
     ylab="Posterior mean of a",
     main="True vs Posterior Mean of a")
abline(0, 1, col="blue", lty="dashed")
plot(exp(logb), exp(pos_summary$param_summary[(n_loc+1):(2*n_loc),"mean"]), 
     xlab="True b",
     ylab="Posterior mean of b",
     main="True vs Posterior Mean of b")
abline(0, 1, col="blue", lty="dashed")
plot(exp(logs), exp(pos_summary$param_summary[(2*n_loc+1):(3*n_loc),"mean"]), 
     xlab="True s",
     ylab="Posterior mean of s",
     main="True vs Posterior Mean of s")
abline(0, 1, col="blue", lty="dashed")


a <- simulatedData$a
logb <- simulatedData$logb
logs <- simulatedData$logs
y <- simulatedData$y
locs <- simulatedData$locs
n_loc <- nrow(locs)

filled.contour(unique(locs$x), unique(locs$y), matrix(a, ncol=sqrt(n_loc)), 
               color.palette = terrain.colors, xlab="Longitude", ylab="Latitude", main="Spatial variation of GEV location parameter a",
               cex.lab=1,cex.axis=1)

filled.contour(unique(locs$x), unique(locs$y), matrix(exp(logb), ncol=sqrt(n_loc)), 
               color.palette = terrain.colors, xlab="Longitude", ylab="Latitude", main="Spatial variation of GEV scale parameter b",
               cex.lab=1,cex.axis=1)


set.seed(123)
n_test <- 100
n_train <- n_loc - n_test
test_ind <- sample(1:n_loc, n_test)
train_ind <- setdiff(1:n_loc, test_ind)
# Obtain coordinate matrices and data lists
locs_test <- locs[test_ind,]
y_test <- y[test_ind]
locs_train <- locs[-test_ind,]
y_train <- y[-test_ind]


# Fit the GEV-GP model to the training set
train_fit_s <- spatialGEV_fit(y = y_train, locs = locs_train, random = "ab",
                              init_param = list(a = rep(0, n_train),
                                                log_b = rep(0,n_train), s = 0,
                                                beta_a = 0, beta_b = 0,
                                                log_sigma_a = 1, log_kappa_a = -2,
                                                log_sigma_b = 1, log_kappa_b = -2),
                              reparam_s = "positive", kernel="spde", silent = T)

train_fit_s

pred_s <- spatialGEV_predict(model = train_fit_s, locs_new = locs_test, n_draw = 2000)
pred_s
pred_summary <- summary(pred_s)
pred_summary[1:5,]


par(mfrow=c(1,2))
plot(sapply(y_test, mean), pred_summary[,"mean"], 
     xlab="Test statistic from test data",
     ylab="Test statistic from predictive distribution",
     main="Test statistic = mean")
abline(0, 1, col="blue", lty="dashed")
plot(sapply(y_test, function(x){quantile(x, probs=0.975)}), pred_summary[,"97.5%"], 
     xlab="Test statistic from test data",
     ylab="Test statistic from predictive distribution",
     main="Test statistic = 97.5% quantile")
abline(0, 1, col="blue", lty="dashed")


#-----------Case study--------------#
library(dplyr)

