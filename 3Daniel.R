#Sheet 3: Project high-dimensional regression: Monte Carlo comparison

#####################################################################################################################


#task 1:
#--(a) Implement a function in R which computes the ridge regression estimator
#--given the observations of the target and design values and
#--with a general tuning parameter ?? for the penalization to be chosen as an input.
#--You shall not use packages that readily include the ridge estimator,
#--but you can use the function optim, or similar ones,
#--to solve optimization problems (similar as for least squares in the course).

# task 1 ------------------------------------------------------------------


ridge_regression<-function(X, y, lambda=1, scale=0){

   #help function: objective
   rr_obj_function <- function(w, X, y, lambda, scale=0) {

    #scaling with 1/n?
    if (scale==1) {
      y<-y/sqrt(length(y))
      X<-X/sqrt(nrow(X))
      }
    return(crossprod(y - X%*%w)+ lambda * sqrt(crossprod(w)))
   }

    result=optim(rep(0, ncol(X)), rr_obj_function, X=X, y=y, lambda=lambda, scale=scale)
    return(result$par)
}

# #test:
# y<-1:10
# X<-diag(10)
# ridge_regression(X,y,lambda=1, scale=1)



########################################################################################################

# task 2 ------------------------------------------------------------------


#task 2:
#--(b) Install the package glmnet and read the description to use it to compute the lasso estimator.


#install.packages("glmnet")
library(glmnet)
#?glmnet

#function to compute lasso
#crossval_lamda asks if the lamda is choosen by crossvalidation
lasso_regression<-function(X, y, lambda=1, scale=0, crossval_lamda=0, intercept = T){

  #scaling with 1/n?
    if (scale==1) {
      y<-y/sqrt(length(y))
      X<-X/sqrt(nrow(X))
    }
  # cross validation lambda?
    if (crossval_lamda==1) {
      lambda <- cv.glmnet(X,y,intercept = intercept)$lambda.min
    }


  result=glmnet(X,y,intercept = intercept, nlambda=1,lambda=lambda)
  return(as.numeric(result$beta))
}


# #test:
# y<-1:10
# X<-diag(10)
# lasso_regression(X,y,lambda=1, scale=0, crossval_lamda=0,intercept=F)


########################################################################################################

# task 3 ------------------------------------------------------------------


#task 3:
#--(c) Simulate a model Y = X?? +??, where ?? = (0.15, ???0.33, 0.25)???, with n = 500, p = 3,
#--and with a design matrix X generated with standard normally distributed entries and scaled by X = scale (X).
#--For ??, generate i.i.d. normally distributed errors with expectation zero and variances 1/4.
#--Compare the estimates in this model across different tuning parameters.
#--Thereto, compute for each ??, averages of the estimates over M iterations.
#--Visualize the distances between least squares and ridge estimates along a sequence of tuning parameters
#--going towards zero.
#--Compare lasso and least squares estimates analogously.
#--Start with M = 10 first and use M = 1000 then.
#--Compare the average estimation errors of the three different methods for ?? = 0.1.


#Simulate model
n<-500
p<-3
beta_true <- c(0.15, -0.33, 0.25)
#beta_true

#once simulate as preparation
design<-scale(matrix(rnorm(n * p), nrow = n, ncol = p))
#scale=standardizing its columns (making the mean 0 and the standard deviation 1)
#get attributes away
attr(design, "scaled:center") <- NULL
attr(design, "scaled:scale") <- NULL
#design

epsilon <- rnorm(n, mean = 0, sd = sqrt(1/4))
#epsilon

y_observed <- design%*%beta_true + epsilon
#y_observed


#now serious:
#estimates

M<-1000
lambda_values <- c(1, 1/10, 1/100, 1/1000, 1/10000, 1/100000, 1/1000000)
estimators_beta_dif_lambda <- matrix(NA, nrow = length(lambda_values), ncol = (3*p))


for (j in 1:length(lambda_values)) {
  lambda_now <- lambda_values[j]
  mean_estimator_ls <- rep(0,p)
  mean_estimator_rr <- rep(0,p)
  mean_estimator_lr <- rep(0,p)



  for (i in 1:M) {

    #model
    design<-scale(matrix(rnorm(n * p), nrow = n, ncol = p))
    attr(design, "scaled:center") <- NULL
    attr(design, "scaled:scale") <- NULL
    epsilon <- rnorm(n, mean = 0, sd = sqrt(1/4))
    y_observed <- design%*%beta_true + epsilon

    #estimators
    #ls
    mean_estimator_ls <-  mean_estimator_ls + ridge_regression(X=design, y=y_observed, lambda=0, scale=0)
    #rr
    mean_estimator_rr <-  mean_estimator_rr + ridge_regression(X=design, y=y_observed, lambda=lambda_now, scale=0)
    #lr  without intercept: F
    mean_estimator_lr <-  mean_estimator_lr + lasso_regression(X=design, y=y_observed, lambda=lambda_now, scale=0, crossval_lamda=0, intercept = F)


  }

  estimators_beta_dif_lambda[j, 1:p] <- mean_estimator_ls/M  #mean until yet only the sum, now the mean
  estimators_beta_dif_lambda[j, (p+1):(2*p)] <- mean_estimator_rr/M
  estimators_beta_dif_lambda[j, (2*p+1):(3*p)] <- mean_estimator_lr/M


}

estimators_beta_dif_lambda

#--Visualize the distances between least squares and ridge estimates along a sequence of tuning parameters
#--going towards zero.

distances_ls_rr <- rep(0, length(lambda_values))
for (i in 1:length(lambda_values)) {
  distances_ls_rr[i] <- sqrt(crossprod(estimators_beta_dif_lambda[i, 1:p]-estimators_beta_dif_lambda[i, (p+1):(2*p)]))
}

plot(lambda_values, distances_ls_rr, log = "y", main = "distances between least squares and ridge estimates")

#################################

# task 3_2 ----------------------------------------------------------------


#--Compare lasso and least squares estimates analogously.
distances_ls_lr <- rep(0, length(lambda_values))
for (i in 1:length(lambda_values)) {
  distances_ls_lr[i] <- sqrt(crossprod(estimators_beta_dif_lambda[i, 1:p]-estimators_beta_dif_lambda[i, (2*p+1):(3*p)]))
}

plot(lambda_values, distances_ls_lr, log = "y", main = "distances between least squares and lasso estimates")


#####################################################

# task 3_3 ----------------------------------------------------------------


#--Compare the average estimation errors of the three different methods for ?? = 0.1.

lambda_now <- 0.1
M<-1000
errors <- matrix(NA, nrow = M, ncol = 3)

for (i in 1:M) {

  #model
  design<-scale(matrix(rnorm(n * p), nrow = n, ncol = p))
  attr(design, "scaled:center") <- NULL
  attr(design, "scaled:scale") <- NULL
  epsilon <- rnorm(n, mean = 0, sd = sqrt(1/4))
  y_observed <- design%*%beta_true + epsilon

  #estimators und error is squared distance of betas
  #ls
  estimator_ls <-  ridge_regression(X=design, y=y_observed, lambda=0, scale=0)
  errors[i,1] <-  crossprod(estimator_ls - beta_true)
  #rr
  estimator_rr <-  ridge_regression(X=design, y=y_observed, lambda=lambda_now, scale=0)
  errors[i,2] <-  crossprod(estimator_rr - beta_true)
  #lr  without intercept: F
  estimator_lr <-  lasso_regression(X=design, y=y_observed, lambda=lambda_now, scale=0, crossval_lamda=0, intercept = F)
  errors[i,3] <-  crossprod(estimator_lr - beta_true)

}

average_est_error <- apply(errors, 2, mean)
average_est_error

########################################################################################################

#task 4:
#--(d) We modify the model from (c) to dimension p = 20, with
#--beta = c (0.15,???0.33,0.25,???0.25, rep (0,5),???0.25,0.12,???0.125, rep (0, 8)).
#--Compute averages of the estimators based on a Monte Carlo simulation with M = 1000 iterations.
#--Consider ?? = 0.1, ?? = 0.125 and ?? = 0.15 for the ridge and the lasso estimator.
#--Present average estimates and empirical L0, L1 and L2 distances to the true ?? (L0 means dimension of non-zero entries).

# task 4 ------------------------------------------------------------------


#Simulate model
n<-500
p<-20
beta_true <- c(0.15, -0.33, 0.25, -0.25, rep (0,5), -0.25,0.12, -0.125, rep (0, 8))
#beta_true

M<-1000  #reduced number
lambda_values <- c(0.1, 0.125, 0.15)
estimators_beta_dif_lambda <- matrix(NA, nrow = length(lambda_values), ncol = (3*p))  #attention, named like before!
empirical_distances<- matrix(NA, nrow = length(lambda_values), ncol = (3*3)) #column 1 is L0, column 2 is L1, column 3 is L2, and first is ls then rr then lr

for (j in 1:length(lambda_values)) {
  lambda_now <- lambda_values[j]

  mean_estimator_ls <- rep(0,p) #better name would be sum
  mean_estimator_rr <- rep(0,p)
  mean_estimator_lr <- rep(0,p)

  ls_sum_distances_L0 <- 0
  ls_sum_distances_L1 <- 0
  ls_sum_distances_L2 <- 0

  rr_sum_distances_L0 <- 0
  rr_sum_distances_L1 <- 0
  rr_sum_distances_L2 <- 0

  lr_sum_distances_L0 <- 0
  lr_sum_distances_L1 <- 0
  lr_sum_distances_L2 <- 0




  for (i in 1:M) {

    #model
    design<-scale(matrix(rnorm(n * p), nrow = n, ncol = p))
    attr(design, "scaled:center") <- NULL
    attr(design, "scaled:scale") <- NULL
    epsilon <- rnorm(n, mean = 0, sd = sqrt(1/4))
    y_observed <- design%*%beta_true + epsilon

    #estimators
    #ls
    est<-ridge_regression(X=design, y=y_observed, lambda=0, scale=0)
    mean_estimator_ls <-  mean_estimator_ls + est
    ls_sum_distances_L0 <- ls_sum_distances_L0 + sum(est != beta_true)
    ls_sum_distances_L1 <- ls_sum_distances_L1 + sum(abs(est-beta_true))
    ls_sum_distances_L2 <- ls_sum_distances_L2 + sqrt(crossprod(est-beta_true))

    #rr
    est <- ridge_regression(X=design, y=y_observed, lambda=lambda_now, scale=0)
    mean_estimator_rr <-  mean_estimator_rr + est
    rr_sum_distances_L0 <- rr_sum_distances_L0 + sum(est != beta_true)
    rr_sum_distances_L1 <- rr_sum_distances_L1 + sum(abs(est-beta_true))
    rr_sum_distances_L2 <- rr_sum_distances_L2 + sqrt(crossprod(est-beta_true))


    #lr  without intercept: F
    est <- lasso_regression(X=design, y=y_observed, lambda=lambda_now, scale=0, crossval_lamda=0, intercept = F)
    mean_estimator_lr <-  mean_estimator_lr + est
    lr_sum_distances_L0 <- lr_sum_distances_L0 + sum(est != beta_true)
    lr_sum_distances_L1 <- lr_sum_distances_L1 + sum(abs(est-beta_true))
    lr_sum_distances_L2 <- lr_sum_distances_L2 + sqrt(crossprod(est-beta_true))
  

  }

  estimators_beta_dif_lambda[j, 1:p] <- mean_estimator_ls/M  #mean until yet only the sum, now the mean
  estimators_beta_dif_lambda[j, (p+1):(2*p)] <- mean_estimator_rr/M
  estimators_beta_dif_lambda[j, (2*p+1):(3*p)] <- mean_estimator_lr/M

  empirical_distances[j,1] <- ls_sum_distances_L0/M
  empirical_distances[j,2] <- ls_sum_distances_L1/M
  empirical_distances[j,3] <- ls_sum_distances_L2/M

  empirical_distances[j,4] <- rr_sum_distances_L0/M
  empirical_distances[j,5] <- rr_sum_distances_L1/M
  empirical_distances[j,6] <- rr_sum_distances_L2/M

  empirical_distances[j,7] <- lr_sum_distances_L0/M
  empirical_distances[j,8] <- lr_sum_distances_L1/M
  empirical_distances[j,9] <- lr_sum_distances_L2/M
  
  print(j)

}

View(estimators_beta_dif_lambda)
View(empirical_distances)
########################################################################################################

# task 5 ------------------------------------------------------------------
#--(e) Consider for increasing dimension p = 20+l ???10, l = 0, . . . , 15,
#--the linear model as above adding more zero entries to ??, i.e.
#--all entries are zero except the ones which are equal as in (d).
#--Run M = 100 iterations to compute estimates.
#--Illustrate the empirical mean absolute deviations of the three different estimators depending on the dimension.
n<-500
M<-100 #reduced number
p<-20

beta_true <- c(0.15, -0.33, 0.25, -0.25, rep (0,5), -0.25,0.12, -0.125, rep (0, 8))
number_of_dimensions <- 16

lambda_now <-0.1

#MAD is mean L1 distance between true beta and our estimates
MAD <- matrix(data = NA, nrow= number_of_dimensions , ncol = 4) #dimension, then ls, then rr, then lr

#!!!: needs so much time to run!
for (l in 0:15) {

  dim<-p+l*10
  MAD[l+1,1] <-dim

  beta_l <- c(beta_true, rep(0,(l*10)))

  ls_sum_absolute_deviation <- 0
  rr_sum_absolute_deviation <- 0
  lr_sum_absolute_deviation <- 0

  for (i in 0:M){

    #model
    entries <-(n * dim)
    design<-scale(matrix(rnorm(entries), nrow = n, ncol = dim))
    attr(design, "scaled:center") <- NULL
    attr(design, "scaled:scale") <- NULL
    epsilon <- rnorm(n, mean = 0, sd = sqrt(1/4))
    y_observed <- design%*%beta_l + epsilon

    #estimators
    #ls
    est<-ridge_regression(X=design, y=y_observed, lambda=0, scale=0)
    ls_sum_absolute_deviation <-  ls_sum_absolute_deviation + sum(abs(est-beta_l))

    #rr
    est <- ridge_regression(X=design, y=y_observed, lambda=lambda_now, scale=0)
    rr_sum_absolute_deviation <-  rr_sum_absolute_deviation + sum(abs(est-beta_l))

    #lr  without intercept: F
    est <- lasso_regression(X=design, y=y_observed, lambda=lambda_now, scale=0, crossval_lamda=0, intercept = F)
    lr_sum_absolute_deviation <-  lr_sum_absolute_deviation + sum(abs(est-beta_l))



  }

  MAD[l+1,2] <- ls_sum_absolute_deviation/M
  MAD[l+1,3] <- rr_sum_absolute_deviation/M
  MAD[l+1,4] <- lr_sum_absolute_deviation/M

  print(l)

}


View(MAD)


#Now:--Illustrate the empirical mean absolute deviations of the three different estimators depending on the dimension.
install.packages("ggplot2")
library(ggplot2)

data <- data.frame(
  x = c(MAD[,1],MAD[,1],MAD[,1]),
  y = c(MAD[,2],MAD[,3],MAD[,4]),
  method = factor(c(rep("ls", length(MAD[,1])),
                   rep("rr", length(MAD[,1])),
                   rep("lr", length(MAD[,1]))))
)
View(data)

# Create the plot with ggplot
ggplot(data, aes(x = x, y = y, color = method)) +
  geom_point(size = 4) +
  labs(title = "empirical mean absolute deviations of the three different estimators depending on the dimension",
       x = "dimension",
       y = "MAD") +
  theme_minimal()

################################

