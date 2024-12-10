#AG_Stat
#Sheet 2
#Project high-dimensional regression: Overfitting

############################
#task 1:
#--(a) Implement a function in R which computes the ridge regression estimator 
#--given the observations of the target and design values and 
#--with a general tuning parameter λ for the penalization to be chosen as an input. 
#--You shall not use packages that readily include the ridge estimator, 
#--but you can use the function optim, or similar ones, to solve optimization problems (similar as for least squares in the course).


ridge_regression_sheet_2 <- function( Y, X, lambda, beta_0) {
  
  ridge_regression_objective <- function(beta) {
    return(sum((Y-X %*% beta)^2) + lambda * sqrt(sum(beta^2) ))
    ## hier war nur sum(beta^2), aber es müsste sqrt(sum(beta^2)) sein
  }
  
  ridge_regression_estimator <- optim(par=beta_0,fn=ridge_regression_objective)
  return(ridge_regression_estimator)
  
}

#testing
Y <- c(1,2,3)
X <- diag(3)  #identity matrix
lambda <- 1
beta_0 <- c(0,0,0)

result <- ridge_regression_sheet_2( Y, X, lambda, beta_0)
ridge_regression_estimator <- result$par
test <- ridge_regression_estimator
rm(X, Y, lambda, beta_0, result, ridge_regression_estimator, test)

############################
#task 2:
#--(b) Install the package glmnet and read the description to use it to compute the lasso estimator.

#install.packages("glmnet")
library(glmnet)

############################
#task 3:
#--(c) We aim to illustrate the overfitting phenomenon for least squares and how regularization can help. 
#--Generate data from the polynomial regression model
#--    yi = β0 + β1xi + β2x^2_i + . . . + βpx^p_i + ϵi   , 1 ⩽ i ⩽ n ,
#--with a “sparse” parameter vector β = (0, 0, 2, 0, . . . , 0)⊤, 
#--i.e. the model parameters are set to zero except for the quadratic term. 

#--Simulate for the design xi ∼ U([−1, 1]), i.e. uniformly distributed on [−1, 1], and ϵi ∼ N (0, 1) i.i.d. standard normal with n = 100.
#--Observations of the response are hence noisy observations of the square function, but this is unknown to the statistician beforehand. 
#--Estimate the parameters for p = 20 with least squares, the ridge and the lasso estimator, respectively. 
#--You can use cv.glmnet to choose a suitable λ for the lasso and simply take the same λ for the ridge estimator. 
#--Compare in a plot the true function and the noisy observations with fits based on the estimated parameters with all three methods. 
#--Iterating the procedure shows that the plots can vary from time to time. 
#--Therefore, run a Monte Carlo simulation and plot Monte-Carlo averages (and quantiles). 
#--You can keep the design fix during all iterations, but re-simulate the noise in each step. Interpret the results.

################
#The design and errors
p <- 20
n<-100

#Simulation of x_i
set.seed(27); x_vector <- runif(n,min=-1, max=1) #the xi 
#summary(x_vector)

#simulation of epsilon_i
set.seed(27); epsilon_vector <- rnorm(n, mean=0, sd=1)
#summary(epsilon_vector)

print(2 * x_vector^2)
#calculation of yi
y_vector <- 2 * x_vector^2 + epsilon_vector
#summary(y_vector)


#calculation of the design
exponents <- 0:p
X <- matrix(0,nrow=n, ncol=p+1)
for (i in 1:n) {
  X[i,] <- x_vector[i] ^exponents
}

#2nd calc of yi
beta_true <- rep(0,p+1)
beta_true[3] = 2
y_2 <- X %*% beta_true + epsilon_vector
#y_2-y_vector #is 0 vector
##-> it works




################
#The estimation

Y <- y_vector
#X <- X
beta_0 <-rep(0,p+1)

?glmnet
#lasso estimator
set.seed(27); lasso_cv <- cv.glmnet(X,Y, alpha=1)
lambda_cv <- lasso_cv$lambda.min
lasso_model <- glmnet(X,Y, alpha=1 , lambda=lambda_cv)  #alpha 1 is lasso
lasso_beta <- as.vector(lasso_model$beta)
print(lasso_beta)

#least squares 
lambda <-0
result_least_squares <- ridge_regression_sheet_2( Y, X, lambda, beta_0)
least_squares_estimator <- result_least_squares$par
least_squares_estimator

#ridge
lambda <-lambda_cv 
result_ridge_regression <- ridge_regression_sheet_2( Y, X, lambda, beta_0)
ridge_regression_estimator <- result_ridge_regression$par


###################################################################
# Plotting a function and data points in R


# Create the plotting area
plot(NA, NA, xlim = c(-2, 2), ylim = c(-10, 10), 
     xlab = "x", ylab = "y", main = "Plot of the true function and the estimators with the noisy data")

# Plot the true function: y = x^2
curve(2* x^2, from = -2, to = 2, col = "blue", lwd = 2, add = TRUE)


# Plot the lasso estimator: 
exponents <- 0:p
lasso_estimator_function_h <- function(x) {
  x_val <- x^ exponents
  return(sum(lasso_beta* x_val))
}
lasso_estimator_function <- Vectorize(lasso_estimator_function_h) # now it could handle vectors as inputs
curve(lasso_estimator_function, from = -2, to = 2, col = "green", lwd = 2, add = TRUE)



# Plot the least squares estimator: 
exponents <- 0:p
least_squares_estimator_function_h <- function(x) {
  x_val <- x^ exponents
  return(sum((least_squares_estimator* x_val)))
}
least_squares_estimator_function <- Vectorize(least_squares_estimator_function_h) # now it could handle vectors as inputs
curve(least_squares_estimator_function, from = -2, to = 2, col = "red", lwd = 2, add = TRUE)




# Plot the ridge_regression estimator: 
exponents <- 0:p
ridge_regression_estimator_function_h <- function(x) {
  x_val <- x^ exponents
  return(sum((ridge_regression_estimator* x_val)))
}
ridge_regression_estimator_function <- Vectorize(ridge_regression_estimator_function_h) # now it could handle vectors as inputs
curve(ridge_regression_estimator_function, from = -2, to = 2, col = "yellow", lwd = 2, add = TRUE)




# Data points
data_points <- data.frame(x = x_vector, y = y_vector)  
# Add the data points
points(data_points$x, data_points$y, col = "purple", pch = 19)


# Optional: Add a legend
legend("topright", legend = c("y = 2 * x^2", "lasso", "least squares", "ridge", "noisy data points"),
       col = c("blue", "green", "red", "yellow", "purple"), lwd = c(2, 2, 2, 2, NA), pch = c(NA, NA, NA, NA, 19))


###########################################################

#--Iterating the procedure shows that the plots can vary from time to time. 
#--Therefore, run a Monte Carlo simulation and plot Monte-Carlo averages (and quantiles). 
#--You can keep the design fix during all iterations, but re-simulate the noise in each step. Interpret the results.

#Monte Carlo

n_simulations <- 100

coefficients_beta <- matrix(0,nrow=n_simulations, ncol=3*(p+1))

#The design and errors
p <- 20
n<-100

#Simulation of x_i
x_vector <- runif(n,min=-1, max=1) #the xi 
#summary(x_vector)

#calculation of the design
exponents <- 0:p
X <- matrix(0,nrow=n, ncol=p+1)
for (i in 1:n) {
  X[i,] <- x_vector[i] ^exponents
}
#X

#simulation of epsilon_i
epsilon_vector <- rnorm(n*n_simulations, mean=0, sd=1)
#summary(epsilon_vector)

for (i in 1: n_simulations) {
  
  
  #calculation of yi
  y_vector <- 2 * x_vector^2 + epsilon_vector[(1+n*(i-1)):(n*i)]
  #summary(y_vector)
  
  
  ################
  #The estimation
  
  Y <- y_vector
  #X <- X
  beta_0 <-rep(0,p+1)
  
  
  #lasso estimator
  lasso_cv <- cv.glmnet(X,Y, alpha=1)
  lambda_cv <- lasso_cv$lambda.min
  lasso_model <- glmnet(X,Y, alpha=1 , lambda=lambda_cv)  #alpha 1 is lasso
  lasso_beta <- as.vector(lasso_model$beta)
  coefficients_beta[i,1:(p+1)] <- lasso_beta 
  
  #least squares 
  lambda <-0
  result_least_squares <- ridge_regression_sheet_2( Y, X, lambda, beta_0)
  least_squares_estimator <- result_least_squares$par
  coefficients_beta[i,((p+1)+1):(2*(p+1))] <- least_squares_estimator
  
  #ridge
  lambda <-lambda_cv 
  result_ridge_regression <- ridge_regression_sheet_2( Y, X, lambda, beta_0)
  ridge_regression_estimator <- result_ridge_regression$par
  coefficients_beta[i,(2*(p+1)+1):(3*(p+1))] <- ridge_regression_estimator
  
  
} 

#calculate means and quantiles:
means_coefficients_beta <-colMeans(coefficients_beta)

quantiles_list <- c(0.025, 0.5, 0.975)

quantiles_coefficients_beta <- apply(coefficients_beta, 2, function(x) {
  quantile(x, probs = quantiles_list)})

#mean makes bad plots, rather median, much better!


#plot
# Create the plotting area
plot(NA, NA, xlim = c(-2, 2), ylim = c(-10, 10), 
     xlab = "x", ylab = "y", main = "Plot of the true function and the Monte Carlo mean estimators and quantiles")

# Plot the true function: y = x^2
curve(2 * x^2, from = -2, to = 2, col = "blue", lwd = 2, add = TRUE)

##################################################
# Plot the lasso estimator: 
exponents <- 0:p
mean_lasso_estimator_function_h <- function(x) {
  x_val <- x^ exponents
  return(sum((quantiles_coefficients_beta[2,1:(p+1)]* x_val)))
}
mean_lasso_estimator_function <- Vectorize(mean_lasso_estimator_function_h) # now it could handle vectors as inputs
curve(mean_lasso_estimator_function, from = -2, to = 2, col = "green", lwd = 2, add = TRUE)

## hier haben wir quantile über die y-werte gebildet und nicht über die beta parameter
q_lasso_estimator_function <- function(x){
  x_val <- x^ exponents
  werte <- coefficients_beta[,(1:(p+1))] %*% x_val
  ergebnis <- quantile(werte, probs = 0.975)
  return(ergebnis)
}
q_lasso_estimator_function <- Vectorize(q_lasso_estimator_function)
curve(q_lasso_estimator_function, from = -2, to = 2, col = "red", lwd = 2, add = TRUE)
## ende neuer teil


q1_mean_lasso_estimator_function_h <- function(x) {
  x_val <- x^ exponents
  return(sum((quantiles_coefficients_beta[1,1:p+1]* x_val)))
}
q1_mean_lasso_estimator_function <- Vectorize(q1_mean_lasso_estimator_function_h) # now it could handle vectors as inputs
curve(q1_mean_lasso_estimator_function, from = -2, to = 2, col = "green", lwd = 2, lty=2, add = TRUE)

q3_mean_lasso_estimator_function_h <- function(x) {
  x_val <- x^ exponents
  return(sum((quantiles_coefficients_beta[3,1:p+1]* x_val)))
}
q3_mean_lasso_estimator_function <- Vectorize(q3_mean_lasso_estimator_function_h) # now it could handle vectors as inputs
curve(q3_mean_lasso_estimator_function, from = -2, to = 2, col = "green", lwd = 2, lty=2, add = TRUE)

##################################################

# Plot the least squares estimator: 
exponents <- 0:p
mean_least_squares_estimator_function_h <- function(x) {
  x_val <- x^ exponents
  return(sum((quantiles_coefficients_beta[2,((p+1)+1):(2*(p+1))]* x_val)))
}
mean_least_squares_estimator_function <- Vectorize(mean_least_squares_estimator_function_h) # now it could handle vectors as inputs
curve(mean_least_squares_estimator_function, from = -2, to = 2, col = "red", lwd = 2, add = TRUE)




# Plot the ridge_regression estimator: 
exponents <- 0:p
mean_ridge_regression_estimator_function_h <- function(x) {
  x_val <- x^ exponents
  return(sum((quantiles_coefficients_beta[2,(2*(p+1)+1):(3*(p+1))]* x_val)))
}
mean_ridge_regression_estimator_function <- Vectorize(mean_ridge_regression_estimator_function_h) # now it could handle vectors as inputs
curve(mean_ridge_regression_estimator_function, from = -2, to = 2, col = "yellow", lwd = 2, add = TRUE)


# Optional: Add a legend
legend("topright", legend = c("y = x^2", "median_lasso", "q1_lasso", "q3_lasso","median_least squares", "median_ridge"),
       col = c("blue", "green", "green", "green", "red", "yellow"), lwd = c(2, 2, 2,2,2, 2), pch = c(NA, NA, NA, NA, NA, NA))
