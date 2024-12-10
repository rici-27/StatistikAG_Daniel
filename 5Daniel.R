#sheet 5
#
#Project covariance estimation: shrinkage
#


#installations
#install.packages("ggplot2") #plotting
#install.packages("MASS") #functions and datasets to support Venables and Ripley, "Modern Applied Statistics with S" 
#install.packages("gridExtra") #subplots
library(ggplot2)
library(MASS)

###############################################################################

#task 0:

#--Implement a Monte Carlo simulation with n = 100 observations of multivariate normal
#--random vectors X ??? R^5, which have expectation 0, and variance-covariance matrix
#--Sigma<???m a t r i x ( c ( 1 , r e p ( . 1 , 4 ) , . 1 , 1 , r e p ( . 1 , 3 ) , . 1 , . 1 , 1 , . 1 , . 1 ,
#--                         r e p ( . 1 , 3 ) , 1 , . 1 , r e p ( . 1 , 4 ) , 1 ) , 5 , 5 ) .

n <- 100
dim <- 5
expectation_task_0 <- rep(0,dim)
#expectation_task_0
sigma_task_0 <- matrix(c(1, rep(0.1,4),
                         0.1,1, rep(0.1,3),
                         0.1,0.1,1,0.1,0.1,
                         rep(0.1,3),1,0.1,
                         rep(0.1,4),1)
                       ,dim,dim)
#sigma_task_0
#variance 1 and covariance 0.1

# Generate multivariate normal data points, X^T in our notation
#return:an n by length(mu) matrix with one sample in each row.
X_matrix_transposed <- mvrnorm(n = n, mu = expectation_task_0, Sigma = sigma_task_0)
#str(X_matrix_transposed)
#View(X_matrix_transposed)

###############################################################################

#task 1:

#--(a) Based on the Monte Carlo simulation with M = 1,000 iterations, 
#--determine the sorted eigenvalues of the sample covariance matrix and 
#--compare them to the true eigenvalues of ?? in a plot. 
#--Compare the sum of the eigenvalues of the sample covariance matrices to the
#--sum of the eigenvalues of ?? by looking at the mean and the empirical distribution in a histogram.
#--Hint: For some suitable matrix A, the command eigen(A)$values gives its eigenvalues.

M<-1000

# Define a 3-dimensional array to store 1000 matrices of size 5x5
sample_covariance_matrices_task_1 <- array(NA, dim = c(dim, dim, M))  # 3x3 matrices, 1000 matrices
#str(sample_covariance_matrices_task_1)

sorted_eig_val_sample_covariance_matrices_task_1 <- array(NA, dim = c(M,dim))

#a vector containing the p eigenvalues of sigma, sorted in decreasing order
sorted_eig_val_sigma <- eigen(x=sigma_task_0, symmetric=TRUE, only.values=TRUE)$values
str(sorted_eig_val_sigma)
sorted_eig_val_sigma

for (i in 1:M) {
  X_transposed <- mvrnorm(n = n, mu = expectation_task_0, Sigma = sigma_task_0)
  S_n <- (t(X_transposed) %*% X_transposed)/n
  sample_covariance_matrices_task_1[,,i] <- S_n
  sorted_eig_val_sample_covariance_matrices_task_1[i,] <- eigen(x=S_n, symmetric=TRUE, only.values=TRUE)$values
  #print(i)
}

#View(sample_covariance_matrices_task_1[,,42])
#sorted_eig_val_sample_covariance_matrices_task_1[42,]

#Plot
boxplot(sorted_eig_val_sample_covariance_matrices_task_1, 
        main = "Boxplots of sample eigenvalues with true eigenvalues", ylab = "eigenvalues", 
        names = c("lambda 5", "lambda 4", "lambda 3", "lambda 2", "lambda 1"))

# Add horizontal lines for the true values
colors_task_1 <- c("green", rep("blue",4))
for (i in 1:length(sorted_eig_val_sigma)) {
  abline(h = sorted_eig_val_sigma[i], col = colors_task_1[i], lwd = 2, lty = 2) #lwd=line width , lty =line type
}

legend("topright", legend = c("true lambda 1 : 1.4", "true lambda 2 to 5: 0.9"), col = c("green", "blue"), lty = 2, lwd = 2)

#--Compare the sum of the eigenvalues of the sample covariance matrices to the
#--sum of the eigenvalues of ?? by looking at the mean and the empirical distribution in a histogram.

sum_eig_val_sigma_task_1 <- sum(sorted_eig_val_sigma)
sum_eig_val_sample_cov_matrices <- rowSums(sorted_eig_val_sample_covariance_matrices_task_1)
#str(sum_eig_val_sample_cov_matrices)

#mean(sum_eig_val_sample_cov_matrices)
#sum_eig_val_sigma_task_1

#mean sum of 5.0026 or 4.9866 almost hits true sum of 5, 
#see lecture notes that expectation of eigenvalues in sum = expected trace is equal


hist(sum_eig_val_sample_cov_matrices, main = "Histogram of eigenvalue sums", xlab = "sum", col = "lightblue", border = "black")
abline(v = mean(sum_eig_val_sample_cov_matrices), col = "purple", lwd = 2)  # 'v' specifies the x-coordinate of the line
abline(v = sum_eig_val_sigma_task_1, col = "red", lwd = 2)

# Add a legend to the plot
legend("topright", legend = c(paste("true sum of eigenvalues = ",sum_eig_val_sigma_task_1 ),
                              paste("mean sum of sample eigenvalues = ",mean(sum_eig_val_sample_cov_matrices) )), col = c("red","purple"), lwd = 2)
#paste converts its arguments (via as.character) to character strings, and concatenates them

###############################################################################

#task 2:

#--(b) Analyse for weights w = l/100, 0 ??? l ??? 100, 
#--the squared bias, variance and mean squared error, defined with the Frobenius norm as in the lecture, 
#--of shrinkage estimators w · ?? · I5 + (1 ??? w) · S ,
#--with S the sample covariance matrix and ?? · I5 the optimal shrinkage target using the true ??, with M = 1,000 iterations. 
#--Illustrate the dependence on w in a plot and determine an empirically optimal w.
#--What do you expect how the plot changes for different n? 
#--Rerun the simulation with n = 1,000 observations.

w <- c(0:100)/100
#str(w)

gamma <- sum(diag(sigma_task_0 %*% diag(dim)))/dim # diag takes diagonal of matrix as vector, normalized form is: /dim
#gamma # is 1?


mean_errors_task_2 <- array(NA,dim = c(length(w), 4))  #1st col: squared bias, 2nd: variance, 3rd: mse, 4th: w

for (i in 1:length(w)){
  weight <- w[i]
  sum_bias <- 0
  sum_var <- 0
  sum_mse <- 0
  for (j in 1:M){
    X_transposed <- mvrnorm(n = n, mu = expectation_task_0, Sigma = sigma_task_0)
    S_n <- (t(X_transposed) %*% X_transposed)/n
    shrink_est <- (1-weight)* S_n + weight*gamma*diag(dim)
    expectation_of_shrink_est <- (1-weight)* sigma_task_0 + weight*gamma*diag(dim)
    
    A <- expectation_of_shrink_est - sigma_task_0  #same for all j in M, not dependent on Sn!!
    #! might be calculated ouit of for loop
    sum_bias <- sum_bias + sum(diag(A %*% t(A)))
    
    B <- shrink_est - expectation_of_shrink_est #depends on j and Sn
    sum_var <- sum_var + sum(diag(B %*% t(B)))
    
    C <- shrink_est - sigma_task_0 #depends on j and Sn
    sum_mse <- sum_mse + sum(diag(C %*% t(C)))
  }
  
  #estimate outer expectation of mse and var with mean (Monte Carlo)
  
  mean_errors_task_2[i,1] <- sum_bias/M 
  mean_errors_task_2[i,2] <- sum_var/M
  mean_errors_task_2[i,3] <- sum_mse/M
  mean_errors_task_2[i,4] <- weight   
  print(weight)
}

#str(mean_errors_task_2)
View(mean_errors_task_2)

#--Illustrate the dependence on w in a plot and determine an empirically optimal w.
# Plot using ggplot

# Convert the matrix into a data frame
data <- data.frame(
  x = rep(mean_errors_task_2[, 4], 3),  # x-values need to be taken 3 times
  y = c(mean_errors_task_2[, 1], mean_errors_task_2[, 2], mean_errors_task_2[, 3]),  # y-values for each function
  label = rep(c("bias^2", "variance", "mse"), each = length(w))   # Labels each 101 times
)

# Plot using ggplot
ggplot(data, aes(x = x, y = y, color = label)) +
  geom_point(size = 1) +  # Add points for each error
  labs(title = "Plot of of bias^2, var and mse of shrinkage estimators depending on w", x = "w", y = "error") +
  scale_color_manual(values = c("red", "blue", "green")) +  # Set custom colors for each function
  theme_minimal() +  # Use a minimal theme
  theme(legend.title = element_blank())  # Remove legend title for clarity

#determine an empirically optimal w
w_optimal_100 <- mean_errors_task_2[which.min(mean_errors_task_2[,3]),4]
w_optimal_100

#we take the w with minimal mse empirically, 
#cause mse is the relevant measure for quality of estimator, as it is sum of bias^2 and var

#--What do you expect how the plot changes for different n? 
#--Rerun the simulation with n = 1,000 observations.

#we expect: 
#sample covariance already better estimator (law of large numbers and Tschebyscheff...)
#this means variance lower and w optimal also lower, and more importance of sample covariance and less of shrinkage

n<- 1000

w <- c(0:100)/100
#str(w)

gamma <- sum(diag(sigma_task_0 %*% diag(dim)))/dim # diag takes diagonal of matrix as vector, normalized form is: /dim
#gamma # is 1?


mean_errors_task_2 <- array(NA,dim = c(length(w), 4))  #1st col: squared bias, 2nd: variance, 3rd: mse, 4th: w

for (i in 1:length(w)){
  weight <- w[i]
  sum_bias <- 0
  sum_var <- 0
  sum_mse <- 0
  for (j in 1:M){
    X_transposed <- mvrnorm(n = n, mu = expectation_task_0, Sigma = sigma_task_0)
    S_n <- (t(X_transposed) %*% X_transposed)/n
    shrink_est <- (1-weight)* S_n + weight*gamma*diag(dim)
    expectation_of_shrink_est <- (1-weight)* sigma_task_0 + weight*gamma*diag(dim)
    
    A <- expectation_of_shrink_est - sigma_task_0  #same for all j in M, not dependent on Sn!!
    #! might be calculated ouit of for loop
    sum_bias <- sum_bias + sum(diag(A %*% t(A)))
    
    B <- shrink_est - expectation_of_shrink_est #depends on j and Sn
    sum_var <- sum_var + sum(diag(B %*% t(B)))
    
    C <- shrink_est - sigma_task_0 #depends on j and Sn
    sum_mse <- sum_mse + sum(diag(C %*% t(C)))
  }
  
  #estimate outer expectation of mse and var with mean (Monte Carlo)
  
  mean_errors_task_2[i,1] <- sum_bias/M 
  mean_errors_task_2[i,2] <- sum_var/M
  mean_errors_task_2[i,3] <- sum_mse/M
  mean_errors_task_2[i,4] <- weight   
  print(weight)
}

#str(mean_errors_task_2)
#View(mean_errors_task_2)

#--Illustrate the dependence on w in a plot and determine an empirically optimal w.
# Plot using ggplot

# Convert the matrix into a data frame
data <- data.frame(
  x = rep(mean_errors_task_2[, 4], 3),  # x-values need to be taken 3 times
  y = c(mean_errors_task_2[, 1], mean_errors_task_2[, 2], mean_errors_task_2[, 3]),  # y-values for each function
  label = rep(c("bias^2", "variance", "mse"), each = length(w))   # Labels each 101 times
)

# Plot using ggplot
ggplot(data, aes(x = x, y = y, color = label)) +
  geom_point(size = 1) +  # Add points for each error
  labs(title = "Plot of of bias^2, var and mse of shrinkage estimators depending on w", x = "w", y = "error") +
  scale_color_manual(values = c("red", "blue", "green")) +  # Set custom colors for each function
  theme_minimal() +  # Use a minimal theme
  theme(legend.title = element_blank())  # Remove legend title for clarity

#determine an empirically optimal w
w_optimal_1000 <- mean_errors_task_2[which.min(mean_errors_task_2[,3]),4]
w_optimal_1000

#now at 0.12 optimal, before at 0.62...as expected :-)



###############################################################################

#task 3:

#--Implement a feasible (non-oracle) shrinkage estimator 
#--based on estimates of the optimal weight and ?? and perform a Monte Carlo simulation with M = 1,000 and n = 100.
#--Determine the eigenvalues of this estimator and compare the sorted eigenvalues 
#--to that of the sample covariance matrix and the true eigenvalues of ?? in a plot. 
#--Compare the estimated optimal weight to your results in (b).

M <- 1000
n<- 100
dim <- 5

sorted_eig_val_sample_covariance_matrices_task_3 <- array(NA, dim = c(M,dim))
sorted_eig_val_feasible_shrinkage_matrices_task_3 <- array(NA, dim = c(M,dim))
optimal_weight_task_3 <- rep(NA, M)


for (j in 1:M){
  #sample cov
  X_transposed <- mvrnorm(n = n, mu = expectation_task_0, Sigma = sigma_task_0)
  S_n <- (t(X_transposed) %*% X_transposed)/n
  sorted_eig_val_sample_covariance_matrices_task_3[j,] <- eigen(x=S_n, symmetric=TRUE, only.values=TRUE)$values
  
  #now feasible shrinkage
  gamma_est <- sum(diag(S_n %*% diag(dim)))/dim # diag takes diagonal of matrix as vector, here identity matrix with dimension 5,
  #normalized form is: /dim
  
  help_delta <- S_n - gamma_est * diag(dim)
  delta_est <- sum(diag(help_delta %*% t(help_delta)))/dim
  
  beta_est_help <- rep(0, n)
  for (i in 1:n){
    help <- outer(t(X_transposed)[,i] , X_transposed[i,])-S_n 
    # outer is matrix product of 2 arrays of same length n to get n times n matrix
    beta_est_help[i] <-sum(diag((help %*% t(help))))/dim
    
  }
  
  beta_est <- min(sum(beta_est_help)/n^2, delta_est)
  
  w_op <- (beta_est/delta_est)
  optimal_weight_task_3[j] <- w_op
  
  #Theorem 3.4
  sigma_opt_feasible <- w_op *gamma_est *diag(dim) + (1-w_op) * S_n
  sorted_eig_val_feasible_shrinkage_matrices_task_3[j,] <- eigen(x=sigma_opt_feasible, symmetric=TRUE, only.values=TRUE)$values
  
  if (j %% 10 == 0) {
    print(j)
  }
}

#plot
# Convert the data into a data frame
data <- data.frame(
  x = rep(1:M, dim*3),  # x-values need to be taken 3 times
  y = c(sorted_eig_val_feasible_shrinkage_matrices_task_3[, 1], sorted_eig_val_feasible_shrinkage_matrices_task_3[, 2], 
        sorted_eig_val_feasible_shrinkage_matrices_task_3[, 3], sorted_eig_val_feasible_shrinkage_matrices_task_3[, 4], 
        sorted_eig_val_feasible_shrinkage_matrices_task_3[, 5],
        sorted_eig_val_sample_covariance_matrices_task_3[, 1], sorted_eig_val_sample_covariance_matrices_task_3[, 2], 
        sorted_eig_val_sample_covariance_matrices_task_3[, 3], sorted_eig_val_sample_covariance_matrices_task_3[, 4], 
        sorted_eig_val_sample_covariance_matrices_task_3[, 5],
        rep(sorted_eig_val_sigma[1],M),  rep(sorted_eig_val_sigma[2],M), 
        rep(sorted_eig_val_sigma[3],M),  rep(sorted_eig_val_sigma[4],M), 
        rep(sorted_eig_val_sigma[5],M) 
  ),  # y-values for each function
  label = rep(c("shrink_5", "shrink_4", "shrink_3","shrink_2", "shrink_1",
                "samp_5", "samp_4", "samp_3","samp_2", "samp_1",
                "true_5", "true_4", "true_3","true_2", "true_1"), each = M)   # Labels each 101 times
)

# Plot using ggplot
ggplot(data, aes(x = x, y = y, color = label)) +
  geom_line(size = 1) +  # Add points for each eig_val
  labs(title = "Plot of the eigenvalues of the feasible shrinkage and the sample covariance estimators 
       as well as the eigenvalues of the true sigma", x = "round j", y = "eigenvalue") +
  #scale_color_manual(values = c(rep("red",5), rep("blue",5), rep("green",5))) +  # Set custom colors for each function
  scale_color_manual(values = c("purple", "brown","blue", "black", "red",  "orange", "pink", "cyan", "magenta", "yellow",rep("green",5)) )+ 
  theme_minimal() +  # Use a minimal theme
  theme(legend.title = element_blank())  # Remove legend title for clarity


#better idea:
library(gridExtra)

x <- 1:M

View(sorted_eig_val_feasible_shrinkage_matrices_task_3)
# Create individual plots
#1
y <- c(sorted_eig_val_feasible_shrinkage_matrices_task_3[, 1], 
       sorted_eig_val_sample_covariance_matrices_task_3[, 1], rep(sorted_eig_val_sigma[1],M))

label = rep(c("shrink_1", "samp_1","true_1"), each = M)

p1 <- ggplot(data.frame(x, y, label), aes(x, y, color = label)) + 
  geom_point(size=1) + 
  labs(title="lambda 1",x = "round j", y = "eigenvalue") +
  scale_color_manual(values = c("pink", "cyan", "black")) +  # Set custom colors for each function
  theme_minimal() +  # Use a minimal theme
  theme(legend.title = element_blank())  # Remove legend title for clarity

#2
y <- c(sorted_eig_val_feasible_shrinkage_matrices_task_3[, 2], 
       sorted_eig_val_sample_covariance_matrices_task_3[, 2], rep(sorted_eig_val_sigma[2],M))

label = rep(c("shrink_2", "samp_2","true_2"), each = M)


p2 <- ggplot(data.frame(x, y, label), aes(x, y, color = label)) + 
  geom_point(size=1) + 
  labs(title="lambda 2",x = "round j", y = "eigenvalue") +
  scale_color_manual(values = c("pink", "cyan", "black")) +  # Set custom colors for each function
  theme_minimal() +  # Use a minimal theme
  theme(legend.title = element_blank())  # Remove legend title for clarity

#3
y <- c(sorted_eig_val_feasible_shrinkage_matrices_task_3[, 3], 
       sorted_eig_val_sample_covariance_matrices_task_3[, 3], rep(sorted_eig_val_sigma[3],M))

label = rep(c("shrink_3", "samp_3","true_3"), each = M)

p3 <- ggplot(data.frame(x, y, label), aes(x, y, color = label)) + 
  geom_point(size=1) + 
  labs(title="lambda 3",x = "round j", y = "eigenvalue") +
  scale_color_manual(values = c("pink", "cyan", "black")) +  # Set custom colors for each function
  theme_minimal() +  # Use a minimal theme
  theme(legend.title = element_blank())  # Remove legend title for clarity

#4
y <- c(sorted_eig_val_feasible_shrinkage_matrices_task_3[, 4], 
       sorted_eig_val_sample_covariance_matrices_task_3[, 4], rep(sorted_eig_val_sigma[4],M))

label = rep(c("shrink_4", "samp_4","true_4"), each = M)


p4 <- ggplot(data.frame(x, y, label), aes(x, y, color = label)) + 
  geom_point(size=1) + 
  labs(title="lambda 4",x = "round j", y = "eigenvalue") +
  scale_color_manual(values = c("pink", "cyan", "black")) +  # Set custom colors for each function
  theme_minimal() +  # Use a minimal theme
  theme(legend.title = element_blank())  # Remove legend title for clarity

#5
y <- c(sorted_eig_val_feasible_shrinkage_matrices_task_3[, 5], 
       sorted_eig_val_sample_covariance_matrices_task_3[, 5], rep(sorted_eig_val_sigma[5],M))

label = rep(c("shrink_5", "samp_5","true_5"), each = M)

p5 <- ggplot(data.frame(x, y, label), aes(x, y, color = label)) + 
  geom_point(size=1) + 
  labs(title="lambda 5",x = "round j", y = "eigenvalue") +
  scale_color_manual(values = c("pink", "cyan", "black")) +  # Set custom colors for each function
  theme_minimal() +  # Use a minimal theme
  theme(legend.title = element_blank())  # Remove legend title for clarity





# Arrange the plots in a grid (2 rows, 3 columns)
grid.arrange(p1, p2, p3, p4, p5, nrow = 3, ncol = 2, 
             top = "The decreasing eigenvalues of the feasible shrinkage estimators, the sample covariance estimators and the true sigma")

#note: the min and max lambda is better with shrinkage, but the lambdas in the middle also? rather not...e.g. lambda 3!
#shrink in direction of mean eigenvalue makes shrinkage, which is around 1, for 2 to 5 to large!

#--Compare the estimated optimal weight to your results in (b).

x <- 1:M

mean_shrink_weight_task_3 <- mean(optimal_weight_task_3)
y <- c(optimal_weight_task_3, 
       rep(w_optimal_100,M), rep(w_optimal_1000,M), rep(mean_shrink_weight_task_3,M))

label = rep(c("shrink_est", "opt_n_100","opt_n_1000", "mean_shrink"), each = M)
# sind die 1er weights auffällig?
ggplot(data.frame(x, y, label), aes(x, y, color = label)) + 
  geom_point(size=1) + 
  labs(title="optimal weight",x = "round j", y = "eigenvalue") +
  scale_color_manual(values = c("pink", "cyan", "magenta" ,"black")) +  # Set custom colors for each function
  theme_minimal() +  # Use a minimal theme
  theme(legend.title = element_blank())

#now optimal weight in mean if w = 0.67, more than opt_100 = 0.57, and a lot more then the optimal for 1000=n
#too much deviation around optimal value? interpretation?

# vielleicht mal M hochdrehen?
# zu viele gewichte gleich 1?

###############################################################################

#task 4:

#--(d) Consider ?? with entries 1 on the diagonal and all entries off the diagonal equal to 0.1.
#--Compare for dimension p = l/10, 1 ??? l ??? 10, the squared norms of estimation errors of ?????1, 
#--estimated with the inverted sample covariance matrix and the inverted shrinkage estimator from (c). 
#--Illustrate the ratios of these squared norms of estimation errors in a plot.

#? p not divide by 10???
p_vector <- 1:10

n<-100

#sigma_max
sigma_max <- array(NA, dim = c(p_vector[length(p_vector)],p_vector[length(p_vector)]))
for (i in 1:p_vector[length(p_vector)]){
  for (j in 1:p_vector[length(p_vector)]){
    if(j == i) {
      sigma_max[i,j]<-1 }
    else {
      sigma_max[i,j]<-0.1
    }
    
  }
}

#View(sigma_max)

expectation_task_4 <- rep(0, p_vector[length(p_vector)])
#expectation_task_4

#the observations
X_transposed_max <- mvrnorm(n = n, mu = expectation_task_4, Sigma = sigma_max)



#squared Frobenius norm of error
estimation_error_task_4 <- array(NA, dim = c(length(p_vector),2))  # first column is sample cov, second shrinkage estimation error

counter<-0

for (p in p_vector){
  
  counter <- counter + 1
  
  sigma_p <- sigma_max[1:p,1:p]
  sigma_p_inv <- solve(sigma_p) #calculate inverse
  
  X_transposed <- matrix(X_transposed_max[,1:p], ncol=p)  # for p=1 we still need a matriux, not a vector
  
  #sample cov
  S_n <- (t(X_transposed) %*% X_transposed)/n
  S_n_inv <- solve(S_n) #calculate inverse
  #distance:
  help <- S_n_inv - sigma_p_inv
  estimation_error_task_4[counter,1] <-sum(diag(help %*% t(help)))/p 
  
  
  #now feasible shrinkage
  gamma_est <- sum(diag(S_n %*% diag(p)))/p # diag takes diagonal of matrix as vector, here identity matrix with dimension p
  
  help_delta <- S_n - gamma_est * diag(p)
  delta_est <- sum(diag(help_delta %*% t(help_delta)))/p
  
  beta_est_help <- rep(0, n)
  for (i in 1:n){
    help <- outer(t(X_transposed)[,i] , X_transposed[i,])-S_n 
    # outer is matrix product of 2 arrays of same length n to get n times n matrix
    beta_est_help[i] <-sum(diag((help %*% t(help))))/p
    
  }
  
  beta_est <- min(sum(beta_est_help)/n^2, delta_est)
  
  if (beta_est==delta_est){    #solution to problem downwards
    w_op=1
  }
  else{w_op <- (beta_est/delta_est)}
  
  #Theorem 3.4
  sigma_opt_feasible <- w_op *gamma_est *diag(p) + (1-w_op) * S_n
  
  help_shrink <- solve(sigma_opt_feasible) - S_n_inv
  estimation_error_task_4[counter,2] <-sum(diag(help_shrink %*% t(help_shrink)))/p 
  
  
  
  print(p)
}

#estimation_error_task_4

#problem with shrink and p=1 : NA value as S_n = gamma_est * diag(p) and delta = 0, so beta/delta =1 but not defined, 
#as solution do we make if condition

x <- p_vector

y <- estimation_error_task_4[,2]/estimation_error_task_4[,1]
label = rep("ratio error shrinkage/error sampel cov", M)

ggplot(data.frame(x, y, label), aes(x, y, color = label)) + 
  geom_point(size=1) + 
  labs(title="error comparison for estimating sigma inverse for different dimensions of p",x = "dimension p", y = "ratio") +
  scale_color_manual(values = c("black")) +  # Set custom colors for each function
  theme_minimal() +  # Use a minimal theme
  theme(legend.title = element_blank())


#interpretation: larger p, shrinkage is better, but here misses Monte Carlo to get away random effects of choosen data sample X


###############################################################################

#task 5:

#--(e) Download the data Stock Bond 2004 to 2006.csv from the WueCampus site 
#--and extract the returns, i.e. the differences of log-prices for each of the 8 stocks. 
#--We analyse the plugin portfolio built upon the optimal weights and estimating ?????1 with the sample covariance matrix 
#--and compare to a portfolio for that we estimate ?? ???1 with the feasible shrinkage estimator. 
#--Moreover, we compare to a benchmark portfolio with constantly equal weight 1/8 on each asset. 
#--Use the data over 336 days backwards as a training sample to estimate ?????1 , assuming it is constant and i.i.d. observations. 
#--We can then compute and compare the three portfolios for 336 days. 
#--Compare the overall returns of the three portfolios.

data_task_5 <- read.csv("C:\\Users\\Daniel\\Desktop\\R_Statistik AG\\Sheet_5\\Stock_Bond_2004_to_2006.csv")
#View(head(data_task_5))
#str(data_task_5)
#colnames(data_task_5)

#? why 8 stocks?

#  These 11 columns likely represent prices or similar metrics for individual stocks:
#GM_AC: General Motors stock price (AC could indicate adjusted closing price or similar).
#F_AC: Ford stock price.
#UTX_AC: United Technologies stock price.
#CAT_AC: Caterpillar stock price.
#MRK_AC: Merck stock price.
#PFE_AC: Pfizer stock price.
#IBM_AC: IBM stock price.
#MSFT_AC: Microsoft stock price.
#C_AC: Citigroup stock price.
#XOM_AC: ExxonMobil stock price.
#SP_AC: S&P 500 index or a related metric.

data_needed_task_5 <-data_task_5[, c("DATE" , "GM_AC" , "F_AC", "UTX_AC"  , "CAT_AC", 
                                     "MRK_AC", "PFE_AC", "IBM_AC" ,  "MSFT_AC" , "C_AC" , "XOM_AC", "SP_AC" )]
data_needed_task_5$DATE <- as.Date(data_needed_task_5$DATE, format = "%m/%d/%Y")
#View(head(data_needed_task_5))
#dim(data_needed_task_5) # 673  12

#define returns
# Neue Spalte berechnen (Differenz der natürlichen Logarithmen)
data_needed_task_5$return_1<- c(NA, diff(log(data_needed_task_5$GM_AC))) #diff subtracts the values in an array with index distance 1
#View(head(data_needed_task_5))

#?works out: return at time t+1 is ln(P+1t-ln(Pt)???, else NA at the end...
data_needed_task_5$return_2<- c(NA, diff(log(data_needed_task_5$F_AC)))
data_needed_task_5$return_3<- c(NA, diff(log(data_needed_task_5$UTX_AC)))
data_needed_task_5$return_4<- c(NA, diff(log(data_needed_task_5$CAT_AC)))
data_needed_task_5$return_5<- c(NA, diff(log(data_needed_task_5$MRK_AC)))
data_needed_task_5$return_6<- c(NA, diff(log(data_needed_task_5$PFE_AC)))
data_needed_task_5$return_7<- c(NA, diff(log(data_needed_task_5$IBM_AC)))
data_needed_task_5$return_8<- c(NA, diff(log(data_needed_task_5$MSFT_AC)))
data_needed_task_5$return_9<- c(NA, diff(log(data_needed_task_5$C_AC)))
data_needed_task_5$return_10<- c(NA, diff(log(data_needed_task_5$XOM_AC)))
data_needed_task_5$return_11<- c(NA, diff(log(data_needed_task_5$SP_AC)))

almost_final_data_task_5 <- data_needed_task_5[,c("DATE","return_1","return_2","return_3","return_4","return_5",
                                                  "return_6","return_7","return_8","return_9", "return_10","return_11")]
final_data_task_5 <- almost_final_data_task_5[-1,]#delete NA row
#View(final_data_task_5)

final_matrix_task_5 <- as.matrix(final_data_task_5[, -1]) #delte DATE column
#View(final_matrix_task_5)
#dim(final_matrix_task_5) # 672  11

days<-336
p <- 11
n <- days

sigma_inv_task_5 <- array(NA, dim = c(2, days, p, p))
w_opt_plugin_task_5 <- array(NA, dim = c(3, days, p))

ones <- rep(1, p)

for (i in 1:days){
  
  X_transposed <- final_matrix_task_5[i:(days+i-1),]  
  
  #sample cov
  S_n <- (t(X_transposed) %*% X_transposed)/n
  S_n_inv <- solve(S_n) #calculate inverse
  sigma_inv_task_5[1,i, , ] <- S_n_inv
  w_opt_plugin_task_5[1,i, ] <- (S_n_inv %*% ones) / ((ones %*% (S_n_inv %*% ones))[1,1])   
  #interpret the denominator not as matrix, but scalar!
  
  #now feasible shrinkage
  gamma_est <- sum(diag(S_n %*% diag(p)))/p # diag takes diagonal of matrix as vector, here identity matrix with dimension p
  
  help_delta <- S_n - gamma_est * diag(p)
  delta_est <- sum(diag(help_delta %*% t(help_delta)))/p
  
  beta_est_help <- rep(0, n)
  for (j in 1:n){
    help <- outer(t(X_transposed)[,j] , X_transposed[j,])-S_n 
    # outer is matrix product of 2 arrays of same length n to get n times n matrix
    beta_est_help[j] <-sum(diag((help %*% t(help))))/p
    
  }
  
  beta_est <- min(sum(beta_est_help)/n^2, delta_est)
  
  if (beta_est==delta_est){    #solution to problem downwards
    w_op=1
  }
  else{w_op <- (beta_est/delta_est)}
  
  #Theorem 3.4
  sigma_opt_feasible <- w_op *gamma_est *diag(p) + (1-w_op) * S_n
  est_sigma_inv <- solve(sigma_opt_feasible)
  sigma_inv_task_5[2,i, , ] <- est_sigma_inv
  w_opt_plugin_task_5[2,i, ] <- (est_sigma_inv %*% ones) / ((ones %*% (est_sigma_inv %*% ones))[1,1])
  
  #a benchmark portfolio with constantly equal weight
  w_opt_plugin_task_5[3,i, ] <- ones/sum(ones)
  
  #print(i)
  
}

#w_opt_plugin_task_5[,1,] #portfolios for the first day : days+1 (row)

#--Compare the overall returns of the three portfolios.

overall_return_sample_cov <- 0
overall_return_shrinkage <- 0
overall_return_benchmark <- 0

for (i in 1:days){
  
  overall_return_sample_cov <- overall_return_sample_cov + sum(w_opt_plugin_task_5[1,i, ] * final_matrix_task_5[(days+i),] )
  overall_return_shrinkage <- overall_return_shrinkage + sum(w_opt_plugin_task_5[2,i, ] * final_matrix_task_5[(days+i),] )
  overall_return_benchmark <- overall_return_benchmark + sum(w_opt_plugin_task_5[3,i, ] * final_matrix_task_5[(days+i),] )
  
  
}

print(paste("Overall return for sample cov: " ,  as.character(overall_return_sample_cov) , " ."))
print(paste("Overall return for shrinkage: " ,  as.character(overall_return_shrinkage) , " ."))
print(paste("Overall return for benchmark: " ,  as.character(overall_return_benchmark) , " ."))

#? is overall return just sum of return w_i * return_i?
#or more to do?
#? normal that benchmark best, then shrink, then sample cov?
#"Overall return for sample cov:  0.107906556212479  ."
#"Overall return for shrinkage:  0.114821809262127  ."
#"Overall return for benchmark:  0.120819651785536  ."


#################################################################################
#not working
data_plot <- data.frame(
  x = mean_errors_task_2[,4], 
  y = c(mean_errors_task_2[,1], mean_errors_task_2[,2], mean_errors_task_2[,3]),
  function_label = rep(c("bias ^2 ", "variance", "mse"), each = length(x))
)

data_plot <- data.frame(mean_errors_task_2)
data_plot


ggplot(data_plot, aes(x = X4, y = c(X1,X2,X3))) +
  geom_line() +  # Add the lines for each function
  labs(title = "Plot of bias, var and mse of shrinkage estimators depending on w", x = "w", y = "error") +
  scale_color_manual(values = c("red", "blue", "green")) +  # Set custom colors for each function
  theme_minimal()   # Use a minimal theme