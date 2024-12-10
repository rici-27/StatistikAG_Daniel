#Sheet 4
#Project regularized regression: Credit data

#########################################
#task 1:
#--(a) Install the packages glmnet and ISLR2 and consider the credit data set str ( Credit ). 
#--The response is balance (average credit card debt for each individual). 
#--Asses(=einschaetzen) first heuristically the relation between each explanatory variable and the response separately.

install.packages("ggplot2")
install.packages("glmnet")
install.packages("ISLR2")
install.packages("gridExtra")
install.packages("tidyr")


library(ggplot2)
library(glmnet)
library(ISLR2)
library(gridExtra) # make subplots
library(tidyr) # Priskas nices tool


#?Credit  # already installed data set
#?glmnet  # do lasso and so on
#?ISLR2  #from statistics book
#?tidyr


# Load the Credit dataset and have a look
original_data <-Credit
str(original_data)

View(original_data)
head(original_data)

#plots of relations response vs explanatory variables

# Gather data into a long format, with balance as the dependent variable 
#long means each datapoint splitted in 10
df_long <- original_data %>%
  gather(key = "variable", value = "value", -Balance)

# Create the plots
ggplot(df_long, aes(x = value, y = Balance)) +
  geom_point() +
  facet_wrap(~ variable, scales = "free", ncol = 3) +   #create a separate plot (facet) for each value of the variable "variable", each column of old data set
  labs(title = "All plots of balance against other variables",
       x = "Value of other variable", y = "average credit card debt") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) 


### separation of factor and numeric
df_long_num <- original_data %>%
  gather(key = "variable", value = "value", -Balance, -Married, -Own, -Student, -Region)

df_long_fac <- original_data %>%
  gather(key = "variable", value = "value", Married, Own, Student, Region)

ggplot(df_long_num, aes(x = value, y = Balance)) +
  geom_point() +
  facet_wrap(~ variable, scales = "free", ncol = 3) +   #create a separate plot (facet) for each value of the variable "variable", each column of old data set
  labs(title = "All plots of balance against numeric variables",
       x = "Value of other variable", y = "average credit card debt") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1)) 

#result: only the variables income, limit and rating seem to have information for balance

ggplot(df_long_fac, aes(x = value, y = Balance)) +
  geom_point() +
  facet_wrap(~ variable, scales = "free", ncol = 2) +   #create a separate plot (facet) for each value of the variable "variable", each column of old data set
  labs(title = "All plots of balance against factor variables",
       x = "Value of other variable", y = "average credit card debt") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1)) # Change x axis tick labels only

#result: all factor variables have no information for balance


########################################################################################
#task 2:
#--(b) Consider a regression with p = 3 and the numeric explanatory variables income, limit and rating. 
#--Standardize these three variables with their empirical standard deviations and then
#--combine them in a design matrix. 
#--Compute ridge regression estimates for this regression based on glmnet (..., alpha = 0, intercept =T,lambda =...) 
#--across different penalty parameters g r i d <??? r e v ( 1 0 ^ seq ( 4 , ???2 , l e n g t h = 1 0 0 0 0 ) ) .
#--Illustrate the estimates jointly in one plot. 
#--Compute the lasso estimates for the same regression with the same grid of penalty parameters. 
#--Create an analogous plot.
#--Compare the ratios of the L2-norms of the ridge regression and the least squares estimates
#--across the grid of penalty parameters. 
#--Perform the comparison first with and then without the intercept parameter. 
#--Illustrate the results in plots and discuss the results.
#--Hint: The command coef( ) can be useful.

p<-3
y_task_2 <- as.vector(original_data$Balance)
str(y_task_2)
data_task_2 <- original_data[, c("Income", "Limit", "Rating")]
str(data_task_2)

sd(data_task_2$Income)
mean(data_task_2$Income)

data_task_2$Income <- (data_task_2$Income) / sd(data_task_2$Income)
# mean abziehen?
sd(data_task_2$Income)
mean(data_task_2$Income)

data_task_2$Limit <- (data_task_2$Limit) / sd(data_task_2$Limit)
data_task_2$Rating <- (data_task_2$Rating ) / sd(data_task_2$Rating)

str(data_task_2)

design_task_2<-as.matrix(data_task_2)
str(design_task_2)


##
#--Compute ridge regression estimates for this regression based on glmnet (..., alpha = 0, intercept =T,lambda =...) 
#--across different penalty parameters g r i d <??? r e v ( 1 0 ^ seq ( 4 , ???2 , l e n g t h = 1 0 0 0 0 ) ) .

vector_grid <- 10^seq(4,-2,length=10000) #10^4 to 10^-2 in 10000 steps
str(vector_grid)
grid <- rev(vector_grid) #change the ordering, first entry now last
str(grid)

rr_task_2 <- glmnet(x=design_task_2, y=y_task_2, alpha=0, intercept=T, lambda=grid)

intercepts_task_2<-rr_task_2$a0
str(intercepts_task_2)
View(intercepts_task_2[10000])
betas_task_2<-rr_task_2$beta
str(betas_task_2)

coefficients_task_2 <- coef(rr_task_2)
str(coefficients_task_2)
head(coefficients_task_2)
#str(coefficients_task_2[2,]) #here we get the estimators for the income parameter

coefficients_task_2_matrix<-t(as.matrix(coefficients_task_2[1:4,])) 
#t means transpose, each of 4 columns now a estimator of one parameter
str(coefficients_task_2_matrix)

#--Illustrate the estimates jointly (=ZUSAMMEN) in one plot. 
boxplot(coefficients_task_2_matrix, 
        main = "Boxplots of Distributions of the parameter estimates", 
        ylab = "Value", 
        col = c("lightblue", "lightgreen", "lightcoral", "orange"), 
        names = c("Intercept", "Income", "Limit", "Rating"))

#?????? Question: ????? why income negativ beta?

#--Compute the lasso estimates for the same regression with the same grid of penalty parameters. 
#--Create an analogous plot.

lr_task_2 <- glmnet(x=design_task_2, y=y_task_2, alpha=1, intercept=T, lambda=grid)

coefficients_task_2_lasso <- coef(lr_task_2)

lasso_coefficients_task_2_matrix<-t(as.matrix(coefficients_task_2_lasso[1:4,])) 
#t means transpose, each of 4 columns now a estimator of one parameter

#--Illustrate the estimates jointly (=ZUSAMMEN) in one plot. 
boxplot(lasso_coefficients_task_2_matrix, 
        main = "Boxplots of Distributions of the lasso parameter estimates", 
        ylab = "Value", 
        col = c("lightblue", "lightgreen", "lightcoral", "orange"), 
        names = c("Intercept", "Income", "Limit", "Rating"))
##
#--Compare the ratios of the L2-norms of the ridge regression and the least squares estimates
#--across the grid of penalty parameters. 
#--Perform the comparison first with and then without the intercept parameter. 

#least squares parameters
x_1 <-design_task_2[,1]
x_2 <-design_task_2[,2]
x_3 <-design_task_2[,3]

least_squares_task_2 <- lm(y_task_2 ~ x_1 + x_2 +x_3)
least_squares_task_2
ls_est_task_2 <- least_squares_task_2$coefficients
#ls_est_task_2

norm_ls_task_2 <- sqrt(sum(ls_est_task_2^2))
#norm_ls_task_2

rr_norms_est_task_2 <- sqrt(rowSums(coefficients_task_2_matrix ^ 2))
#head(rr_norms_est_task_2)

ratios_L2_task_2 <- rr_norms_est_task_2/norm_ls_task_2
head(ratios_L2_task_2)

##
#Illustrate the results in plots and discuss the results.

# Create the scatterplot
data_rr_norm <- data.frame(ratios_L2_task_2 ,rr_task_2$lambda) #rr_task_2$lambda is grid in descending order, so = rev(grid)
ggplot(data_rr_norm, aes(x = rr_task_2$lambda , y = ratios_L2_task_2 )) +
  geom_point(size = .1) +
  labs(title = "ratios of the L2-norms of the ridge regression and the least squares estimates", x = "lambdas", y = "ratios")


#
#now without intercept:  
#

ls_est_task_2_without_intercept <- ls_est_task_2[2:length(ls_est_task_2)]
#ls_est_task_2_without_intercept

norm_ls_task_2_without_intercept <- sqrt(sum(ls_est_task_2_without_intercept^2))
#norm_ls_task_2_without_intercept

coefficients_task_2_matrix_without_intercept <- coefficients_task_2_matrix[,2:ncol(coefficients_task_2_matrix)]
#head(coefficients_task_2_matrix_without_intercept)

rr_norms_est_task_2_without_intercept <- sqrt(rowSums(coefficients_task_2_matrix_without_intercept ^ 2))
#head(rr_norms_est_task_2)

ratios_L2_task_2_without_intercept <- rr_norms_est_task_2_without_intercept/norm_ls_task_2_without_intercept
#head(ratios_L2_task_2_without_intercept)

##
#Illustrate the results in plots and discuss the results.

# Create the scatterplot
data_rr_norm <- data.frame(ratios_L2_task_2_without_intercept ,rr_task_2$lambda) #rr_task_2$lambda is grid in descending order, so = rev(grid)
ggplot(data_rr_norm, aes(x = rr_task_2$lambda , y = ratios_L2_task_2_without_intercept )) +
  geom_point(size = .1) +
  labs(title = "ratios of the L2-norms of the rr and the ls estimates without intercept", x = "lambdas", y = "ratios")

#Discussion: the greater lambda, the smaller the L2 norm of the rr estimator, because this 
#is minimized in the penalty. But, the intercept is always the same (see boxplots).
#this i why in the plot with intercept, we only converge falling to 0.7 not 0 L2 norm.

########################################################################################
#task 3:
#--(c) One could reasonably have measured income in dollars instead of thousands of dollars,
#--which would result in a rescaling of the observed values of income by a factor of 1,000.
#--Consider explanatory variables income, limit and rating, here without standardization.
#--How should rescaling influence least squares 
#--and the regularized regression results theoretically? 
#--Compute the estimates before and after rescaling, the lasso and ridge regression
#--once more with glmnet and with ?? = 100, and interpret the outcomes.

#Theoretically influence of rescaling: 
#ls: 
#no influence for estimation, just the scaled variable induces a beta_i 
#scaled with 1 over the scaling m?,
#cause then just multiply beta form left with  
#diagonal matrix that scales beta_i with m?, or X from right:  X*M*beta.
#X*M scales the ith column of X with m?...

#rr or lasso: 
#you like to minimize the norm of b, means that beta_i is already small as divide by 1000.
#get the other values small and do not approximate this value beta_i so well?


#Compute it:
lambda_3 <-100


y_task_3 <- as.vector(original_data$Balance)
#str(y_task_3)
data_task_3 <- original_data[, c("Income", "Limit", "Rating")]
#str(data_task_3)

design_task_3<-as.matrix(data_task_3)
#str(design_task_3)


#rescaled:
m <-1000
design_task_3_rescaled <-design_task_3
# Scale the first column by 1000
design_task_3_rescaled[, 1] <- design_task_3_rescaled[, 1] * m
#head(design_task_3_rescaled)

#estimate

#ls:
#least squares parameters
x_1 <-design_task_3[,1]
x_2 <-design_task_3[,2]
x_3 <-design_task_3[,3]

least_squares_task_3 <- lm(y_task_3 ~ x_1 + x_2 +x_3)
ls_est_task_3 <- least_squares_task_3$coefficients
ls_est_task_3

#scaled
x_1 <-design_task_3_rescaled[,1]
x_2 <-design_task_3_rescaled[,2]
x_3 <-design_task_3_rescaled[,3]

least_squares_task_3_rescaled <- lm(y_task_3 ~ x_1 + x_2 +x_3)
ls_est_task_3_rescaled <- least_squares_task_3_rescaled$coefficients
ls_est_task_3_rescaled

#we see everything same, just beta_income divided by 1000, just as we thought in theory, nice!
#different intercept: not standardized X!


#
#rr:
#

rr_task_3 <- glmnet(x=design_task_3, y=y_task_3, alpha=0, intercept=T, lambda=lambda_3)

coefficients_task_3 <- coef(rr_task_3)
str(coefficients_task_3)
head(coefficients_task_3)

coefficients_task_3_matrix<-t(as.matrix(coefficients_task_3[1:4,])) 
#t means transpose, each of 4 columns now a estimator of one parameter
str(coefficients_task_3_matrix)

#scaled:
rr_task_3_rescaled <- glmnet(x=design_task_3_rescaled, y=y_task_3, alpha=0, intercept=T, lambda=lambda_3)

coefficients_task_3_rescaled <- coef(rr_task_3_rescaled)
#str(coefficients_task_3_rescaled)
#head(coefficients_task_3_rescaled)

coefficients_task_3_matrix_rescaled<-t(as.matrix(coefficients_task_3_rescaled[1:4,])) 
#t means transpose, each of 4 columns now a estimator of one parameter
str(coefficients_task_3_matrix_rescaled)


#
#lasso lr:
#

lr_task_3 <- glmnet(x=design_task_3, y=y_task_3, alpha=1, intercept=T, lambda=lambda_3)

lr_coefficients_task_3 <- coef(lr_task_3)
#str(lr_coefficients_task_3)
#head(lr_coefficients_task_3)

lr_coefficients_task_3_matrix<-t(as.matrix(lr_coefficients_task_3[1:4,])) 
#t means transpose, each of 4 columns now a estimator of one parameter
#str(lr_coefficients_task_3_matrix)

#scaled:
lr_task_3_rescaled <- glmnet(x=design_task_3_rescaled, y=y_task_3, alpha=1, intercept=T, lambda=lambda_3)

lr_coefficients_task_3_rescaled <- coef(lr_task_3_rescaled)
#str(lr_coefficients_task_3_rescaled)
#head(lr_coefficients_task_3_rescaled)

lr_coefficients_task_3_matrix_rescaled<-t(as.matrix(lr_coefficients_task_3_rescaled[1:4,])) 
#t means transpose, each of 4 columns now a estimator of one parameter
#str(lr_coefficients_task_3_matrix_rescaled)


#show together
# Create a matrix with 3 rows and 4 columns
task_3_summary_matrix <- matrix(, nrow = 4, ncol = 6)
task_3_summary_matrix[,1]<-t(ls_est_task_3)
task_3_summary_matrix[,2]<-t(ls_est_task_3_rescaled)
task_3_summary_matrix[,3]<-t(coefficients_task_3_matrix)
task_3_summary_matrix[,4]<-t(coefficients_task_3_matrix_rescaled)
task_3_summary_matrix[,5]<-t(lr_coefficients_task_3_matrix)
task_3_summary_matrix[,6]<-t(lr_coefficients_task_3_matrix_rescaled)

# Set column names
colnames(task_3_summary_matrix) <- c("ls", "ls_scaled", "rr", "rr_scaled","lr", "lr_scaled")

# Print the matrix with column names
print(task_3_summary_matrix)

#for ridge: just as at ls, beta_income divided by 1000, 
#for lasso: income no influence, set to 0, so rescaling does not change anything, 
#beta is 0 is reasonable as income on first plot had no clear influence...

#Has Rici same results??????

#Chat tells us:
#In regularized regression (Ridge, Lasso), rescaling is crucial 
#because the regularization penalty depends on the scale of the variables. 
#Without rescaling, the model may treat variables of different scales unequally, 
#which can distort the shrinkage effect and feature selection. 
#Standardization ensures that the regularization term applies uniformly 
#to all features, improving fairness and interpretability.

########################################################################################
#task 4:
#--(d) Perform the regressions with the complete designs given in the data. 
#--Use the command predict ( ) to compare predictions based on the three 
#--different regression estimates, least squares, ridge regression and lasso. 
#--Make predictions for the explanatory variables from the 400th individual 
#--in the sample as an input. 
#--Use data driven choices of ?? here with cv . glmnet ( ) $lambda . min .
#--Hint: The command model.matrix( ) can be useful.

#The model.matrix() function in R is used 
#to create a design matrix from a formula and a data frame.
#e.g. X <- model.matrix(response ~ numeric_var1 + numeric_var2 + factor_var, data = df)

library(dplyr)
#str(original_data)

#true value
y_400_true <- original_data[400,11]
y_400_true


#least squares with lm:
cut_data_task_4 <- original_data[-nrow(original_data),]
str(cut_data_task_4)
least_squares_task_4 <- lm(Balance~Income+Limit+Rating+Cards+Age
                           +Education+Own+Student+Married+Region, 
                           data = cut_data_task_4)
ls_est_task_4 <- least_squares_task_4$coefficients
ls_est_task_4

ls_prediction <- 
  predict(object=least_squares_task_4, newdata=(original_data[400,])[-length(original_data[400,])])
#cut away balance value
ls_prediction

#distance to true value
ls_task_4_dist <- abs(y_400_true-ls_prediction)
ls_task_4_dist

#
#now rr and lasso
#

#make design matrix out of the data, especially dummy variables for categorical stuff
design_matrix_x_data_task_4 <- model.matrix(Balance~Income+Limit+Rating+Cards+Age
                                     +Education+Own+Student+Married+Region, 
                                     data = original_data)
#str(design_matrix_x_data_task_4)
print(design_matrix_x_data_task_4[400,])


y_data_task_4 <- original_data %>% 
  select(Balance)
#str(y_data_task_4)


#in order to predict the last value, we exclude it from the data we make the regressions on
# Remove the last row using indexing
design_matrix_x_data_task_4_final <- 
  design_matrix_x_data_task_4[-nrow(design_matrix_x_data_task_4), ]
y_data_task_4_final <- y_data_task_4[-nrow(y_data_task_4),]
#str(design_matrix_x_data_task_4_final)
#str(y_data_task_4_final)

#combined_data_task_4 <- cbind(design_matrix_x_data_task_4_final, y_data_task_4_final )
#str(combined_data_task_4)
#head(combined_data_task_4)


#
#train models
#

#rr:
#???Intercept True or False? Is in design already, do we need it with true?
lambda_rr_task_4 <- cv.glmnet(x=design_matrix_x_data_task_4_final,y=y_data_task_4_final, 
                              alpha=0, intercept=FALSE)$lambda.min
#lambda_rr_task_4 

rr_task_4 <- glmnet(x=design_matrix_x_data_task_4_final,y=y_data_task_4_final, 
                    alpha=0, intercept=FALSE, lambda=lambda_rr_task_4 )
rr_est_task_4 <- rr_task_4$beta
rr_est_task_4

#lr:
#???Intercept True or False? Is in design already, do we need it with true?
lambda_lr_task_4 <- cv.glmnet(x=design_matrix_x_data_task_4_final,y=y_data_task_4_final, 
                              alpha=1, intercept=FALSE)$lambda.min
#lambda_lr_task_4 

lr_task_4 <- glmnet(x=design_matrix_x_data_task_4_final,y=y_data_task_4_final, 
                    alpha=1, intercept=FALSE, lambda=lambda_lr_task_4 )
lr_est_task_4 <- lr_task_4$beta
lr_est_task_4

?cv.glmnet
#
#predictions
#

rr_prediction <- 
  predict(rr_task_4, newx=design_matrix_x_data_task_4[400,])
rr_prediction
#distance to true value
rr_task_4_dist <- abs(y_400_true-rr_prediction)
rr_task_4_dist

lr_prediction <- 
  predict(lr_task_4, newx=design_matrix_x_data_task_4[400,])
lr_prediction
#distance to true value
lr_task_4_dist <- abs(y_400_true-lr_prediction)
lr_task_4_dist

#rr pretty bad, but because of large lambda, around 225, vs at lasso: 1
#nevertheless, zero for rating strange at lasso parameters...

########################################################################################
#task 5:
#--(e) Select ?? for the lasso within the grid (0:100)/10 by splitting the data 
#--randomly in one half of training and another half of test sample, 
#--the random decompositions iterated M = 1,000 times, and 
#--minimizing the overall empirical L2-prediction error of the forecasts 
#--for the test samples based on the regressions with the training samples. 
#--Illustrate the dependence of the empirical L2-prediction error on ?? in a plot. 
#--Compare to analogous empirical L2-prediction errors based on least squares estimates 
#--and of simply fitting a model with just an intercept.

grid_lambda <- (0:100)/10
#grid_lambda


data_task_5 <- original_data
#str(data_task_5)

#make design matrix out of the data, especially dummy variables for categorical stuff
design_matrix_x_data_task_5 <- model.matrix(Balance~Income+Limit+Rating+Cards+Age
                                            +Education+Own+Student+Married+Region, 
                                            data = data_task_5)
#str(design_matrix_x_data_task_4)
#design_matrix_x_data_task_4[1,]


y_data_task_5 <- data_task_5 %>% 
  select(Balance)
#str(y_data_task_5)


#
#split the data randomly
#

lambda_of_iteration <- 1
rows_task_5 <- 1:400

indices_resampled <- sample(x=rows_task_5,size=length(rows_task_5), replace=FALSE) #shuffle the indices
training_design <- design_matrix_x_data_task_5[indices_resampled[1:200],]
test_design <- design_matrix_x_data_task_5[indices_resampled[201:400],]


training_y <- y_data_task_5[indices_resampled[1:200],1]
test_y <- y_data_task_5[indices_resampled[201:400],1]

#
#train on training data:
#

#???Intercept True or False? Is in design already, do we need it with true?

lr_task_5 <- glmnet(x=training_design,y=training_y, 
                    alpha=1, intercept=FALSE, lambda=lambda_of_iteration )
#lr_est_task_5 <- lr_task_5$beta
#lr_est_task_5

#
#predictions
#

lr_prediction <- predict(lr_task_5, newx=test_design)
#L2-distance to true Balances , in a mean value divide by 200
lr_task_5_dist <- sqrt(sum((test_y-lr_prediction)^2))/(400/2)
lr_task_5_dist

#
#now in a loop:
#

M<-1000
counter <- 0 #tells us in which iteration stepo we are

L2_empirical_prediction_error_df <- 
  data.frame(Lambdas = grid_lambda, L2_empirical_prediction_error = rep(NA, length(grid_lambda)))

#str(L2_empirical_prediction_error_df)

for (lambda_of_iteration in grid_lambda){
  
  sum_L2_empirical_prediction_error <- 0
  counter <- counter+1
  
  for(iteration in 1:M){
    
    indices_resampled <- sample(x=rows_task_5,size=length(rows_task_5), replace=FALSE) #shuffle the indices
    
    training_design <- design_matrix_x_data_task_5[indices_resampled[1:200],]
    test_design <- design_matrix_x_data_task_5[indices_resampled[201:400],]
    training_y <- y_data_task_5[indices_resampled[1:200],1]
    test_y <- y_data_task_5[indices_resampled[201:400],1]
    
    #???Intercept True or False? Is in design already, do we need it with true?
    
    lr_task_5 <- glmnet(x=training_design,y=training_y, 
                        alpha=1, intercept=FALSE, lambda=lambda_of_iteration )
    
    lr_prediction <- predict(lr_task_5, newx=test_design)
    #L2-distance to true Balances , in a mean value divide by 200
    lr_task_5_dist <- sqrt(sum((test_y-lr_prediction)^2))/(400/2)
    sum_L2_empirical_prediction_error <- sum_L2_empirical_prediction_error + lr_task_5_dist
    
    
  }
  
  L2_empirical_prediction_error_df[counter,2] <- sum_L2_empirical_prediction_error/M
  #again divide by M for distance per data point
  
  #control of progress
  print(lambda_of_iteration)
}

#results:
str(L2_empirical_prediction_error_df)
#head(L2_empirical_prediction_error_df)

#plot
# Plot the data with ggplot and minimal theme
ggplot(L2_empirical_prediction_error_df, aes(x = grid_lambda, y = L2_empirical_prediction_error)) +
  geom_point() +  # Add points to the plot
  labs(title = "Lasso - L2_empirical_prediction_error", x = "lambda", y = "error") +
  theme_minimal() # Apply the minimal theme

#Select lambda
min_error_lambda <- 
  L2_empirical_prediction_error_df[which.min(L2_empirical_prediction_error_df$L2_empirical_prediction_error), ]
min_error_lambda

#--> lambda = 0.1 is optimal, which is strange, as larger lambda should prevent overfitting...

#
#
#Now: least squares
#
#

ls_sum_L2_empirical_prediction_error <- 0

for(iteration in 1:M){
  
  indices_resampled <- sample(x=rows_task_5,size=length(rows_task_5), replace=FALSE) #shuffle the indices
  
  training_data <- data_task_5[indices_resampled[1:200],]
  test_data <- subset(data_task_5[indices_resampled[201:400],], select=-Balance) #cut Balance away
  
  test_y <- y_data_task_5[indices_resampled[201:400],1]
  
  least_squares_task_5 <- lm(Balance~Income+Limit+Rating+Cards+Age
                             +Education+Own+Student+Married+Region, 
                             data = training_data)
  #ls_est_task_5 <- least_squares_task_5$coefficients
  #ls_est_task_5
  
  ls_prediction <- 
    predict(object=least_squares_task_5, newdata= test_data) 
  #ls_prediction
  
  #L2-distance to true Balances , in a mean value divide by 200
  ls_task_5_dist <- sqrt(sum((test_y-ls_prediction)^2))/(400/2)
  ls_sum_L2_empirical_prediction_error <- ls_sum_L2_empirical_prediction_error + ls_task_5_dist
  
  
}

ls_L2_empirical_prediction_error <- ls_sum_L2_empirical_prediction_error/M
#again divide by M for distance per data point


ls_L2_empirical_prediction_error

#7.20 was now smaller than for best lambda 0.1 the value 8.67



#
#
#Now: just intercept
#
#

intercept_sum_L2_empirical_prediction_error <- 0

for(iteration in 1:M){
  
  indices_resampled <- sample(x=rows_task_5,size=length(rows_task_5), replace=FALSE) #shuffle the indices
  
  training_data <- as.data.frame(data_task_5[indices_resampled[1:200],]$Balance) #only need y or Balance vector
  names(training_data)[1] <- "Balance"
  
  test_y <- y_data_task_5[indices_resampled[201:400],1]
  
  intercept_task_5 <- lm(Balance~1, 
                             data = training_data)
  intercept_est_task_5 <- intercept_task_5$coefficients
  intercept_est_task_5
  
  intercept_prediction <- 
    predict(object=intercept_task_5, newdata= as.data.frame(rep(0,200))) # newdata is just a length...
  #intercept_prediction
  
  #L2-distance to true Balances , in a mean value divide by 200
  intercept_task_5_dist <- sqrt(sum((test_y-intercept_prediction)^2))/(400/2)
  intercept_sum_L2_empirical_prediction_error <- intercept_sum_L2_empirical_prediction_error + intercept_task_5_dist
  
  
}

intercept_L2_empirical_prediction_error <- intercept_sum_L2_empirical_prediction_error/M
#again divide by M for distance per data point


intercept_L2_empirical_prediction_error

#32.56 was now not smaller than for best lambda 0.1 the value 8.67, only intercept is bad...


#
#Alternatively: intercept is just mean of the given balances, as mean minimizes least sqaures to a constant:
#

new_intercept_sum_L2_empirical_prediction_error <- 0

for(iteration in 1:M){
  
  indices_resampled <- sample(x=rows_task_5,size=length(rows_task_5), replace=FALSE) #shuffle the indices
  
  training_data <- as.data.frame(data_task_5[indices_resampled[1:200],]$Balance) #only need y or Balance vector
  names(training_data)[1] <- "Balance"
  
  test_y <- y_data_task_5[indices_resampled[201:400],1]
  
  mean_Balance <-  mean(training_data$Balance)
  
  intercept_prediction <- rep(mean_Balance,200)
  #intercept_prediction
  
  #L2-distance to true Balances , in a mean value divide by 200
  intercept_task_5_dist <- sqrt(sum((test_y-intercept_prediction)^2))/(400/2)
  new_intercept_sum_L2_empirical_prediction_error <- new_intercept_sum_L2_empirical_prediction_error + intercept_task_5_dist
  
  
}

new_intercept_L2_empirical_prediction_error <- new_intercept_sum_L2_empirical_prediction_error/M
#again divide by M for distance per data point


new_intercept_L2_empirical_prediction_error

#32.64 was now not smaller than for best lambda 0.1 the value 8.67, only intercept is bad...

########################################################################################
#trash that did not work

plots<-list()
for (i in 1:10){
  plots<- append(plots,ggplot(original_data) + 
    geom_point(aes(x = original_data[[i]], y = Balance), size=3) + 
    labs(title = paste("Scatter Plot: Balance for",names(original_data)[i] ),
         x = "explanatory variable", y = "average credit card debt") +
    theme_minimal())
  
}

#ggplot(original_data) + 
#  geom_point(aes(x = Income, y = Balance), size=3) + http://127.0.0.1:19433/graphics/c9f81be7-3131-4e63-b2ab-96cd0b80396c.png
#  labs(title = "Scatter Plot: Balance for Income",
#       x = "Income in 1000$", y = "average credit card debt")


grid.arrange(grobs=plots)
#grid.arrange(plots[1], plots[1], plots[3], plots[4],plots[5], plots[6],
 #                                          plots[7], plots[8],plots[9], plots[10],ncol = 5)

########################################################################################
 