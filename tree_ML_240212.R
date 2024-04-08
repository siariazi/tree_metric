# extracting trees from BEAST2 output 
library(ape)
library(treespace)
#library(phangorn)
#library(TreeDist)
library(tracerer)
library(ggplot2)
library(randomForest)
library(caret)
library(reshape2)
library(RColorBrewer)

directory <- "/Users/siavashriazi/SFU/tree metrics/tree metrics codes"


# set working directory 
setwd(directory)

file_name = "RSV2"

# reading trees directly from .trees file using read.nexus() from ape package
beast_trees <- read.nexus(paste(file_name,".trees",sep=""),force.multi = TRUE)
# parsing log file
file_name_log <- paste(file_name,".log",sep = "")  
# Get the path of the 'inst/extdata' directory in the tracerer package
extdata_path <- system.file("extdata", package = "tracerer")

# Set the path of the log file within the 'inst/extdata' directory
new_file_path_log <- file.path(extdata_path, file_name_log)

# Copy the log file to 'inst/extdata' directory
file.copy(from = file.path(directory, file_name_log), to = new_file_path_log, overwrite = TRUE)

# Obtain the full path of the log file within the tracerer package
filename_log <- get_tracerer_path(file_name_log)
# Parse that log file
beast_log <- parse_beast_tracelog_file(filename_log)

calc_matrices <- function(initial,maxn){

  kcmat = sprmat = rfmat = wrfmat = kfmat = pathmat = jrfmat = msdmat = nymat = irfmat = tdmat = likemat = matrix(NA,maxn,maxn)
  
  for (i in seq(1,maxn)){
    for (j in seq(1,maxn)){
      #print(c(i,j))
      kcmat[i,j] = treeDist(beast_trees[[i+initial]],beast_trees[[j+initial]])
      sprmat[i,j] = SPR.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])
      rfmat[i,j] = RF.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])
      wrfmat[i,j] = wRF.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])
      kfmat[i,j] = KF.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])
      pathmat[i,j] = path.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])
      jrfmat[i,j] = JaccardRobinsonFoulds(beast_trees[[i+initial]],beast_trees[[j+initial]])
      msdmat[i,j] = MatchingSplitDistance(beast_trees[[i+initial]],beast_trees[[j+initial]])
      nymat[i,j] = NyeSimilarity(beast_trees[[i+initial]],beast_trees[[j+initial]])
      irfmat[i,j] = InfoRobinsonFoulds(beast_trees[[i+initial]],beast_trees[[j+initial]])
      tdmat[i,j] = TreeDistance(beast_trees[[i+initial]],beast_trees[[j+initial]])
      
      likemat[i,j] = abs(beast_log$likelihood[i+initial] - beast_log$likelihood[j+initial])
    }
  }
  
  # Prepare the data
  data <- data.frame(like = c(likemat), kc = c(kcmat), spr = c(sprmat), 
                     rf = c(rfmat), wrf = c(wrfmat), irf = c(irfmat), 
                     jrf = c(jrfmat), kf = c(kfmat), path = c(pathmat), 
                     ms = c(msdmat), ny = c(nymat), td = c(tdmat))

  
  return(data)
}
initial = 50
maxn = 200
model_data <- calc_matrices(initial,maxn)
plot(model_data$like~model_data$wrf,xlab="d(metrics)",ylab = "d(Likelihoods)",main=paste("wrf: ",initial,"_",initial+maxn,sep = ""))
plot(model_data$like~model_data$jrf,ylab = "d(Likelihoods)",main=paste(initial,"_",initial+maxn,sep = ""))

plot(beast_log$posterior[1:800],xlab='Samples',ylab = 'Posterior',main=paste("Samples: ",initial,"_",initial+maxn,sep = ""))
abline(v=initial,col='red')
abline(v=maxn+initial,col='red')

initial2 = 200
maxn2 = 200
model_data2 <- calc_matrices(initial2,maxn2)
plot(model_data2$like~model_data2$kf,xlab="d(metrics)",ylab = "d(Likelihoods)",main=paste("wrf: ",initial2,"_",initial2+maxn2,sep = ""))
cor(model_data2$wrf,model_data2$like,method = "pearson")
plot(model_data2$like~model_data2$kf,ylab = "d(Likelihoods)",main=paste(initial2,"_",initial2+maxn2,sep = ""))


# Assuming model_data2 is your dataframe
model_data_long <- pivot_longer(model_data2[,c("like","kc","wrf","kf")], 
                                cols = c(kc, kf, wrf), 
                                names_to = "Metric", 
                                values_to = "Value")


# Assuming model_data2 is your dataframe
model_data_long <- pivot_longer(model_data2[,c("like","wrf","path","ms")], 
                                cols = c(wrf,path, ms), 
                                names_to = "Metric", 
                                values_to = "Value")


ggplot(model_data_long, aes(x = Value, y = like)) + 
  geom_point(shape = 21, fill = NA,  alpha = 0.2) +  # Create scatter plots
  facet_wrap(~ Metric, scales = "free", nrow = 1) +  # Arrange plots in a row
  labs(x = "d(Metrics)", y = "d(Likelihoods)") +  # General labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve readability of x-axis labels


plot(beast_log$likelihood[1:800],xlab='Trees',ylab = 'Likelihood',main=paste("Trees: ",initial2,"_",initial2+maxn2,sep = ""))
abline(v=initial2,col='red')
abline(v=maxn2+initial2,col='red')

# Use ggplot to create the plot
ggplot(beast_log[1:800,], aes(x = Sample/500, y = likelihood)) +
  geom_line() +  # Plot the likelihood as a line
  geom_vline(xintercept = c(initial2, initial2 + maxn2), col = 'red') +  # Add vertical lines at initial2 and initial2 + maxn2
  labs(x = 'Trees', y = 'Likelihood', 
       title = paste("Trees: ", initial2, "-", initial2 + maxn2, sep = "")) +
  theme_minimal()


initial3 = 400
maxn3 = 200
# kf is also good
model_data3 <- calc_matrices(initial3,maxn3)
plot(model_data3$like~model_data3$wrf,xlab="d(metrics)",ylab = "d(Likelihoods)",main=paste("wrf: ",initial3,"_",initial3+maxn3,sep = ""))

initial4 = 600
maxn4 = 200
# kf is also good
model_data4 <- calc_matrices(initial4,maxn4)
plot(model_data4$like~model_data4$wrf,xlab="d(metrics)",ylab = "d(Likelihoods)",main=paste("wrf: ",initial4,"_",initial4+maxn4,sep = ""))

initial5 = 200
maxn5 = 600
# kf is also good
#model_data5 <- calc_matrices(initial5,maxn5)
plot(model_data5$like~model_data5$wrf,xlab="d(metrics)",ylab = "d(Likelihoods)",main=paste("wrf: ",initial5,"_",initial5+maxn5,sep = ""))

plot(beast_log$posterior[1:800],xlab='Samples',ylab = 'Posterior',main=paste("Samples: ",initial5,"_",initial5+maxn5,sep = ""))
abline(v=initial5,col='red')
abline(v=maxn5+initial5,col='red')


######## PCA plot
# Assuming 'data' is your dataframe
# Perform PCA on the metrics, excluding the 'like' variable
model_data22 <- model_data2[,-c(11)] # removing ny 
model_data23 <- model_data22[,c(1,5,9,10)] # picking wrf, kf and kc as they showed the highest correlation

pca_data <- model_data23
metrics <- pca_data[, -1]  # Exclude 'like'

# Assuming your dataframe of metrics is named 'metrics'
melted_metrics <- melt(metrics, variable.name = "Metric", value.name = "Value")

pca_result <- prcomp(metrics, scale. = TRUE)

# Get PCA scores
scores <- as.data.frame(pca_result$x)

# Add an identifier for each observation to merge with melted metrics
scores$ID <- rownames(scores)

melted_metrics$ID <- rep(rownames(metrics), each = ncol(metrics))

# Merge melted metrics with PCA scores based on the ID
plot_data <- merge(melted_metrics, scores, by = "ID")

# Define a color palette
color_palette <- brewer.pal(n = length(unique(plot_data$Metric)), name = "Set3")

# Plot PCA with the chosen color palette
ggplot(plot_data, aes(x = PC1, y = PC2, color = Metric)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = color_palette) +   # Apply the color palette
   labs(title = "PCA of Metrics Colored by Metric Name",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

# Plot PCA with outlined points (borders only, no fill)
ggplot(plot_data, aes(x = PC1, y = PC2, color = Metric)) +
  geom_point(shape = 21, fill = NA, size = 2) +  # Shape 21 with no fill
  labs(title = "PCA of Metrics Colored by Metric Name",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

############################



# Summary of PCA to see variance explained
summary(pca_result)

# Plot the PCA results
biplot(pca_result, scale = 0)

# Calculate correlations of each metric with 'like'
correlations <- sapply(pca_data[-1], function(x) cor(x, pca_data$like))

# Sort and visualize the correlations
sorted_cor <- sort(correlations, decreasing = TRUE)

# Simple bar plot of correlations
barplot(sorted_cor, las = 2, main = "Correlation with Likelihood", ylab = "Correlation Coefficient", col = 'blue')


# Assuming 'metrics' contains your data excluding the 'like' column
pca_result <- prcomp(metrics, scale. = TRUE)


###########


# Split data into features and target
model_features <- model_data2[, -1]  # Exclude the 'like' column
target <- model_data2$like

# Train a Random Forest model
rf_model <- randomForest(x = model_features, y = target, ntree = 500)

# Check the model's importance scores to see which metrics are most predictive
importance(rf_model)

# Predict and evaluate using your preferred method; here, we use the model on the training data for simplicity
predictions <- predict(rf_model, model_features)

# For correlation, though not the primary goal, it's still informative
cor(target, predictions, method = "spearman")

# Assuming 'predictions' are your model's predictions and 'data$like' are the actual values
model_data2$predictions <- predictions  # Add predictions to your dataframe for plotting

ggplot(model_data2, aes(x = like, y = predictions)) +
  geom_point(alpha = 0.5) +  # Plot actual vs predicted values
  geom_abline(color = "red", linetype = "dashed") +  # Add a y=x reference line
  labs(x = "Actual Values", y = "Predicted Values", title = "Actual vs. Predicted") +
  theme_minimal()

# Plot variable importance
varImpPlot(rf_model, type = 2, main = "Feature Importance")

# cross validation the results
test_data <- calc_matrices(initial2+maxn2+1,400)

# Prepare test features and actual values
test_features <- test_data[, -1]  # Exclude the 'like' column
test_actual <- test_data$like

# Predict on test data
test_predictions <- predict(rf_model, test_features)

# Calculate Spearman's rank correlation for the test set
test_correlation <- cor(test_actual, test_predictions, method = "spearman")

# Output the correlation
print(test_correlation)
###########################
# Assuming the randomForest package is being used
model_features2 <- model_features[, c(1,4,7)]  # Exclude the 'like' column
test_features2 <- test_features[,c(1,4,7)]
rf_model_adjusted <- randomForest(x = model_features2, y = target, ntree = 500, 
                                  maxnodes = 1000,  # Limit the depth of trees
                                  mtry = round(sqrt(ncol(model_features2))))  # Use a subset of features

# Re-evaluate on training data
predictions_adjusted <- predict(rf_model_adjusted, model_features2)
cor(target, predictions_adjusted, method = "spearman")

# Re-test on test data
test_predictions_adjusted <- predict(rf_model_adjusted, test_features2)
test_correlation_adjusted <- cor(test_actual, test_predictions_adjusted, method = "spearman")
print(test_correlation_adjusted)
##########################
library(xgboost)

# Assuming model_data2 is your training set and model_data5 is your test set
data_matrix <- xgb.DMatrix(data = as.matrix(model_data2[, c(-1,-13)]), label = model_data2$like)
test_matrix <- xgb.DMatrix(data = as.matrix(model_data5[, -1]), label = model_data5$like)

# Parameters for the xgboost model
params <- list(booster = "gbtree", objective = "reg:squarederror", eta = 0.1, max_depth = 6, subsample = 0.5, colsample_bytree = 0.5)

# Training the model
model_xgb <- xgb.train(params = params, data = data_matrix, nrounds = 100, verbose = 0)

# Predicting on test data
predictions_xgb <- predict(model_xgb, test_matrix)

# Calculate Spearman's rank correlation for the test set
test_correlation_xgb <- cor(model_data5$like, predictions_xgb, method = "spearman")

# Output the correlation
print(test_correlation_xgb)


#########################
# Define a parameter grid (example)
params_grid <- list(
  eta = c(0.01, 0.05, 0.1),
  max_depth = c(3, 6, 9),
  subsample = c(0.5, 0.7, 1),
  colsample_bytree = c(0.5, 0.7, 1)
)

# Use a validation set or cross-validation to evaluate parameter combinations
# Example: This step is manual but should be automated with loops or using caret/tune/gridExtra packages

# Adjust based on results
best_params <- list(
  eta = 0.05,
  max_depth = 6,
  subsample = 0.7,
  colsample_bytree = 0.7,
  objective = "reg:squarederror"
)

# Train with best parameters found
final_model <- xgb.train(params = best_params, data = data_matrix, nrounds = 1000, watchlist = list(val = test_matrix), early_stopping_rounds = 10)

new_predictions <- predict(final_model, test_matrix)

# Calculate Spearman's rank correlation for the test set
test_correlation_xgb <- cor(model_data5$like, new_predictions, method = "spearman")

# Output the correlation
print(test_correlation_xgb)
########################## caret
# Define control parameters for cross-validation
train_control <- trainControl(method = "boot", number = 5, search = "random")

# Train the model with cross-validation
model_cv <- train(x = features, y = target,
                  method = "rf",  # Specify Random Forest
                  trControl = train_control,
                  tuneLength = 3)  # Number of tuning parameters to try

# Output results
print(model_cv)