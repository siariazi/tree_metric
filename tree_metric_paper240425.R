# R script to do scatter plot, optimization, correlations and pca 
library(ape)
library(treespace)
library(phangorn)
library(TreeDist)
library(tracerer)
library(ggplot2)
library(randomForest)
library(caret)
library(reshape2)
library(RColorBrewer)
library(distory)
library(tidyr)

library(patchwork)

# A data frame to keep all calculated metric distances and likelihood 
all_model_data <- data.frame()

# A dataframe to keep all the weights from optimization
all_weights <- data.frame()

# A data frame to keep all correlations
all_corr_data <- data.frame()

# Data frames to keep the PCA results
all_loadings <- data.frame()
all_loadings_long <- data.frame()
all_scores <- data.frame()
# File input directory
directory <- "/Users/siavashriazi/SFU/tree metrics/tree metrics codes"

# A directory for saving the plots
save_dir <- "/Users/siavashriazi/Dropbox/Apps/Overleaf/what-is-a-good-tree-metric/Figures"

# set working directory 
setwd(directory)


file_names <- c("RSV2","hcv_coal","Ireland_alpha","Ireland_delta")
#file_names <- c("RSV2","hcv_coal")

for (file_name in file_names){
  #file_name = "Ireland_delta"
  
  # reading trees directly from .trees file using read.nexus() from ape package
  beast_trees <- read.nexus(paste(file_name,".trees",sep=""),force.multi = TRUE)
  
  if (file_name=="Ireland_alpha" || file_name=="Ireland_delta"){
    beast_trees = lapply(beast_trees, 
                         function(tree) { 
                           tree$tip.label = paste("tip", order(tree$tip.label), sep="")
                           return(tree)})
  }
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
    
    kcmat = sprmat = rfmat = wrfmat = kfmat = pathmat = jrfmat = msdmat = nsmat = bhvmat = irfmat = msifmat = nnimat = likemat = matrix(NA,maxn,maxn)
    
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
        nsmat[i,j] = NyeSimilarity(beast_trees[[i+initial]],beast_trees[[j+initial]])
        bhvmat[i,j] = dist.multiPhylo(list(beast_trees[[i+initial]],beast_trees[[j+initial]]))
        irfmat[i,j] = InfoRobinsonFoulds(beast_trees[[i+initial]],beast_trees[[j+initial]])
        msifmat[i,j] = MatchingSplitInfoDistance(beast_trees[[i+initial]],beast_trees[[j+initial]])
        nnimat[i,j] = NNIDist(beast_trees[[i+initial]],beast_trees[[j+initial]])[2]
        likemat[i,j] = abs(beast_log$likelihood[i+initial] - beast_log$likelihood[j+initial])
      }
    }
    
    # Prepare the data
    data <- data.frame(Like = c(likemat), KC = c(kcmat), SPR = c(sprmat), 
                       RF = c(rfmat), WRF = c(wrfmat), IRF = c(irfmat), 
                       JRF = c(jrfmat), KF = c(kfmat), PD = c(pathmat), 
                       MSD = c(msdmat), NS = c(nsmat), BHV = c(bhvmat), MSI = c(msifmat),
                       NNI = c(nnimat)) 

    return(data)
  }

  initial = 200
  # this maxn should be 200 
  maxn = 200
  model_data <- calc_matrices(initial,maxn)
  plot(model_data$Like~model_data$KF,xlab="d(metrics)",ylab = "d(Likelihoods)",main=paste("wrf: ",initial,"_",initial+maxn,sep = ""))
  cor(model_data$WRF,model_data$Like,method = "pearson")
  plot(model_data$Like~model_data$KF,ylab = "d(Likelihoods)",main=paste(initial,"_",initial+maxn,sep = ""))
  
  # Assuming model_data is your dataframe
  model_data_long <- pivot_longer(model_data, 
                                  cols = colnames(model_data)[2:dim(model_data)[2]], 
                                  names_to = "Metric", 
                                  values_to = "Value")
  
  
  
  plot(beast_log$likelihood[1:800],xlab='Trees',ylab = 'Likelihood',main=paste("Trees: ",initial,"_",initial+maxn,sep = ""))
  abline(v=initial,col='red')
  abline(v=maxn+initial,col='red')
  
  # Use ggplot to create the plot
  ggplot(beast_log[1:800,], aes(x = Sample/500, y = likelihood)) +
    geom_line() +  # Plot the likelihood as a line
    geom_vline(xintercept = c(initial, initial + maxn), col = 'red') +  # Add vertical lines at initial2 and initial2 + maxn2
    labs(x = 'Trees', y = 'Likelihood', 
         title = paste("Trees: ", initial, "-", initial + maxn, sep = "")) +
    theme_minimal()
  
  ##############################
  # Optimization
  
  predict_likelihood_difference <- function(weights, metrics) {
    # Calculate the weighted sum of metrics
    prediction <- as.matrix(metrics) %*% weights
    return(prediction)
  }
  
  # Define the optimization function
  optimize_correlation <- function(weights) {
    metrics <- model_data[, -1] # Exclude the 'like' column
    predictions <- predict_likelihood_difference(weights, metrics)
    # We can use Pearson or Spearman correlation 
    corr <- cor(model_data$Like, predictions, method = "pearson")
    # Return negative correlation because optim() minimizes
    return(-corr)
  }
  
  # Initial weights
  initial_weights <- rep(1, ncol(model_data) - 1) # Assuming starting with equal weights
  
  optimize_correlation(initial_weights)
  
  # Optimize
  result <- optim(initial_weights, optimize_correlation, method = "L-BFGS-B")
  
  # Optimized weights
  optimized_weights <- result$par
  
  # Check the results
  optimized_weights
  
  weight_data <- data.frame(Metric = colnames(model_data[-1]),Weights = round(optimized_weights,2),Data=file_name)
  all_weights <- rbind(all_weights,weight_data)
  # Calculate predicted values
  fitmat <- as.matrix(model_data[,-1]) %*% optimized_weights
  
  optim_result <- cor(model_data$Like,c(fitmat))
  
  plot(model_data$Like~c(fitmat),xlab="Optimized Tree Distance",ylab = "Likelihood")
  
  model_data2 <- model_data
  model_data2$Optim <- c(fitmat)
  
  model_data2_long <- pivot_longer(model_data2, 
                                   cols = colnames(model_data2)[2:dim(model_data2)[2]], 
                                   names_to = "Metric", 
                                   values_to = "Value")
  model_data2_long$Data <- file_name
  all_model_data <- rbind(all_model_data,model_data2_long)
  ##############################
  # Calculate correlations of each metric with 'like'
  correlations <- sapply(model_data2[-1], function(x) cor(x, model_data2$Like))
  
  # Sort and visualize the correlations
  sorted_cor <- sort(correlations, decreasing = TRUE)
  
  # Simple bar plot of correlations
  barplot(sorted_cor, las = 2, main = "Correlation with Likelihood", ylab = "Correlation Coefficient", col = 'blue')
  
  
  # Create a data frame
  corr_data <- data.frame(Metric = names(correlations), Correlation = correlations)
  
  corr_data <- corr_data[order(-corr_data$Correlation),]
  
  
  corr_data$Metric <- factor(corr_data$Metric, levels = corr_data$Metric)
  corr_data$Data <- file_name
  all_corr_data <- rbind(all_corr_data,corr_data)
  
  ######################### PCA plot
  # Assuming 'data' is your dataframe
  # Perform PCA on the metrics, excluding the 'like' variable
  
  metrics <- model_data2[, -1]  # Exclude 'like'
  
  pca_result <- prcomp(metrics, scale. = TRUE)
  
  
  # Summary of PCA to see variance explained
  summary(pca_result)
  
  # Plot the PCA results
  biplot(pca_result, scale = 1)
  # plotting factors with autoplot
  #autoplot(pca_result,
  #         loadings = TRUE, loadings.label = TRUE, loadings.label.size  = 5)
  
  
  # Alternative way of plotting factors 
  
  # Extract Scores
  scores_df <- as.data.frame(pca_result$x)
  scores_df$Data <- file_name
  all_scores <- rbind(all_scores,scores_df)
  
  # Extract Loadings
  loadings_df <- as.data.frame(pca_result$rotation)
  loadings_df$Metric <- rownames(loadings_df)
  loadings_df$Data <- file_name
  all_loadings <- rbind(all_loadings,loadings_df)
  
  # Get loadings (rotation matrix)
  loadings <- pca_result$rotation
  # Convert the loadings to a long format for plotting
  loadings_long <- reshape2::melt(loadings)
  loadings_long$Data <- file_name
  all_loadings_long <- rbind(all_loadings_long,loadings_long)
  
  # closing the loop
}
########################################
# Function to modify names in data
modify_data <- function(df) {
  df$Data[df$Data == "RSV2"] <- "RSV"
  df$Data[df$Data == "hcv_coal"] <- "HCV"
  df$Data[df$Data == "Ireland_alpha"] <- "Alpha"
  df$Data[df$Data == "Ireland_delta"] <- "Delta"
  return(df)
}
all_corr_data <- modify_data(all_corr_data)
all_loadings <- modify_data(all_loadings)
all_loadings_long <- modify_data(all_loadings_long)
all_scores <- modify_data(all_scores)
all_model_data <- modify_data(all_model_data)
all_weights <- modify_data(all_weights)

wide_weights <- pivot_wider(all_weights, names_from = Data, values_from = Weights)
write.csv(wide_weights,file.path(directory,"all_weights.csv"),row.names = FALSE)
#########################################

plot_scatter_all <- ggplot(all_model_data, aes(x = Value, y = Like, color = Data)) + 
  geom_point(shape = 21, fill = NA,  alpha = 0.5, size = 0.5) +  # Create scatter plots
  facet_grid(Metric ~ Data, scales = "free") +  # Arrange plots in a row
  labs(x = "Metric Distance", y = "|log(L(T1) - log(L(T2)|") +  # General labels
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve readability of x-axis labels

print(plot_scatter_all)

# The order of metrics according to their correlation
desired_order <- c("Optim", "WRF", "BHV", "PD", "MSD", "KF", "KC", "MSI", "JRF", "NNI","RF", "IRF", "SPR", "NS")
# Convert 'Metric' to a factor with specified order
all_model_data2 <- all_model_data
all_model_data2$Metric <- factor(all_model_data2$Metric, levels = desired_order)

RSV_data <- subset(all_model_data2,Data == "RSV")
HCV_data <- subset(all_model_data2,Data == "HCV")
Alpha_data <- subset(all_model_data2,Data == "Alpha")
Delta_data <- subset(all_model_data2,Data == "Delta")

custom_palette <- c("RSV" = "cornflowerblue", "HCV" = "chartreuse2", "Alpha" = "firebrick3", "Delta" = "purple")

dot_size <- 0.05
dot_alpha <- 0.05
plot_scatter_rsv <- ggplot(RSV_data, aes(x = Value, y = Like, color = Data)) + 
  geom_point(shape = 21, fill = NA,  alpha = dot_alpha, size = dot_size) +  # Create scatter plots
  facet_wrap(~ Metric, scales = "free") +  # Arrange plots in a row
  labs(x = "Metric Distance", y = "|log(L(T1) - log(L(T2)|") +  # General labels
  theme_bw() +  scale_color_manual(values = custom_palette) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")  # Improve readability of x-axis labels
#print(plot_scatter_rsv)

plot_scatter_hcv <- ggplot(HCV_data, aes(x = Value, y = Like, color = Data)) + 
  geom_point(shape = 21, fill = NA, alpha = dot_alpha, size = dot_size) +  # Create scatter plots
  facet_wrap(~ Metric, scales = "free") +  # Arrange plots in a row
  labs(x = "Metric Distance", y = "|log(L(T1) - log(L(T2)|") +  # General labels
  theme_bw() + scale_color_manual(values = custom_palette) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")  # Improve readability of x-axis labels
#print(plot_scatter_hcv)

plot_scatter_alpha <- ggplot(Alpha_data, aes(x = Value, y = Like, color = Data)) + 
  geom_point(shape = 21, fill = NA, alpha = dot_alpha, size = dot_size) +  # Create scatter plots
  facet_wrap(~ Metric, scales = "free") +  # Arrange plots in a row
  labs(x = "Metric Distance", y = "|log(L(T1) - log(L(T2)|") +  # General labels
  theme_bw() + scale_color_manual(values = custom_palette) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")  # Improve readability of x-axis labels
#print(plot_scatter_alpha)

plot_scatter_delta <- ggplot(Delta_data, aes(x = Value, y = Like, color = Data)) + 
  geom_point(shape = 21, fill = NA, alpha = dot_alpha, size = dot_size) +  # Create scatter plots
  facet_wrap(~ Metric, scales = "free") +  # Arrange plots in a row
  labs(x = "Metric Distance", y = "|log(L(T1) - log(L(T2)|") +  # General labels
  theme_bw() + scale_color_manual(values = custom_palette) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")  # Improve readability of x-axis labels
#print(plot_scatter_delta)

# Combine the plots
combined_plot <- plot_scatter_rsv + plot_scatter_hcv + plot_scatter_alpha + plot_scatter_delta +
  plot_layout(guides = 'collect') & theme(legend.position = "none")

#print(combined_plot)

############
# Legend plot
# Data for creating the custom legend
legend_data <- data.frame(Data = c("RSV", "HCV", "Alpha", "Delta"))

# Create a custom legend plot
legend_plot <- ggplot(legend_data, aes(x = Data, y = 1, color = Data)) +
  geom_point(size = 5) +  # Large points for visibility
  scale_color_manual(values = custom_palette) +
  theme_void() +  # Remove all unnecessary plot elements
  theme(legend.position = "bottom") +
  labs(color = "Data") +  # Custom legend title
  guides(color = guide_legend(override.aes = list(size = 6)))

# Print just the legend
#print(legend_plot)

# Combine plots with the custom legend plot, 900-1000 is good
final_plot <- combined_plot / legend_plot + 
  plot_layout(heights = c(1, 0.1))  # Adjust the height ratio between plots and legend

# Print the final plot with custom legend
print(final_plot)
#ggsave(final_plot, filename = "plot_scatter_all.png", path = save_dir, width = 8, height = 5, dpi = 300)
#############################
# Create the histogram
plot_corr_all <- ggplot(all_corr_data, aes(x = Metric, y = Correlation, fill = Data)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +  # Use position_dodge to separate bars by Data within each Metric
  labs(x = "Metric", y = "Correlation") +
  scale_fill_manual(values = custom_palette)  +  # Optional: use a color palette that is distinct
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve x-axis label readability

# Display the plot
print(plot_corr_all)
ggsave(plot_corr_all, filename = "plot_corr_all.png", path = save_dir, width = 5, height = 4, dpi = 300)
################
plot_weights_all <- ggplot(all_weights, aes(x = Metric, y = Weights, fill = Data)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +  # Use position_dodge to separate bars by Data within each Metric
  labs(x = "Metric", y = "Weights") +
  scale_fill_manual(values = custom_palette)  +  # Optional: use a color palette that is distinct
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve x-axis label readability

plot_weights_all
#################
# Create base plot
plot_pca <- ggplot(data = all_scores, aes(x = PC1, y = PC2)) +
  geom_point(shape=1, alpha = 0.1,size = 0.5,color='grey') +  # Plot scores
  theme_bw() + facet_wrap(~ Data, scales = "free") +
  labs(x = "Principal Component 1", y = "Principal Component 2")

# Add loadings as longer, colored arrows
plot_pca <- plot_pca + geom_segment(data = all_loadings, aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5, color = Metric),
                                    arrow = arrow(type = "closed", length = unit(0.05, "inches")))

plot_pca
ggsave(plot_pca, filename = "plot_pca_all.png", path = save_dir, width = 5, height = 4, dpi = 300)

##############
# Create a heatmap: Contribution of Metrics to PCs
plot_loadings <- ggplot(all_loadings_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
  labs(x = "Metric", y = "Principal Component") +
  theme_bw()  + facet_wrap(~ Data, scales = "free") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Increase size and adjust angle
        axis.text.y = element_text(size = 10))  # Increase size of y-axis labels if needed

plot_loadings

ggsave(plot_loadings, filename = "plot_loadings_all.png", path = save_dir, width = 8, height = 5, dpi = 300)
#################################


