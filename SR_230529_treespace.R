# extracting trees from BEAST2 output 
library(ape)
library(treespace)
library(phangorn)
library(TreeDist)
library(tracerer)
library(ggplot2)
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

beast_log2 <- beast_log[1:200,]
beast_log2$"time" <- rep(seq(1,20),each=10)

# Define the color palette
color_palette <- brewer.pal(9, "Blues")

beast_log2$time <- as.factor(beast_log2$time)
# Create the plot
ggplot(beast_log2, aes(x = Sample, y = posterior, color = time)) +
  geom_point()  + labs(x = "Sample", y = "Posterior", title = "Sample vs. Posterior") + theme_bw() +
  geom_vline(xintercept = 24500, linetype = "dashed", color = "red")


plot(beast_log$posterior[1:200])
abline(v=50,col='red')

plot(beast_log$posterior[1:200]~beast_log$Sample[1:200],xlab='Sample',ylab = 'Posterior')
abline(v=24500,col='red')

initial = 50
maxn = 200
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

plot(c(likemat)~c(kcmat),xlab="Kendal-Colijn",ylab = "Likelihood")
cor(c(likemat),c(kcmat))
plot(c(likemat)~c(sprmat),xlab="Subtree Prune Regraft",ylab = "Likelihood")
plot(c(likemat)~c(rfmat),xlab="Robinson Foulds",ylab = "Likelihood")
plot(c(likemat)~c(wrfmat),xlab="Robinson Foulds Weighted",ylab = "Likelihood")
plot(c(likemat)~c(irfmat),xlab="Info Robinson Foulds",ylab = "Likelihood")
plot(c(likemat)~c(jrfmat),xlab="Jaccard Robinson Foulds",ylab = "Likelihood")
plot(c(likemat)~c(kfmat),xlab="Kuhner Felsenstien",ylab = "Likelihood")
plot(c(likemat)~c(pathmat),xlab="Path Difference",ylab = "Likelihood")
plot(c(likemat)~c(msdmat),xlab="Matching Split Distance",ylab = "Likelihood")
plot(c(likemat)~c(nymat),xlab="Nye Similarity",ylab = "Likelihood")
plot(c(likemat)~c(tdmat),xlab="Tree Distance",ylab = "Likelihood")
cor(c(likemat),c(tdmat))


###############
# putting all the metrics into a dataframe 
metrics <- c("Kendal-Colijn", "Subtree Prune Regraft", "Robinson Foulds",
             "Robinson Foulds Weighted", "Info Robinson Foulds", "Jaccard Robinson Foulds",
             "Kuhner Felsenstien", "Path Difference", "Matching Split Distance",
             "Nye Similarity", "Tree Distance")

# Define a vector of the variables representing your data for each metric
data <- list(kcmat, sprmat, rfmat, wrfmat, irfmat, jrfmat, kfmat, pathmat, msdmat, nymat, tdmat)

# Calculate correlation values

correlations <- sapply(data, function(d) cor(c(likemat), c(d)))

# Create a data frame
df <- data.frame(Metric = metrics, Correlation = correlations)

# Set rownames
rownames(df) <- df$Metric

# Remove the Metric column
#df$Metric <- NULL

# Write data frame to .csv file
write.csv(df, file = paste0(file_name, "_correlations.csv"), row.names = TRUE)
###############
# optimization

minfunc <- function(parms){
  # run comments if you want lasso regression  
  #lambda = 10
  total_sum = sum((kcmat*parms[1]-likemat)^2) + sum((sprmat*parms[2]-likemat)^2) + sum((rfmat*parms[3]-likemat)^2) +
  sum((wrfmat*parms[4]-likemat)^2) + sum((irfmat*parms[5]-likemat)^2) + sum((jrfmat*parms[6]-likemat)^2) +
  sum((kfmat*parms[7]-likemat)^2) + sum((pathmat*parms[8]-likemat)^2) + sum((msdmat*parms[9]-likemat)^2) +
  sum((nymat*parms[10]-likemat)^2) + sum((tdmat*parms[11]-likemat)^2) + parms[12]
  # regularzation part (LASSO regression)
  #lambda*sum(abs(parms[1])+abs(parms[2])+abs(parms[3])+abs(parms[4])+abs(parms[5])
  #   +abs(parms[6])+abs(parms[7])+abs(parms[8])+abs(parms[9])+abs(parms[10])+abs(parms[11]))
  return(total_sum)
  
}

par =  c(rep(0.1,12))
minfunc(par)
opt = optim(par, minfunc, method = "Nelder-Mead") # BFGS had the best performance 

df3 <- data.frame(Metric = c(metrics,"intercept"), Coefficients = opt$par)
write.csv(df3, file = paste0(file_name, "_optim2.csv"), row.names = TRUE)

optimat = matrix(NA,maxn,maxn)

for (i in seq(1,maxn)){
  for (j in seq(1,maxn)){
    #print(c(i,j))
    optimat[i,j] = treeDist(beast_trees[[i+initial]],beast_trees[[j+initial]])*opt$par[1] + 
                   SPR.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])*opt$par[2] +
                   RF.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])*opt$par[3] +
                   wRF.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])*opt$par[4] +
                   KF.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])*opt$par[5] +
                   path.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])*opt$par[6] +
                   JaccardRobinsonFoulds(beast_trees[[i+initial]],beast_trees[[j+initial]])*opt$par[7] +
                   MatchingSplitDistance(beast_trees[[i+initial]],beast_trees[[j+initial]])*opt$par[8] +
                   NyeSimilarity(beast_trees[[i+initial]],beast_trees[[j+initial]])*opt$par[9] +
                   InfoRobinsonFoulds(beast_trees[[i+initial]],beast_trees[[j+initial]])*opt$par[10] +
                   TreeDistance(beast_trees[[i+initial]],beast_trees[[j+initial]])*opt$par[11] + opt$par[12]

  }
}

plot(c(likemat)~c(optimat),xlab="Optimized Tree Distance",ylab = "Likelihood")
cor(c(likemat),c(optimat))

# using lm function
data <- data.frame(like = c(likemat), kc = c(kcmat), spr = c(sprmat), 
                   rf = c(rfmat), wrf = c(wrfmat), irf = c(irfmat), 
                   jrf = c(jrfmat), kf = c(kfmat), path = c(pathmat), 
                   ms = c(msdmat), ny = c(nymat), td = c(tdmat))
model <- lm(like ~ kc + spr + rf + wrf + irf + jrf + kf + path + ms + ny + td, data = data)
model <- lm(like ~ kc + spr + rf + wrf + irf + jrf + kf + path + ms, data = data)

summary(model)

# plotting the result of linear regression 
fitmat = matrix(NA,maxn,maxn)

for (i in seq(1,maxn)){
  for (j in seq(1,maxn)){
    #print(c(i,j))
      fitmat[i,j] =  
      treeDist(beast_trees[[i+initial]],beast_trees[[j+initial]])*as.numeric(model$coefficients[2]) +
      SPR.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])*as.numeric(model$coefficients[3]) +
      RF.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])*as.numeric(model$coefficients[4]) +
      wRF.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])*as.numeric(model$coefficients[5]) +
      InfoRobinsonFoulds(beast_trees[[i+initial]],beast_trees[[j+initial]])*as.numeric(model$coefficients[6]) 
      JaccardRobinsonFoulds(beast_trees[[i+initial]],beast_trees[[j+initial]])*as.numeric(model$coefficients[7]) +
      KF.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])*as.numeric(model$coefficients[8]) +
      path.dist(beast_trees[[i+initial]],beast_trees[[j+initial]])*as.numeric(model$coefficients[9]) +
      MatchingSplitDistance(beast_trees[[i+initial]],beast_trees[[j+initial]])*as.numeric(model$coefficients[10]) +
      as.numeric(model$coefficients[1])
      #NyeSimilarity(beast_trees[[i+initial]],beast_trees[[j+initial]])*as.numeric(model$coefficients[12]) +
      #TreeDistance(beast_trees[[i+initial]],beast_trees[[j+initial]])*as.numeric(model$coefficients[13])
       
  }
}

metrics2 <- c("intercept",metrics[0:9])
df2 <- data.frame(Metric = metrics2, coefficients = unname(model$coefficients))
# Set rownames
rownames(df2) <- df2$Metric
# Remove the Metric column
#df2$Metric <- NULL
# Write data frame to .csv file
write.csv(df2, file = paste0(file_name, "_coefficients_lm.csv"), row.names = TRUE)

plot(c(likemat)~c(fitmat),xlab="Optimized Tree Distance (lm)",ylab = "Likelihood")
cor(c(likemat),c(fitmat))

heatmap(likemat,Colv = NA, Rowv = NA)
heatmap(kcmat,Colv = NA, Rowv = NA)
heatmap(sprmat,Colv = NA, Rowv = NA)
heatmap(rfmat,Colv = NA, Rowv = NA)
heatmap(wrfmat,Colv = NA, Rowv = NA)
heatmap(irfmat,Colv = NA, Rowv = NA)
heatmap(jrfmat,Colv = NA, Rowv = NA)
heatmap(kfmat,Colv = NA, Rowv = NA)
heatmap(pathmat,Colv = NA, Rowv = NA)
heatmap(msdmat,Colv = NA, Rowv = NA)
heatmap(nymat,Colv = NA, Rowv = NA)
heatmap(tdmat,Colv = NA, Rowv = NA)

#https://thibautjombart.github.io/treespace/articles/DengueVignette.html
#######################################

beast_trees2 <- beast_trees[1:200]
BEASTscape <- treespace(beast_trees2,method="patristic", nf=2, lambda=0, return.tree.vectors=FALSE, processors=1)
plotGrovesD3(BEASTscape$pco)

#####################
# find clusters or 'groves':
BEASTGroves <- findGroves(BEASTscape, nclust=4, clustering = "single")

# find median tree(s) per cluster:
BEASTMeds <- medTree(beast_trees2, groups=BEASTGroves$groups)
# for each cluster, select a single median tree to represent it:
BEASTMedTrees <- c(BEASTMeds$`1`$trees[[1]],
                   BEASTMeds$`2`$trees[[1]],
                   BEASTMeds$`3`$trees[[1]],
                   BEASTMeds$`4`$trees[[1]])


# extract the numbers from the tree list 'BEASTtrees' which correspond to the median trees: 
BEASTMedTreeNums <-c(which(BEASTGroves$groups==1)[[BEASTMeds$`1`$treenumbers[[1]]]],
                     which(BEASTGroves$groups==2)[[BEASTMeds$`2`$treenumbers[[1]]]],
                     which(BEASTGroves$groups==3)[[BEASTMeds$`3`$treenumbers[[1]]]],
                     which(BEASTGroves$groups==4)[[BEASTMeds$`4`$treenumbers[[1]]]])
# prepare a vector to highlight median and MCC trees
highlightTrees <- rep(1,201)
highlightTrees[[201]] <- 2
highlightTrees[BEASTMedTreeNums] <- 2
# prepare colours:
BEASTcols <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3")

# plot:
plotGrovesD3(BEASTscape$pco,
             groups=as.vector(BEASTGroves$groups),
             colors=BEASTcols,
             col_lab="Cluster")

####################################
# find clusters or 'groves':
BEASTGroves <- findGroves(BEASTscape, nclust=3, clustering = "single")

# find median tree(s) per cluster:
BEASTMeds <- medTree(beast_trees2, groups=BEASTGroves$groups)
# for each cluster, select a single median tree to represent it:
BEASTMedTrees <- c(BEASTMeds$`1`$trees[[1]],
                   BEASTMeds$`2`$trees[[1]],
                   BEASTMeds$`3`$trees[[1]])


# extract the numbers from the tree list 'BEASTtrees' which correspond to the median trees: 
BEASTMedTreeNums <-c(which(BEASTGroves$groups==1)[[BEASTMeds$`1`$treenumbers[[1]]]],
                     which(BEASTGroves$groups==2)[[BEASTMeds$`2`$treenumbers[[1]]]],
                     which(BEASTGroves$groups==3)[[BEASTMeds$`3`$treenumbers[[1]]]])
# prepare a vector to highlight median and MCC trees
highlightTrees <- rep(1,201)
highlightTrees[[201]] <- 2
highlightTrees[BEASTMedTreeNums] <- 2
# prepare colours:
BEASTcols <- c("#66c2a5","#fc8d62","#8da0cb")

# plot:
plotGrovesD3(BEASTscape$pco,
             groups=as.vector(BEASTGroves$groups),
             colors=BEASTcols,
             col_lab="Cluster")
#############################
# find clusters or 'groves':
BEASTGroves <- findGroves(BEASTscape, nclust=2, clustering = "single")

# find median tree(s) per cluster:
BEASTMeds <- medTree(beast_trees2, groups=BEASTGroves$groups)
# for each cluster, select a single median tree to represent it:
BEASTMedTrees <- c(BEASTMeds$`1`$trees[[1]],
                   BEASTMeds$`2`$trees[[1]])


# extract the numbers from the tree list 'BEASTtrees' which correspond to the median trees: 
BEASTMedTreeNums <-c(which(BEASTGroves$groups==1)[[BEASTMeds$`1`$treenumbers[[1]]]],
                     which(BEASTGroves$groups==2)[[BEASTMeds$`2`$treenumbers[[1]]]])
# prepare a vector to highlight median and MCC trees
highlightTrees <- rep(1,201)
highlightTrees[[201]] <- 2
highlightTrees[BEASTMedTreeNums] <- 2
# prepare colours:
BEASTcols <- c("#66c2a5","#8da0cb")

# plot:
plotGrovesD3(BEASTscape$pco,
             groups=as.vector(BEASTGroves$groups),
             colors=BEASTcols,
             col_lab="Cluster")
############################# Mahsa's plot
res<-treespace(beast_trees2, method = "patristic" ,nf=3)
table.paint(as.matrix(res$D))
scatter(res$pco)
plot(res$pco$li[,1],res$pco$li[,2])
resdf <- res$pco$li
resdf$"time" <- rep(seq(1,20),each=10)
if(require(ggplot2)){
  resplot <- ggplot(resdf, aes(x=A1, y=A2)) # create plot
  resplot + geom_density2d(colour="gray80") + # contour lines
    geom_point(size=6, shape=1, colour="gray50") + # grey edges
    geom_point(size=6, alpha=0.2, colour="navy") + # transparent blue points
    xlab("") + ylab("") + theme_bw(base_family="") # remove axis labels and grey background
}

resdf$"time" <- rep(seq(1,20),each=10)
resdf$time <- as.factor(resdf$time)

if(require(ggplot2)){
  resplot <- ggplot(resdf, aes(x = A1, y = A2)) +
    geom_density2d(aes(colour = time)) +  # Specify the colour aesthetic within aes()
    geom_point(size = 6, shape = 1, aes(colour = time)) +
    scale_colour_manual(values = rainbow(length(unique(resdf$time)))) +  # Specify color palette
    xlab("") + ylab("") + theme_bw(base_family = "") 
  print(resplot)
}

if (require(ggplot2)) {
      resdf$time <- factor(resdf$time, levels = unique(resdf$time))  # Specify the order of levels
       
         resplot <- ggplot(resdf, aes(x = A1, y = A2, colour = time)) +
         geom_point(size = 6, shape = 1) +
         scale_colour_manual(values = rainbow(length(unique(resdf$time))),
                                                                 breaks = levels(resdf$time),
                                                                 labels = levels(resdf$time)) +
         xlab("") + ylab("") + theme_bw(base_family = "")
        
        print(resplot)
}


# plotting subset
resdf_sub <- resdf[50:nrow(resdf),]

if (require(ggplot2)) {
  resplot <- ggplot(resdf_sub, aes(x = A1, y = A2, colour = time)) +
    geom_point(size = 6, shape = 1) +
    scale_colour_manual(values = rainbow(length(unique(resdf$time))),
                        breaks = levels(resdf$time),
                        labels = levels(resdf$time)) +
    xlab("") + ylab("") + theme_bw(base_family = "")
  
  print(resplot)
}

if (require(ggplot2)) {
  resdf_sub$time <- factor(resdf_sub$time, levels = unique(resdf_sub$time))  # Specify the order of levels
  
  resplot <- ggplot(resdf_sub, aes(x = A1, y = A2, colour = time)) +
    geom_point(size = 6, shape = 1) +
    scale_colour_manual(values = rainbow(length(unique(resdf$time))),
                        breaks = levels(resdf$time),
                        labels = levels(resdf$time)) +
    xlab("") + ylab("") + theme_bw(base_family = "")
  
  print(resplot)
}
## Not run:
if(require(rgl)){
  plot3d(resdf[,1], resdf[,2], resdf[,3], type="s", size=1.5,
         col="navy", alpha=0.5, xlab="", ylab="", zlab="")
}

############################# file conversion
file_name_log_csv <- paste(file_name,"_log.csv",sep = "")

#write.csv(beast_log,file=paste('/Users/siavashriazi/Desktop/SFU/BEAST files/',file_name_log_csv,sep = ""),na='')
fileName <- file.path(directory, file_name_log_csv)
write.csv(beast_log,file=fileName)

beast_log2 <- read.csv(file = file_name_log_csv,header = TRUE)
beast_log2 <- beast_log2[,-1]
############################################# 
# log file has more rows that trees so I filter log and keep only values that are present in trees
library(stringr)

# extracing all state values
state_values <- sapply(beast_trees, function(x) attr(x, "STATE"))

# converting state_values to a list
list_state_values <- list(state_values)

# converting the list of state_values to a string so I can process the string with stringr package
string_values <- as.character(list_state_values)

# Extract numbers from states
numbers <- str_extract_all(string_values, "(?<=STATE_)[0-9]+(?=\\s*=)")

# converting numbers to integers
numbers2 <- as.numeric(unlist(numbers[1]))

# filtering log file to only keep the values which we have tree for those
beast_log_filtered <- beast_log[beast_log$Sample %in% numbers2, ]

# example, if you want to extract likelihood of the first tree 
i=1

# unname discard the header of dataframe, column 3 is likelihood
unname(beast_log_filtered[i,])[3]

################### heatmap
hmcol <- colorRampPalette(brewer.pal(11, "RdBu"))(20)
library(RColorBrewer)
library(gplots)
heatmap.2(rfmat,col=hmcol,
          scale="none", trace="none", dendrogram = "none",symbreaks = T, 
          breaks = seq(-100,100,length.out = 21), 
          keysize=2, srtCol=45, Colv=FALSE,labRow=F,cexCol = 1)

# metrics from phangorn
treeDist(i,j) # Kendal-Colijn distance
#treedist(i, j) # 
SPR.dist(i, j)
#sprdist(i,j)
RF.dist(i,j)
wRF.dist(i,j)
KF.dist(i,j)
path.dist(i,j)

# metrics from TreeDist
PathDist(i,j) # same as path.dist(i,j)
JaccardRobinsonFoulds(i,j)
MatchingSplitDistance(i,j)
NNIDist(i,j)
MASTSize(i,j,rooted=TRUE)
NyeSimilarity(i,j)
InfoRobinsonFoulds(i,j)
RobinsonFoulds(i,j) # same as RF.dist()
#RobinsonFouldsMatching(i,j)
#RobinsonFouldsSplits(i,j)
SPRDist(i,j) # same as spr.dist()
TreeDistance(i,j)
