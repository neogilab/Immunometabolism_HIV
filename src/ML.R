library(RColorBrewer)
library(pROC)
library(scales)
library(randomForest)
library(caret)
library(Boruta)
library(ROCR)
library(gplots)

#' Save randomForest algorithm and associated features selection
#'
#'@importFrom randomForest randomForest
#'
#' @param matrix matrix containing data and first column conditions
#' @param title title of analysis
#' @return random forest object
#' @export
RF_FS <- function(matrix, title){
  dir.create(paste0("results/ML/", title))
  # select number of trees
  model_trees(matrix, title)
  nb_trees <- readline(prompt="Enter Number of trees: ")
  nb_trees <- as.integer(nb_trees) # convert character into integer
  print(paste("You have selected ", nb_trees, " trees for this analysis.", sep = ""))
  # select important features with botura package
  imp_features <- B_feature_select(nb_trees, matrix, title)
  print(paste(" The number of important features for this analysis is ", length(imp_features), ".", sep = ""))
  # Select only important features expression data
  final_data = matrix[, c(imp_features, "Condition")]
  #  Data partition for training and test of model to be built
  Index = createDataPartition(final_data$Condition, p=0.90, list = FALSE, times = 1)
  Train = final_data[Index, ]
  print(paste("Number of samples for training is ", nrow(Train), ".", sep = ""))
  Test = final_data[-Index, ]
  print(paste("Number of samples for test is ", nrow(Test), ".", sep = ""))
  # select mtry parameter
  mtry_selection(Train, nb_trees, title)
  nb_mtry <- readline(prompt="Enter Number of mtry: ")
  nb_mtry <- as.integer(nb_mtry) # convert character into integer
  print(paste("You have selected ", nb_mtry, " as parameter mtry for this analysis.", sep = ""))
  # build final model
  Test$Condition <- factor(Test$Condition)
  Train$Condition <- factor(Train$Condition)
  RFModel.Final = randomForest(Condition~.,data = Train, mtry= nb_mtry,ntree=nb_trees)
  RFModel.Final
  path_rf_model <- paste0("results/ML/", title, "/RF_model_classification_", title, ".RData")
  save(RFModel.Final, file= path_rf_model)
  path_features <- paste0("results/ML/", title, "/Important_feature_", title, ".pdf")
  pdf(path_features, width=8,height=6,paper='special')
  randomForest::varImpPlot(RFModel.Final)
  dev.off()
  # Fit the model on the test dataset to predict Condition
  print("Now we are going to test the model.")
  test.forest = predict(RFModel.Final, type = "class", newdata = Test)
  # show confusion matrix
  cfm <- confusionMatrix(Test$Condition,test.forest)
  print(cfm)
  path_cfm <- paste0("results/ML/", title, "/Confusion_matrix_", title, ".pdf")
  ggplotConfusionMatrix(cfm)
  ggsave(path_cfm)
  return(RFModel.Final)
}


#' Save randomForest algorithm and associated features selection
#'
#'@importFrom randomForest randomForest
#'
#' @param matrix matrix containing data and first column conditions
#' @param title title of analysis
#' @return random forest object
#' @export
feature_select_Boruta <- function(matrix, title){
  dir.create(paste0("results/ML/", title))
  # select number of trees
  model_trees(matrix, title)
  nb_trees <- readline(prompt="Enter Number of trees: ")
  nb_trees <- as.integer(nb_trees) # convert character into integer
  print(paste("You have selected ", nb_trees, " trees for this analysis.", sep = ""))
  # select important features with botura package
  imp_features <- B_feature_select(nb_trees, matrix, title)
  print(paste(" The number of important features for this analysis is ", length(imp_features), ".", sep = ""))
  # Select only important features expression data
  final_data = matrix[, c(imp_features, "Condition")]
  return(final_data)
}


#' Save randomForest algorithm and associated features selection
#'
#'@importFrom randomForest randomForest
#'
#' @param matrix matrix containing data and first column conditions
#' @param title title of analysis
#' @return random forest object
#' @export
feature_select_Boruta_2 <- function(matrix, title){
  dir.create(paste0("results/ML/", title))
  # select number of trees
  nb_trees <- 1000
  nb_trees <- as.integer(nb_trees) # convert character into integer
  # select important features with botura package
  imp_features <- B_feature_select(nb_trees, matrix, title)
  print(paste(" The number of important features for this analysis is ", length(imp_features), ".", sep = ""))
  # Select only important features expression data
  final_data = matrix[, c(imp_features, "Condition")]
  return(final_data)
}



#' Print and save model to select number of trees
#'
#' @param matrix matrix containing data and first column conditions
#' @param title title of analysis
#' @export
model_trees <- function(matrix, title){
  RFModel1 = randomForest(Condition~., data = matrix, ntree=1000, importance = TRUE)
  path_model <- paste0("results/ML/", title, "/Model_trees_", title, ".pdf")
  plot(RFModel1)
  dev.copy(x11)
  dev.print(pdf, path_model)
}

#' Feature selection with Boruta
#'
#'#'@import Boruta
#'
#' @param nb_trees Number of trees selected for the analysis
#' @param matrix matrix containing data and first column conditions
#' @param title title of analysis
#' @return List of features selected in order to differentiate this two conditions
#' @export
B_feature_select <- function(nb_trees, matrix, title){
  BT = Boruta(Condition~., ntree = nb_trees, data = matrix)
  SelectedAttributes = getSelectedAttributes(BT)
  path_features <- paste0("results/ML/", title, "/SelectedAttributes_", title, ".csv")
  write.csv(SelectedAttributes, file = path_features)
  a=read.csv(path_features)
  a=a[, -1]
  a=as.character(a)
  return(a)
}

#' Select mtry according to OOB error curve
#'
#' @param Train matrix containing training data set
#' @param nb_trees Number of trees selected 
#' @param title title of analysis
#' @export
mtry_selection <- function(Train, nb_trees, title){
  set.seed(234)
  mtry=tuneRF(Train[,-3],
              Train[,3],
              mtryStart = 5,
              ntreeTry = nb_trees,
              stepFactor = 1.5, 
              improve = 0.001, 
              trace=TRUE, 
              plot = TRUE,
              doBest = TRUE, 
              importance = TRUE)
  path_OOB <- paste0("results/ML/", title, "/OOB_error_", title, ".pdf")
  dev.copy(x11)
  dev.print(pdf, path_OOB)
}

#' Produce heatmap for random forest
#'
#' @param subset subset of n genes defining the model
#' @param group 2 column dataframe with name of sample and condition
#' @param sample_name name of the sample
#' @param nb_features number of features selected to define the model
#' @param control control condition for color (HC or VP or EC)
#' @export
my_heatmap <- function(subset, group, sample_name, nb_features, control){
    path_heatmap <- paste("results/figures/", sample_name, "_", nb_features, "_heatmap.pdf", sep ="")
    pdf(path_heatmap, width=6, height=4, paper='special')
    hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256)
    colnames(subset) <- group # This will label the heatmap columns
    csc <- rep(hmcol[50],6)
    csc[group == control] <- hmcol[200]   # column side color will be purple for T and orange for B
    heatmap(as.matrix(subset),scale="row", col=hmcol,ColSideColors=csc)
    dev.off()
}

#' Benchmarking randomForest
#'
#' @param subset_1 dataframe with proteomics expression 1
#' @param subset_2 dataframe with proteomics expression 2
#' @param subset_3 dataframe with proteomics expression 3
#' @param name_analysis name of analysis
#' @export
benchmarking <- function(subset_1, subset_2, subset_3, name_analysis, condition_1, condition_2){
    a <- subset_random_forest(subset_1, condition_1, condition_2)
    b <- subset_random_forest(subset_2, condition_1, condition_2)
    c <- subset_random_forest(subset_3, condition_1, condition_2)
    path_ROC <- paste("results/ML/RF_CD4/", "ROC_", name_analysis, ".pdf", sep ="")
    pdf(path_ROC,width=8,height=6,paper='special')
    par(pty="s")
    plot(a, col = "#8e44ad")
    lines(b, col = "#2471a3")
    lines(c, col =  "#cb4335")
    legend("bottomright", 
           legend = c("100 features", "50 features", "30 features"), 
           col = c("#8e44ad", "#2471a3", "#cb4335"), 
           pch = 20, 
           bty = "n", 
           pt.cex = 2, 
           cex = 1.2, 
           text.col = "black", 
           horiz = F , 
           inset = c(0.1, 0.1))
    dev.off()
    print(auc(a))
    print(auc(b))
    print(auc(c))
}




#' RandomForest for benchmarking
#'
#' @param subset dataframe with proteomics expression (subset)
#' @param condition_1 conditions to compare
#' @param condition_2 conditions to compare
#' @return rf.roc
#' @export
subset_random_forest <- function(subset, condition_1, condition_2){
    rf_frame <- data.frame(t(as.matrix(subset, rownames = TRUE)),
                           group = c(rep(condition_1,6),rep(condition_2,6)), check.names=TRUE)
    rf_frame$group <-as.factor(as.character(rf_frame$group))
    set.seed(71)
    rf<-randomForest(group ~ ., data= rf_frame,ntree=10)
    rf.roc<-roc(rf_frame$group, rf$votes[,2])
    return(rf.roc)
}

#' Produce confusion matrix for random forest object
#'
#' @import ggplot2
#' @import pca3d
#' @import gplots
#' @import scales
#'
#' @param m random forest object
#' @return p
#' @export
ggplotConfusionMatrix <- function(m){
  mytitle <- paste("Accuracy", percent_format()(m$overall[1]),
                   "Kappa", percent_format()(m$overall[2]))
  p <-
    ggplot(data = as.data.frame(m$table) ,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue", limits=c(0, 10)) +
    geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
    theme(legend.position = "none") +
    ggtitle(mytitle)
  return(p)
}
