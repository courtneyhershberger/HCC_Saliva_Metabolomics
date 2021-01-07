

#The chosen metabolites
best_met <- oob[118,"metabolites"]
best_met <- unlist(strsplit(best_met, ","))

df <- proc.data.log.scaled
df <- df[,best_met]
df$mrn <- rownames(df)
df <- merge(df, all.meta[,c("mrn","diagnosis")], by="mrn")
df$mrn <- NULL
df <- df[which(df$diagnosis %in% c("Healthy","Cirrhosis","HCC")),]
df$diagnosis <- factor(df$diagnosis, levels = c("Healthy", "Cirrhosis", "HCC"))


######################################LOOCV
#makes a table of predictions and truths for leave one out method
rf_loop <- function(x, df){
  #make model with LOO
  df2 <- df[-x,]
  rf_classifier = randomForest(diagnosis ~ ., data=df2, ntree=myntree, importance=TRUE)
  
  #use model to predict missing
  rfpred <- predict(rf_classifier, df[x,])
  myrow <- data.frame(run=x, prediction=rfpred, truth=df[x,"diagnosis"],stringsAsFactors = F)
}

pred.truth <- do.call(rbind,lapply(1:nrow(df), rf_loop, df=df))

#Extracts accuracy and other numbers from a table of predictions and truths
truthtable <- function(group, pred.truth){
  TP <- nrow(pred.truth[which(pred.truth$prediction %in% group & pred.truth$truth %in% group),])
  TN <- nrow(pred.truth[which(!pred.truth$prediction %in% group & !pred.truth$truth %in% group),])
  FP <- nrow(pred.truth[which(pred.truth$prediction %in% group & !pred.truth$truth %in% group),])
  FN <- nrow(pred.truth[which(!pred.truth$prediction %in% group & pred.truth$truth %in% group),])
  total <- nrow(pred.truth)
  
  myrow <- data.frame(
    accuracy= (TP+TN)/total,
    balanced_accuracy=((TP/(TP+FN))+(TN/(FP+TN)))/2,
    misclassification = (FP+FN)/total,
    sensitivity= TP/(TP+FN),
    specificity =TN/(TN+FP),
    PPV = TP/total,
    NPV = TN/total,
    stringsAsFactors = F
  )
  rownames(myrow) <- group
  myrow <- data.frame(t(myrow),stringsAsFactors = F)
  myrow <- apply(myrow, 2, function(x) round(x*100,digits = 2))
  return(myrow)
}

rf8_accuracy_metrics <- do.call(cbind, lapply(c("Healthy","HCC","Cirrhosis"), truthtable, pred.truth=pred.truth))




