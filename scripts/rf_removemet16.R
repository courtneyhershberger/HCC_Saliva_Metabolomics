#repeat all models without metabolite_16

#Top 8
best_met <- oob[118,"metabolites"]
best_met <- unlist(strsplit(best_met, ","))
best_met_no16 <- best_met[which(!best_met %in% "metabolite_16")]
best_met_form <- paste(best_met_no16, sep=' + ', collapse=" + ")

df <- proc.data.log.scaled
#remove metabolite 16
df <- df[,!colnames(df) %in% "metabolite_16"]
df$mrn <- rownames(df)
df <- merge(df, all.meta[,c("mrn","diagnosis")], by="mrn")
df$mrn <- NULL
df <- df[which(df$diagnosis %in% c("Healthy","Cirrhosis","HCC")),]
df$diagnosis <- factor(df$diagnosis, levels = c("Healthy", "Cirrhosis", "HCC"))


######################################Random Forest
#makes a table of predictions and truths for leave one out method
rf_loop <- function(x, df, variables){
  df <- df[,c(variables,"diagnosis")]
  #make model with LOO
  df2 <- df[-x,]
  rf_classifier = randomForest(diagnosis ~ ., data=df2, ntree=100, importance=TRUE)
  
  #use model to predict missing
  rfpred <- predict(rf_classifier, df[x,])
  myrow <- data.frame(run=x, prediction=rfpred, truth=df[x,"diagnosis"],stringsAsFactors = F)
}

pred.truth.7 <- do.call(rbind,lapply(1:110, rf_loop, df=df, variables=best_met_no16))
pred.truth.124 <- do.call(rbind,lapply(1:110, rf_loop, df=df, variables=colnames(df)))

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

rf_7_accuracy <- data.frame(do.call(cbind, lapply(c("Healthy","HCC","Cirrhosis"), truthtable, pred.truth=pred.truth.7)))
rf_124_accuracy <- data.frame(do.call(cbind, lapply(c("Healthy","HCC","Cirrhosis"), truthtable, pred.truth=pred.truth.124)))

##########Decision Tree
dt_loop <- function(x, df, my_ntree, my_mtry){
  #make model with LOO
  df2 <- df[-x,]
  fit <- rpart(as.formula(paste0("diagnosis ~", best_met_form)), 
               method="class",
               data=df)
  #use model to predict missing
  rfpred <- predict(fit, df[x,], type = "class")
  myrow <- data.frame(run=x, prediction=rfpred, truth=df[x,"diagnosis"],stringsAsFactors = F)
}

pred.truth <- do.call(rbind,lapply(1:110, dt_loop, df=df))

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

dt_accuracy <- data.frame(do.call(cbind, lapply(c("Healthy","HCC","Cirrhosis"), truthtable, pred.truth=pred.truth)))

rf_124_accuracy$Model <- "iRF125"
rf_124_accuracy$metric <- rownames(rf_124_accuracy)
rf_7_accuracy$Model <- "iRF8"
rf_7_accuracy$metric <- rownames(rf_7_accuracy)
dt_accuracy$Model <- "DT"
dt_accuracy$metric <- rownames(dt_accuracy)

no_16 <- rbind(rf_124_accuracy, rf_7_accuracy, dt_accuracy)
no_16 <- gather(no_16, Disease, Percent, Healthy:Cirrhosis, factor_key=TRUE)
no_16$metric <- factor(no_16$metric, ordered = T, levels=c("sensitivity","specificity","accuracy","balanced_accuracy","misclassification","PPV","NPV"))
no_16$metric <- capitalize(sub("_"," ",no_16$metric))
no_16$Model <- factor(no_16$Model, ordered = T, level=c("iRF125","iRF8","DT"))


no_16_accuracy_barplot <- ggplot(no_16, aes(Disease, Percent, fill=Disease, alpha=Model, color=Disease)) +
    geom_col(position=position_dodge(width=0.8), width=.6) +
    facet_grid(metric~.) +
    theme_classic() +
    my_scale_fill("main3") +
    my_scale_color("main3") + 
    ylim(0,100)


