#update: decided to let mtry values float. 

#Random forest on metabolites
set.seed(42)

#prepare dataframe
#df <- proc.data.log.scaled #(this part is in workflow)
df$mrn <- rownames(df)
df <- merge(all.meta[,c("mrn","diagnosis")], df, by="mrn")
df$mrn <- NULL
df <- df[which(df$diagnosis %in% c("Healthy","Cirrhosis","HCC")),]
df$diagnosis <- factor(df$diagnosis, levels = c("Healthy", "Cirrhosis", "HCC"))

#select multiple parameters
ntree <- c(50,100,150,200)

#my functions
#makes a table of predictions and truths for leave one out method
rf_loop <- function(x, df, my_ntree){
  #make model with LOO
  df2 <- df[-x,]
  rf_classifier = randomForest(diagnosis ~ ., data=df2, ntree=my_ntree, importance=TRUE)
  
  #use model to predict missing
  rfpred <- predict(rf_classifier, df[x,])
  myrow <- data.frame(run=x, prediction=rfpred, truth=df[x,"diagnosis"],stringsAsFactors = F)
}


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

#Test parameters ntree
testparams <- function(my_ntree, df){

  pred.truth <- do.call(rbind,lapply(1:111, rf_loop, df=df, my_ntree=my_ntree))
  
  summary <- do.call(cbind, lapply(c("Healthy","HCC","Cirrhosis"), truthtable, pred.truth=pred.truth))
  #print(str(summary))
  return(summary)
}


summary_list <- lapply(ntree, testparams, df=df)
#now there is a list of summary tables.

#Function for making figure to show parameter selection
maketable <- function(summary, summary_list){
  mytable <- do.call(rbind, lapply(summary_list, function(x) x[summary,]))
  mytable <- data.frame(Healthy=mytable[,"Healthy"], HCC=mytable[,"HCC"], Cirrhosis=mytable[,"Cirrhosis"], Average = round(rowMeans(mytable[,c('Healthy', 'HCC','Cirrhosis')], na.rm=TRUE),digits = 2))
  mytable <- cbind(mytable, ntree)
  rownames(mytable) <- NULL
  mytable$metric_type <- capitalize(summary)
  mytable$metric_type <- sub("_"," ",mytable$metric_type)

  
  return(mytable)
}

mytable <- do.call("rbind",lapply(rownames(summary_list[[1]]), maketable, summary_list=summary_list))

mytable_long <- gather(mytable, diagnosis, metric, Healthy,Cirrhosis,HCC, factor_key=TRUE)
mytable_long$best <- ifelse(mytable_long$ntree==150, "Yes","No")

accuracy_metrics_plot <- ggplot(data=mytable_long, aes(x=ntree, y=metric, fill=diagnosis, color=diagnosis)) +
  geom_line(aes(color=diagnosis)) +
  geom_point(aes(color=diagnosis, fill=diagnosis), pch=21) +
  geom_point(data=mytable_long[which(mytable_long$best %in% "Yes"),], aes(x=ntree, y=metric, fill=diagnosis), color="black", pch=22, size=1.7) +
  my_scale_fill(palette = "main3")+
  my_scale_color(palette = "main3")+
  facet_grid(metric_type ~ ., scales="free") +
  xlab("Number of Trees") +
  ylab("") +
  theme_bw() +
  theme(legend.position="bottom") +
  guides(color=guide_legend("Disease State"), fill = FALSE)

print(accuracy_metrics_plot)

