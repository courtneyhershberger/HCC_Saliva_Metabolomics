#The chosen metabolites
best_met <- oob[114,"metabolites"]
best_met <- unlist(strsplit(best_met, ","))
best_met_form <- paste(best_met, sep=' + ', collapse=" + ")

df <- proc.data.log.scaled
df <- df[,best_met]
df$mrn <- rownames(df)
df <- merge(df, all.meta[,c("mrn","diagnosis")], by="mrn")
rownames(df) <- df$mrn
df$mrn <- NULL
df <- df[which(df$diagnosis %in% c("Healthy","Cirrhosis","HCC")),]
df$diagnosis <- factor(df$diagnosis, levels = c("Healthy", "Cirrhosis", "HCC"))


#grow tree

fit <- rpart(as.formula(paste("diagnosis ~", best_met_form)), 
             method="class",
             data=df)

#display results
printcp(fit)
plotcp(fit)
x <- summary(fit)

rules <- rpart.plot::rpart.rules(fit)


#plot tree
plot(fit, uniform=TRUE, main="Classification Tree for Diagnosis of Cirrhosis or HCC")
text(fit, use.n=TRUE, all=TRUE, cex=.8)
print(fancyRpartPlot(fit, caption = NULL))

#Prediction accuracy
df$where <- fit$where
leafnodes <- data.frame(with(df, table(diagnosis, where)))
colnames(leafnodes) <- c("disease_status","leaf","Freq")
prediction <- data.frame(leafnodes %>% group_by(leaf) %>% top_n(n=1), stringsAsFactors = F)
colnames(prediction) <- c("prediction","leaf","Freq")
leafnodes <- merge(leafnodes, prediction[,c("prediction","leaf")])
totals <- leafnodes %>% group_by(leaf) %>% summarise(total=sum(Freq))
leafnodes <- merge(leafnodes, totals)
leafnodes$prop <- leafnodes$Freq/leafnodes$total

#leafnodes <- aggregate(leafnodes$Freq, by=list(disease_status=leafnodes$disease_status, prediction=leafnodes$prediction), FUN=sum, na.rm=T)

 makebarplots <- function(leaf, leafnodes){
   leafnodes <- leafnodes[which(leafnodes$leaf %in% leaf),]
    leafbarplot <- ggplot(leafnodes, aes(fill=disease_status, y=prop, x=disease_status)) + 
    geom_col() +
    #facet_grid(.~leaf) +
    my_scale_fill(palette="main3") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          strip.text.x = element_blank())+
    labs(x="",y="") +
    ylim(0,1)
    print(leafbarplot)
    }

leafbarplots <- lapply(unique(leafnodes$leaf), makebarplots, leafnodes=leafnodes)

#table
spread(leafnodes, prediction, Freq)

######################################LOOCV
#makes a table of predictions and truths for leave one out method
rf_loop <- function(x, df){
  #make model with LOO
  df2 <- df[-x,]
  fit <- rpart(as.formula(paste0("diagnosis ~", best_met_form)), 
               method="class",
               data=df)
  #use model to predict missing
  rfpred <- predict(fit, df[x,], type = "class")
  myrow <- data.frame(run=x, prediction=rfpred, truth=df[x,"diagnosis"],stringsAsFactors = F)
}

pred.truth <- do.call(rbind,lapply(1:110, rf_loop, df=df))

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

summary <- do.call(cbind, lapply(c("Healthy","HCC","Cirrhosis"), truthtable, pred.truth=pred.truth))




