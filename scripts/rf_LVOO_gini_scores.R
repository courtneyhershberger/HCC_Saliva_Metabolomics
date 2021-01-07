#Random forest on metabolites
set.seed(42)

#prepare dataframe
df <- proc.data.log.scaled #(this part is in workflow)
df$mrn <- rownames(df)
df <- merge(all.meta[,c("mrn","diagnosis")], df, by="mrn")
df$mrn <- NULL
df <- df[which(df$diagnosis %in% c("Healthy","Cirrhosis","HCC")),]
df$diagnosis <- factor(df$diagnosis, levels = c("Healthy", "Cirrhosis", "HCC"))


#collect gini scores from each leave one out loop
rf_loop <- function(x, df){
  #make model with LOO
  df2 <- df[-x,]
  rf_classifier = randomForest(diagnosis ~ ., data=df2, ntree=myntree, importance=TRUE)
  
  #get gini
  gini <- rf_classifier$importance[,"MeanDecreaseGini",drop=F]
  colnames(gini) <- x
  return(gini)
}

gini.scores <- do.call(cbind,lapply(1:nrow(df), rf_loop, df=df))
gini.scores <- as.matrix(gini.scores)

gini.sum <- data.frame(Metabolite=rownames(gini.scores), Mean=rowMeans(gini.scores),Min=rowMins(gini.scores), Max=rowMaxs(gini.scores))
myorder <- gini.sum[order(gini.sum$Mean, decreasing = F),"Metabolite"]
gini.sum$Metabolite <- factor(gini.sum$Metabolite, levels=myorder)
gini.sum <- gini.sum[order(gini.sum$Metabolite, decreasing=T),]

p = ggplot(data=gini.sum,
           aes(x = Metabolite,y = Mean, ymin = Min, ymax = Max ))+
  geom_pointrange() +
  xlab('Metabolite')+ ylab("Range of Gini Scores")+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip() +
  theme(axis.text.x=element_text(size=15, vjust=0.5, color = 'black'),  # x-axis labels
        axis.text.y=element_text(size=15, vjust=0.5, color = 'black'),  # y-axis labels
        axis.title.x=element_text(size=17.5, vjust=0.1),                # x-title justification  
        axis.title.y=element_text(size=17.5, vjust=1.5),                 # y-title justification
        legend.title=element_blank(),
        legend.position= "bottom",       
        legend.text = element_text(size = 14.5)) +
  theme_classic()
  
print(p)
