#Random forest on metabolites
set.seed(42)

#prepare dataframe
df <- proc.data.log.scaled #(this part is in workflow)
df$mrn <- rownames(df)
df <- merge(all.meta[,c("mrn","diagnosis")], df, by="mrn")
df$mrn <- NULL
df <- df[which(df$diagnosis %in% c("Healthy","Cirrhosis","HCC")),]
df$diagnosis <- factor(df$diagnosis, levels = c("Healthy", "Cirrhosis", "HCC"))


oob <- data.frame(run=as.numeric(),
                  Mean=as.numeric(),
                 Min=as.numeric(), 
                 Max=as.numeric(),
                 metabolites <- as.character(),
                 stringsAsFactors=FALSE) 

#collect gini scores and out of bag error from each leave one out loop
rf_loop <- function(x, df){
  #make model with LOO
  df2 <- df[-x,]
  rf_classifier = randomForest(diagnosis ~ ., data=df2, ntree=myntree, importance=TRUE)
  oob <- rf_classifier$err.rate[,1][[99]]
  oob_df <- data.frame(MeanDecreaseGini=oob)
  rownames(oob_df) <- "oob"
  
  #get gini
  gini <- rf_classifier$importance[,"MeanDecreaseGini",drop=F]
  gini <- rbind(gini, oob_df)
  colnames(gini) <- x
  return(gini)
}

for(i in 1:125){
  gini.scores <- do.call(cbind,lapply(1:nrow(df), rf_loop, df=df))
  gini.scores <- as.matrix(gini.scores)
  gini.sum <- data.frame(Metabolite=rownames(gini.scores), Mean=rowMeans(gini.scores),Min=rowMins(gini.scores), Max=rowMaxs(gini.scores))
  
  #get oob table
  myrun <- data.frame(gini.sum["oob",c("Mean","Min","Max")])
  myrun$run <- i
  myrun$metabolites <- paste(gini.sum$Metabolite[!gini.sum$Metabolite %in% "oob"], collapse=",")
  oob <- rbind(oob, myrun)

  #drop metabolite with lowest gini score
  gini.sum <- gini.sum[!(row.names(gini.sum) %in% "oob"), ]
  gini.sum <- gini.sum[order(gini.sum$Mean),]
  lowest <- gini.sum[1,"Metabolite"]
  df <- df[ , -which(names(df) %in% lowest)]
  
}

#plot the change in out of bag error as we remove one metabolite at a time
t2.rect1 <- data.frame(xmin=112.5, xmax=125.5, ymin=-Inf, ymax=Inf)
oob$metabolite_num <- 126-oob$run

drop_metabolites <- ggplot(data=oob,
           aes(x = metabolite_num,y = Mean, ymin = Min, ymax = Max))+
  scale_x_reverse() +
  xlab('Number of Selected Metabolites')+ ylab("LOOCV OOB Error")+
  #geom_rect(data=t2.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="yellow", alpha=0.4, inherit.aes = FALSE) +
  geom_pointrange(size=.07) +
  geom_pointrange(data=oob[oob$metabolite_num==4,], aes(x=metabolite_num,y=Mean,ymin=Min,ymax=Max),color=my_colors["dark red"],size=.07) +
  geom_pointrange(data=oob[oob$metabolite_num==12,], aes(x=metabolite_num,y=Mean,ymin=Min,ymax=Max),color=my_colors["dark teal"],size=.07) +
  geom_pointrange(data=oob[oob$metabolite_num==125,], aes(x=metabolite_num,y=Mean,ymin=Min,ymax=Max),color=my_colors["dark yellow"],size=.07) +
  theme_classic() +
  theme(axis.text.x=element_text(size=10, vjust=0.5, color = 'black'),  # x-axis labels
        axis.text.y=element_text(size=10, vjust=0.5, color = 'black'),  # y-axis labels
        axis.title.x=element_text(size=11, vjust=0.1),                # x-title justification  
        axis.title.y=element_text(size=11, vjust=1.5),                 # y-title justification
        legend.title=element_blank(),
        legend.position= "bottom",       
        legend.text = element_text(size = 12)) 
print(drop_metabolites)


#Top metabolites
irf4 <- oob[122,"metabolites"]
irf4  <- unlist(strsplit(irf4 , ","))
irf12 <- oob[114,"metabolites"]
irf12  <- unlist(strsplit(irf12 , ","))

##### New gini plots with additional data
LOOCV.gini.sum$irf4 <- ifelse(LOOCV.gini.sum$Metabolite %in% irf4, "Yes","No")
LOOCV.gini.sum$irf12 <- ifelse(LOOCV.gini.sum$Metabolite %in% irf12, "Yes","No")

#format titles
LOOCV.gini.sum$Metabolite <- replace_ids(LOOCV.gini.sum$Metabolite, met_key[,c("ID","metabolite_name")])
LOOCV.gini.sum$Metabolite <- capitalize(LOOCV.gini.sum$Metabolite)
LOOCV.gini.sum$Metabolite <- factor(LOOCV.gini.sum$Metabolite, ordered = T, levels=LOOCV.gini.sum[order(LOOCV.gini.sum$Mean),"Metabolite"])
LOOCV.gini.sum <- LOOCV.gini.sum[order(LOOCV.gini.sum$Metabolite),]
t2.rect2 <- data.frame(xmin=min(which(LOOCV.gini.sum$best %in% "Yes"))-1, xmax=length(LOOCV.gini.sum$Metabolite)+.5, ymin=-Inf, ymax=Inf)


giniplot_1 = ggplot(data=LOOCV.gini.sum,
                    aes(x = Metabolite,y = Mean, ymin = Min, ymax = Max)) +
  geom_pointrange(size=.3) +
  geom_pointrange(data=LOOCV.gini.sum[which(LOOCV.gini.sum$irf12=="Yes"),], aes(x=Metabolite,y=Mean,ymin=Min,ymax=Max),color=my_colors["dark teal"],size=.3) +
  geom_pointrange(data=LOOCV.gini.sum[which(LOOCV.gini.sum$irf4=="Yes"),], aes(x=Metabolite,y=Mean,ymin=Min,ymax=Max),color=my_colors["dark red"],size=.3) +
  xlab('All Metabolites (iRF125)')+ ylab("Range of Gini Scores")+
  coord_flip() +
  theme_classic() +
  theme(axis.text.x=element_text(size=12, vjust=0.5, color = 'black'),  # x-axis labels
        axis.text.y=element_blank(),  # y-axis labels
        axis.title.x=element_text(size=13, vjust=0.1),                # x-title justification  
        axis.title.y=element_text(size=13, vjust=1.5),                 # y-title justification,
        legend.position= "none") 
          
print(giniplot_1)

        
 ###########Plot2
giniplot_2 = ggplot(data=LOOCV.gini.sum[which(LOOCV.gini.sum$irf12 %in% "Yes"),],
                    aes(x = Metabolite, y = Mean, ymin = Min, ymax = Max ))+
  geom_pointrange(color=my_colors["dark teal"], size=.3) +
  xlab('Selected Metabolites (irf12)')+ ylab("Range of Gini Scores")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip() +
  theme_classic()

print(giniplot_2)

############Plot3
   
giniplot_3 = ggplot(data=LOOCV.gini.sum[which(LOOCV.gini.sum$irf4 %in% "Yes"),],
                    aes(x = Metabolite, y = Mean, ymin = Min, ymax = Max ))+
  geom_pointrange(color=my_colors["dark red"], size=.3) +
  xlab('iRF4')+ ylab("Range of Gini Scores")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip() +
  theme_classic()

print(giniplot_3)        
