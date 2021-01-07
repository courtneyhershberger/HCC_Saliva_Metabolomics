#################
## Saliva Metabolome Workflow ##
## PCA ##
#################


#### Perform PCA ####
saverownames <- rownames(pca.data)
pca.data <- data.frame(sapply(pca.data, as.character))
pca.data <- data.frame(sapply(pca.data, as.numeric))
rownames(pca.data) <- saverownames

my.meta <- col.meta.data[c("label","comment","batch")]
my.meta[grep("pool", my.meta$comment),"comment"] <- "pool"
my.meta <- merge(my.meta, all.meta, by="label", all.x = T)
my.meta[is.na(my.meta$diagnosis),"diagnosis"] <- "pool"

#which ones are duplicated?
IDs <- my.meta[duplicated(my.meta$mrn)&!is.na(my.meta$mrn),"mrn"]
my.meta$duplicated <- ifelse(my.meta$mrn %in% IDs, "duplicated","unique")
my.meta$diagnosis <- factor(my.meta$diagnosis, ordered = T, levels=c("Healthy","Cirrhosis","HCC","Mixed HCC-ICC","pool"))
my.meta$batch <- as.factor(my.meta$batch)
myIDs <- data.frame("mrn"=unique(my.meta$mrn), "ID"=paste0("P_",1:length(unique(my.meta$mrn))))
my.meta <- merge(my.meta, myIDs, by="mrn")




pcafunction <- function(batch, title){
  #subset
  mysamples <- my.meta[which(my.meta$batch %in% batch),"label"]
  pca.data <- pca.data[which(row.names(pca.data) %in% mysamples),]
  
  pca <- prcomp(pca.data,scale=T)
  #print("ScreePlot")
  #print(screeplot(pca))
  #print(summary(pca))
  
  pca_df <- as.data.frame(pca$x)
  pca_df <- pca_df[,1:2]
  colnames(pca_df) <- c("PrinComp1","PrinComp2")
  pca_df$label <- rownames(pca_df)
  pca_df <- merge(pca_df, my.meta,by="label")

  p <- ggplot(pca_df, aes(x=PrinComp1, y=PrinComp2,color=diagnosis,shape=batch)) +
  geom_point() +
    ggtitle(title) +
    theme_bw()
  print(p)
  return(pca_df)
  
}

pca_df1 <- pcafunction(batch="1",title="PCA_Experiment_1")
pca_df2 <- pcafunction(batch="2",title="PCA_Experiment_2")
pca_df_both <- pcafunction(batch=c("1","2"),title="PCA_all_samples")

#check the duplicates
print("Within experiment one, there are 3 duplicates of the HCC Diagnosis, they are labeled below")
p <- ggplot(pca_df1, aes(PrinComp1, PrinComp2, color=diagnosis)) +
 geom_point(alpha=0.5) +
  theme_bw()
p <- p +  geom_text(data=subset(pca_df1, duplicated %in% "duplicated" & diagnosis %in% "HCC"),
            aes(PrinComp1,PrinComp2,label=mrn)) +
  ggforce::facet_zoom(x = duplicated == 'duplicated', y = duplicated == 'duplicated') + 
  theme_bw() +
  ggtitle("PCA Experiment 1: Technical Replicates")
print(p)


print("Between the two experimental batches, 6 healthy samples are replicated, they are displayed below")

p <- ggplot(pca_df_both[which(pca_df_both$diagnosis %in% c("HCC","Healthy","Cirrhosis")),], aes(PrinComp1, PrinComp2, color=diagnosis, shape=batch)) +
  geom_point() +
  labs(col="Disease Status",shape="Batch") +
  scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) +
  xlab("PC1") +
  ylab("PC2")+
  theme_bw()
  
p <- p +  geom_label_repel(data=subset(pca_df_both, duplicated %in% "duplicated" & diagnosis %in% "Healthy"),
                    aes(PrinComp1,PrinComp2,label=ID)) +
  theme_bw() +
  ggtitle("PCA")
print(p)

print(p + ggforce::facet_zoom(x = duplicated == 'duplicated', y = duplicated == 'duplicated'))


