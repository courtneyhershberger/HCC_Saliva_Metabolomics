#compare models

formattable <- function(table, type){
  table[,"Disease"] <- rownames(table)
  table[,"model"] <- type
  return(table)
  
}
df1 <- formattable(table = decision_tree_metrics, type = "Decision Tree")
df2 <- formattable(table = random_forest_metrics, type = "iRF125")
df3 <- formattable(table = final_rf_accuracy_metrics, type = "iRF4")
df4 <- formattable(table = rf12_accuracy_metrics, type = "irf12")

df_wide <- do.call("rbind", list(df1,df2,df3,df4))
df_long <- gather(df_wide, metric_type, Percent, accuracy:NPV, factor_key=TRUE)
df_long$Model<- factor(df_long$model, ordered = T, levels=c("iRF125","irf12","iRF4","Decision Tree"))
df_long$Disease <- factor(df_long$Disease, ordered = T, levels=c("Healthy","Cirrhosis","HCC"))
df_long$colors <- paste(df_long$Disease, df_long$model, sep="_")
df_long$colors <- factor(df_long$colors, ordered=T, levels=c("Healthy_iRF125","Healthy_irf12","Healthy_iRF4","Healthy_Decision Tree",
                                                             "Cirrhosis_iRF125","Cirrhosis_irf12","Cirrhosis_iRF4","Cirrhosis_Decision Tree",
                                                             "HCC_iRF125","HCC_irf12","HCC_iRF4","HCC_Decision Tree"))


makeplot <- function(mymetric){
  mytitle <- sub("_"," ",as.character(mymetric))
  mytitle <- capitalize(mytitle)
  p <- ggplot(df_long[which(df_long$metric_type %in% mymetric),], aes(Disease, Percent, fill=colors, color=Disease)) +
    geom_col(position=position_dodge(width=0.8), width=.6) +
    #facet_grid(metric_type~., scales="free") +
    theme_classic() +
    my_scale_fill("gradient") +
    my_scale_color("darkest") + 
    #ggtitle(mytitle) +
    labs(x="",y=mytitle) +
    theme(legend.position = "None") +
    ylim(0,100)
  print(p)
  return(p)
}

accuracy_metrics_plotlist <- lapply(c("sensitivity","specificity","balanced_accuracy","misclassification"), makeplot)





