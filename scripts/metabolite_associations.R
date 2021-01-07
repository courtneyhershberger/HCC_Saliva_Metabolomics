#Associations for each metabolite
metabolites <- colnames(df)

df$mrn <- rownames(df)
df <- merge(all.meta, df, by="mrn")
df <- df[which(df$diagnosis %in% c("Healthy","Cirrhosis","HCC")),]
df$diagnosis <- factor(df$diagnosis, levels = c("Healthy", "Cirrhosis", "HCC"))


#Logistic Regression
glmfit <- function(metabolites, control, condition, covar, df){
  df <- df[which(df$diagnosis %in% c(control,condition)),]
  df$outcome <- ifelse(df$diagnosis %in% control, 0, 1)
  
  fit <- glm(df[,"outcome"]~df[,metabolites]+df[,"age_days"]+df[,"sex"],data=df, family = "binomial")
  model <- summary(fit)
  p.value <- model$coefficients[2,4]
  beta <- round(model$coefficients[2,1], digits=3)
  std.err <- round(model$coefficients[2,2], digits=3)
  
  #get fold change of median
  group1_med <- median(exp(df[which(df$diagnosis %in% control),metabolites]))
  group2_med <- median(exp(df[which(df$diagnosis %in% condition),metabolites]))
  FC <- group2_med/group1_med
  log2FC <- log2(FC)
  
  myrow <- data.frame(Variable=metabolites, Group1=control, Group2=condition, Pvalue=p.value, Beta=beta, StdErr=std.err, group1_med=group1_med, group2_med=group2_med, FC=FC, log2FC=log2FC)
  return(myrow)
}

healthy.v.cirr <- do.call(rbind, lapply(metabolites, glmfit, control="Healthy", condition="Cirrhosis", df=df))
healthy.v.hcc <- do.call(rbind, lapply(metabolites, glmfit, control="Healthy", condition="HCC", df=df))
cirr.v.hcc <- do.call(rbind, lapply(metabolites, glmfit, control="Cirrhosis", condition="HCC", df=df))
met_pvalues <- do.call(rbind, list(healthy.v.cirr,healthy.v.hcc,cirr.v.hcc))

#adjust for multiple hypothesis testing, format pvalue and FDR for tables and graphing
met_pvalues$FDR <- p.adjust(met_pvalues$Pvalue, method="fdr") 
met_pvalues$neg_log_10FDR <- -log10(met_pvalues$FDR)
met_pvalues$FDR <- format.pval(met_pvalues$FDR,  eps=.001, digits=3) 
met_pvalues$Pvalue <- format.pval(met_pvalues$Pvalue,  eps=.001, digits=3)

#significant
print(sig_met <- met_pvalues[which(met_pvalues$neg_log_10FDR > 0.69897),])

#boxplots

makeboxplot <- function(met, met_key){
  print(met)
  met_name <- met_key[which(met_key$ID %in% met),"metabolite_name"]
  met_pvalues_sub <- met_pvalues[which(met_pvalues$Variable %in% met),]
  names(met_pvalues_sub)[names(met_pvalues_sub) == 'Variable'] <- ".y."
  names(met_pvalues_sub)[names(met_pvalues_sub) == 'Group1'] <- "group1"
  names(met_pvalues_sub)[names(met_pvalues_sub) == 'Group2'] <- "group2"
  df <- df[,c(met,"diagnosis")]
  colnames(df) <- c("metabolite","Diagnosis")
  
  ###setting up pvalues
  #myincrease <- (max(df$metabolite)-min(df$metabolite))*.06
  ypositions <- get_y_position(df,metabolite~Diagnosis, step.increase = .3)
  met_pvalues_sub <- merge(met_pvalues_sub, ypositions)
  xcoord <- data.frame(group=levels(df$Diagnosis),xmin=c(1:length(levels(df$Diagnosis))),xmax=c(1:length(levels(df$Diagnosis))))
  met_pvalues_sub <- merge(met_pvalues_sub, xcoord[,c("group","xmin")], by.x="group1", by.y="group")
  met_pvalues_sub <- merge(met_pvalues_sub, xcoord[,c("group","xmax")], by.x="group2", by.y="group")
  met_pvalues_sub <- met_pvalues_sub[,c(".y.","group1","group2","FDR","y.position","groups","xmin","xmax")]
  met_pvalues_sub[which(met_pvalues_sub$FDR>.2),"FDR"] <- "n.s."
  

  p <- ggboxplot(df, x = "Diagnosis", y = "metabolite",
                 color = "Diagnosis", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                 add = "jitter") +
    stat_pvalue_manual(met_pvalues_sub, label="FDR", tip.length= 0.01) +
    ggtitle(capitalize(met_name)) +
    ylab("Relative Abundance") +
    theme(legend.position = "none") +
    theme(axis.title.x = element_blank(), plot.title = element_text(size=12, face="bold"))
  print(p)
  
}

metabolite_boxplot_list <- lapply(unique(sig_met$Variable), makeboxplot, met_key=met_key)

############Add Volcano plots for each comparison
#met_pvalues (Beta vs FDR)
df <- met_pvalues
df <- merge(df, met_key, by.x="Variable", by.y="ID")
df$comparison <- paste(df$Group1, "vs", df$Group2)
df$FDR
df$significant <- ifelse(df$neg_log_10FDR>0.69897, "Is_significant","Not_significant")

volcanos <- ggplot(data=df, aes(x=log2FC, y=neg_log_10FDR, color=significant)) +
  geom_hline(yintercept = 0.69, linetype="dashed", col = my_colors["dark red"]) +
  geom_vline(xintercept = 0, col = my_colors["dark gray"]) +
  geom_point(size=1) +
  geom_text_repel(aes(label=ifelse(neg_log_10FDR>.69897,as.character(capitalize(metabolite_name)),'')),size=2.7, seed=42, force=2) +
  my_scale_color(palette="two") +
  ylab(expression(Log[10](FDR))) +
  facet_grid(.~comparison) +
  theme_bw() +
  theme(legend.position = "none")
  
print(volcanos)
