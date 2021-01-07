

#Table 1: Summary statistics
maketable <- function(type){
  df <- all.meta[all.meta$diagnosis %in% type,]
  
  #total
  total <- nrow(df)
  
  #ages
  ages <- df$age_days
  mean_age <- round(mean(ages)/365, digits = 1)
  min_age <- round(min(ages)/365)
  max_age <- round(max(ages)/365)
  ages <- data.frame(age = paste0(mean_age, " (", min_age, "-", max_age,")"))
  
  #Male
  male <- length(which(df$sex %in% "Male"))
  male_percent <- round((male/total)*100)
  male <- paste0(male, " (", male_percent, "%)")
  
  #Female
  female <- length(which(df$sex %in% "Female"))
  female_percent <- round((female/total)*100)
  female <- paste0(female, " (", female_percent, "%)")
  
  total <- as.character(total)
  
  mycolumn = c(total, ages, "", male, female)
  
}
table1 <- data.frame(do.call("cbind", lapply(c("Healthy","Cirrhosis","HCC"), maketable)))
colnames(table1) <- c("Healthy","Cirrhosis","HCC")
table1$Characteristic <- c("Total (n)", "Mean age (min-max)", "Sex", "   Male (%)", "   Female (%)")
table1 <- table1[, c("Characteristic", "Healthy", "Cirrhosis", "HCC")]

#Table 1 Extended
clinical_features <- clinical_features[which(clinical_features$MRN %in% all.meta$mrn),]
totals <- clinical_features %>% dplyr::count(Diagnosis)

binarycounts <- function(x){
  df <- aggregate(clinical_features[,x], by=list(Category=clinical_features$Diagnosis), FUN=sum)
  df <- merge(df, totals, by.x="Category",by.y="Diagnosis")
  df$percent <- round((df[,x]/df[,"n"])*100)
  df$final <- paste0(df[,x], "(",df[,"percent"],"%)")
  myrow <- data.frame(Characteristic=x, Healthy=df[which(df$Category=="hernia clinic"),"final"], Cirrhosis=df[which(df$Category=="Cirrhosis"),"final"],HCC=df[which(df$Category=="HCC"),"final"])
  return(myrow)
}
table1.2 <- do.call("rbind",lapply(c("Diabetes mellitus type 2","Hypertension","Coronary artery disease","Hyperlipidemia","Psychiatric disorder","COPD/Asthma/OSA","Other cancer history","Thyroid","Other PMH","Ascites","Encephalopathy"), binarycounts))

mean_sem <- function(x){
  print(x)
  st.err <- function(x) {
    sd(x)/sqrt(length(x))
  }
  df <- clinical_features[,c("Diagnosis",x)]
  df <- df[complete.cases(df),]
  df1 <- aggregate(df[,x], by=list(Category=df$Diagnosis), FUN=mean)
  df1$sterr <- aggregate(df[,x], by=list(Category=df$Diagnosis), FUN=st.err)[,x]
  df1$final <- paste0(round(df1[,x],digits = 1)," (",round(df1[,"sterr"],digits = 1),")")
  myrow <- data.frame(Characteristic=paste0("Mean ",x," (std.err)"), Healthy=df1[which(df1$Category=="hernia clinic"),"final"], Cirrhosis=df1[which(df1$Category=="Cirrhosis"),"final"],HCC=df1[which(df1$Category=="HCC"),"final"])
  return(myrow)
}
table1.3 <- do.call("rbind",lapply(c("Hemoglobin (g/dl)","Platelets (k/uL)","AST (U/L)","ALT (U/L)","ALP (U/L)","Bilirubin, Total (mg/dL)","Albumin (g/dL)","PT-INR","Glucose (mg/dL)","Creatinine (mg/dL)"), mean_sem))

table1 <- do.call("rbind", list(table1, table1.2, table1.3))

#Table 2: Most significant metabolites
assoc_tables <- function(table, variable, variable_name){
  table <- table[order(table$Variable),]
  table$Beta_SE <- paste0(table[,"Beta"], " (",table[,"StdErr"], ")")
  table <- table[,c("Group1", "Group2","Variable","Beta_SE","Pvalue","FDR")]
  colnames(table) <- c("Group1 (Reference)", "Group2", variable_name, "Coefficient (SE)", "P Value", "FDR P")
  if(variable=="n"){
    table <- table[,c("Group1 (Reference)", "Group2", "Coefficient (SE)", "P Value", "FDR P")]
  }
  return(table)
}
table2 <- assoc_tables(sig_met, variable="y", variable_name="Metabolite")
table2$Metabolite <- capitalize(replace_ids(table2$Metabolite, met_key[,c("ID","metabolite_name")]))
table2 <- table2[order(table2$`Group1 (Reference)`, table2$Group2),]
supptable1 <- assoc_tables(age_pvalues, variable="n", variable_name="Age")
supptable2 <- assoc_tables(sex_pvalues, variable="n", variable_name="Sex")
supptable3 <- assoc_tables(smoker_pvalues, variable="n", variable_name="Smoking_Status")

#Table 3: Accuracy Metrics
formattable <- function(table, type){
  table <- round(table, digits=1)
  table <- rbind(table, round(colMeans(table),digits=1))
  table[,"Disease"] <- rownames(table)
  table[,"Model"] <- type
  table[4,"Disease"] <- "Average"
  return(table)
  
}
df4 <- formattable(table = decision_tree_metrics, type = "Decision Tree")
df1 <- formattable(table = random_forest_metrics, type = "iRF125")
df3 <- formattable(table = final_rf_accuracy_metrics, type = "iRF4")
df2 <- formattable(table = rf12_accuracy_metrics, type = "iRF12")

table3 <- do.call("rbind", list(df1,df2,df3,df4))
table3 <- table3[,c("Disease","Model","sensitivity","specificity","balanced_accuracy","misclassification","PPV","NPV")]
colnames(table3) <- capitalize(colnames(table3))
colnames(table3) <- sub("_"," ",colnames(table3))

#Supp Table 3: All metabolites
supptable4 <- met_pvalues
supptable4 <- merge(supptable4, met_key, by.x="Variable","ID")
supptable4 <- supptable4[order(supptable4$Variable),]
supptable4$metabolite_name <- capitalize(supptable4$metabolite_name)
supptable4$Comparison <- paste(supptable4$Group1, "vs", supptable4$Group2)
supptable4 <- supptable4[,c("Comparison","metabolite_name", "Beta","StdErr","Pvalue","FDR")]
supptable4$Beta_SE <- paste0(supptable4$Beta," (",supptable4$StdErr,")")
supptable4 <- supptable4[,c("Comparison","metabolite_name","Beta_SE","Pvalue","FDR")]
colnames(supptable4) <- c("Comparison","Metabolite","Coefficient (SE)", "P Value", "FDR P")
supptable4.a <- supptable4[which(supptable4$Comparison %in% "Healthy vs Cirrhosis"),-1]
supptable4.c <- supptable4[which(supptable4$Comparison %in% "Healthy vs HCC"),-1]
supptable4.b <- supptable4[which(supptable4$Comparison %in% "Cirrhosis vs HCC"),-1]

################################################Main Tables
rtffile <- RTF("tables/manuscript_tables.doc")

#Table 1: Summary Statistics of Cohort

addParagraph(rtffile, "Summary statistics for study cohort\n")
addTable(rtffile, data.frame(table1, row.names = NULL), col.justify="C", header.col.justify="C")
addParagraph(rtffile, "\n")

addParagraph(rtffile, "Significant disease status associations with metabolite abudnance\n")
addTable(rtffile, data.frame(table2, row.names = NULL), col.justify="C", header.col.justify="C")
addParagraph(rtffile, "\n")

addParagraph(rtffile, "Accuracy metrics for predicting disease status\n")
addTable(rtffile, data.frame(table3, row.names = NULL), col.justify="C", header.col.justify="C")
addParagraph(rtffile, "\n")

done(rtffile)

################################################Supplementary Tables
rtffile <- RTF("tables/supplementary_tables.doc")

addParagraph(rtffile, "Supplementary Table 1: Disease status associations with age (days)\n")
addTable(rtffile, data.frame(supptable1[,-3], row.names = NULL), col.justify="C", header.col.justify="C")
addParagraph(rtffile, "\n")

addParagraph(rtffile, "Supplementary Table 2: Disease status associations with sex (male)\n")
addTable(rtffile, data.frame(supptable2[,-3], row.names = NULL), col.justify="C", header.col.justify="C")
addParagraph(rtffile, "\n")

addParagraph(rtffile, "Supplementary Table 3: Disease status associations with Smoking Status\n")
addTable(rtffile, data.frame(supptable3[,-3], row.names = NULL), col.justify="C", header.col.justify="C")
addParagraph(rtffile, "\n")

addParagraph(rtffile, "Supplementary Table 4: Associations of metabolites in Healthy (reference) versus Cirrhosis\n")
addTable(rtffile, data.frame(supptable4.a, row.names = NULL), col.justify="C", header.col.justify="C")
addParagraph(rtffile, "\n")

addParagraph(rtffile, "Supplementary Table 5: Associations of metabolites in Cirrhosis (reference) vs HCC\n")
addTable(rtffile, data.frame(supptable4.b, row.names = NULL), col.justify="C", header.col.justify="C")
addParagraph(rtffile, "\n")

addParagraph(rtffile, "Supplementary Table 6: Associations of metabolites in Healthy (reference) vs HCC\n")
addTable(rtffile, data.frame(supptable4.c, row.names = NULL), col.justify="C", header.col.justify="C")
addParagraph(rtffile, "\n")

done(rtffile)
