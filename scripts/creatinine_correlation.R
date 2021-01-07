
creatinine <- met_key[which(met_key$metabolite_name=="creatinine"),"ID"]

clinical_creat <- clinical_features[,c("MRN","Creatinine (mg/dL)")]
rownames(proc.data.log.scaled) 

df <- proc.data.log.scaled[,creatinine,drop=F]
df$label <- rownames(df)
df <- merge(df, all.meta)
df <- df[,c(creatinine,"mrn")]

df <- merge(df, clinical_creat, by.x="mrn", by.y="MRN")
colnames(df) <- c("mrn","metabolite_creatinine","clinical_creatinine")

print(ggplot(df, aes(x=metabolite_creatinine, y=clinical_creatinine)) +
  geom_point() +
  theme_classic())
