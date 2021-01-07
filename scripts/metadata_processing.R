#joining and combining unblinded files with redcap and EMR data

#fake mrns are removed during the merge because they don't exist in redcap or EMR
#batch 1 gets covariates from combination of redcap and EMR
all.meta.1 <- merge(unblind1, covar1)
all.meta.1 <- all.meta.1[,c("mrn", "WCMC Sample ID", "Gender", "Age", "diagnosis")]
colnames(all.meta.1) <- c("mrn","label","sex","age_days","diagnosis")

#batch 2 gets some covariates from redcap
all.meta.2 <- merge(unblind2, covar2, by.x="mrn",by.y="MRN")
all.meta.2 <- all.meta.2[,c("mrn","WCMC Sample ID", "Gender", "Date  saliva pre-treatment _age_days", "diagnosis")]
colnames(all.meta.2) <- c("mrn","label","sex","age_days","diagnosis")

#add missing values from redcap with information from EMR
all.meta.2.na <- all.meta.2[!complete.cases(all.meta.2),]
all.meta.2.na <- merge(all.meta.2.na, covar1)
all.meta.2.na <- all.meta.2.na[,c("mrn", "label", "Gender", "Age", "diagnosis")]
colnames(all.meta.2.na) <- c("mrn","label","sex","age_days","diagnosis")

all.meta.2 <- all.meta.2[complete.cases(all.meta.2),]

#add smoking data
all.meta <- do.call(rbind, list(all.meta.1, all.meta.2, all.meta.2.na))
all.meta <- merge(all.meta, covar3[,c("MRN","Current smoker at pre treatment")], by.x="mrn", by.y="MRN")
colnames(all.meta) <- c("mrn","label","sex","age_days","diagnosis","smoker")

write.csv(all.meta, "processed_data/saliva_etabolome_GCTOF_MS_meta_processed.csv",quote = F, row.names = F)


#clean up clinical data for table 1
col_names <- colnames(clinical_features)
clinical_features <- clinical_features[which(clinical_features$Diagnosis %in% c("hernia clinic","HCC","Cirrhosis")),]
clinical_features <- data.frame(clinical_features)

#Impute <x with half the value
impute_less <- function(col){
  print(col)
  x <- clinical_features[,col]
  if(any(grepl("<",x))){
    print(str(x))
    mynum <- x[grep("<",x)]
    mynum <- as.numeric(sub("<","",mynum))
    mynum <- mynum/2
    x[grep("<",x)] <- mynum
    x <- as.numeric(x)
    return(x)
  }
  else{
    return(x)
  }
}

newcols <- data.frame(do.call("cbind",lapply(colnames(clinical_features)[21:30], impute_less)))
colnames(newcols) <- colnames(clinical_features)[21:30]
clinical_features <- cbind(clinical_features[1:20], newcols)
clinical_features <- tibble(clinical_features)
colnames(clinical_features) <- col_names
