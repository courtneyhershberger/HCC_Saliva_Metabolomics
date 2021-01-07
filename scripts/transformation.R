#this script takes the processed relative abundance data and assesses distribution of random 25 metabolites
#The distributions are not normal, the relative abundances are log transformed and assessed a second time

#wide to long format
proc.data.long <- proc.data
proc.data.long$ID <- rownames(proc.data.long)
proc.data.long <- gather(proc.data.long, metabolite, relative_abundance, -ID, factor_key=TRUE, na.rm=T)
proc.data.long <- proc.data.long[complete.cases(proc.data.long),]
proc.data.long$relative_abundance <- as.numeric(proc.data.long$relative_abundance)

#distribution of relative abundance
print("normalized relative abundance")
print(ggplot(data=proc.data.long, aes(x=relative_abundance))+ 
  geom_histogram(bins=10) + 
  facet_wrap(~metabolite, scales="free", ncol=5) + 
  theme_classic())

#log transform & re-plot
print("log transformed normalized relative abundance")
proc.data.long$relative_abundance <- log(proc.data.long$relative_abundance)
print(ggplot(data=proc.data.long, aes(x=relative_abundance))+ 
        geom_histogram(bins=10) + 
        facet_wrap(~metabolite, scales="free",ncol=5) + 
        theme_classic())

#log transform proc.data
proc.data.log <- data.frame(sapply(proc.data, as.numeric))
proc.data.log <- log(proc.data.log)
rownames(proc.data.log) <- rownames(proc.data)
