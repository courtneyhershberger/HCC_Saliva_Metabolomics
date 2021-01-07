#################
## Saliva Metabolome Workflow ##
## 1. Data Preprocessing
#################

file.prefix <- "saliva_metabolome_GCTOF_MS_"
output.dir <- "processed_data/"

#### Collect sample meta data ####
datalist <- data

#pool samples do not have unique identifier
fixpools <- function(data){
  data <- data.frame(data, stringsAsFactors = F)
  data$order <- 1:nrow(data)
  myrow <- grep("label",data[,8])
  poollabel <- as.character(data[myrow,])
  poollabel2 <- make.unique(poollabel,sep="_")
  data[myrow,] <- poollabel2
  data$order <- NULL
  return(data)
}
datalist <- lapply(datalist, fixpools)

getmeta <- function(data, run){
  myrow <- grep("BinBase",data[,1])
  col.meta.data <- as.data.frame(data[1:myrow,8:ncol(data)])
  col.meta.data <- t(col.meta.data)
  colnames(col.meta.data) <- col.meta.data[1,]
  col.meta.data <- col.meta.data[-1,]
  col.meta.data <- data.frame(col.meta.data, stringsAsFactors = F)
  col.meta.data[,"file_id"] <- rownames(col.meta.data)
  col.meta.data <- col.meta.data[,c("mx.class","mx.sample","label","comment","species","organ","treatment","file_id")]  
  col.meta.data$batch <- run
  col.meta.data <- data.frame(lapply(col.meta.data, as.character))
  return(col.meta.data)
}

col.meta.data <- do.call("rbind",mapply(getmeta, data=datalist, run=c("1","2"), SIMPLIFY = F))
col.meta.data <- col.meta.data[which(col.meta.data$label %in% all.meta$label| grepl("pool",col.meta.data$label)),]

write.csv(col.meta.data, file = file.path(output.dir, paste0(file.prefix, "sample_metadata.csv")))
print("beginning of sample_metadata.csv")
print(head(col.meta.data))


#### Collect metabolite meta data ####
getmetarow <- function(data){
  myrow <- grep("BinBase",data[,1]) #identify column headers
  myrow2 <- myrow + 1 #identify first row of data
  row.meta.data <- data[myrow2:dim(data)[1], 1:8]
  colnames(row.meta.data) <- data[myrow,1:8]
  row.meta.data[,"ret.index"] <- as.numeric(row.meta.data[,"ret.index"])
  row.meta.data[,"quant mz"] <- as.numeric(row.meta.data[,"quant mz"])
  row.meta.data[,"BB id"] <- as.numeric(row.meta.data[,"BB id"])
  row.meta.data[,"PubChem"] <- as.numeric(row.meta.data[,"PubChem"])
  row.meta.data[,"mass spec"] <- NULL
  row.meta.data <- row.meta.data[!is.na(row.meta.data$PubChem),]
  return(row.meta.data)
}

row.meta.data <- do.call("rbind",lapply(datalist, getmetarow))
row.meta.data <- unique(row.meta.data) 
row.meta.data$metabolite_name <- row.meta.data$`BinBase name`
row.meta.data$`BinBase name` <- NULL

write.csv(row.meta.data, file = file.path(output.dir, paste0(file.prefix, "metabolite_metadata.csv")))
print("beginning of metabolite_metadata.csv")
print(head(row.meta.data))


#### Process Data ####
getdata <- function(data){
  
  #get important row numbers
  myrow_binbase <- grep("BinBase",data[,1]) #identify column headers
  myrow_data <- myrow_binbase + 1 #identify first row of data
  myrowlabel <- grep("label",data[,8])
  
  proc.data <- data[myrow_data:nrow(data), 9:ncol(data)]
  colnames(proc.data) <- data[myrowlabel,9:ncol(data)]
  rownames(proc.data) <- data[myrow_data:nrow(data),1]

  return(proc.data)
}
proc.data.list <- lapply(datalist, getdata)

#bring two input file together and transpose
proc.data <- merge(proc.data.list[[1]], proc.data.list[[2]], by='row.names')
rownames(proc.data) <- proc.data$Row.names
proc.data$Row.names <- NULL
proc.data <- proc.data[which(rownames(proc.data) %in% row.meta.data$metabolite_name),]
proc.data <- proc.data[complete.cases(proc.data),]

#get new IDS
met_key <- data.frame("metabolite_name" = rownames(proc.data))
met_key$ID <- paste("metabolite",1:nrow(met_key), sep="_")
rownames(proc.data) <- replace_ids(myIDs=rownames(proc.data), mykey=met_key[,c("metabolite_name","ID")])
proc.data <- data.frame(t(proc.data), stringsAsFactors = F)

#remove samples with fake mrn
proc.data <- proc.data[which(rownames(proc.data) %in% all.meta$label | grepl("pool",rownames(proc.data))),]

###Missing data
#remove 0s by imputing half the minimum for each metabolite.
x <- colSums(as.matrix(proc.data) == 0) #only 3 missing values

saverownames <- rownames(proc.data)
min_impute <- function(mycol){
  mycol <- as.numeric(mycol)
  mycol[mycol==0] <- min(mycol[mycol>0],na.rm=TRUE)/2
  return(mycol)
}

proc.data <- data.frame(apply(proc.data, 2, min_impute))
rownames(proc.data) <- saverownames

#drop patient who was mis-diagnosed. Not HCC but cholangiocarcinoma 75378629
drop <- all.meta[which(all.meta$mrn==75378629),"label"]
proc.data <- proc.data[which(rownames(proc.data)!=drop),]
all.meta <- all.meta[which(all.meta$mrn!=75378629),]
all.meta <- all.meta[which(all.meta$diagnosis %in% c("Healthy","HCC","Cirrhosis")),]


#write to file
write.csv(proc.data, file = file.path(output.dir, paste0(file.prefix, "data_processed.csv")))
print("sample of data_processed.csv")
print(proc.data[1:5,1:5])



