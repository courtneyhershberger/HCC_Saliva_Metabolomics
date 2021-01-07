library(waffle)

#I have a list of rules
#column 4-6 divides the first node. 
rules$leaf <- 1:nrow(rules)
tiers <- (ncol(rules)-3)/4
leaf <- data.frame(fit$where)
leaf$node <- as.numeric(as.factor(leaf$fit.where))
leaf <- merge(leaf, all.meta[,c("mrn","diagnosis","sex")], by.x="row.names",by.y="mrn")

#reformat rules
combine <- function(tier, rules){
  start <- 4*tier
  middle <- start + 1
  end <- start + 2
  columnname <- paste0("rule_",tier)
  rule <- data.frame(x = paste(rules[,start],rules[,middle],rules[,end]))
  colnames(rule) <- columnname
  return(rule)
}
rules2 <- do.call("cbind",lapply(1:tiers, combine, rules=rules))
rules2$leaf <- 1:nrow(rules2)
rules2[rules2=="  "] <- NA

#
getnodes <- function(tier, rules2){
  if(tier==0){
    leaf=paste(rules2$leaf, collapse="_")
    df2 <- data.frame(tier=tier, leaf=leaf)
  }else{
  df <- rules2[,c(paste0("rule_",tier), "leaf")]
  colnames(df) <- c("rule","leaf")
  df <- df[which(df$rule!=" "),]
  df2 <- df %>% group_by(rule) %>% summarise(leaf = paste(unique(leaf), collapse = '_'))
  df2$tier <- tier
  df2 <- df2[,c("tier","leaf")]
  }
  return(df2)
  }
nodes <- data.frame(do.call("rbind",lapply(tiers:0, getnodes, rules2=rules2)))

#proportions in each node
propnode <- function(myrow, leaf, outcome){
  groups <- unlist(strsplit(myrow["leaf"], "_"))
  df <- leaf[which(leaf$node %in% groups),]
  df2 <- data.frame(table(df[,outcome]))
  df2$total <- sum(df2$Freq)
  #df2$tier=myrow["tier"]
  df2$leaf=myrow["leaf"]
  df2$outcome <- outcome
  return(df2)
}

mytable_diagnosis <- do.call("rbind",apply(nodes, 1, propnode, leaf=leaf, outcome="diagnosis"))
mytable_diagnosis <- spread(mytable_diagnosis, Var1, Freq)
mytable_diagnosis[is.na(mytable_diagnosis)] <- 0

mytable_sex <- do.call("rbind",apply(nodes, 1, propnode, leaf=leaf, outcome="sex"))
mytable_sex <- spread(mytable_sex, Var1, Freq)




makebreakfast <- function(myrow, dim, mycolors){
  filename <- paste0("figures/waffleplot_group", myrow["leaf"],"_", myrow["outcome"], ".pdf")
  print(filename)
  
  if(dim=="square"){
    height <- round(sqrt(as.numeric(myrow["total"])))
    width <- height
  }
  if(dim=="rectangle"){
    width <- trunc(sqrt(1.5*as.numeric(myrow["total"])))
    height <- round(as.numeric(myrow["total"])/width)
    
  }
  if(dim=="custom"){
    width <- 14
    height <- round(as.numeric(myrow["total"])/width)
  }

  waffle_input <- within(as.list(myrow), rm(total, tier, leaf, outcome))
  waffle_input[is.na(waffle_input)] <- 0
  print(waffle_input)
  w <- waffle(as.numeric(waffle_input), 
              colors=mycolors, 
              legend_pos = "none", 
              rows = height,
              size=1)
  print(w)

  pdf(filename, height = height/4, width = width/4)
  print(w)
  dev.off()
  return(w)
}

colors_diagnosis <- c(my_colors["dark yellow"],my_colors["dark red"],my_colors["dark teal"], "white")
colors_sex <- c("red", "blue", "white")

my_waffles <- apply(mytable_diagnosis[which(mytable_diagnosis$leaf %in% c("1","2","3","4","5","3_4")),], 1, makebreakfast, dim="square", mycolors=colors_diagnosis)
my_waffles <- apply(mytable_diagnosis[which(mytable_diagnosis$leaf %in% c("1_2","3_4_5")),], 1, makebreakfast, dim="rectangle", mycolors=colors_diagnosis)
my_waffles <- apply(mytable_diagnosis[which(mytable_diagnosis$leaf == "1_2_3_4_5"),], 1, makebreakfast, dim="custom", mycolors=colors_diagnosis)

my_waffles <- apply(mytable_sex[which(mytable_sex$leaf %in% c("1","2","3","4","5","3_4")),], 1, makebreakfast, dim="square", mycolors=colors_sex)
my_waffles <- apply(mytable_sex[which(mytable_sex$leaf %in% c("1_2","3_4_5")),], 1, makebreakfast, dim="rectangle", mycolors=colors_sex)
my_waffles <- apply(mytable_sex[which(mytable_sex$leaf == "1_2_3_4_5"),], 1, makebreakfast, dim="custom", mycolors=colors_sex)



