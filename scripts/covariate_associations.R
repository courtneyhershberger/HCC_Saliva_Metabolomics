#GLM
#My glm function
glmfunction <- function(control, condition, var, df){
  df <- df[which(df$diagnosis %in% c(control,condition)),]
  df$outcome <- ifelse(df$diagnosis %in% control, 0, 1)
  fit <- glm(df[,"outcome"]~df[,var], family = "binomial")
  model <- summary(fit)
  print(model)
  p.value <- model$coefficients[2,4]
  beta <- model$coefficients[2,1]
  sterr <- model$coefficients[2,2]
  myrow <- data.frame(Group1=condition, Group2=control, Variable=var, Pvalue=p.value, Beta=beta, StdErr=sterr)
  return(myrow)
}

#Age and outcome
x <- glmfunction(control="Healthy", condition="Cirrhosis", var="age_days", df=all.meta)
y <- glmfunction(control="Healthy", condition="HCC", var="age_days", df=all.meta)
z <- glmfunction(control="Cirrhosis", condition="HCC", var="age_days", df=all.meta)
age_pvalues <- do.call(rbind, list(x,y,z))
age_pvalues$FDR <- p.adjust(age_pvalues$Pvalue, method="fdr")
age_pvalues$FDR <- format.pval(age_pvalues$FDR,  eps=.001, digits=3) 
age_pvalues$Pvalue <- format.pval(age_pvalues$Pvalue, esp=.001, digits=3)
age_pvalues$Beta <- format.pval(age_pvalues$Beta, esp=.001)
age_pvalues$StdErr <- format.pval(age_pvalues$StdErr, esp=0.001)

print(ggplot(all.meta, aes(x=diagnosis, y=age_days)) + 
        geom_boxplot() + theme_classic())
print(age_pvalues)

#Sex and outcome
x <- glmfunction(control="Healthy", condition="Cirrhosis", var="sex", df=all.meta)
y <- glmfunction(control="Healthy", condition="HCC", var="sex", df=all.meta)
z <- glmfunction(control="Cirrhosis", condition="HCC", var="sex", df=all.meta)
sex_pvalues <- do.call(rbind, list(x,y,z))
sex_pvalues$FDR <- p.adjust(sex_pvalues$Pvalue, method = "fdr")
sex_pvalues$FDR <- format.pval(sex_pvalues$FDR,  eps=.001, digits=3) 
sex_pvalues$Pvalue <- format.pval(sex_pvalues$Pvalue, esp=.001, digits=3)
sex_pvalues$Beta <- format.pval(sex_pvalues$Beta, esp=.001)
sex_pvalues$StdErr <- format.pval(sex_pvalues$StdErr, esp=0.001)

print(table(all.meta$sex, all.meta$diagnosis))

print(sex_pvalues)

#smoker and outcome
x <- glmfunction(control="Healthy", condition="Cirrhosis", var="smoker", df=all.meta)
y <- glmfunction(control="Healthy", condition="HCC", var="smoker", df=all.meta)
z <- glmfunction(control="Cirrhosis", condition="HCC", var="smoker", df=all.meta)
smoker_pvalues <- do.call(rbind, list(x,y,z))
smoker_pvalues$FDR <- p.adjust(smoker_pvalues$Pvalue, method = "fdr")
smoker_pvalues$FDR <- format.pval(smoker_pvalues$FDR,  eps=.001, digits=3) 
smoker_pvalues$Pvalue <- format.pval(smoker_pvalues$Pvalue, esp=.001, digits=3)
smoker_pvalues$Beta <- format.pval(smoker_pvalues$Beta, esp=.001)
smoker_pvalues$StdErr <- format.pval(smoker_pvalues$StdErr, esp=0.001)

print(table(all.meta$smoker, all.meta$diagnosis))

print(smoker_pvalues)

