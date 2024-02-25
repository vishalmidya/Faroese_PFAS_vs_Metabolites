library(foreach)
library(doParallel)
library(car)
library(MASS)
library(ggplot2)
library(data.table)
library(iterators)
library(parallel)
library(carData)
library(tidyverse)

### calculate FI

cores=detectCores()
cl <- makeCluster(5) 
registerDoParallel(cl)
# /sc/arion/projects/Faroese/FI
# C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese
all_sig_hits<- read.csv("/sc/arion/projects/Faroese/pfas_met/FI/sig_metabolites_longitudinal.csv")
FI<- data.frame(PFAS_age=rep(NA,37),
                Age=rep(NA,37),
                PFAS=rep(NA, 37),
                Met_id=rep(NA,37),
                Mode=rep(NA,37),
                lower_FI=rep(NA,37),
                higher_FI=rep(NA,37))

######################## PFAS at 0 #######################
### pantothenic acid
#### met at 7 (HILIC)
matched_data_pfos_at_0_met_at_7<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/matched_data_pfos_at_0_met_at_7.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 7, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met1269"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_7
  data_permuted[, "Met1269"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met1269"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met1269 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_7 <- suppressMessages(lm(Met1269 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[1, ]<- c("0", "7", "PFOS", "Met1269", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 14 (HILIC)
matched_data_pfos_at_0_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/matched_data_pfos_at_0_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 14, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met1269"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_14
  data_permuted[, "Met1269"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met1269"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met1269 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_7 <- suppressMessages(lm(Met1269 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[2, ]<- c("0", "14", "PFOS", "Met1269", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 22 (HILIC)
matched_data_pfos_at_0_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/matched_data_pfos_at_0_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 22, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met1269"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_22
  data_permuted[, "Met1269"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met1269"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met1269 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_7 <- suppressMessages(lm(Met1269 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[3, ]<- c("0", "22", "PFOS", "Met1269", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (HILIC)
matched_data_pfos_at_0_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/matched_data_pfos_at_0_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 28, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met1269"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_28
  data_permuted[, "Met1269"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met1269"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met1269 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_7 <- suppressMessages(lm(Met1269 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[4, ]<- c("0", "28", "PFOS", "Met1269", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


######################################################

### Heptanoic acid
#### met at 7 (HILIC)
matched_data_pfna_at_0_met_at_7<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_0/minerva_data_pfna_0_metabolites_7/matched_data_pfna_at_0_met_at_7.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 7, PFAS == "PFNA", Mode == "C18", Met_id == "Met87"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_0_met_at_7
  data_permuted[, "Met87"][data_permuted[, "cpfna0"] == 1] <- (data_permuted[, "Met87"][data_permuted[, "cpfna0"] == 1] + values[j] - ts)
  s <- summary(lm(Met87 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna0"] <- sample(data_permuted[, "cpfna0"])
    lmer_7 <- suppressMessages(lm(Met87 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[5, ]<- c("0", "7", "PFNA", "Met87", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 14 (C18)
matched_data_pfna_at_0_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_0/minerva_data_pfna_0_metabolites_14/matched_data_pfna_at_0_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 14, PFAS == "PFNA", Mode == "C18", Met_id == "Met87"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_0_met_at_14
  data_permuted[, "Met87"][data_permuted[, "cpfna0"] == 1] <- (data_permuted[, "Met87"][data_permuted[, "cpfna0"] == 1] + values[j] - ts)
  s <- summary(lm(Met87 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna0"] <- sample(data_permuted[, "cpfna0"])
    lmer_7 <- suppressMessages(lm(Met87 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[6, ]<- c("0", "14", "PFNA", "Met87", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 22 (C18)
matched_data_pfna_at_0_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_0/minerva_data_pfna_0_metabolites_22/matched_data_pfna_at_0_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 22, PFAS == "PFNA", Mode == "C18", Met_id == "Met87"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_0_met_at_22
  data_permuted[, "Met87"][data_permuted[, "cpfna0"] == 1] <- (data_permuted[, "Met87"][data_permuted[, "cpfna0"] == 1] + values[j] - ts)
  s <- summary(lm(Met87 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna0"] <- sample(data_permuted[, "cpfna0"])
    lmer_7 <- suppressMessages(lm(Met87 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[7, ]<- c("0", "22", "PFNA", "Met87", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (C18)
matched_data_pfna_at_0_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_0/minerva_data_pfna_0_metabolites_28/matched_data_pfna_at_0_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 28, PFAS == "PFNA", Mode == "C18", Met_id == "Met87"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_0_met_at_28
  data_permuted[, "Met87"][data_permuted[, "cpfna0"] == 1] <- (data_permuted[, "Met87"][data_permuted[, "cpfna0"] == 1] + values[j] - ts)
  s <- summary(lm(Met87 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna0"] <- sample(data_permuted[, "cpfna0"])
    lmer_7 <- suppressMessages(lm(Met87 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[8, ]<- c("0", "28", "PFNA", "Met87", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



######################################################

### DL-Glutamate

#### met at 7 (C18)
matched_data_pfos_at_0_met_at_7<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/matched_data_pfos_at_0_met_at_7.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 7, PFAS == "PFOS", Mode == "C18", Met_id == "Met253"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_7
  data_permuted[, "Met253"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met253"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met253 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_7 <- suppressMessages(lm(Met253 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[9, ]<- c("0", "7", "PFOS", "Met253", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 14 (C18)
matched_data_pfos_at_0_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/matched_data_pfos_at_0_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 14, PFAS == "PFOS", Mode == "C18", Met_id == "Met253"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_14
  data_permuted[, "Met253"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met253"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met253 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_7 <- suppressMessages(lm(Met253 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[10, ]<- c("0", "14", "PFOS", "Met253", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 22 (C18)
matched_data_pfos_at_0_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/matched_data_pfos_at_0_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 22, PFAS == "PFOS", Mode == "C18", Met_id == "Met253"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_22
  data_permuted[, "Met253"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met253"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met253 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_7 <- suppressMessages(lm(Met253 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[11, ]<- c("0", "22", "PFOS", "Met253", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (C18)
matched_data_pfos_at_0_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/matched_data_pfos_at_0_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 28, PFAS == "PFOS", Mode == "C18", Met_id == "Met253"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_28
  data_permuted[, "Met253"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met253"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met253 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_7 <- suppressMessages(lm(Met253 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[12, ]<- c("0", "28", "PFOS", "Met253", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))







######################## PFAS at 7 

### Asymmetric dimethylarginine/Symmetric dimethylarginine (Dimethyl-Arginine)

#### met at 14 (HILIC)
matched_data_pfoa_at_7_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_14/matched_data_pfoa_at_7_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 14, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met1055"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_7_met_at_14
  data_permuted[, "Met1055"][data_permuted[, "cpfoa7"] == 1] <- (data_permuted[, "Met1055"][data_permuted[, "cpfoa7"] == 1] + values[j] - ts)
  s <- summary(lm(Met1055 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa7"] <- sample(data_permuted[, "cpfoa7"])
    lmer_7 <- suppressMessages(lm(Met1055 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[13, ]<- c("7", "14", "PFOA", "Met1055", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 22 (HILIC)
matched_data_pfoa_at_7_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_22/matched_data_pfoa_at_7_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 22, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met1055"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_7_met_at_22
  data_permuted[, "Met1055"][data_permuted[, "cpfoa7"] == 1] <- (data_permuted[, "Met1055"][data_permuted[, "cpfoa7"] == 1] + values[j] - ts)
  s <- summary(lm(Met1055 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa7"] <- sample(data_permuted[, "cpfoa7"])
    lmer_7 <- suppressMessages(lm(Met1055 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[14, ]<- c("7", "22", "PFOA", "Met1055", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (HILIC)
matched_data_pfoa_at_7_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_28/matched_data_pfoa_at_7_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 28, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met1055"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_7_met_at_28
  data_permuted[, "Met1055"][data_permuted[, "cpfoa7"] == 1] <- (data_permuted[, "Met1055"][data_permuted[, "cpfoa7"] == 1] + values[j] - ts)
  s <- summary(lm(Met1055 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa7"] <- sample(data_permuted[, "cpfoa7"])
    lmer_7 <- suppressMessages(lm(Met1055 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[15, ]<- c("7", "28", "PFOA", "Met1055", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


######################################################

# N-Acetylneuraminic acid/N-Acetyl-a-neuraminic acid

#### met at 14 (C18)
matched_data_pfhxs_at_7_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_14/matched_data_pfhxs_at_7_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 14, PFAS == "PFHxS", Mode == "C18", Met_id == "Met707"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_7_met_at_14
  data_permuted[, "Met707"][data_permuted[, "cpfhxs7"] == 1] <- (data_permuted[, "Met707"][data_permuted[, "cpfhxs7"] == 1] + values[j] - ts)
  s <- summary(lm(Met707 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs7"] <- sample(data_permuted[, "cpfhxs7"])
    lmer_7 <- suppressMessages(lm(Met707 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[16, ]<- c("7", "14", "PFHxS", "Met707", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 22 (C18)
matched_data_pfhxs_at_7_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_22/matched_data_pfhxs_at_7_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 22, PFAS == "PFHxS", Mode == "C18", Met_id == "Met707"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_7_met_at_22
  data_permuted[, "Met707"][data_permuted[, "cpfhxs7"] == 1] <- (data_permuted[, "Met707"][data_permuted[, "cpfhxs7"] == 1] + values[j] - ts)
  s <- summary(lm(Met707 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs7"] <- sample(data_permuted[, "cpfhxs7"])
    lmer_7 <- suppressMessages(lm(Met707 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[17, ]<- c("7", "22", "PFHxS", "Met707", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (C18)
matched_data_pfhxs_at_7_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_28/matched_data_pfhxs_at_7_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 28, PFAS == "PFHxS", Mode == "C18", Met_id == "Met707"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_7_met_at_28
  data_permuted[, "Met707"][data_permuted[, "cpfhxs7"] == 1] <- (data_permuted[, "Met707"][data_permuted[, "cpfhxs7"] == 1] + values[j] - ts)
  s <- summary(lm(Met707 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs7"] <- sample(data_permuted[, "cpfhxs7"])
    lmer_7 <- suppressMessages(lm(Met707 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[18, ]<- c("7", "28", "PFHxS", "Met707", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



######################################################

# N-Acetylneuraminic acid/N-Acetyl-a-neuraminic acid

#### met at 14 (C18)
matched_data_pfhxs_at_7_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_14/matched_data_pfhxs_at_7_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 14, PFAS == "PFHxS", Mode == "C18", Met_id == "Met695"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_7_met_at_14
  data_permuted[, "Met695"][data_permuted[, "cpfhxs7"] == 1] <- (data_permuted[, "Met695"][data_permuted[, "cpfhxs7"] == 1] + values[j] - ts)
  s <- summary(lm(Met695 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs7"] <- sample(data_permuted[, "cpfhxs7"])
    lmer_7 <- suppressMessages(lm(Met695 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[19, ]<- c("7", "14", "PFHxS", "Met695", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 22 (C18)
matched_data_pfhxs_at_7_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_22/matched_data_pfhxs_at_7_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 22, PFAS == "PFHxS", Mode == "C18", Met_id == "Met695"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_7_met_at_22
  data_permuted[, "Met695"][data_permuted[, "cpfhxs7"] == 1] <- (data_permuted[, "Met695"][data_permuted[, "cpfhxs7"] == 1] + values[j] - ts)
  s <- summary(lm(Met695 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs7"] <- sample(data_permuted[, "cpfhxs7"])
    lmer_7 <- suppressMessages(lm(Met695 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[20, ]<- c("7", "22", "PFHxS", "Met695", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (C18)
matched_data_pfhxs_at_7_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_28/matched_data_pfhxs_at_7_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 28, PFAS == "PFHxS", Mode == "C18", Met_id == "Met695"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_7_met_at_28
  data_permuted[, "Met695"][data_permuted[, "cpfhxs7"] == 1] <- (data_permuted[, "Met695"][data_permuted[, "cpfhxs7"] == 1] + values[j] - ts)
  s <- summary(lm(Met695 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs7"] <- sample(data_permuted[, "cpfhxs7"])
    lmer_7 <- suppressMessages(lm(Met695 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[21, ]<- c("7", "28", "PFHxS", "Met695", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


########################################### 
# PFAS at 14


# Adenine

#### met at 22 (HILIC)
matched_data_pfda_at_14_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_14/minerva_data_pfda_14_metabolites_22/matched_data_pfda_at_14_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 22, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met459"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_14_met_at_22
  data_permuted[, "Met459"][data_permuted[, "cpfda14"] == 1] <- (data_permuted[, "Met459"][data_permuted[, "cpfda14"] == 1] + values[j] - ts)
  s <- summary(lm(Met459 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda14"] <- sample(data_permuted[, "cpfda14"])
    lmer_7 <- suppressMessages(lm(Met459 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[22, ]<- c("14", "22", "PFDA", "Met459", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (HILIC)
matched_data_pfda_at_14_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_14/minerva_data_pfda_14_metabolites_28/matched_data_pfda_at_14_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met459"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_14_met_at_28
  data_permuted[, "Met459"][data_permuted[, "cpfda14"] == 1] <- (data_permuted[, "Met459"][data_permuted[, "cpfda14"] == 1] + values[j] - ts)
  s <- summary(lm(Met459 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda14"] <- sample(data_permuted[, "cpfda14"])
    lmer_7 <- suppressMessages(lm(Met459 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[23, ]<- c("14", "28", "PFDA", "Met459", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


########################################### 
# Iminodiacetic acid

#### met at 22 (HILIC)
matched_data_pfhxs_at_14_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_22/matched_data_pfhxs_at_14_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 22, PFAS == "PFHxS", Mode == "HILIC", Met_id == "Met407"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_14_met_at_22
  data_permuted[, "Met407"][data_permuted[, "cpfhxs14"] == 1] <- (data_permuted[, "Met407"][data_permuted[, "cpfhxs14"] == 1] + values[j] - ts)
  s <- summary(lm(Met407 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs14"] <- sample(data_permuted[, "cpfhxs14"])
    lmer_7 <- suppressMessages(lm(Met407 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[24, ]<- c("14", "22", "PFHxS", "Met407", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (HILIC)
matched_data_pfhxs_at_14_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_28/matched_data_pfhxs_at_14_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 28, PFAS == "PFHxS", Mode == "HILIC", Met_id == "Met407"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_14_met_at_28
  data_permuted[, "Met407"][data_permuted[, "cpfhxs14"] == 1] <- (data_permuted[, "Met407"][data_permuted[, "cpfhxs14"] == 1] + values[j] - ts)
  s <- summary(lm(Met407 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs14"] <- sample(data_permuted[, "cpfhxs14"])
    lmer_7 <- suppressMessages(lm(Met407 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[25, ]<- c("14", "28", "PFHxS", "Met407", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


########################################### 
# Iminodiacetic acid

#### met at 22 (HILIC)
matched_data_pfhxs_at_14_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_22/matched_data_pfhxs_at_14_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 22, PFAS == "PFHxS", Mode == "HILIC", Met_id == "Met363"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_14_met_at_22
  data_permuted[, "Met363"][data_permuted[, "cpfhxs14"] == 1] <- (data_permuted[, "Met363"][data_permuted[, "cpfhxs14"] == 1] + values[j] - ts)
  s <- summary(lm(Met363 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs14"] <- sample(data_permuted[, "cpfhxs14"])
    lmer_7 <- suppressMessages(lm(Met363 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[26, ]<- c("14", "22", "PFHxS", "Met363", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (HILIC)
matched_data_pfhxs_at_14_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_28/matched_data_pfhxs_at_14_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 28, PFAS == "PFHxS", Mode == "HILIC", Met_id == "Met363"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_14_met_at_28
  data_permuted[, "Met363"][data_permuted[, "cpfhxs14"] == 1] <- (data_permuted[, "Met363"][data_permuted[, "cpfhxs14"] == 1] + values[j] - ts)
  s <- summary(lm(Met363 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs14"] <- sample(data_permuted[, "cpfhxs14"])
    lmer_7 <- suppressMessages(lm(Met363 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[27, ]<- c("14", "28", "PFHxS", "Met363", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))

########################################### 
# Alpha-Tocopherol

#### met at 22 (HILIC)
matched_data_pfoa_at_14_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/matched_data_pfoa_at_14_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 22, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met1936"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_14_met_at_22
  data_permuted[, "Met1936"][data_permuted[, "cpfoa14"] == 1] <- (data_permuted[, "Met1936"][data_permuted[, "cpfoa14"] == 1] + values[j] - ts)
  s <- summary(lm(Met1936 ~ cpfoa14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa14"] <- sample(data_permuted[, "cpfoa14"])
    lmer_7 <- suppressMessages(lm(Met1936 ~ cpfoa14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[28, ]<- c("14", "22", "PFOA", "Met1936", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (HILIC)
matched_data_pfoa_at_14_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/matched_data_pfoa_at_14_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 28, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met1936"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_14_met_at_28
  data_permuted[, "Met1936"][data_permuted[, "cpfoa14"] == 1] <- (data_permuted[, "Met1936"][data_permuted[, "cpfoa14"] == 1] + values[j] - ts)
  s <- summary(lm(Met1936 ~ cpfoa14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa14"] <- sample(data_permuted[, "cpfoa14"])
    lmer_7 <- suppressMessages(lm(Met1936 ~ cpfoa14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[29, ]<- c("14", "28", "PFOA", "Met1936", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


########################################### 
# L-Octanoylcarnitine

#### met at 22 (HILIC)
matched_data_pfda_at_14_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_14/minerva_data_pfda_14_metabolites_22/matched_data_pfda_at_14_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 22, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met1571"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_14_met_at_22
  data_permuted[, "Met1571"][data_permuted[, "cpfda14"] == 1] <- (data_permuted[, "Met1571"][data_permuted[, "cpfda14"] == 1] + values[j] - ts)
  s <- summary(lm(Met1571 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda14"] <- sample(data_permuted[, "cpfda14"])
    lmer_7 <- suppressMessages(lm(Met1571 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[30, ]<- c("14", "22", "PFDA", "Met1571", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (HILIC)
matched_data_pfda_at_14_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_14/minerva_data_pfda_14_metabolites_28/matched_data_pfda_at_14_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met1571"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_14_met_at_28
  data_permuted[, "Met1571"][data_permuted[, "cpfda14"] == 1] <- (data_permuted[, "Met1571"][data_permuted[, "cpfda14"] == 1] + values[j] - ts)
  s <- summary(lm(Met1571 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda14"] <- sample(data_permuted[, "cpfda14"])
    lmer_7 <- suppressMessages(lm(Met1571 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[31, ]<- c("14", "28", "PFDA", "Met1571", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


########################################### 
# 81114-Eicosatrienoic acid (linolenic acid) omega-6

#### met at 22 (C18)
matched_data_pfna_at_14_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_14/minerva_data_pfna_14_metabolites_22/matched_data_pfna_at_14_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 22, PFAS == "PFNA", Mode == "C18", Met_id == "Met676"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_14_met_at_22
  data_permuted[, "Met676"][data_permuted[, "cpfna14"] == 1] <- (data_permuted[, "Met676"][data_permuted[, "cpfna14"] == 1] + values[j] - ts)
  s <- summary(lm(Met676 ~ cpfna14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna14"] <- sample(data_permuted[, "cpfna14"])
    lmer_7 <- suppressMessages(lm(Met676 ~ cpfna14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[32, ]<- c("14", "22", "PFNA", "Met676", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




#### met at 28 (C18)
matched_data_pfna_at_14_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_14/minerva_data_pfna_14_metabolites_28/matched_data_pfna_at_14_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 28, PFAS == "PFNA", Mode == "C18", Met_id == "Met676"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_14_met_at_28
  data_permuted[, "Met676"][data_permuted[, "cpfna14"] == 1] <- (data_permuted[, "Met676"][data_permuted[, "cpfna14"] == 1] + values[j] - ts)
  s <- summary(lm(Met676 ~ cpfna14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna14"] <- sample(data_permuted[, "cpfna14"])
    lmer_7 <- suppressMessages(lm(Met676 ~ cpfna14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[33, ]<- c("14", "28", "PFNA", "Met676", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))

########################################### 
# PFAS at 22

# L-Aspartic acid/D-Aspartic acid/Iminodiacetic acid

#### met at 28 (HILIC)
matched_data_pfoa_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/matched_data_pfoa_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met398"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_22_met_at_28
  data_permuted[, "Met398"][data_permuted[, "cpfoa22"] == 1] <- (data_permuted[, "Met398"][data_permuted[, "cpfoa22"] == 1] + values[j] - ts)
  s <- summary(lm(Met398 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa22"] <- sample(data_permuted[, "cpfoa22"])
    lmer_7 <- suppressMessages(lm(Met398 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[34, ]<- c("22", "28", "PFOA", "Met398", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


########################################### 
# Ureidopropionic acid

#### met at 28 (HILIC)
matched_data_pfda_at_22_met_at_28<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met389"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met389"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met389"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met389 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_7 <- suppressMessages(lm(Met389 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[35, ]<- c("22", "28", "PFDA", "Met389", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


########################################### 
# Midazolam

#### met at 28 (HILIC)
matched_data_pfhxs_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/hilic/pfhxs/pfhxs_22/minerva_data_pfhxs_22_metabolites_28/matched_data_pfhxs_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFHxS", Mode == "HILIC", Met_id == "Met1771"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_22_met_at_28
  data_permuted[, "Met1771"][data_permuted[, "cpfhxs22"] == 1] <- (data_permuted[, "Met1771"][data_permuted[, "cpfhxs22"] == 1] + values[j] - ts)
  s <- summary(lm(Met1771 ~ cpfhxs22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs22"] <- sample(data_permuted[, "cpfhxs22"])
    lmer_7 <- suppressMessages(lm(Met1771 ~ cpfhxs22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[36, ]<- c("22", "28", "PFHxS", "Met1771", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))

########################################### 
# Betaine

#### met at 28 (HILIC)
matched_data_pfna_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/hilic/pfna/pfna_22/minerva_data_pfna_22_metabolites_28/matched_data_pfna_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFNA", Mode == "HILIC", Met_id == "Met1473"))$beta
# 
# func_val <- function(ts){
#   tol = 0.15
#   if(ts >0){
#     return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
#   }
#   else
#     return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
# }

values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))
for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_22_met_at_28
  data_permuted[, "Met1473"][data_permuted[, "cpfna22"] == 1] <- (data_permuted[, "Met1473"][data_permuted[, "cpfna22"] == 1] + values[j] - ts)
  s <- summary(lm(Met1473 ~ cpfna22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna22"] <- sample(data_permuted[, "cpfna22"])
    lmer_7 <- suppressMessages(lm(Met1473 ~ cpfna22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}


FI[37, ]<- c("22", "28", "PFNA", "Met1473", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


write.table(FI,"/sc/arion/projects/Faroese/FI/FI.txt", row.names = FALSE)



