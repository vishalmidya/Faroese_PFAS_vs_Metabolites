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
# /sc/arion/projects/Faroese/pfas_met/FI/sig_metabolites_longitudinal.csv


all_sig_hits<- read.csv("/sc/arion/projects/Faroese/pfas_met/FI/sig_metabolites_longitudinal.csv")

FI<- data.frame(PFAS_age=rep(NA,49),
                Age=rep(NA,49),
                PFAS=rep(NA, 49),
                Met_id=rep(NA,49),
                Mode=rep(NA,49),
                lower_FI=rep(NA,49),
                higher_FI=rep(NA,49))


################ function ################
# FI<- function(data, Met_id, PFAS_time,  Met_id_no, PFAS_time_no, age){
# for(j in 1:length(values)){
#   pp <- NA_real_
#   data_permuted <- data
#   data_permuted[, Met_id][data_permuted[, PFAS_time] == 1] <- (data_permuted[, Met_id][data_permuted[, PFAS_time] == 1] + values[j] - ts)
#   s <- summary(lm(Met_id_no ~ PFAS_time_no + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age, data=data_permuted))
#   pp <- foreach(i=1:M, .combine='c') %dopar% {
#     set.seed(runif(1, 0, 1e4))
#     data_permuted[, PFAS_time] <- sample(data_permuted[, PFAS_time])
#     lmer_7 <- suppressMessages(lm(Met_id_no ~ PFAS_time_no + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age, data=data_permuted))
#     s_lmer_7 <- summary(lmer_7)
#     pp <- s_lmer_7$coefficients[2,"Estimate"]
#   }
#   
#   pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
# }
# }

######################## PFAS at 0 #######################
### Cimetidine desaturated VS PFDA ###
### Met at 7 
matched_data_pfda_at_0_met_at_7<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_0/minerva_data_pfda_0_metabolites_7/matched_data_pfda_at_0_met_at_7.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
      filter(PFAS_age == 0, Age == 7, PFAS == "PFDA", Mode == "C18", Met_id == "Met300"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_0_met_at_7
  data_permuted[, "Met300"][data_permuted[, "cpfda0"] == 1] <- (data_permuted[, "Met300"][data_permuted[, "cpfda0"] == 1] + values[j] - ts)
  s <- summary(lm(Met300 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda0"] <- sample(data_permuted[, "cpfda0"])
    lmer_7 <- suppressMessages(lm(Met300 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[1, ]<- c("0", "7", "PFDA", "Met300", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### met314 VS PFDA ###
### Met at 7
matched_data_pfda_at_0_met_at_7<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_0/minerva_data_pfda_0_metabolites_7/matched_data_pfda_at_0_met_at_7.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 7, PFAS == "PFDA", Mode == "C18", Met_id == "Met314"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_0_met_at_7
  data_permuted[, "Met314"][data_permuted[, "cpfda0"] == 1] <- (data_permuted[, "Met314"][data_permuted[, "cpfda0"] == 1] + values[j] - ts)
  s <- summary(lm(Met314 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda0"] <- sample(data_permuted[, "cpfda0"])
    lmer_7 <- suppressMessages(lm(Met314 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[2, ]<- c("0", "7", "PFDA", "Met314", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met258 VS PFOS ###
### Met at 7 
matched_data_pfos_at_0_met_at_7<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/matched_data_pfos_at_0_met_at_7.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 7, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met258"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_7
  data_permuted[, "Met258"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met258"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met258 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_7 <- suppressMessages(lm(Met258 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[3, ]<- c("0", "7", "PFOS", "Met258", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


### Met867 VS PFOA ###
### Met at 7 
matched_data_pfoa_at_0_met_at_7<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_7/matched_data_pfoa_at_0_met_at_7.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 7, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met867"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_0_met_at_7
  data_permuted[, "Met867"][data_permuted[, "cpfoa0"] == 1] <- (data_permuted[, "Met867"][data_permuted[, "cpfoa0"] == 1] + values[j] - ts)
  s <- summary(lm(Met867 ~ cpfoa0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa0"] <- sample(data_permuted[, "cpfoa0"])
    lmer_7 <- suppressMessages(lm(Met867 ~ cpfoa0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[4, ]<- c("0", "7", "PFOA", "Met867", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))





### Met580 VS PFNA ###
### Met at 7 
matched_data_pfna_at_0_met_at_7<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_0/minerva_data_pfna_0_metabolites_7/matched_data_pfna_at_0_met_at_7.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 7, PFAS == "PFNA", Mode == "C18", Met_id == "Met580"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_0_met_at_7
  data_permuted[, "Met580"][data_permuted[, "cpfna0"] == 1] <- (data_permuted[, "Met580"][data_permuted[, "cpfna0"] == 1] + values[j] - ts)
  s <- summary(lm(Met580 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna0"] <- sample(data_permuted[, "cpfna0"])
    lmer_7 <- suppressMessages(lm(Met580 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7 ,data=data_permuted))
    s_lmer_7 <- summary(lmer_7)
    pp <- s_lmer_7$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[5, ]<- c("0", "7", "PFNA", "Met580", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met258 VS PFOS ###
### Met at 14 
matched_data_pfos_at_0_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/matched_data_pfos_at_0_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 14, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met258"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_14
  data_permuted[, "Met258"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met258"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met258 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_14 <- suppressMessages(lm(Met258 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_14 <- summary(lmer_14)
    pp <- s_lmer_14$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}



FI[6, ]<- c("0", "14", "PFOS", "Met258", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met867 VS PFOA ###
### Met at 14 
matched_data_pfoa_at_0_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_14/matched_data_pfoa_at_0_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 14, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met867"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_0_met_at_14
  data_permuted[, "Met867"][data_permuted[, "cpfoa0"] == 1] <- (data_permuted[, "Met867"][data_permuted[, "cpfoa0"] == 1] + values[j] - ts)
  s <- summary(lm(Met867 ~ cpfoa0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa0"] <- sample(data_permuted[, "cpfoa0"])
    lmer_14 <- suppressMessages(lm(Met867 ~ cpfoa0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_14 <- summary(lmer_14)
    pp <- s_lmer_14$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[7, ]<- c("0", "14", "PFOA", "Met867", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met300 VS PFDA ###
### Met at 14 
matched_data_pfda_at_0_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_0/minerva_data_pfda_0_metabolites_14/matched_data_pfda_at_0_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 14, PFAS == "PFDA", Mode == "C18", Met_id == "Met300"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_0_met_at_14
  data_permuted[, "Met300"][data_permuted[, "cpfda0"] == 1] <- (data_permuted[, "Met300"][data_permuted[, "cpfda0"] == 1] + values[j] - ts)
  s <- summary(lm(Met300 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda0"] <- sample(data_permuted[, "cpfda0"])
    lmer_14 <- suppressMessages(lm(Met300 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_14 <- summary(lmer_14)
    pp <- s_lmer_14$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[8, ]<- c("0", "14", "PFDA", "Met300", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met314 VS PFDA ###
### Met at 14 
matched_data_pfda_at_0_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_0/minerva_data_pfda_0_metabolites_14/matched_data_pfda_at_0_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 14, PFAS == "PFDA", Mode == "C18", Met_id == "Met314"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_0_met_at_14
  data_permuted[, "Met314"][data_permuted[, "cpfda0"] == 1] <- (data_permuted[, "Met314"][data_permuted[, "cpfda0"] == 1] + values[j] - ts)
  s <- summary(lm(Met314 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda0"] <- sample(data_permuted[, "cpfda0"])
    lmer_14 <- suppressMessages(lm(Met314 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_14 <- summary(lmer_14)
    pp <- s_lmer_14$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[9, ]<- c("0", "14", "PFDA", "Met314", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))





### Met580 VS PFNA ###
### Met at 14 
matched_data_pfna_at_0_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_0/minerva_data_pfna_0_metabolites_14/matched_data_pfna_at_0_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 14, PFAS == "PFNA", Mode == "C18", Met_id == "Met580"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_0_met_at_14
  data_permuted[, "Met580"][data_permuted[, "cpfna0"] == 1] <- (data_permuted[, "Met580"][data_permuted[, "cpfna0"] == 1] + values[j] - ts)
  s <- summary(lm(Met580 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna0"] <- sample(data_permuted[, "cpfna0"])
    lmer_14 <- suppressMessages(lm(Met580 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_14 <- summary(lmer_14)
    pp <- s_lmer_14$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[10, ]<- c("0", "14", "PFNA", "Met580", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))





### Met258 VS PFOS ###
### Met at 22
matched_data_pfos_at_0_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/matched_data_pfos_at_0_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 22, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met258"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_22
  data_permuted[, "Met258"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met258"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met258 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_22 <- suppressMessages(lm(Met258 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_22 <- summary(lmer_22)
    pp <- s_lmer_22$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[11, ]<- c("0", "22", "PFOS", "Met258", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))






### Met867 VS PFOA ###
### Met at 22
matched_data_pfoa_at_0_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_22/matched_data_pfoa_at_0_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 22, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met867"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_0_met_at_22
  data_permuted[, "Met867"][data_permuted[, "cpfoa0"] == 1] <- (data_permuted[, "Met867"][data_permuted[, "cpfoa0"] == 1] + values[j] - ts)
  s <- summary(lm(Met867 ~ cpfoa0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa0"] <- sample(data_permuted[, "cpfoa0"])
    lmer_22 <- suppressMessages(lm(Met867 ~ cpfoa0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_22 <- summary(lmer_22)
    pp <- s_lmer_22$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[12, ]<- c("0", "22", "PFOA", "Met867", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))







### Met300 VS PFDA ###
### Met at 22
matched_data_pfda_at_0_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_0/minerva_data_pfda_0_metabolites_22/matched_data_pfda_at_0_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 22, PFAS == "PFDA", Mode == "C18", Met_id == "Met300"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_0_met_at_22
  data_permuted[, "Met300"][data_permuted[, "cpfda0"] == 1] <- (data_permuted[, "Met300"][data_permuted[, "cpfda0"] == 1] + values[j] - ts)
  s <- summary(lm(Met300 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda0"] <- sample(data_permuted[, "cpfda0"])
    lmer_22 <- suppressMessages(lm(Met300 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_22 <- summary(lmer_22)
    pp <- s_lmer_22$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[13, ]<- c("0", "22", "PFDA", "Met300", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))








### Met314 VS PFDA ###
### Met at 22
matched_data_pfda_at_0_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_0/minerva_data_pfda_0_metabolites_22/matched_data_pfda_at_0_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 22, PFAS == "PFDA", Mode == "C18", Met_id == "Met314"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_0_met_at_22
  data_permuted[, "Met314"][data_permuted[, "cpfda0"] == 1] <- (data_permuted[, "Met314"][data_permuted[, "cpfda0"] == 1] + values[j] - ts)
  s <- summary(lm(Met314 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda0"] <- sample(data_permuted[, "cpfda0"])
    lmer_22 <- suppressMessages(lm(Met314 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_22 <- summary(lmer_22)
    pp <- s_lmer_22$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[14, ]<- c("0", "22", "PFDA", "Met314", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met580 VS PFNA ###
### Met at 22
matched_data_pfna_at_0_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_0/minerva_data_pfna_0_metabolites_22/matched_data_pfna_at_0_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 22, PFAS == "PFNA", Mode == "C18", Met_id == "Met580"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_0_met_at_22
  data_permuted[, "Met580"][data_permuted[, "cpfna0"] == 1] <- (data_permuted[, "Met580"][data_permuted[, "cpfna0"] == 1] + values[j] - ts)
  s <- summary(lm(Met580 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna0"] <- sample(data_permuted[, "cpfna0"])
    lmer_22 <- suppressMessages(lm(Met580 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_22 <- summary(lmer_22)
    pp <- s_lmer_22$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[15, ]<- c("0", "22", "PFNA", "Met580", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))





### Met258 VS PFOS ###
### Met at 28
matched_data_pfos_at_0_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/matched_data_pfos_at_0_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 28, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met258"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_0_met_at_28
  data_permuted[, "Met258"][data_permuted[, "cpfos0"] == 1] <- (data_permuted[, "Met258"][data_permuted[, "cpfos0"] == 1] + values[j] - ts)
  s <- summary(lm(Met258 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos0"] <- sample(data_permuted[, "cpfos0"])
    lmer_28 <- suppressMessages(lm(Met258 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[16, ]<- c("0", "28", "PFOS", "Met258", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met867 VS PFOA ###
### Met at 28
matched_data_pfoa_at_0_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_28/matched_data_pfoa_at_0_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 28, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met867"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_0_met_at_28
  data_permuted[, "Met867"][data_permuted[, "cpfoa0"] == 1] <- (data_permuted[, "Met867"][data_permuted[, "cpfoa0"] == 1] + values[j] - ts)
  s <- summary(lm(Met867 ~ cpfoa0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa0"] <- sample(data_permuted[, "cpfoa0"])
    lmer_28 <- suppressMessages(lm(Met867 ~ cpfoa0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[17, ]<- c("0", "28", "PFOA", "Met867", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met300 VS PFDA ###
### Met at 28
matched_data_pfda_at_0_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_0/minerva_data_pfda_0_metabolites_28/matched_data_pfda_at_0_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 28, PFAS == "PFDA", Mode == "C18", Met_id == "Met300"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_0_met_at_28
  data_permuted[, "Met300"][data_permuted[, "cpfda0"] == 1] <- (data_permuted[, "Met300"][data_permuted[, "cpfda0"] == 1] + values[j] - ts)
  s <- summary(lm(Met300 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda0"] <- sample(data_permuted[, "cpfda0"])
    lmer_28 <- suppressMessages(lm(Met300 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[18, ]<- c("0", "28", "PFDA", "Met300", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))





### Met314 VS PFDA ###
### Met at 28
matched_data_pfda_at_0_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_0/minerva_data_pfda_0_metabolites_28/matched_data_pfda_at_0_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 28, PFAS == "PFDA", Mode == "C18", Met_id == "Met314"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_0_met_at_28
  data_permuted[, "Met314"][data_permuted[, "cpfda0"] == 1] <- (data_permuted[, "Met314"][data_permuted[, "cpfda0"] == 1] + values[j] - ts)
  s <- summary(lm(Met314 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda0"] <- sample(data_permuted[, "cpfda0"])
    lmer_28 <- suppressMessages(lm(Met314 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[19, ]<- c("0", "28", "PFDA", "Met314", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met580 VS PFNA ###
### Met at 28
matched_data_pfna_at_0_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_0/minerva_data_pfna_0_metabolites_28/matched_data_pfna_at_0_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 0, Age == 28, PFAS == "PFNA", Mode == "C18", Met_id == "Met580"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_0_met_at_28
  data_permuted[, "Met580"][data_permuted[, "cpfna0"] == 1] <- (data_permuted[, "Met580"][data_permuted[, "cpfna0"] == 1] + values[j] - ts)
  s <- summary(lm(Met580 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna0"] <- sample(data_permuted[, "cpfna0"])
    lmer_28 <- suppressMessages(lm(Met580 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[20, ]<- c("0", "28", "PFNA", "Met580", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### PFAS Age at 7 ###

### Met192 VS PFOS ###
### Met at 14
matched_data_pfos_at_7_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_7/minerva_data_pfos_7_metabolites_14/matched_data_pfos_at_7_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 14, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met192"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_7_met_at_14
  data_permuted[, "Met192"][data_permuted[, "cpfos7"] == 1] <- (data_permuted[, "Met192"][data_permuted[, "cpfos7"] == 1] + values[j] - ts)
  s <- summary(lm(Met192 ~ cpfos7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos7"] <- sample(data_permuted[, "cpfos7"])
    lmer_14 <- suppressMessages(lm(Met192 ~ cpfos7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_14 <- summary(lmer_14)
    pp <- s_lmer_14$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[21, ]<- c("7", "14", "PFOS", "Met192", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))

### Met389 VS PFOA ###
### Met at 14
matched_data_pfoa_at_7_met_at_14<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_14/matched_data_pfoa_at_7_met_at_14.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 14, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met389"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_7_met_at_14
  data_permuted[, "Met389"][data_permuted[, "cpfoa7"] == 1] <- (data_permuted[, "Met389"][data_permuted[, "cpfoa7"] == 1] + values[j] - ts)
  s <- summary(lm(Met389 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa7"] <- sample(data_permuted[, "cpfoa7"])
    lmer_14 <- suppressMessages(lm(Met389 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14 ,data=data_permuted))
    s_lmer_14 <- summary(lmer_14)
    pp <- s_lmer_14$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[22, ]<- c("7", "14", "PFOA", "Met389", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


### Met192 VS PFOS ###
### Met at 22
matched_data_pfos_at_7_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_7/minerva_data_pfos_7_metabolites_22/matched_data_pfos_at_7_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 22, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met192"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_7_met_at_22
  data_permuted[, "Met192"][data_permuted[, "cpfos7"] == 1] <- (data_permuted[, "Met192"][data_permuted[, "cpfos7"] == 1] + values[j] - ts)
  s <- summary(lm(Met192 ~ cpfos7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos7"] <- sample(data_permuted[, "cpfos7"])
    lmer_22 <- suppressMessages(lm(Met192 ~ cpfos7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_22 <- summary(lmer_22)
    pp <- s_lmer_22$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[23, ]<- c("7", "22", "PFOS", "Met192", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met389 VS PFOA ###
### Met at 22
matched_data_pfoa_at_7_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_22/matched_data_pfoa_at_7_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 22, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met389"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_7_met_at_22
  data_permuted[, "Met389"][data_permuted[, "cpfoa7"] == 1] <- (data_permuted[, "Met389"][data_permuted[, "cpfoa7"] == 1] + values[j] - ts)
  s <- summary(lm(Met389 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa7"] <- sample(data_permuted[, "cpfoa7"])
    lmer_22 <- suppressMessages(lm(Met389 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_22 <- summary(lmer_22)
    pp <- s_lmer_22$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[24, ]<- c("7", "22", "PFOA", "Met389", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met192 VS PFOS ###
### Met at 28
matched_data_pfos_at_7_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_7/minerva_data_pfos_7_metabolites_28/matched_data_pfos_at_7_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 28, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met192"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_7_met_at_28
  data_permuted[, "Met192"][data_permuted[, "cpfos7"] == 1] <- (data_permuted[, "Met192"][data_permuted[, "cpfos7"] == 1] + values[j] - ts)
  s <- summary(lm(Met192 ~ cpfos7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos7"] <- sample(data_permuted[, "cpfos7"])
    lmer_28 <- suppressMessages(lm(Met192 ~ cpfos7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[25, ]<- c("7", "28", "PFOS", "Met192", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met389 VS PFOA ###
### Met at 28
matched_data_pfoa_at_7_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_28/matched_data_pfoa_at_7_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 7, Age == 28, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met389"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_7_met_at_28
  data_permuted[, "Met389"][data_permuted[, "cpfoa7"] == 1] <- (data_permuted[, "Met389"][data_permuted[, "cpfoa7"] == 1] + values[j] - ts)
  s <- summary(lm(Met389 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa7"] <- sample(data_permuted[, "cpfoa7"])
    lmer_28 <- suppressMessages(lm(Met389 ~ cpfoa7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[26, ]<- c("7", "28", "PFOA", "Met389", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


### PFAS Age = 14 ###

### Met519 VS PFHxS ###
### Met at 22
matched_data_pfhxs_at_14_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_22/matched_data_pfhxs_at_14_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 22, PFAS == "PFHxS", Mode == "HILIC", Met_id == "Met519"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_14_met_at_22
  data_permuted[, "Met519"][data_permuted[, "cpfhxs14"] == 1] <- (data_permuted[, "Met519"][data_permuted[, "cpfhxs14"] == 1] + values[j] - ts)
  s <- summary(lm(Met519 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs14"] <- sample(data_permuted[, "cpfhxs14"])
    lmer_22 <- suppressMessages(lm(Met519 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_22 <- summary(lmer_22)
    pp <- s_lmer_22$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[27, ]<- c("14", "22", "PFHxS", "Met519", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))





### Met793 VS PFDA ###
### Met at 22
matched_data_pfda_at_14_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_14/minerva_data_pfda_14_metabolites_22/matched_data_pfda_at_14_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 22, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met793"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_14_met_at_22
  data_permuted[, "Met793"][data_permuted[, "cpfda14"] == 1] <- (data_permuted[, "Met793"][data_permuted[, "cpfda14"] == 1] + values[j] - ts)
  s <- summary(lm(Met793 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda14"] <- sample(data_permuted[, "cpfda14"])
    lmer_22 <- suppressMessages(lm(Met793 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_22 <- summary(lmer_22)
    pp <- s_lmer_22$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[28, ]<- c("14", "22", "PFDA", "Met793", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met438 VS PFNA ###
### Met at 22
matched_data_pfna_at_14_met_at_22<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_14/minerva_data_pfna_14_metabolites_22/matched_data_pfna_at_14_met_at_22.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 22, PFAS == "PFNA", Mode == "C18", Met_id == "Met438"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_14_met_at_22
  data_permuted[, "Met438"][data_permuted[, "cpfna14"] == 1] <- (data_permuted[, "Met438"][data_permuted[, "cpfna14"] == 1] + values[j] - ts)
  s <- summary(lm(Met438 ~ cpfna14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna14"] <- sample(data_permuted[, "cpfna14"])
    lmer_22 <- suppressMessages(lm(Met438 ~ cpfna14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22 ,data=data_permuted))
    s_lmer_22 <- summary(lmer_22)
    pp <- s_lmer_22$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[29, ]<- c("14", "22", "PFNA", "Met438", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met519 VS PFHxS ###
### Met at 28
matched_data_pfhxs_at_14_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_28/matched_data_pfhxs_at_14_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 28, PFAS == "PFHxS", Mode == "HILIC", Met_id == "Met519"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfhxs_at_14_met_at_28
  data_permuted[, "Met519"][data_permuted[, "cpfhxs14"] == 1] <- (data_permuted[, "Met519"][data_permuted[, "cpfhxs14"] == 1] + values[j] - ts)
  s <- summary(lm(Met519 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfhxs14"] <- sample(data_permuted[, "cpfhxs14"])
    lmer_28 <- suppressMessages(lm(Met519 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[30, ]<- c("14", "28", "PFHxS", "Met519", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met793 VS PFDA ###
### Met at 28
matched_data_pfda_at_14_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_14/minerva_data_pfda_14_metabolites_28/matched_data_pfda_at_14_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met793"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_14_met_at_28
  data_permuted[, "Met793"][data_permuted[, "cpfda14"] == 1] <- (data_permuted[, "Met793"][data_permuted[, "cpfda14"] == 1] + values[j] - ts)
  s <- summary(lm(Met793 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda14"] <- sample(data_permuted[, "cpfda14"])
    lmer_28 <- suppressMessages(lm(Met793 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[31, ]<- c("14", "28", "PFDA", "Met793", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))





### Met438 VS PFNA ###
### Met at 28
matched_data_pfna_at_14_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_14/minerva_data_pfna_14_metabolites_28/matched_data_pfna_at_14_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 14, Age == 28, PFAS == "PFNA", Mode == "C18", Met_id == "Met438"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfna_at_14_met_at_28
  data_permuted[, "Met438"][data_permuted[, "cpfna14"] == 1] <- (data_permuted[, "Met438"][data_permuted[, "cpfna14"] == 1] + values[j] - ts)
  s <- summary(lm(Met438 ~ cpfna14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfna14"] <- sample(data_permuted[, "cpfna14"])
    lmer_28 <- suppressMessages(lm(Met438 ~ cpfna14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[32, ]<- c("14", "28", "PFNA", "Met438", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


### PFAS Age = 22 ###

### Met120 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met120"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met120"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met120"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met120 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_28 <- suppressMessages(lm(Met120 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[33, ]<- c("22", "28", "PFDA", "Met120", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met190 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met190"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met190"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met190"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met190 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_28 <- suppressMessages(lm(Met190 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[34, ]<- c("22", "28", "PFDA", "Met190", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met206 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met206"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met206"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met206"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met206 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_28 <- suppressMessages(lm(Met206 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[35, ]<- c("22", "28", "PFDA", "Met206", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))






### Met267 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met267"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met267"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met267"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met267 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_28 <- suppressMessages(lm(Met267 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[36, ]<- c("22", "28", "PFDA", "Met267", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met348 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met348"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met348"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met348"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met348 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_28 <- suppressMessages(lm(Met348 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[37, ]<- c("22", "28", "PFDA", "Met348", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met452 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met452"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met452"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met452"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met452 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_28 <- suppressMessages(lm(Met452 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[38, ]<- c("22", "28", "PFDA", "Met452", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met583 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met583"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met583"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met583"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met583 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_28 <- suppressMessages(lm(Met583 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[39, ]<- c("22", "28", "PFDA", "Met583", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met676 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met676"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met676"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met676"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met676 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_28 <- suppressMessages(lm(Met676 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[40, ]<- c("22", "28", "PFDA", "Met676", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met677 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met677"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met677"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met677"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met677 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_28 <- suppressMessages(lm(Met677 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[41, ]<- c("22", "28", "PFDA", "Met677", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met728 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met728"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met728"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met728"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met728 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_28 <- suppressMessages(lm(Met728 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[42, ]<- c("22", "28", "PFDA", "Met728", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met780 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "HILIC", Met_id == "Met780"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfda_at_22_met_at_28
  data_permuted[, "Met780"][data_permuted[, "cpfda22"] == 1] <- (data_permuted[, "Met780"][data_permuted[, "cpfda22"] == 1] + values[j] - ts)
  s <- summary(lm(Met780 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfda22"] <- sample(data_permuted[, "cpfda22"])
    lmer_28 <- suppressMessages(lm(Met780 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[43, ]<- c("22", "28", "PFDA", "Met780", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met610 VS PFOS ###
### Met at 28
matched_data_pfos_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_22/minerva_data_pfos_22_metabolites_28/matched_data_pfos_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met610"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_22_met_at_28
  data_permuted[, "Met610"][data_permuted[, "cpfos22"] == 1] <- (data_permuted[, "Met610"][data_permuted[, "cpfos22"] == 1] + values[j] - ts)
  s <- summary(lm(Met610 ~ cpfos22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos22"] <- sample(data_permuted[, "cpfos22"])
    lmer_28 <- suppressMessages(lm(Met610 ~ cpfos22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[44, ]<- c("22", "28", "PFOS", "Met610", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met742 VS PFOS ###
### Met at 28
matched_data_pfos_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_22/minerva_data_pfos_22_metabolites_28/matched_data_pfos_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFOS", Mode == "HILIC", Met_id == "Met742"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfos_at_22_met_at_28
  data_permuted[, "Met742"][data_permuted[, "cpfos22"] == 1] <- (data_permuted[, "Met742"][data_permuted[, "cpfos22"] == 1] + values[j] - ts)
  s <- summary(lm(Met742 ~ cpfos22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfos22"] <- sample(data_permuted[, "cpfos22"])
    lmer_28 <- suppressMessages(lm(Met742 ~ cpfos22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[45, ]<- c("22", "28", "PFOS", "Met742", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))



### Met79 VS PFOA ###
### Met at 28
matched_data_pfoa_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/matched_data_pfoa_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met79"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_22_met_at_28
  data_permuted[, "Met79"][data_permuted[, "cpfoa22"] == 1] <- (data_permuted[, "Met79"][data_permuted[, "cpfoa22"] == 1] + values[j] - ts)
  s <- summary(lm(Met79 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa22"] <- sample(data_permuted[, "cpfoa22"])
    lmer_28 <- suppressMessages(lm(Met79 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[46, ]<- c("22", "28", "PFOA", "Met79", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




### Met86 VS PFOA ###
### Met at 28
matched_data_pfoa_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/matched_data_pfoa_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met86"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_22_met_at_28
  data_permuted[, "Met86"][data_permuted[, "cpfoa22"] == 1] <- (data_permuted[, "Met86"][data_permuted[, "cpfoa22"] == 1] + values[j] - ts)
  s <- summary(lm(Met86 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa22"] <- sample(data_permuted[, "cpfoa22"])
    lmer_28 <- suppressMessages(lm(Met86 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[47, ]<- c("22", "28", "PFOA", "Met86", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


### Met886 VS PFOA ###
### Met at 28
matched_data_pfoa_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/matched_data_pfoa_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met886"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_22_met_at_28
  data_permuted[, "Met886"][data_permuted[, "cpfoa22"] == 1] <- (data_permuted[, "Met886"][data_permuted[, "cpfoa22"] == 1] + values[j] - ts)
  s <- summary(lm(Met886 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa22"] <- sample(data_permuted[, "cpfoa22"])
    lmer_28 <- suppressMessages(lm(Met886 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[48, ]<- c("22", "28", "PFOA", "Met886", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))


### Met977 VS PFOA ###
### Met at 28
matched_data_pfoa_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/matched_data_pfoa_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFOA", Mode == "HILIC", Met_id == "Met977"))$beta


values <- seq(ts-0.30, ts+0.30, 0.005)
# values <- func_val(ts)
pval_FI <- rep(NA_real_,length(values))
beta_FI <- rep(NA_real_,length(values))


for(j in 1:length(values)){
  pp <- NA_real_
  data_permuted <- matched_data_pfoa_at_22_met_at_28
  data_permuted[, "Met977"][data_permuted[, "cpfoa22"] == 1] <- (data_permuted[, "Met977"][data_permuted[, "cpfoa22"] == 1] + values[j] - ts)
  s <- summary(lm(Met977 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted))
  pp <- foreach(i=1:M, .combine='c') %dopar% {
    set.seed(runif(1, 0, 1e4))
    data_permuted[, "cpfoa22"] <- sample(data_permuted[, "cpfoa22"])
    lmer_28 <- suppressMessages(lm(Met977 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28 ,data=data_permuted))
    s_lmer_28 <- summary(lmer_28)
    pp <- s_lmer_28$coefficients[2,"Estimate"]
  }
  
  pval_FI[j] <- mean(abs(pp) > abs(s$coefficients[2,"Estimate"] - ts))
}

FI[49, ]<- c("22", "28", "PFOA", "Met977", "HILIC", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))




## export results

write.table(FI,"/sc/arion/projects/Faroese/pfas_met/FI/FI.txt", row.names = FALSE)





stopCluster(cl)











