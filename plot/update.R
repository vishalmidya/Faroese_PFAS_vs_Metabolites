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


all_sig_hits<- read.csv("/sc/arion/projects/Faroese/pfas_met/FI/sig_metabolites_longitudinal.csv")



FI<- data.frame(PFAS_age=rep(NA,1),
                Age=rep(NA,1),
                PFAS=rep(NA, 1),
                Met_id=rep(NA,1),
                Mode=rep(NA,1),
                lower_FI=rep(NA,1),
                higher_FI=rep(NA,1))


### Met452 VS PFDA ###
### Met at 28
matched_data_pfda_at_22_met_at_28<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

M <- 5e3

ts = (all_sig_hits %>% 
        filter(PFAS_age == 22, Age == 28, PFAS == "PFDA", Mode == "C18", Met_id == "Met452"))$beta


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

FI[1, ]<- c("22", "28", "PFDA", "Met452", "C18", round(min(values[pval_FI > 0.05]),4), round(max(values[pval_FI > 0.05]),4))





write.table(FI,"/sc/arion/projects/Faroese/pfas_met/FI/update.txt", row.names = FALSE)

