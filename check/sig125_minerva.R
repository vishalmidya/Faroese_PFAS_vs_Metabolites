library(readr)
library(dplyr)
library(knitr)
library(xtable)
library(glmnet)
library(corrplot)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(GGally)
library(gplots)
library(tidyr)
library(tidyverse)   # clustering algorithms
library(Hmisc)
library(data.table)




### STANDARD ERROR
std<- data.frame(PFAS_age=rep(NA,22),
                 Age=rep(NA,22),
                 PFAS=rep(NA, 22),
                 Met_id=rep(NA,22),
                 Mode=rep(NA,22),
                 std.err.125=rep(NA, 22))


#### 1.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_7/matched_data_pfoa_at_0_met_at_7.csv", check.names = F)

model <- (lm(Met372 ~ cpfoa0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[1, ]<- c("0", "7", "pfoa", "Met372", "HILIC", s)

#### 2.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/c18/pfna/pfna_0/minerva_data_pfna_0_metabolites_7/matched_data_pfna_at_0_met_at_7.csv", check.names = F)

model <- (lm(Met36 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[2, ]<- c("0", "7", "pfna", "Met36", "C18", s)

#### 3.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_14/matched_data_pfhxs_at_0_met_at_14.csv", check.names = F)

model <- (lm(Met196 ~ cpfhxs0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[3, ]<- c("0", "14", "pfhxs", "Met196", "HILIC", s)


#### 4.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/c18/pfda/pfda_0/minerva_data_pfda_0_metabolites_14/matched_data_pfda_at_0_met_at_14.csv", check.names = F)

model <- (lm(Met315 ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[4, ]<- c("0", "14", "pfda", "Met315", "C18", s)

#### 5.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/c18/pfna/pfna_0/minerva_data_pfna_0_metabolites_14/matched_data_pfna_at_0_met_at_14.csv", check.names = F)

model <- (lm(Met580 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[5, ]<- c("0", "14", "pfna", "Met580", "C18", s)


#### 6.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_22/matched_data_pfhxs_at_7_met_at_22.csv", check.names = F)

model <- (lm(Met930 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[6, ]<- c("7", "22", "pfhxs", "Met930", "HILIC", s)


#### 7.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/c18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_22/matched_data_pfhxs_at_7_met_at_22.csv", check.names = F)

model <- (lm(Met193 ~ cpfhxs7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[7, ]<- c("7", "22", "pfhxs", "Met193", "C18", s)


#### 8.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/matched_data_pfoa_at_14_met_at_22.csv", check.names = F)

model <- (lm(Met565 ~ cpfoa14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[8, ]<- c("14", "22", "pfoa", "Met565", "HILIC", s)


#### 9.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/c18/pfna/pfna_14/minerva_data_pfna_14_metabolites_22/matched_data_pfna_at_14_met_at_22.csv", check.names = F)

model <- (lm(Met149 ~ cpfna14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[9, ]<- c("14", "22", "pfna", "Met149", "C18", s)


#### 10.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_28/matched_data_pfhxs_at_14_met_at_28.csv", check.names = F)

model <- (lm(Met582 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[10, ]<- c("14", "28", "pfhxs", "Met582", "HILIC", s)


#### 11.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_28/matched_data_pfhxs_at_14_met_at_28.csv", check.names = F)

model <- (lm(Met793 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[11, ]<- c("14", "28", "pfhxs", "Met793", "HILIC", s)


#### 12.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfda/pfda_14/minerva_data_pfda_14_metabolites_28/matched_data_pfda_at_14_met_at_28.csv", check.names = F)

model <- (lm(Met793 ~ cpfda14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[12, ]<- c("14", "28", "pfda", "Met793", "HILIC", s)


#### 13.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/c18/pfna/pfna_14/minerva_data_pfna_14_metabolites_28/matched_data_pfna_at_14_met_at_28.csv", check.names = F)

model <- (lm(Met438 ~ cpfna14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[13, ]<- c("14", "28", "pfna", "Met438", "C18", s)


#### 14.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

model <- (lm(Met75 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[14, ]<- c("22", "28", "pfda", "Met75", "HILIC", s)


#### 15.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

model <- (lm(Met780 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[15, ]<- c("22", "28", "pfda", "Met780", "HILIC", s)


#### 16.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfos/pfos_22/minerva_data_pfos_22_metabolites_28/matched_data_pfos_at_22_met_at_28.csv", check.names = F)

model <- (lm(Met138 ~ cpfos22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[16, ]<- c("22", "28", "pfos", "Met138", "HILIC", s)


#### 17.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/matched_data_pfoa_at_22_met_at_28.csv", check.names = F)

model <- (lm(Met79 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[17, ]<- c("22", "28", "pfoa", "Met79", "HILIC", s)


#### 18.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/matched_data_pfoa_at_22_met_at_28.csv", check.names = F)

model <- (lm(Met86 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[18, ]<- c("22", "28", "pfoa", "Met86", "HILIC", s)


#### 19.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/matched_data_pfoa_at_22_met_at_28.csv", check.names = F)

model <- (lm(Met240 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[19, ]<- c("22", "28", "pfoa", "Met240", "HILIC", s)


#### 20.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/matched_data_pfoa_at_22_met_at_28.csv", check.names = F)

model <- (lm(Met886 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[20, ]<- c("22", "28", "pfoa", "Met886", "HILIC", s)

#### 21.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/c18/pfos/pfos_22/minerva_data_pfos_22_metabolites_28/matched_data_pfos_at_22_met_at_28.csv", check.names = F)

model <- (lm(Met24 ~ cpfos22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[21, ]<- c("22", "28", "pfos", "Met24", "C18", s)


#### 22.
data<- read.csv("/sc/arion/projects/Faroese/pfas_met_125/pfas_met/c18/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/matched_data_pfoa_at_22_met_at_28.csv", check.names = F)

model <- (lm(Met219 ~ cpfoa22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[22, ]<- c("22", "28", "pfoa", "Met219", "C18", s)


std$PFAS<- factor(std$PFAS,
                  levels = c("pfos", "pfoa", "pfna", "pfda", "pfhxs"),
                  labels = c("PFOS", "PFOA", "PFNA", "PFDA", "PFHxS"))


write.table(std,"/sc/arion/projects/Faroese/pfas_met_125/pfas_met/check/std125_sig125.txt", row.names = FALSE)








