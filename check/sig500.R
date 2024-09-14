
library(car)
library(readr)
library(dplyr)
library(knitr)
library(xtable)
library(glmnet)
library(corrplot)
library(ggpubr)
library("EnvStats")
library(lmerTest)
library("merTools")
library(reshape2)
library(ggplot2)
library(GGally)
library(mgcv)
library(gplots)
library(tidyr)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(Hmisc)
library(gtsummary)
library(data.table)





### STANDARD ERROR
std<- data.frame(PFAS_age=rep(NA,9),
                 Age=rep(NA,9),
                 PFAS=rep(NA, 9),
                 Met_id=rep(NA,9),
                 Mode=rep(NA,9),
                 std.err.500=rep(NA, 9))


#### 1.
data<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_0/minerva_data_pfna_0_metabolites_7/matched_data_pfna_at_0_met_at_7.csv", check.names = F)

model <- (lm(Met36 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[1, ]<- c("0", "7", "pfna", "Met36", "C18", s)


#### 2.
data<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/matched_data_pfos_at_0_met_at_14.csv", check.names = F)

model <- (lm(Met258 ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[2, ]<- c("0", "14", "pfos", "Met258", "HILIC", s)

#### 3.
data<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_0/minerva_data_pfna_0_metabolites_22/matched_data_pfna_at_0_met_at_22.csv", check.names = F)

model <- (lm(Met8 ~ cpfna0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[3, ]<- c("0", "22", "pfna", "Met8", "C18", s)

#### 4.
data<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/matched_data_pfna_at_7_met_at_28.csv", check.names = F)

model <- (lm(Met322 ~ cpfna7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[4, ]<- c("7", "28", "pfna", "Met322", "HILIC", s)

#### 5.
data<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_7/minerva_data_pfos_7_metabolites_28/matched_data_pfos_at_7_met_at_28.csv", check.names = F)

model <- (lm(Met5 ~ cpfos7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[5, ]<- c("7", "28", "pfos", "Met5", "HILIC", s)

#### 6.
data<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_7/minerva_data_pfos_7_metabolites_28/matched_data_pfos_at_7_met_at_28.csv", check.names = F)

model <- (lm(Met259 ~ cpfos7 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[6, ]<- c("7", "28", "pfos", "Met259", "HILIC", s)


#### 7.
data<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_22/matched_data_pfhxs_at_14_met_at_22.csv", check.names = F)

model <- (lm(Met895 ~ cpfhxs14 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[7, ]<- c("14", "22", "pfhxs", "Met895", "HILIC", s)


#### 8.
data<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/matched_data_pfda_at_22_met_at_28.csv", check.names = F)

model <- (lm(Met780 ~ cpfda22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[8, ]<- c("22", "28", "pfda", "Met780", "HILIC", s)


#### 9.
data<- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_22/minerva_data_pfhxs_22_metabolites_28/matched_data_pfhxs_at_22_met_at_28.csv", check.names = F)

model <- (lm(Met434 ~ cpfhxs22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
s <- summary(model)$coefficients[2, "Std. Error"]


std[9, ]<- c("22", "28", "pfhxs", "Met434", "C18", s)






std$PFAS<- factor(std$PFAS,
                  levels = c("pfos", "pfoa", "pfna", "pfda", "pfhxs"),
                  labels = c("PFOS", "PFOA", "PFNA", "PFDA", "PFHxS"))



write.table(std,"~/Projects/Faroese_Minerva/check/std500_sig500.txt", row.names = FALSE)


