# Full matching with estimated propensity score with caliper adjustment - if some of the covariates are already matched, they are taken out one by one and the reamining ones are matched

library(optmatch)
library(MatchIt)
library(readxl)
library(cobalt)
library(tidyverse)

## import data
merged_omics <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/merged_omics_hilic.csv", check.names = F)


## all covariates: sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7

## matching
#---------------------------  at age 14
m.out1.pfhxs14_age14 <- matchit(cpfhxs14 ~  sex + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14,
                              data = merged_omics[merged_omics$Year == 14,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.2)

summary(m.out1.pfhxs14_age14)
f <- love.plot(m.out1.pfhxs14_age14)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


#---------------------------  at age 22
m.out1.pfhxs14_age22 <- matchit(cpfhxs14 ~  sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, 
                              data = merged_omics[merged_omics$Year == 22,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.55)

summary(m.out1.pfhxs14_age22)
f <- love.plot(m.out1.pfhxs14_age22)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


#---------------------------  at age 28
m.out1.pfhxs14_age28 <- matchit(cpfhxs14 ~ sex  + smokepreg_2 + cmatfishpreg  + cparity + age28, 
                              data = merged_omics[merged_omics$Year == 28,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.5)

summary(m.out1.pfhxs14_age28)
f <- love.plot(m.out1.pfhxs14_age28)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


## export matched data
m.out1.pfhxs14_age14.matched <- match.data(m.out1.pfhxs14_age14)
write.csv(m.out1.pfhxs14_age14.matched, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_14/matched_data_pfhxs_at_14_met_at_14.csv",
          row.names = F)

m.out1.pfhxs14_age22.matched <- match.data(m.out1.pfhxs14_age22)
write.csv(m.out1.pfhxs14_age22.matched, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_22/matched_data_pfhxs_at_14_met_at_22.csv",
          row.names = F)

m.out1.pfhxs14_age28.matched <- match.data(m.out1.pfhxs14_age28)
write.csv(m.out1.pfhxs14_age28.matched, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_28/matched_data_pfhxs_at_14_met_at_28.csv",
          row.names = F)




