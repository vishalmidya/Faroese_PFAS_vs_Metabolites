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
#---------------------------  at age 7
m.out1.pfoa7_age7 <- matchit(cpfoa7 ~   mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity, 
                             data = merged_omics[merged_omics$Year == 7,], discard = "both", method = "full",
                             distance = "glm", caliper = 0.2)

summary(m.out1.pfoa7_age7)
f <- love.plot(m.out1.pfoa7_age7)
f + labs(x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)
# plot(m.out1.pfoa7_age7, type = "jitter", interactive = FALSE)


#---------------------------  at age 14
m.out1.pfoa7_age14 <- matchit(cpfoa7 ~  sex + mage + cmatfishpreg + age14,
                              data = merged_omics[merged_omics$Year == 14,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.1)

summary(m.out1.pfoa7_age14)
f <- love.plot(m.out1.pfoa7_age14)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


#---------------------------  at age 22
m.out1.pfoa7_age22 <- matchit(cpfoa7 ~   sex + mage + mbmi + cmatfishpreg, 
                              data = merged_omics[merged_omics$Year == 22,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.3)

summary(m.out1.pfoa7_age22)
f <- love.plot(m.out1.pfoa7_age22)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


#---------------------------  at age 28
m.out1.pfoa7_age28 <- matchit(cpfoa7 ~ mage + mbmi  + smokepreg_2  + cparity + age28, 
                              data = merged_omics[merged_omics$Year == 28,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.5)

summary(m.out1.pfoa7_age28)
f <- love.plot(m.out1.pfoa7_age28)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


## export matched data
m.out1.pfoa7_age7.matched <- match.data(m.out1.pfoa7_age7)
write.csv(m.out1.pfoa7_age7.matched, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_7/matched_data_pfoa_at_7_met_at_7.csv",
          row.names = F)

m.out1.pfoa7_age14.matched <- match.data(m.out1.pfoa7_age14)
write.csv(m.out1.pfoa7_age14.matched, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_14/matched_data_pfoa_at_7_met_at_14.csv",
          row.names = F)

m.out1.pfoa7_age22.matched <- match.data(m.out1.pfoa7_age22)
write.csv(m.out1.pfoa7_age22.matched, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_22/matched_data_pfoa_at_7_met_at_22.csv",
          row.names = F)

m.out1.pfoa7_age28.matched <- match.data(m.out1.pfoa7_age28)
write.csv(m.out1.pfoa7_age28.matched, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_28/matched_data_pfoa_at_7_met_at_28.csv",
          row.names = F)




