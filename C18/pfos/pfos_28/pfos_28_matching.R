# Full matching with estimated propensity score with caliper adjustment - if some of the covariates are already matched, they are taken out one by one and the reamining ones are matched

library(optmatch)
library(MatchIt)
library(readxl)
library(cobalt)
library(tidyverse)

## import data
merged_omics <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/merged_omics_c18.csv", check.names = F)

## all covariates: sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28

## matching
#---------------------------  at age 28
m.out1.pfos28_age28 <- matchit(cpfos28 ~  mage + mbmi  + cmatfishpreg  + cparity, 
                               data = merged_omics[merged_omics$Year == 28,], discard = "both", method = "full", 
                               distance = "glm", caliper = 0.3)

summary(m.out1.pfos28_age28)
f <- love.plot(m.out1.pfos28_age28)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


## export matched data
m.out1.pfos28_age28.matched <- match.data(m.out1.pfos28_age28)
write.csv(m.out1.pfos28_age28.matched, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_28/minerva_data_pfos_28_metabolite_28/matched_data_pfos_at_28_met_at_28.csv",
          row.names = F)


