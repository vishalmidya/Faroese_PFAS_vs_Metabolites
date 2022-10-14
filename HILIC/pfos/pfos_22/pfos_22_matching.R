# Full matching with estimated propensity score with caliper adjustment - if some of the covariates are already matched, they are taken out one by one and the reamining ones are matched

library(optmatch)
library(MatchIt)

sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7

merged_omics <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/New_faroese/HILIC/merged_omics_hilic.csv", check.names = F)


m.out1.pfos22_age22 <- matchit(cpfos22 ~  sex    + smokepreg_2 + cmatfishpreg  + cparity + age22, 
                              data = merged_omics[merged_omics$Year == 22,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.05)

summary(m.out1.pfos22_age22)
f <- love.plot(m.out1.pfos22_age22)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


m.out1.pfos22_age28 <- matchit(cpfos22 ~  sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, 
                              data = merged_omics[merged_omics$Year == 28,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.5)

summary(m.out1.pfos22_age28)
f <- love.plot(m.out1.pfos22_age28)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


m.out1.pfos22_age22.matched <- match.data(m.out1.pfos22_age22)
write.csv(m.out1.pfos22_age22.matched, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 22/minerva_data_pfos_22_metabolite_22/matched_data_pfos_at_22_met_at_22.csv",
          row.names = F)


m.out1.pfos22_age28.matched <- match.data(m.out1.pfos22_age28)
write.csv(m.out1.pfos22_age28.matched, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 22/minerva_data_pfos_22_metabolite_28/matched_data_pfos_at_22_met_at_28.csv",
          row.names = F)



