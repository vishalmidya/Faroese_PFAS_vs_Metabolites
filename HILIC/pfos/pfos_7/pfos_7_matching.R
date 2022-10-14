# Full matching with estimated propensity score with caliper adjustment - if some of the covariates are already matched, they are taken out one by one and the reamining ones are matched

library(optmatch)
library(MatchIt)

sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7

pfos AT AGE 7 

# pfos7

# sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7

m.out1.pfos7_age7 <- matchit(cpfos7 ~ sex + mage   + smokepreg_2 + cmatfishpreg  + cparity  , 
                             data = merged_omics[merged_omics$Year == 7,], discard = "both", method = "full",
                             distance = "glm", caliper = 0.5)

summary(m.out1.pfos7_age7)
f <- love.plot(m.out1.pfos7_age7)
f + labs(title  = "caliper = 0.5\n Metabolites at age 7 vs. pfos ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


m.out1.pfos7_age14 <- matchit(cpfos7 ~  sex + mage   + smokepreg_2 + cmatfishpreg   + age14 ,
                              data = merged_omics[merged_omics$Year == 14,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.5)

summary(m.out1.pfos7_age14)
f <- love.plot(m.out1.pfos7_age14)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


m.out1.pfos7_age22 <- matchit(cpfos7 ~  sex + mage + mbmi  + smokepreg_2   + age22, 
                              data = merged_omics[merged_omics$Year == 22,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.1)

summary(m.out1.pfos7_age22)
f <- love.plot(m.out1.pfos7_age22)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


m.out1.pfos7_age28 <- matchit(cpfos7 ~  sex + mage   + smokepreg_2 + cmatfishpreg   + age28, 
                              data = merged_omics[merged_omics$Year == 28,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.5)

summary(m.out1.pfos7_age28)
f <- love.plot(m.out1.pfos7_age28)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)

m.out1.pfos7_age7.matched <- match.data(m.out1.pfos7_age7)
write.csv(m.out1.pfos7_age7.matched, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 7/minerva_data_pfos_7_metabolite_7/matched_data_pfos_at_7_met_at_7.csv",
          row.names = F)

m.out1.pfos7_age14.matched <- match.data(m.out1.pfos7_age14)
write.csv(m.out1.pfos7_age14.matched, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 7/minerva_data_pfos_7_metabolite_14/matched_data_pfos_at_7_met_at_14.csv",
          row.names = F)

m.out1.pfos7_age22.matched <- match.data(m.out1.pfos7_age22)
write.csv(m.out1.pfos7_age22.matched, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 7/minerva_data_pfos_7_metabolite_22/matched_data_pfos_at_7_met_at_22.csv",
          row.names = F)


m.out1.pfos7_age28.matched <- match.data(m.out1.pfos7_age28)
write.csv(m.out1.pfos7_age28.matched, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 7/minerva_data_pfos_7_metabolite_28/matched_data_pfos_at_7_met_at_28.csv",
          row.names = F)


