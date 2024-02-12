library(car)
library(readr)
library(lattice)
library(nlme)
library(ggplot2)
library(GGally)
library(nnet)
library(foreign)
library(biotools)
library(glmmML)
library(MASS)
library(lme4)
library(multcomp)
library(dplyr)
library(knitr)
library(xtable)
library(kableExtra)
library(DT)
library(glmnet)
library(corrplot)
library(ggpubr)
library(lmerTest)
library(reshape2)
library(ggplot2)
library(GGally)
library(mgcv)
library(gplots)
library(tidyr)
library(cluster)
library(factoextra)
library(psych)
library(bkmr)
library(mi)
library(stargazer)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) 
library(spatstat)
library(Hmisc)
library(gtsummary)
library(blme)
library("labelled")
library(lavaan)
library(igraph)
library(grpreg)
library(robustHD)
library(gWQS)
library(gridExtra)
library(ggcorrplot)
library(BWQS)
library(qwraps2)
library(MatchIt)
library(data.table)
library(DescTools)
library(mice)
library(ggrepel)

# PFAS exposure at age 22
## HILIC +ve

### pfoa
#### met at 28

data.pfoa.22.met_at_28 <- read.table("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfoa\\pfoa_22\\minerva_data_pfoa_22_metabolites_28\\pfoa_22_met_28_beta_fisher_hilic.txt", header=TRUE)
data.pfoa.22.met_at_28$PFAS <-  rep("pfoa", nrow(data.pfoa.22.met_at_28))
data.pfoa.22.met_at_28$Age = rep("28", nrow(data.pfoa.22.met_at_28))
pfoa_met <- rbind( data.pfoa.22.met_at_28)


### pfos
#### met at 28

data.pfos.22.met_at_28 <- read.table("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfos\\pfos_22\\minerva_data_pfos_22_metabolites_28\\pfos_22_met_28_beta_fisher_hilic.txt", header=TRUE)
data.pfos.22.met_at_28$PFAS <-  rep("pfos", nrow(data.pfos.22.met_at_28))
data.pfos.22.met_at_28$Age = rep("28", nrow(data.pfos.22.met_at_28))
pfos_met <- rbind( data.pfos.22.met_at_28)

### pfda
#### met at 28

data.pfda.22.met_at_28 <- read.table("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfda\\pfda_22\\minerva_data_pfda_22_metabolites_28\\pfda_22_met_28_beta_fisher_hilic.txt", header=TRUE)
data.pfda.22.met_at_28$PFAS <-  rep("pfda", nrow(data.pfda.22.met_at_28))
data.pfda.22.met_at_28$Age = rep("28", nrow(data.pfda.22.met_at_28))
pfda_met <- rbind( data.pfda.22.met_at_28)

### pfna
#### met at 28

data.pfna.22.met_at_28 <- read.table("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfna\\pfna_22\\minerva_data_pfna_22_metabolites_28\\pfna_22_met_28_beta_fisher_hilic.txt", header=TRUE)
data.pfna.22.met_at_28$PFAS <-  rep("pfna", nrow(data.pfna.22.met_at_28))
data.pfna.22.met_at_28$Age = rep("28", nrow(data.pfna.22.met_at_28))
pfna_met <- rbind( data.pfna.22.met_at_28)

### pfhxs
#### met at 28

data.pfhxs.22.met_at_28 <- read.table("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfhxs\\pfhxs_22\\minerva_data_pfhxs_22_metabolites_28\\pfhxs_22_met_28_beta_fisher_hilic.txt", header=TRUE)
data.pfhxs.22.met_at_28$PFAS <-  rep("pfhxs", nrow(data.pfhxs.22.met_at_28))
data.pfhxs.22.met_at_28$Age = rep("28", nrow(data.pfhxs.22.met_at_28))
pfhxs_met <- rbind( data.pfhxs.22.met_at_28)



pfas_met_tab_hilic <- rbind(pfhxs_met,pfda_met,pfna_met,pfos_met,pfoa_met)
pfas_met_tab_hilic$Mode <- rep("HILIC", nrow(pfas_met_tab_hilic))
write.table(pfas_met_tab_hilic, "C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfas_22_metabolites_all_HILIC.txt", row.names = F)



