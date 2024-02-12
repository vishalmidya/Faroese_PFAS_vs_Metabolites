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

# PFAS exposure at age 0
## C18 +ve

### pfoa age 0 
#### met at 7

data.pfoa.0.met_at_7 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_7\\pfoa_0_met_7_beta_fisher_C18.txt", header=TRUE)
data.pfoa.0.met_at_7$PFAS <-  rep("pfoa", nrow(data.pfoa.0.met_at_7))
data.pfoa.0.met_at_7$Age = rep("7", nrow(data.pfoa.0.met_at_7))

#### met at 14

data.pfoa.0.met_at_14 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_14\\pfoa_0_met_14_beta_fisher_C18.txt", header=TRUE)
data.pfoa.0.met_at_14$PFAS <-  rep("pfoa", nrow(data.pfoa.0.met_at_14))
data.pfoa.0.met_at_14$Age = rep("14", nrow(data.pfoa.0.met_at_14))

#### met at 22

data.pfoa.0.met_at_22 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_22\\pfoa_0_met_22_beta_fisher_C18.txt", header=TRUE)
data.pfoa.0.met_at_22$PFAS <-  rep("pfoa", nrow(data.pfoa.0.met_at_22))
data.pfoa.0.met_at_22$Age = rep("22", nrow(data.pfoa.0.met_at_22))

#### met at 28

data.pfoa.0.met_at_28 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_28\\pfoa_0_met_28_beta_fisher_C18.txt", header=TRUE)
data.pfoa.0.met_at_28$PFAS <-  rep("pfoa", nrow(data.pfoa.0.met_at_28))
data.pfoa.0.met_at_28$Age = rep("28", nrow(data.pfoa.0.met_at_28))

pfoa_met <- rbind(data.pfoa.0.met_at_7,data.pfoa.0.met_at_14, data.pfoa.0.met_at_22, data.pfoa.0.met_at_28)


### pfos age 0 
#### met at 7

data.pfos.0.met_at_7 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_7\\pfos_0_met_7_beta_fisher_C18.txt", header=TRUE)
data.pfos.0.met_at_7$PFAS <-  rep("pfos", nrow(data.pfos.0.met_at_7))
data.pfos.0.met_at_7$Age = rep("7", nrow(data.pfos.0.met_at_7))

#### met at 14

data.pfos.0.met_at_14 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_14\\pfos_0_met_14_beta_fisher_C18.txt", header=TRUE)
data.pfos.0.met_at_14$PFAS <-  rep("pfos", nrow(data.pfos.0.met_at_14))
data.pfos.0.met_at_14$Age = rep("14", nrow(data.pfos.0.met_at_14))

#### met at 22

data.pfos.0.met_at_22 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_22\\pfos_0_met_22_beta_fisher_C18.txt", header=TRUE)
data.pfos.0.met_at_22$PFAS <-  rep("pfos", nrow(data.pfos.0.met_at_22))
data.pfos.0.met_at_22$Age = rep("22", nrow(data.pfos.0.met_at_22))

#### met at 28

data.pfos.0.met_at_28 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_28\\pfos_0_met_28_beta_fisher_C18.txt", header=TRUE)
data.pfos.0.met_at_28$PFAS <-  rep("pfos", nrow(data.pfos.0.met_at_28))
data.pfos.0.met_at_28$Age = rep("28", nrow(data.pfos.0.met_at_28))

pfos_met <- rbind(data.pfos.0.met_at_7,data.pfos.0.met_at_14, data.pfos.0.met_at_22, data.pfos.0.met_at_28)

### pfna age 0 
#### met at 7

data.pfna.0.met_at_7 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_7\\pfna_0_met_7_beta_fisher_C18.txt", header=TRUE)
data.pfna.0.met_at_7$PFAS <-  rep("pfna", nrow(data.pfna.0.met_at_7))
data.pfna.0.met_at_7$Age = rep("7", nrow(data.pfna.0.met_at_7))

#### met at 14

data.pfna.0.met_at_14 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_14\\pfna_0_met_14_beta_fisher_C18.txt", header=TRUE)
data.pfna.0.met_at_14$PFAS <-  rep("pfna", nrow(data.pfna.0.met_at_14))
data.pfna.0.met_at_14$Age = rep("14", nrow(data.pfna.0.met_at_14))

#### met at 22

data.pfna.0.met_at_22 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_22\\pfna_0_met_22_beta_fisher_C18.txt", header=TRUE)
data.pfna.0.met_at_22$PFAS <-  rep("pfna", nrow(data.pfna.0.met_at_22))
data.pfna.0.met_at_22$Age = rep("22", nrow(data.pfna.0.met_at_22))

#### met at 28

data.pfna.0.met_at_28 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_28\\pfna_0_met_28_beta_fisher_C18.txt", header=TRUE)
data.pfna.0.met_at_28$PFAS <-  rep("pfna", nrow(data.pfna.0.met_at_28))
data.pfna.0.met_at_28$Age = rep("28", nrow(data.pfna.0.met_at_28))

pfna_met <- rbind(data.pfna.0.met_at_7,data.pfna.0.met_at_14, data.pfna.0.met_at_22, data.pfna.0.met_at_28)


### pfda age 0 
#### met at 7

data.pfda.0.met_at_7 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_7\\pfda_0_met_7_beta_fisher_C18.txt", header=TRUE)
data.pfda.0.met_at_7$PFAS <-  rep("pfda", nrow(data.pfda.0.met_at_7))
data.pfda.0.met_at_7$Age = rep("7", nrow(data.pfda.0.met_at_7))

#### met at 14

data.pfda.0.met_at_14 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_14\\pfda_0_met_14_beta_fisher_C18.txt", header=TRUE)
data.pfda.0.met_at_14$PFAS <-  rep("pfda", nrow(data.pfda.0.met_at_14))
data.pfda.0.met_at_14$Age = rep("14", nrow(data.pfda.0.met_at_14))

#### met at 22

data.pfda.0.met_at_22 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_22\\pfda_0_met_22_beta_fisher_C18.txt", header=TRUE)
data.pfda.0.met_at_22$PFAS <-  rep("pfda", nrow(data.pfda.0.met_at_22))
data.pfda.0.met_at_22$Age = rep("22", nrow(data.pfda.0.met_at_22))

#### met at 28

data.pfda.0.met_at_28 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_28\\pfda_0_met_28_beta_fisher_C18.txt", header=TRUE)
data.pfda.0.met_at_28$PFAS <-  rep("pfda", nrow(data.pfda.0.met_at_28))
data.pfda.0.met_at_28$Age = rep("28", nrow(data.pfda.0.met_at_28))

pfda_met <- rbind(data.pfda.0.met_at_7,data.pfda.0.met_at_14, data.pfda.0.met_at_22, data.pfda.0.met_at_28)


### pfhxs age 0 
#### met at 7

data.pfhxs.0.met_at_7 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_7\\pfhxs_0_met_7_beta_fisher_C18.txt", header=TRUE)
data.pfhxs.0.met_at_7$PFAS <-  rep("pfhxs", nrow(data.pfhxs.0.met_at_7))
data.pfhxs.0.met_at_7$Age = rep("7", nrow(data.pfhxs.0.met_at_7))

#### met at 14

data.pfhxs.0.met_at_14 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_14\\pfhxs_0_met_14_beta_fisher_C18.txt", header=TRUE)
data.pfhxs.0.met_at_14$PFAS <-  rep("pfhxs", nrow(data.pfhxs.0.met_at_14))
data.pfhxs.0.met_at_14$Age = rep("14", nrow(data.pfhxs.0.met_at_14))

#### met at 22

data.pfhxs.0.met_at_22 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_22\\pfhxs_0_met_22_beta_fisher_C18.txt", header=TRUE)
data.pfhxs.0.met_at_22$PFAS <-  rep("pfhxs", nrow(data.pfhxs.0.met_at_22))
data.pfhxs.0.met_at_22$Age = rep("22", nrow(data.pfhxs.0.met_at_22))

#### met at 28

data.pfhxs.0.met_at_28 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_28\\pfhxs_0_met_28_beta_fisher_C18.txt", header=TRUE)
data.pfhxs.0.met_at_28$PFAS <-  rep("pfhxs", nrow(data.pfhxs.0.met_at_28))
data.pfhxs.0.met_at_28$Age = rep("28", nrow(data.pfhxs.0.met_at_28))

pfhxs_met <- rbind(data.pfhxs.0.met_at_7,data.pfhxs.0.met_at_14, data.pfhxs.0.met_at_22, data.pfhxs.0.met_at_28)


pfas_met_tab_C18 <- rbind(pfhxs_met,pfda_met,pfna_met,pfos_met,pfoa_met)
pfas_met_tab_C18$Mode <- rep("C18", nrow(pfas_met_tab_C18))
write.table(pfas_met_tab_C18, "C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfas_0_metabolites_all_C18.txt", row.names = F)

