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

# PFAS exposure at age 7
## C18 

### pfoa
#### met at 14

data.pfoa.7.met_at_14 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_14/pfoa_7_met_14_beta_fisher_C18.txt", header=TRUE)
data.pfoa.7.met_at_14$PFAS <-  rep("pfoa", nrow(data.pfoa.7.met_at_14))
data.pfoa.7.met_at_14$Age = rep("14", nrow(data.pfoa.7.met_at_14))

#### met at 22

data.pfoa.7.met_at_22 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_22/pfoa_7_met_22_beta_fisher_C18.txt", header=TRUE)
data.pfoa.7.met_at_22$PFAS <-  rep("pfoa", nrow(data.pfoa.7.met_at_22))
data.pfoa.7.met_at_22$Age = rep("22", nrow(data.pfoa.7.met_at_22))

#### met at 28

data.pfoa.7.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_28/pfoa_7_met_28_beta_fisher_C18.txt", header=TRUE)
data.pfoa.7.met_at_28$PFAS <-  rep("pfoa", nrow(data.pfoa.7.met_at_28))
data.pfoa.7.met_at_28$Age = rep("28", nrow(data.pfoa.7.met_at_28))

pfoa_met <- rbind(data.pfoa.7.met_at_14, data.pfoa.7.met_at_22, data.pfoa.7.met_at_28)


### pfos
#### met at 14

data.pfos.7.met_at_14 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_7/minerva_data_pfos_7_metabolites_14/pfos_7_met_14_beta_fisher_C18.txt", header=TRUE)
data.pfos.7.met_at_14$PFAS <-  rep("pfos", nrow(data.pfos.7.met_at_14))
data.pfos.7.met_at_14$Age = rep("14", nrow(data.pfos.7.met_at_14))

#### met at 22

data.pfos.7.met_at_22 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_7/minerva_data_pfos_7_metabolites_22/pfos_7_met_22_beta_fisher_C18.txt", header=TRUE)
data.pfos.7.met_at_22$PFAS <-  rep("pfos", nrow(data.pfos.7.met_at_22))
data.pfos.7.met_at_22$Age = rep("22", nrow(data.pfos.7.met_at_22))

#### met at 28

data.pfos.7.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_7/minerva_data_pfos_7_metabolites_28/pfos_7_met_28_beta_fisher_C18.txt", header=TRUE)
data.pfos.7.met_at_28$PFAS <-  rep("pfos", nrow(data.pfos.7.met_at_28))
data.pfos.7.met_at_28$Age = rep("28", nrow(data.pfos.7.met_at_28))

pfos_met <- rbind(data.pfos.7.met_at_14, data.pfos.7.met_at_22, data.pfos.7.met_at_28)



### pfna
#### met at 14

data.pfna.7.met_at_14 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_14/pfna_7_met_14_beta_fisher_C18.txt", header=TRUE)
data.pfna.7.met_at_14$PFAS <-  rep("pfna", nrow(data.pfna.7.met_at_14))
data.pfna.7.met_at_14$Age = rep("14", nrow(data.pfna.7.met_at_14))

#### met at 22

data.pfna.7.met_at_22 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_22/pfna_7_met_22_beta_fisher_C18.txt", header=TRUE)
data.pfna.7.met_at_22$PFAS <-  rep("pfna", nrow(data.pfna.7.met_at_22))
data.pfna.7.met_at_22$Age = rep("22", nrow(data.pfna.7.met_at_22))

#### met at 28

data.pfna.7.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/pfna_7_met_28_beta_fisher_C18.txt", header=TRUE)
data.pfna.7.met_at_28$PFAS <-  rep("pfna", nrow(data.pfna.7.met_at_28))
data.pfna.7.met_at_28$Age = rep("28", nrow(data.pfna.7.met_at_28))

pfna_met <- rbind(data.pfna.7.met_at_14, data.pfna.7.met_at_22, data.pfna.7.met_at_28)


### pfda
#### met at 14

data.pfda.7.met_at_14 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_7/minerva_data_pfda_7_metabolites_14/pfda_7_met_14_beta_fisher_C18.txt", header=TRUE)
data.pfda.7.met_at_14$PFAS <-  rep("pfda", nrow(data.pfda.7.met_at_14))
data.pfda.7.met_at_14$Age = rep("14", nrow(data.pfda.7.met_at_14))

#### met at 22

data.pfda.7.met_at_22 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_7/minerva_data_pfda_7_metabolites_22/pfda_7_met_22_beta_fisher_C18.txt", header=TRUE)
data.pfda.7.met_at_22$PFAS <-  rep("pfda", nrow(data.pfda.7.met_at_22))
data.pfda.7.met_at_22$Age = rep("22", nrow(data.pfda.7.met_at_22))

#### met at 28

data.pfda.7.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_7/minerva_data_pfda_7_metabolites_28/pfda_7_met_28_beta_fisher_C18.txt", header=TRUE)
data.pfda.7.met_at_28$PFAS <-  rep("pfda", nrow(data.pfda.7.met_at_28))
data.pfda.7.met_at_28$Age = rep("28", nrow(data.pfda.7.met_at_28))

pfda_met <- rbind(data.pfda.7.met_at_14, data.pfda.7.met_at_22, data.pfda.7.met_at_28)


### pfhxs
#### met at 14

data.pfhxs.7.met_at_14 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_14/pfhxs_7_met_14_beta_fisher_C18.txt", header=TRUE)
data.pfhxs.7.met_at_14$PFAS <-  rep("pfhxs", nrow(data.pfhxs.7.met_at_14))
data.pfhxs.7.met_at_14$Age = rep("14", nrow(data.pfhxs.7.met_at_14))

#### met at 22

data.pfhxs.7.met_at_22 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_22/pfhxs_7_met_22_beta_fisher_C18.txt", header=TRUE)
data.pfhxs.7.met_at_22$PFAS <-  rep("pfhxs", nrow(data.pfhxs.7.met_at_22))
data.pfhxs.7.met_at_22$Age = rep("22", nrow(data.pfhxs.7.met_at_22))

#### met at 28

data.pfhxs.7.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_28/pfhxs_7_met_28_beta_fisher_C18.txt", header=TRUE)
data.pfhxs.7.met_at_28$PFAS <-  rep("pfhxs", nrow(data.pfhxs.7.met_at_28))
data.pfhxs.7.met_at_28$Age = rep("28", nrow(data.pfhxs.7.met_at_28))

pfhxs_met <- rbind(data.pfhxs.7.met_at_14, data.pfhxs.7.met_at_22, data.pfhxs.7.met_at_28)


pfas_met_tab_C18 <- rbind(pfhxs_met,pfda_met,pfna_met,pfos_met,pfoa_met)
pfas_met_tab_C18$Mode <- rep("C18", nrow(pfas_met_tab_C18))
write.table(pfas_met_tab_C18, "/sc/arion/projects/Faroese/pfas_met/c18/pfas_7_metabolites_all_C18.txt", row.names = F)

