library(foreach)
library(doParallel)
library(car)
library(MASS)
library(ggplot2)
library(data.table)
library(iterators)
library(parallel)
library(carData)


# PFAS exposure at age 14
## C18 

### pfoa
#### met at 22

data.pfoa.14.met_at_22 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/pfoa_14_met_22_beta_fisher_c18.txt", header=TRUE)
data.pfoa.14.met_at_22$PFAS <-  rep("pfoa", nrow(data.pfoa.14.met_at_22))
data.pfoa.14.met_at_22$Age = rep("22", nrow(data.pfoa.14.met_at_22))

#### met at 28

data.pfoa.14.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/pfoa_14_met_28_beta_fisher_c18.txt", header=TRUE)
data.pfoa.14.met_at_28$PFAS <-  rep("pfoa", nrow(data.pfoa.14.met_at_28))
data.pfoa.14.met_at_28$Age = rep("28", nrow(data.pfoa.14.met_at_28))

pfoa_met <- rbind(data.pfoa.14.met_at_22, data.pfoa.14.met_at_28)

### pfos
#### met at 22

data.pfos.14.met_at_22 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_14/minerva_data_pfos_14_metabolites_22/pfos_14_met_22_beta_fisher_c18.txt", header=TRUE)
data.pfos.14.met_at_22$PFAS <-  rep("pfos", nrow(data.pfos.14.met_at_22))
data.pfos.14.met_at_22$Age = rep("22", nrow(data.pfos.14.met_at_22))

#### met at 28

data.pfos.14.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_14/minerva_data_pfos_14_metabolites_28/pfos_14_met_28_beta_fisher_c18.txt", header=TRUE)
data.pfos.14.met_at_28$PFAS <-  rep("pfos", nrow(data.pfos.14.met_at_28))
data.pfos.14.met_at_28$Age = rep("28", nrow(data.pfos.14.met_at_28))

pfos_met <- rbind(data.pfos.14.met_at_22, data.pfos.14.met_at_28)

### pfna
#### met at 22

data.pfna.14.met_at_22 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_14/minerva_data_pfna_14_metabolites_22/pfna_14_met_22_beta_fisher_c18.txt", header=TRUE)
data.pfna.14.met_at_22$PFAS <-  rep("pfna", nrow(data.pfna.14.met_at_22))
data.pfna.14.met_at_22$Age = rep("22", nrow(data.pfna.14.met_at_22))

#### met at 28

data.pfna.14.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_14/minerva_data_pfna_14_metabolites_28/pfna_14_met_28_beta_fisher_c18.txt", header=TRUE)
data.pfna.14.met_at_28$PFAS <-  rep("pfna", nrow(data.pfna.14.met_at_28))
data.pfna.14.met_at_28$Age = rep("28", nrow(data.pfna.14.met_at_28))

pfna_met <- rbind(data.pfna.14.met_at_22, data.pfna.14.met_at_28)

### pfda
#### met at 22

data.pfda.14.met_at_22 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_14/minerva_data_pfda_14_metabolites_22/pfda_14_met_22_beta_fisher_c18.txt", header=TRUE)
data.pfda.14.met_at_22$PFAS <-  rep("pfda", nrow(data.pfda.14.met_at_22))
data.pfda.14.met_at_22$Age = rep("22", nrow(data.pfda.14.met_at_22))

#### met at 28

data.pfda.14.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_14/minerva_data_pfda_14_metabolites_28/pfda_14_met_28_beta_fisher_c18.txt", header=TRUE)
data.pfda.14.met_at_28$PFAS <-  rep("pfda", nrow(data.pfda.14.met_at_28))
data.pfda.14.met_at_28$Age = rep("28", nrow(data.pfda.14.met_at_28))

pfda_met <- rbind(data.pfda.14.met_at_22, data.pfda.14.met_at_28)

### pfhxs
#### met at 22

data.pfhxs.14.met_at_22 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_22/pfhxs_14_met_22_beta_fisher_c18.txt", header=TRUE)
data.pfhxs.14.met_at_22$PFAS <-  rep("pfhxs", nrow(data.pfhxs.14.met_at_22))
data.pfhxs.14.met_at_22$Age = rep("22", nrow(data.pfhxs.14.met_at_22))

#### met at 28

data.pfhxs.14.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_28/pfhxs_14_met_28_beta_fisher_c18.txt", header=TRUE)
data.pfhxs.14.met_at_28$PFAS <-  rep("pfhxs", nrow(data.pfhxs.14.met_at_28))
data.pfhxs.14.met_at_28$Age = rep("28", nrow(data.pfhxs.14.met_at_28))

pfhxs_met <- rbind(data.pfhxs.14.met_at_22, data.pfhxs.14.met_at_28)

pfas_met_tab_c18 <- rbind(pfhxs_met,pfda_met,pfna_met,pfos_met,pfoa_met)
pfas_met_tab_c18$Mode <- rep("C18", nrow(pfas_met_tab_c18))
write.table(pfas_met_tab_c18, "/sc/arion/projects/Faroese/pfas_met/c18/pfas_14_metabolites_all_c18.txt", row.names = F)

