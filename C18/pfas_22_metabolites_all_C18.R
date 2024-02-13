library(foreach)
library(doParallel)
library(car)
library(MASS)
library(ggplot2)
library(data.table)
library(iterators)
library(parallel)
library(carData)


# PFAS exposure at age 22
## C18 +ve

### pfoa
#### met at 28

data.pfoa.22.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_c18.txt", header=TRUE)
data.pfoa.22.met_at_28$PFAS <-  rep("pfoa", nrow(data.pfoa.22.met_at_28))
data.pfoa.22.met_at_28$Age = rep("28", nrow(data.pfoa.22.met_at_28))
pfoa_met <- rbind( data.pfoa.22.met_at_28)


### pfos
#### met at 28

data.pfos.22.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_22/minerva_data_pfos_22_metabolites_28/pfos_22_met_28_beta_fisher_c18.txt", header=TRUE)
data.pfos.22.met_at_28$PFAS <-  rep("pfos", nrow(data.pfos.22.met_at_28))
data.pfos.22.met_at_28$Age = rep("28", nrow(data.pfos.22.met_at_28))
pfos_met <- rbind( data.pfos.22.met_at_28)

### pfda
#### met at 28

data.pfda.22.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/pfda_22_met_28_beta_fisher_c18.txt", header=TRUE)
data.pfda.22.met_at_28$PFAS <-  rep("pfda", nrow(data.pfda.22.met_at_28))
data.pfda.22.met_at_28$Age = rep("28", nrow(data.pfda.22.met_at_28))
pfda_met <- rbind( data.pfda.22.met_at_28)

### pfna
#### met at 28

data.pfna.22.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_22/minerva_data_pfna_22_metabolites_28/pfna_22_met_28_beta_fisher_c18.txt", header=TRUE)
data.pfna.22.met_at_28$PFAS <-  rep("pfna", nrow(data.pfna.22.met_at_28))
data.pfna.22.met_at_28$Age = rep("28", nrow(data.pfna.22.met_at_28))
pfna_met <- rbind( data.pfna.22.met_at_28)

### pfhxs
#### met at 28

data.pfhxs.22.met_at_28 <- read.table("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_22/minerva_data_pfhxs_22_metabolites_28/pfhxs_22_met_28_beta_fisher_c18.txt", header=TRUE)
data.pfhxs.22.met_at_28$PFAS <-  rep("pfhxs", nrow(data.pfhxs.22.met_at_28))
data.pfhxs.22.met_at_28$Age = rep("28", nrow(data.pfhxs.22.met_at_28))
pfhxs_met <- rbind( data.pfhxs.22.met_at_28)



pfas_met_tab_c18 <- rbind(pfhxs_met,pfda_met,pfna_met,pfos_met,pfoa_met)
pfas_met_tab_c18$Mode <- rep("C18", nrow(pfas_met_tab_c18))
write.table(pfas_met_tab_c18, "/sc/arion/projects/Faroese/pfas_met/c18/pfas_22_metabolites_all_c18.txt", row.names = F)



