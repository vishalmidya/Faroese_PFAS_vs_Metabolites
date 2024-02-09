library(car)
library(readr)
library(readxl)
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
library(qwraps2)
library(knitr)
library(xtable)
library(kableExtra)
library(DT)
library(glmnet)
library(corrplot)
library(ggpubr)
library("EnvStats")
library(lmerTest)
library("merTools")
library(reshape2)
library(ggplot2)
library(GGally)
library(mgcv)
library(gplots)
library(tidyr)
library(cluster)
library(factoextra)
library(psych)
library(xgboost)
library(SHAPforxgboost)
library(bkmr)
library(mi)
library(stargazer)
library(factoextra) 
library(spatstat)
library(Hmisc)
library(robustHD)
library(gWQS)
library(gridExtra)
library(ggcorrplot)
library(BWQS)
library(qwraps2)
library(mgcv)
library(ggcorrplot)
library("MatchIt")
library(cobalt)
library(WeightIt)
library(ggalluvial)
library(mice)
library(data.table)
library(omu)
library(renv)

# Set working directory in local repository

# DDD

# setwd("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/git/Faroese_PFAS_vs_Metabolites/HILIC")


# Always start with renv to capture generate workflow

# renv::init()

# Epi data
t<-as.data.frame(read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/17.12.19_December2020.csv"))
d <- as.data.frame(read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/17.12.19_December2020.csv"))
d$id <- as.character(d$id)

# eligible id: 500 ids

d1 <- as.data.frame(read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/eligible_ids.csv"))
d1$id <- as.character(as.numeric(as.character(d1$id)))

# Match the eligible ids with epi data

y <- rep(NA_real_, length(d1$id))
for(i in 1:length(d1$id)){
  y[i] = sum(d$id %in% d1$id[i])
}

dt <- d[which(d$id %in% d1$id),]

# Select the PFAS exposure columns and covariates

d <- dt[,c("id",'pfos_0', 'pfos_7', 'pfos_14', 'pfos_22', 'pfos_28', 'pfoa_0', 'pfoa_7', 'pfoa_14', 'pfoa_22', 'pfoa_28', 
           'pfhxs_0', 'pfhxs_7', 'pfhxs_14', 'pfhxs_22', 'pfhxs_28',  'pfna_0', 'pfna_7', 'pfna_14', 'pfna_22', 'pfna_28', 
           'pfda_0', 'pfda_7', 'pfda_14', 'pfda_22', 'pfda_28', 'bmi_28', 'waistcirc_28', 'diaBP_dxt_28', 'sysBP_dxt_28', 
           'uhdl_28y', 'trig_28y', 'dldl_28y', 'chol_28y','cir_28', 'insulinauc_28y', 'glucoseauc_28y', 'matsuda_28', 
           'homair_28', 'igi_28', "sex","mage","mbmi","parity","smokepreg_2","matfishpreg_cat2","breastfed_tot",
           "age_7","age_14","age22","age28")]

d$id <- as.character(d$id)

# 493 ids

# Impute Data using MICE
sapply(d, function (x) sum(is.na(x)))


init = mice(d, maxit=0) 
meth = init$method
predM = init$predictorMatrix
set.seed(1234)
imputed = mice(d, method="norm", predictorMatrix=predM, m=3, ntree = 100)
imputed <- complete(imputed, action = 2) # this is your final imputed dataset


d_subset <- imputed[,c("id",'pfos_0', 'pfos_7', 'pfos_14', 'pfos_22', 'pfos_28', 'pfoa_0', 'pfoa_7', 'pfoa_14', 'pfoa_22', 
                       'pfoa_28', 'pfhxs_0', 'pfhxs_7', 'pfhxs_14', 'pfhxs_22', 'pfhxs_28',  'pfna_0', 'pfna_7', 'pfna_14', 
                       'pfna_22', 'pfna_28', 'pfda_0', 'pfda_7', 'pfda_14', 'pfda_22', 'pfda_28', 'bmi_28', 'waistcirc_28', 
                       'diaBP_dxt_28', 'sysBP_dxt_28', 'uhdl_28y', 'trig_28y', 'dldl_28y', 'chol_28y','cir_28', 'insulinauc_28y', 
                       'glucoseauc_28y', 'matsuda_28', 'homair_28', 'igi_28', "sex","mage","mbmi","parity","smokepreg_2",
                       "matfishpreg_cat2","breastfed_tot",
                       "age_7","age_14","age22","age28")]

d_long <- reshape(d_subset, idvar = "id", direction = "long", timevar = "Year", 
                  varying = list(c('pfos_0', 'pfos_7', 'pfos_14', 'pfos_22', 'pfos_28'), 
                                 c('pfoa_0', 'pfoa_7', 'pfoa_14', 'pfoa_22', 'pfoa_28'),
                                 c('pfhxs_0', 'pfhxs_7', 'pfhxs_14', 'pfhxs_22', 'pfhxs_28'),
                                 c('pfna_0', 'pfna_7', 'pfna_14', 'pfna_22', 'pfna_28'),
                                 c('pfda_0', 'pfda_7', 'pfda_14', 'pfda_22', 'pfda_28')),
                  v.names = c("PFOS","PFOA","PFHxS","PFNA","PFDA"))

d_long$Year <- (d_long$Year - 1)*7
d_long$Year <- ifelse(d_long$Year == 21, 22, d_long$Year)
d_long$Year <- factor(d_long$Year, levels = c("0","7","14","22","28"))


# Read Metabolomics data - HILIC+

mapping_hilic <- read.delim("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/Perfluoroalkyl_ALL_mapping_hilicpos.txt")

# read median summarized HILIC+ table
median_sum_feature_hilic <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/MetabComBat_HILICpos_Merged_mediansummarized_featuretable_final.csv")

#7109

### remove VT_201114_M456_133_A, VT_201114_M456_139_A, VT_201114_M456_133_B, VT_201114_M456_139_B (Not patient sample - duplicated)

median_sum_feature_hilic <- median_sum_feature_hilic[,!(colnames(median_sum_feature_hilic) %in% c("VT_201114_M456_133_A","VT_201114_M456_139_A",
                                                                                                  "VT_201114_M456_133_B","VT_201114_M456_139_B"))]

# Extract File name ID_YEAR
x <- colnames(median_sum_feature_hilic)[3:ncol(median_sum_feature_hilic)]
y <- rep(NA_character_, length(x))
for(i in 1:length(x)){
  y[i] <- sum(mapping_hilic$File.Name %in% x[i])
}

# median_sum_feature_hilic <- median_sum_feature_hilic[,!(colnames(median_sum_feature_hilic) %in% x[which(y == "0")])]
x <- colnames(median_sum_feature_hilic)[3:ncol(median_sum_feature_hilic)]
y <- rep(NA_character_, length(x))
for(i in 1:length(x)){
  y[i] = mapping_hilic$Sample.ID[which(mapping_hilic$File.Name %in% x[i])] 
}
colnames(median_sum_feature_hilic)[3:ncol(median_sum_feature_hilic)] <- y


# valid columns - columns with year_id
valid_colnames <- NA_character_
for(i in 1:length(colnames(median_sum_feature_hilic))){
  if(is.na(strsplit(colnames(median_sum_feature_hilic)[i],"_")[[1]][2])==F & (is.na((strsplit(strsplit(colnames(median_sum_feature_hilic)[i],"_")[[1]][2],"")[[1]][3] == "y")) == F & (strsplit(strsplit(colnames(median_sum_feature_hilic)[i],"_")[[1]][2],"")[[1]][3] == "y")) ||
     (is.na((strsplit(strsplit(colnames(median_sum_feature_hilic)[i],"_")[[1]][2],"")[[1]][2] == "y")) == F & (strsplit(strsplit(colnames(median_sum_feature_hilic)[i],"_")[[1]][2],"")[[1]][2] == "y"))){
    valid_colnames <- c(valid_colnames, colnames(median_sum_feature_hilic)[i])
  }
}

valid_colnames <- unique(valid_colnames[-1])
length(valid_colnames)/4      # 309  - year_id



# Metabolites data filtering

## 1. confirmed metabolites

### a). import confirmed metabolites, with additional metabolites

conf_met_org <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/confirmed_metabolites.csv", fileEncoding = "Latin1")
conf_met_add <- read_excel("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/confirmed_metabolites_added.xlsx")

conf_met_add<- conf_met_add %>% 
               mutate(Method.RT = case_when(Method == "HILIC/ESI+" ~ "HILIC+",
                                            Method %in% c("C18/ESI-", "C18/ES-") ~ "C18-",
                                            .default = NA))
conf_met_org<- conf_met_org %>%  
               select(Metabolite, m.z..adduct., time, Method.RT)
names(conf_met_org)[3]<- "Time"


##### clean additional metabolites: remove unknown method
conf_met_add<- conf_met_add %>%  
  filter(Method.RT %in% c("HILIC+", "C18-")) %>% 
  select(Metabolite, `m/z (Adduct)`, `RT (sec)`, Method.RT)

names(conf_met_add)[2]<- "m.z..adduct."
names(conf_met_add)[3]<- "Time"
names(conf_met_add)[4]<- "Method.RT"

##### combine metabolites together
conf_met<- rbind(conf_met_org, conf_met_add)

##### clean metabolites: remove unknown m/z or time
conf_met$mz_ratio <- rep(NA_real_, nrow(conf_met))
conf_met$adduct <- rep(NA_character_, nrow(conf_met))

for(i in 1:nrow(conf_met)){
  conf_met$mz_ratio[i] = as.numeric(strsplit(conf_met$m.z..adduct.[i]," ")[[1]][1])
  x <- strsplit(strsplit(conf_met$m.z..adduct.[i]," ")[[1]][2],"")[[1]]
  conf_met$adduct[i] <- as.character(knitr::combine_words(x[x!= "(" & x!= ")"], and = "", sep = ""))
}

for(i in 1:nrow(conf_met)){
  conf_met$time[i] = as.numeric(strsplit(conf_met$Time[i],",")[[1]][1])
}

dim(conf_met)

conf_met<- conf_met %>% 
            filter(is.na(time)==FALSE & is.na(mz_ratio)==FALSE)

dim(conf_met)

##### clean metabolites: restrict to HILIC+
conf_met<- conf_met %>% 
           filter(Method.RT == "HILIC+") %>% 
           select(Metabolite, mz_ratio, time, Method.RT)

dim(conf_met)

# [1] 754   4


## 2. median file
##### for each metabolite, non-zero value > 75%
median_sum_feature_hilic$seq <- paste0("chem_",seq(1, nrow(median_sum_feature_hilic)))
median_sum_feature_hilic<- median_sum_feature_hilic[,c(valid_colnames,"mz","time","seq")]
median_sum_feature_hilic$zero<- rowSums(median_sum_feature_hilic[,valid_colnames]==0) 
median_sum_feature_hilic$zero_per<- median_sum_feature_hilic$zero/length(valid_colnames)*100

median_sum_feature_hilic<- median_sum_feature_hilic[median_sum_feature_hilic$zero_per<25,]
dim(median_sum_feature_hilic)
#5067

##### for each metabolite, CV% > 30%
median_sum_feature_hilic<- median_sum_feature_hilic %>% 
  rowwise() %>% 
  mutate(
    mean = mean(c_across(valid_colnames[1]:valid_colnames[1236])),
    sd = sd(c_across(valid_colnames[1]:valid_colnames[1236])),
    cv = 100*(sd/mean))
median_sum_feature_hilic<- median_sum_feature_hilic[median_sum_feature_hilic$cv>30,]
dim(median_sum_feature_hilic)
#4205




# Characterize metabolites from median_sum_feature_hilic using confirmed metabolites
median_sum_feature_hilic$Metabolite <- rep(NA_character_, nrow(median_sum_feature_hilic))
median_sum_feature_hilic$seq <- paste0("chem_",seq(1, nrow(median_sum_feature_hilic)))


for(i in 1:nrow(median_sum_feature_hilic)){
  
  rrt <- conf_met$Metabolite[which(abs(median_sum_feature_hilic$mz[i] - conf_met$mz) < 3)]  ###### criteria
  a <- unique(rrt[!is.na(rrt)])
  b <- unique(conf_met$Metabolite[conf_met$Metabolite %in% a][abs(median_sum_feature_hilic$time[i] - conf_met$time[conf_met$Metabolite %in% a]) < 5])

  if(length(b)!= 0){
    
    dft <- conf_met[conf_met$Metabolite %in% b,c("Metabolite")]
    median_sum_feature_hilic$Metabolite[i] = knitr::combine_words(unique(dft), and = "", sep = "/")
    
  }
}


length(unique(median_sum_feature_hilic$Metabolite)) 
# [1] 675


length(median_sum_feature_hilic[is.na(median_sum_feature_hilic$Metabolite) == FALSE, ]$Metabolite)
# 1107


# median_sum_feature_hilic %>% select(mz, time, Metabolite)



# export data
data_hilic <- median_sum_feature_hilic[is.na(median_sum_feature_hilic$Metabolite) == FALSE,c(valid_colnames,"mz","time","seq")]
data_hilic$Met_id <- paste0("Met",seq(1:nrow(data_hilic)))
dim(data_hilic) 
# [1] 1107 1240


write.csv(data_hilic, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/data_hilic.csv", row.names = F)


###################################################################################################################

data_met_clinical_hilic <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/data_hilic.csv", check.names = F)

# transpose met data
t_data_met_clinical_hilic <- data.table::transpose(data_met_clinical_hilic)

# get row and colnames in order
Met_name<- data_met_clinical_hilic$Met_id
colnames(t_data_met_clinical_hilic) <- Met_name
t_data_met_clinical_hilic$id_Year <- colnames(data_met_clinical_hilic)
dim(t_data_met_clinical_hilic)
# [1] 1240 1108


t_data_met_clinical_hilic$id <- rep(NA_character_,nrow(t_data_met_clinical_hilic))
t_data_met_clinical_hilic$Year <- rep(NA_character_,nrow(t_data_met_clinical_hilic))

for(i in 1:nrow(t_data_met_clinical_hilic)){
  t_data_met_clinical_hilic$id[i] <- strsplit(t_data_met_clinical_hilic$id_Year[i],"_")[[1]][1]
  t_data_met_clinical_hilic$Year[i] <- as.numeric(knitr::combine_words(strsplit(strsplit(t_data_met_clinical_hilic$id_Year[i],"_")[[1]][2],"y")[[1]][1],and = "", sep =""))
}


# combine with Epi data
merged_met_clinical_hilic <- merge(d_long,t_data_met_clinical_hilic, by = c("id","Year"))
merged_met_clinical_hilic$id <- as.character(merged_met_clinical_hilic$id)
dim(merged_met_clinical_hilic)
# [1] 1236 1140

# check available sample size at different timepoints
nrow(merged_met_clinical_hilic %>% filter(Year == "7"))
nrow(merged_met_clinical_hilic %>% filter(Year == "14"))
nrow(merged_met_clinical_hilic %>% filter(Year == "22")) #493
nrow(merged_met_clinical_hilic %>% filter(Year == "28")) #493

## Subset for only 125 ids
merged_met_clinical_hilic <- merged_met_clinical_hilic[which(merged_met_clinical_hilic$id %in% as.vector(t_data_met_clinical_hilic[t_data_met_clinical_hilic$Year == "7","id"])),]
dim(merged_met_clinical_hilic)
# [1] 500 1140



# Met data standardization: by metabolite and by year
merged_met_clinical_hilic<- merged_met_clinical_hilic %>% 
  group_by(Year) %>% 
  mutate_at(Met_name, ~(as.numeric(.) %>% as.vector))  %>% 
  mutate_at(Met_name, ~(log(. + 1, base = 2) %>% as.vector))%>% 
  mutate_at(Met_name,  ~(scale(.) %>% as.vector)) %>% 
  ungroup()


# combine with dichotomous PFAS
d_subset$cpfoa0 <- ifelse(d_subset$pfoa_0 > median(d_subset$pfoa_0), 1, 0)
d_subset$cpfos0 <- ifelse(d_subset$pfos_0 > median(d_subset$pfos_0), 1, 0)
d_subset$cpfhxs0 <- ifelse(d_subset$pfhxs_0 > median(d_subset$pfhxs_0), 1, 0)
d_subset$cpfna0 <- ifelse(d_subset$pfna_0 > median(d_subset$pfna_0), 1, 0)
d_subset$cpfda0 <- ifelse(d_subset$pfda_0 > median(d_subset$pfda_0), 1, 0)


d_subset$cpfoa7 <- ifelse(d_subset$pfoa_7 > median(d_subset$pfoa_7), 1, 0)
d_subset$cpfos7 <- ifelse(d_subset$pfos_7 > median(d_subset$pfos_7), 1, 0)
d_subset$cpfhxs7 <- ifelse(d_subset$pfhxs_7 > median(d_subset$pfhxs_7), 1, 0)
d_subset$cpfna7 <- ifelse(d_subset$pfna_7 > median(d_subset$pfna_7), 1, 0)
d_subset$cpfda7 <- ifelse(d_subset$pfda_7 > median(d_subset$pfda_7), 1, 0)

d_subset$cpfoa14 <- ifelse(d_subset$pfoa_14 > median(d_subset$pfoa_14), 1, 0)
d_subset$cpfos14 <- ifelse(d_subset$pfos_14 > median(d_subset$pfos_14), 1, 0)
d_subset$cpfhxs14 <- ifelse(d_subset$pfhxs_14 > median(d_subset$pfhxs_14), 1, 0)
d_subset$cpfna14 <- ifelse(d_subset$pfna_14 > median(d_subset$pfna_14), 1, 0)
d_subset$cpfda14 <- ifelse(d_subset$pfda_14 > median(d_subset$pfda_14), 1, 0)

d_subset$cpfoa22 <- ifelse(d_subset$pfoa_22 > median(d_subset$pfoa_22), 1, 0)
d_subset$cpfos22 <- ifelse(d_subset$pfos_22 > median(d_subset$pfos_22), 1, 0)
d_subset$cpfhxs22 <- ifelse(d_subset$pfhxs_22 > median(d_subset$pfhxs_22), 1, 0)
d_subset$cpfna22 <- ifelse(d_subset$pfna_22 > median(d_subset$pfna_22), 1, 0)
d_subset$cpfda22 <- ifelse(d_subset$pfda_22 > median(d_subset$pfda_22), 1, 0)

d_subset$cpfoa28 <- ifelse(d_subset$pfoa_28 > median(d_subset$pfoa_28), 1, 0)
d_subset$cpfos28 <- ifelse(d_subset$pfos_28 > median(d_subset$pfos_28), 1, 0)
d_subset$cpfhxs28 <- ifelse(d_subset$pfhxs_28 > median(d_subset$pfhxs_28), 1, 0)
d_subset$cpfna28 <- ifelse(d_subset$pfna_28 > median(d_subset$pfna_28), 1, 0)
d_subset$cpfda28 <- ifelse(d_subset$pfda_28 > median(d_subset$pfda_28), 1, 0)


merged_met_clinical_hilic <- merge(merged_met_clinical_hilic, d_subset[,c("pfoa_0","pfos_0","pfhxs_0","pfna_0","pfda_0","id",
                                                                          "cpfoa0","cpfos0","cpfhxs0","cpfna0","cpfda0",
                                                                          "cpfoa7","cpfos7","cpfhxs7","cpfna7","cpfda7",
                                                                          "cpfoa14","cpfos14","cpfhxs14","cpfna14","cpfda14",
                                                                          "cpfoa22","cpfos22","cpfhxs22","cpfna22","cpfda22",
                                                                          "cpfoa28","cpfos28","cpfhxs28","cpfna28","cpfda28")],by = "id")



# reformate variables
merged_met_clinical_hilic$Year <- as.character(merged_met_clinical_hilic$Year)
merged_met_clinical_hilic$Year <- factor(merged_met_clinical_hilic$Year, levels  = c("7","14","22","28"))

merged_met_clinical_hilic$cmbmi <- ifelse(merged_met_clinical_hilic$mbmi >= 25, "High","Low")
merged_met_clinical_hilic$cmage <- ifelse(merged_met_clinical_hilic$mage >= as.numeric(quantile(merged_met_clinical_hilic$mage)[4]), "High","Low")
merged_met_clinical_hilic$cparity <- ifelse(merged_met_clinical_hilic$parity >= 1, "1","0")
merged_met_clinical_hilic$csmokepreg <- ifelse(merged_met_clinical_hilic$smokepreg_2 == 1, "1","0")
merged_met_clinical_hilic$cbreastfed_tot <- ifelse(merged_met_clinical_hilic$breastfed_tot >= as.numeric(quantile(merged_met_clinical_hilic$breastfed_tot)[2]), "High","Low")
merged_met_clinical_hilic$cmatfishpreg <- as.numeric(ifelse(merged_met_clinical_hilic$matfishpreg_cat2 >= 0.5, "1","0"))


merged_met_clinical_hilic$age7 <- merged_met_clinical_hilic$age_7
merged_met_clinical_hilic$age14 <- merged_met_clinical_hilic$age_14


dim(merged_met_clinical_hilic[merged_met_clinical_hilic$Year == 7,(colnames(merged_met_clinical_hilic) %in% Met_name)])
dim(merged_met_clinical_hilic[merged_met_clinical_hilic$Year == 14,(colnames(merged_met_clinical_hilic) %in% Met_name)])
dim(merged_met_clinical_hilic[merged_met_clinical_hilic$Year == 22,(colnames(merged_met_clinical_hilic) %in% Met_name)])
dim(merged_met_clinical_hilic[merged_met_clinical_hilic$Year == 28,(colnames(merged_met_clinical_hilic) %in% Met_name)])
#[1]  125 1107

write.csv(merged_met_clinical_hilic, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/merged_omics_hilic.csv", row.names = F)

