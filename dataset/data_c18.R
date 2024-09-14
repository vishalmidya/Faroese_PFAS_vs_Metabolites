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
t<-as.data.frame(read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/17.12.19_December2020.csv"))
d <- as.data.frame(read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/17.12.19_December2020.csv"))
d$id <- as.character(d$id)

# eligible id: 500 ids

d1 <- as.data.frame(read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/eligible_ids.csv"))
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


# Read Metabolomics data - c18+

mapping_c18 <- read.delim("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/Perfluoroalkyl_ALL_mapping_c18_rs.txt")

# read median summarized c18+ table
median_sum_feature_c18 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/MetabComBat_c18neg_merged_mediansummarized_featuretable_final.csv")

median_sum_feature_c18 <- median_sum_feature_c18 %>% 
  dplyr::select(-ends_with("_A"), -ends_with("_B"))

# Extract File name ID_YEAR

x <- colnames(median_sum_feature_c18)[3:ncol(median_sum_feature_c18)]
y <- rep(NA_character_, length(x))
for(i in 1:length(x)){
  y[i] <- sum(mapping_c18$File.Name %in% x[i])
}

# median_sum_feature_c18 <- median_sum_feature_c18[,!(colnames(median_sum_feature_c18) %in% x[which(y == "0")])]
x <- colnames(median_sum_feature_c18)[3:ncol(median_sum_feature_c18)]
y <- rep(NA_character_, length(x))
for(i in 1:length(x)){
  y[i] = mapping_c18$Sample.ID[which(mapping_c18$File.Name %in% x[i])] 
}
colnames(median_sum_feature_c18)[3:ncol(median_sum_feature_c18)] <- y


# valid columns - columns with year_id
valid_colnames <- NA_character_
for(i in 1:length(colnames(median_sum_feature_c18))){
  if(is.na(strsplit(colnames(median_sum_feature_c18)[i],"_")[[1]][2])==F & (is.na((strsplit(strsplit(colnames(median_sum_feature_c18)[i],"_")[[1]][2],"")[[1]][3] == "y")) == F & (strsplit(strsplit(colnames(median_sum_feature_c18)[i],"_")[[1]][2],"")[[1]][3] == "y")) ||
     (is.na((strsplit(strsplit(colnames(median_sum_feature_c18)[i],"_")[[1]][2],"")[[1]][2] == "y")) == F & (strsplit(strsplit(colnames(median_sum_feature_c18)[i],"_")[[1]][2],"")[[1]][2] == "y"))){
    valid_colnames <- c(valid_colnames, colnames(median_sum_feature_c18)[i])
  }
}

valid_colnames <- unique(valid_colnames[-1])
length(valid_colnames)/4      # 309  - year_id



# Metabolites data filtering



################################################################################
## 1. confirmed metabolites

### a). import confirmed metabolites, with additional metabolites


conf_met_org <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/confirmed_metabolites.csv", fileEncoding = "Latin1")
conf_met_add <- read_excel("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/confirmed_metabolites_added.xlsx")

conf_met_add<- conf_met_add %>%
  mutate(Method.RT = case_when(Method == "HILIC/ESI+" ~ "HILIC+",
                               Method %in% c("C18/ESI-", "C18/ES-") ~ "C18-",
                               .default = NA))
conf_met_org<- conf_met_org %>%  
  dplyr::select(Metabolite, m.z..adduct., time, Method.RT)
names(conf_met_org)[3]<- "Time"


#### clean additional metabolites: remove unknown method
conf_met_add<- conf_met_add %>%
  filter(Method.RT %in% c("HILIC+", "C18-")) %>%
  dplyr::select(Metabolite, `m/z (Adduct)`, `RT (sec)`, Method.RT)

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

##### clean metabolites: restrict to c18+
conf_met<- conf_met %>% 
  filter(Method.RT == "C18-") %>% 
  dplyr::select(Metabolite, mz_ratio, time, Method.RT)

dim(conf_met)

# [1] 754   4


## 2. median file
##### for each metabolite, non-zero value > 75%
median_sum_feature_c18$seq <- paste0("chem_",seq(1, nrow(median_sum_feature_c18)))
median_sum_feature_c18<- median_sum_feature_c18[,c(valid_colnames,"mz","time","seq")]
median_sum_feature_c18$zero<- rowSums(median_sum_feature_c18[,valid_colnames]==0) 
median_sum_feature_c18$zero_per<- median_sum_feature_c18$zero/length(valid_colnames)*100

median_sum_feature_c18<- median_sum_feature_c18[median_sum_feature_c18$zero_per<25,]
dim(median_sum_feature_c18)
#5067

##### for each metabolite, CV% > 30%
median_sum_feature_c18<- median_sum_feature_c18 %>% 
  rowwise() %>% 
  mutate(
    mean = mean(c_across(valid_colnames[1]:valid_colnames[1235])),
    sd = sd(c_across(valid_colnames[1]:valid_colnames[1235])),
    cv = 100*(sd/mean))
median_sum_feature_c18<- median_sum_feature_c18[median_sum_feature_c18$cv>30,]
dim(median_sum_feature_c18)
#[1] 3224 1243




# Characterize metabolites from median_sum_feature_c18 using confirmed metabolites
median_sum_feature_c18$Metabolite_multiple <- rep(NA_character_, nrow(median_sum_feature_c18))
median_sum_feature_c18$seq <- paste0("chem_",seq(1, nrow(median_sum_feature_c18)))


for(i in 1:nrow(median_sum_feature_c18)){
  
  rrt <- conf_met$Metabolite[which(abs(median_sum_feature_c18$mz[i] - conf_met$mz) < 3)]  ###### criteria
  a <- unique(rrt[!is.na(rrt)])
  b <- unique(conf_met$Metabolite[conf_met$Metabolite %in% a][abs(median_sum_feature_c18$time[i] - conf_met$time[conf_met$Metabolite %in% a]) < 5])
  
  if(length(b)!= 0){
    
    dft <- conf_met[conf_met$Metabolite %in% b,c("Metabolite")]
    median_sum_feature_c18$Metabolite_multiple[i] = knitr::combine_words(unique(dft), and = "", sep = "/")
    
  }
}


length(unique(median_sum_feature_c18$Metabolite_multiple)) 
# [1] 323


length(median_sum_feature_c18[is.na(median_sum_feature_c18$Metabolite_multiple) == FALSE, ]$Metabolite_multiple)
# 652 


# median_sum_feature_c18 %>% dplyr::select(mz, time, Metabolite_multiple)



# remove Metabolite_multiples unable to annotate
data_c18 <- median_sum_feature_c18[is.na(median_sum_feature_c18$Metabolite_multiple) == FALSE,c(valid_colnames,"mz","time","seq","Metabolite_multiple")]


data_c18$Met_id <- paste0("Met",seq(1:nrow(data_c18)))


################################################################################





################################################################################
# restrict to only original confirmed metabolites and annotate with closest confirmed metabolite!!!!!!!!!!!!!!!!!!!!!
conf_met_org$mz_ratio <- rep(NA_real_, nrow(conf_met_org))
conf_met_org$adduct <- rep(NA_character_, nrow(conf_met_org))

for(i in 1:nrow(conf_met_org)){
  conf_met_org$mz_ratio[i] = as.numeric(strsplit(conf_met_org$m.z..adduct.[i]," ")[[1]][1])
  x <- strsplit(strsplit(conf_met_org$m.z..adduct.[i]," ")[[1]][2],"")[[1]]
  conf_met_org$adduct[i] <- as.character(knitr::combine_words(x[x!= "(" & x!= ")"], and = "", sep = ""))
}

for(i in 1:nrow(conf_met_org)){
  conf_met_org$time[i] = as.numeric(strsplit(conf_met_org$Time[i],",")[[1]][1])
}

dim(conf_met_org)

conf_met_org<- conf_met_org %>% 
  filter(is.na(time)==FALSE & is.na(mz_ratio)==FALSE)

dim(conf_met_org)


conf_met<- conf_met_org %>% 
  filter(Method.RT == "C18-") %>% 
  dplyr::select(Metabolite, mz_ratio, time, Method.RT)

dim(conf_met)
# 156    4



data_c18$Metabolite <- rep(NA_character_, nrow(data_c18))
for(i in 1:nrow(data_c18)){
  
  rrt <- conf_met$Metabolite[which(abs(data_c18$mz[i] - conf_met$mz) < 3)]  ###### criteria
  a <- unique(rrt[!is.na(rrt)]) 
  b <- unique(conf_met$Metabolite[conf_met$Metabolite %in% a][abs(data_c18$time[i] - conf_met$time[conf_met$Metabolite %in% a]) < 5])
  c <- conf_met$Metabolite[conf_met$Metabolite %in% b][which.min(abs(data_c18$time[i] - conf_met$time[conf_met$Metabolite %in% b]) + abs(data_c18$mz[i] - conf_met$mz[conf_met$Metabolite %in% b]))]
  
  if(length(c)!= 0){
    
    dft <- conf_met[conf_met$Metabolite %in% c,c("Metabolite")]
    data_c18$Metabolite[i] = knitr::combine_words(unique(dft), and = "", sep = "/")
    
  }
}



length(data_c18[is.na(data_c18$Metabolite) == FALSE, ]$Metabolite)
# 326

# remove metabolites unable to annotate
data_c18 <- data_c18[is.na(data_c18$Metabolite) == FALSE,]

dim(data_c18) 
# [1] 326 1241

# data_c18 %>% dplyr::select(mz, time, Metabolite,Met_id)

keep_metabolites_c18 = data.frame(Met_id = data_c18$Met_id)

write.csv(data_c18, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/data_c18.csv", row.names = F)

write.csv(keep_metabolites_c18, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/keep_metabolites_c18.csv", row.names = F)

################################################################################################################
data_met_clinical_c18 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/data_c18.csv", check.names = F)

# transpose met data
t_data_met_clinical_c18 <- data.table::transpose(data_met_clinical_c18)

# get row and colnames in order
Met_name<- data_met_clinical_c18$Met_id
colnames(t_data_met_clinical_c18) <- Met_name
t_data_met_clinical_c18$id_Year <- colnames(data_met_clinical_c18)
dim(t_data_met_clinical_c18)
# [1] 1240  327


t_data_met_clinical_c18$id <- rep(NA_character_,nrow(t_data_met_clinical_c18))
t_data_met_clinical_c18$Year <- rep(NA_character_,nrow(t_data_met_clinical_c18))

for(i in 1:nrow(t_data_met_clinical_c18)){
  t_data_met_clinical_c18$id[i] <- strsplit(t_data_met_clinical_c18$id_Year[i],"_")[[1]][1]
  t_data_met_clinical_c18$Year[i] <- as.numeric(knitr::combine_words(strsplit(strsplit(t_data_met_clinical_c18$id_Year[i],"_")[[1]][2],"y")[[1]][1],and = "", sep =""))
}


# combine with Epi data
merged_met_clinical_c18 <- merge(d_long,t_data_met_clinical_c18, by = c("id","Year"))
merged_met_clinical_c18$id <- as.character(merged_met_clinical_c18$id)
dim(merged_met_clinical_c18)
# [1] 1235  359

# check available sample size at different timepoints
nrow(merged_met_clinical_c18 %>% filter(Year == "7"))
nrow(merged_met_clinical_c18 %>% filter(Year == "14"))
nrow(merged_met_clinical_c18 %>% filter(Year == "22")) #492
nrow(merged_met_clinical_c18 %>% filter(Year == "28")) #493

## !!!!!!!!!!!!!!!!!!!! Subset for only 125 ids (also keep another one with all eligible ids)
merged_met_clinical_c18 <- merged_met_clinical_c18[which(merged_met_clinical_c18$id %in% as.vector(t_data_met_clinical_c18[t_data_met_clinical_c18$Year == "7","id"])),]
dim(merged_met_clinical_c18)
# [1] 500 359



# Met data standardization: by metabolite and not by year (group_by(Year) %>% )
merged_met_clinical_c18<- merged_met_clinical_c18 %>% 
  mutate_at(Met_name, ~(as.numeric(.) %>% as.vector))  %>% 
  mutate_at(Met_name, ~(log(. + 1, base = 2) %>% as.vector))%>% 
  mutate_at(Met_name,  ~(scale(.) %>% as.vector)) 

# impute PFAS value below LOD
d_subset$pfoa_0<- ifelse(d_subset$pfoa_0 < 0.03, 0.03/sqrt(2), d_subset$pfoa_0)
d_subset$pfos_0<- ifelse(d_subset$pfos_0 < 0.03, 0.03/sqrt(2), d_subset$pfos_0)
d_subset$pfhxs_0<- ifelse(d_subset$pfhxs_0 < 0.05, 0.05/sqrt(2), d_subset$pfhxs_0)
d_subset$pfna_0<- ifelse(d_subset$pfna_0 < 0.05, 0.05/sqrt(2), d_subset$pfna_0)
d_subset$pfda_0<- ifelse(d_subset$pfoa_0 < 0.03, 0.03/sqrt(2), d_subset$pfda_0)

d_subset$pfoa_7<- ifelse(d_subset$pfoa_7 < 0.03, 0.03/sqrt(2), d_subset$pfoa_7)
d_subset$pfos_7<- ifelse(d_subset$pfos_7 < 0.03, 0.03/sqrt(2), d_subset$pfos_7)
d_subset$pfhxs_7<- ifelse(d_subset$pfhxs_7 < 0.03, 0.03/sqrt(2), d_subset$pfhxs_7)
d_subset$pfna_7<- ifelse(d_subset$pfna_7 < 0.03, 0.03/sqrt(2), d_subset$pfna_7)
d_subset$pfda_7<- ifelse(d_subset$pfoa_7 < 0.03, 0.03/sqrt(2), d_subset$pfda_7)

d_subset$pfoa_14<- ifelse(d_subset$pfoa_14 < 0.03, 0.03/sqrt(2), d_subset$pfoa_14)
d_subset$pfos_14<- ifelse(d_subset$pfos_14 < 0.03, 0.03/sqrt(2), d_subset$pfos_14)
d_subset$pfhxs_14<- ifelse(d_subset$pfhxs_14 < 0.03, 0.03/sqrt(2), d_subset$pfhxs_14)
d_subset$pfna_14<- ifelse(d_subset$pfna_14 < 0.03, 0.03/sqrt(2), d_subset$pfna_14)
d_subset$pfda_14<- ifelse(d_subset$pfoa_14 < 0.03, 0.03/sqrt(2), d_subset$pfda_14)

d_subset$pfoa_22<- ifelse(d_subset$pfoa_22 < 0.03, 0.03/sqrt(2), d_subset$pfoa_22)
d_subset$pfos_22<- ifelse(d_subset$pfos_22 < 0.03, 0.03/sqrt(2), d_subset$pfos_22)
d_subset$pfhxs_22<- ifelse(d_subset$pfhxs_22 < 0.03, 0.03/sqrt(2), d_subset$pfhxs_22)
d_subset$pfna_22<- ifelse(d_subset$pfna_22 < 0.03, 0.03/sqrt(2), d_subset$pfna_22)
d_subset$pfda_22<- ifelse(d_subset$pfoa_22 < 0.03, 0.03/sqrt(2), d_subset$pfda_22)

d_subset$pfoa_28<- ifelse(d_subset$pfoa_28 < 0.03, 0.03/sqrt(2), d_subset$pfoa_28)
d_subset$pfos_28<- ifelse(d_subset$pfos_28 < 0.03, 0.03/sqrt(2), d_subset$pfos_28)
d_subset$pfhxs_28<- ifelse(d_subset$pfhxs_28 < 0.03, 0.03/sqrt(2), d_subset$pfhxs_28)
d_subset$pfna_28<- ifelse(d_subset$pfna_28 < 0.03, 0.03/sqrt(2), d_subset$pfna_28)
d_subset$pfda_28<- ifelse(d_subset$pfoa_28 < 0.03, 0.03/sqrt(2), d_subset$pfda_28)




merged_met_clinical_c18 <- merge(merged_met_clinical_c18, d_subset[,c("pfoa_0","pfos_0","pfhxs_0","pfna_0","pfda_0","id",
                                                                      "pfoa_7","pfos_7","pfhxs_7","pfna_7","pfda_7",
                                                                      "pfoa_14","pfos_14","pfhxs_14","pfna_14","pfda_14",
                                                                      "pfoa_22","pfos_22","pfhxs_22","pfna_22","pfda_22",
                                                                      "pfoa_28","pfos_28","pfhxs_28","pfna_28","pfda_28")],by = "id")


# combine with dichotomous PFAS
merged_met_clinical_c18$cpfoa0 <- ifelse(merged_met_clinical_c18$pfoa_0 > median(merged_met_clinical_c18$pfoa_0), 1, 0)
merged_met_clinical_c18$cpfos0 <- ifelse(merged_met_clinical_c18$pfos_0 > median(merged_met_clinical_c18$pfos_0), 1, 0)
merged_met_clinical_c18$cpfhxs0 <- ifelse(merged_met_clinical_c18$pfhxs_0 > median(merged_met_clinical_c18$pfhxs_0), 1, 0)
merged_met_clinical_c18$cpfna0 <- ifelse(merged_met_clinical_c18$pfna_0 > median(merged_met_clinical_c18$pfna_0), 1, 0)
merged_met_clinical_c18$cpfda0 <- ifelse(merged_met_clinical_c18$pfda_0 > median(merged_met_clinical_c18$pfda_0), 1, 0)


merged_met_clinical_c18$cpfoa7 <- ifelse(merged_met_clinical_c18$pfoa_7 > median(merged_met_clinical_c18$pfoa_7), 1, 0)
merged_met_clinical_c18$cpfos7 <- ifelse(merged_met_clinical_c18$pfos_7 > median(merged_met_clinical_c18$pfos_7), 1, 0)
merged_met_clinical_c18$cpfhxs7 <- ifelse(merged_met_clinical_c18$pfhxs_7 > median(merged_met_clinical_c18$pfhxs_7), 1, 0)
merged_met_clinical_c18$cpfna7 <- ifelse(merged_met_clinical_c18$pfna_7 > median(merged_met_clinical_c18$pfna_7), 1, 0)
merged_met_clinical_c18$cpfda7 <- ifelse(merged_met_clinical_c18$pfda_7 > median(merged_met_clinical_c18$pfda_7), 1, 0)

merged_met_clinical_c18$cpfoa14 <- ifelse(merged_met_clinical_c18$pfoa_14 > median(merged_met_clinical_c18$pfoa_14), 1, 0)
merged_met_clinical_c18$cpfos14 <- ifelse(merged_met_clinical_c18$pfos_14 > median(merged_met_clinical_c18$pfos_14), 1, 0)
merged_met_clinical_c18$cpfhxs14 <- ifelse(merged_met_clinical_c18$pfhxs_14 > median(merged_met_clinical_c18$pfhxs_14), 1, 0)
merged_met_clinical_c18$cpfna14 <- ifelse(merged_met_clinical_c18$pfna_14 > median(merged_met_clinical_c18$pfna_14), 1, 0)
merged_met_clinical_c18$cpfda14 <- ifelse(merged_met_clinical_c18$pfda_14 > median(merged_met_clinical_c18$pfda_14), 1, 0)

merged_met_clinical_c18$cpfoa22 <- ifelse(merged_met_clinical_c18$pfoa_22 > median(merged_met_clinical_c18$pfoa_22), 1, 0)
merged_met_clinical_c18$cpfos22 <- ifelse(merged_met_clinical_c18$pfos_22 > median(merged_met_clinical_c18$pfos_22), 1, 0)
merged_met_clinical_c18$cpfhxs22 <- ifelse(merged_met_clinical_c18$pfhxs_22 > median(merged_met_clinical_c18$pfhxs_22), 1, 0)
merged_met_clinical_c18$cpfna22 <- ifelse(merged_met_clinical_c18$pfna_22 > median(merged_met_clinical_c18$pfna_22), 1, 0)
merged_met_clinical_c18$cpfda22 <- ifelse(merged_met_clinical_c18$pfda_22 > median(merged_met_clinical_c18$pfda_22), 1, 0)

merged_met_clinical_c18$cpfoa28 <- ifelse(merged_met_clinical_c18$pfoa_28 > median(merged_met_clinical_c18$pfoa_28), 1, 0)
merged_met_clinical_c18$cpfos28 <- ifelse(merged_met_clinical_c18$pfos_28 > median(merged_met_clinical_c18$pfos_28), 1, 0)
merged_met_clinical_c18$cpfhxs28 <- ifelse(merged_met_clinical_c18$pfhxs_28 > median(merged_met_clinical_c18$pfhxs_28), 1, 0)
merged_met_clinical_c18$cpfna28 <- ifelse(merged_met_clinical_c18$pfna_28 > median(merged_met_clinical_c18$pfna_28), 1, 0)
merged_met_clinical_c18$cpfda28 <- ifelse(merged_met_clinical_c18$pfda_28 > median(merged_met_clinical_c18$pfda_28), 1, 0)


# reformate variables
merged_met_clinical_c18$Year <- as.character(merged_met_clinical_c18$Year)
merged_met_clinical_c18$Year <- factor(merged_met_clinical_c18$Year, levels  = c("7","14","22","28"))

merged_met_clinical_c18$cmbmi <- ifelse(merged_met_clinical_c18$mbmi >= 25, "High","Low")
merged_met_clinical_c18$cmage <- ifelse(merged_met_clinical_c18$mage >= as.numeric(quantile(merged_met_clinical_c18$mage)[4]), "High","Low")
merged_met_clinical_c18$cparity <- ifelse(merged_met_clinical_c18$parity >= 1, "1","0")
merged_met_clinical_c18$csmokepreg <- ifelse(merged_met_clinical_c18$smokepreg_2 == 1, "1","0")
merged_met_clinical_c18$cbreastfed_tot <- ifelse(merged_met_clinical_c18$breastfed_tot >= as.numeric(quantile(merged_met_clinical_c18$breastfed_tot)[2]), "High","Low")
merged_met_clinical_c18$cmatfishpreg <- as.numeric(ifelse(merged_met_clinical_c18$matfishpreg_cat2 >= 0.5, "1","0"))


merged_met_clinical_c18$age7 <- merged_met_clinical_c18$age_7
merged_met_clinical_c18$age14 <- merged_met_clinical_c18$age_14

dim(merged_met_clinical_c18[merged_met_clinical_c18$Year == 7,(colnames(merged_met_clinical_c18) %in% Met_name)])
dim(merged_met_clinical_c18[merged_met_clinical_c18$Year == 14,(colnames(merged_met_clinical_c18) %in% Met_name)])
dim(merged_met_clinical_c18[merged_met_clinical_c18$Year == 22,(colnames(merged_met_clinical_c18) %in% Met_name)])
dim(merged_met_clinical_c18[merged_met_clinical_c18$Year == 28,(colnames(merged_met_clinical_c18) %in% Met_name)])
#[1]  125 326

write.csv(merged_met_clinical_c18, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/merged_omics_c18.csv", row.names = F)

