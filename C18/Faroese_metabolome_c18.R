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
library(DescTools)

# Set working directory in local repository

# DDD

# setwd("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC")


# Always start with renv to capture generate workflow

# renv::init()

# Epi data

d <- as.data.frame(read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/17.12.19_December2020.csv"))
d$id <- as.character(d$id)

# 500 ids that were eligible 

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



# Impute Data using MICE

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



### remove VT_201114_M456_133_A, VT_201114_M456_139_A, VT_201114_M456_133_B, VT_201114_M456_139_B (Not patient sample - duplicated)

median_sum_feature_c18 <- median_sum_feature_c18[,!(colnames(median_sum_feature_c18) %like any% c("%A", "%B"))]

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


## confirmed metabolites

conf_met <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/confirmed_metabolites.csv")
conf_met <- conf_met[conf_met$Method.RT == "C18-",]
conf_met$mz_ratio <- rep(NA_real_, nrow(conf_met))
conf_met$adduct <- rep(NA_character_, nrow(conf_met))

for(i in 1:nrow(conf_met)){
  conf_met$mz_ratio[i] = as.numeric(strsplit(conf_met$m.z..adduct.[i]," ")[[1]][1])
  x <- strsplit(strsplit(conf_met$m.z..adduct.[i]," ")[[1]][2],"")[[1]]
  conf_met$adduct[i] <- as.character(knitr::combine_words(x[x!= "(" & x!= ")"], and = "", sep = ""))
}
conf_met <- conf_met[is.na(conf_met$mz_ratio)==F,]

## match the confirmed metabolites with stage4 and keep the ones which match wrt Elemental.Composition - chem formula

stage4_c18 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/Stage4.csv")
common_chems_c18 <- intersect(unique(stage4_c18$Formula), conf_met$Elemental.Composition)

ft <- conf_met[conf_met$Elemental.Composition %in% common_chems_c18,c("KEGG.ID","HMDB.ID","Elemental.Composition")]
colnames(ft)[3] <- c("Formula")

## restrict stage4 to the common metabolites only
### restrict by chemical formula

stage4_c18 <- stage4_c18[which(stage4_c18$Formula %in% common_chems_c18),]
stage4_c18 <- merge(stage4_c18, ft, by = "Formula") 

### Match by adduct and common metabolites

rnames <- NA_character_
for(i in 1:nrow(stage4_c18)){
  if(length(intersect(stage4_c18$Adduct[i], unique(conf_met$adduct[conf_met$Elemental.Composition %in% stage4_c18$Formula[i]]))) != 0){
    rnames <- c(rnames,as.character(i)) 
  }
}

stage4_c18 <- stage4_c18[rownames(stage4_c18) %in% rnames,]
stage4_c18$KEGG.ID <- ifelse(stage4_c18$KEGG.ID == "-",NA_character_,stage4_c18$KEGG.ID)
stage4_c18$HMDB.ID <- ifelse(stage4_c18$HMDB.ID == "-",NA_character_,stage4_c18$HMDB.ID)

# dim(stage4_c18)
# head(stage4_c18)
# unique(stage4_c18$KEGG.ID)

# Characterize metabolites from median_sum_feature_c18 using confirmed metabolites

median_sum_feature_c18$KEGG <- rep(NA_character_, nrow(median_sum_feature_c18))
median_sum_feature_c18$seq <- paste0("chem_",seq(1, nrow(median_sum_feature_c18)))


for(i in 1:nrow(median_sum_feature_c18)){
  
  rrt <- stage4_c18$KEGG.ID[which(abs(median_sum_feature_c18$mz[i] - stage4_c18$mz) < 5)]
  a <- unique(rrt[!is.na(rrt)])
  b <- unique(stage4_c18$KEGG.ID[stage4_c18$KEGG.ID %in% a][abs(median_sum_feature_c18$time[i] - stage4_c18$time[stage4_c18$KEGG.ID %in% a]) < 30])
  
  if(length(b)!= 0){
    
    score <- max(stage4_c18$Annotation.confidence.score[stage4_c18$KEGG.ID %in% b])
    dft <- stage4_c18[stage4_c18$KEGG.ID %in% b,c("KEGG.ID","Annotation.confidence.score")]
    median_sum_feature_c18$KEGG[i] = knitr::combine_words(unique(dft$KEGG.ID[dft$Annotation.confidence.score == score]), and = "", sep = "/")
    
  }
  
}


length(unique(median_sum_feature_c18$KEGG)) #61


# Data standardization

median_sum_feature_c18 <- median_sum_feature_c18[!is.na(median_sum_feature_c18$KEGG),] # get rid of rows with missing KEGG ID
median_sum_feature_c18[,valid_colnames] <- log(median_sum_feature_c18[,valid_colnames]+1, base = 2)
median_sum_feature_c18[,valid_colnames] <- scale(median_sum_feature_c18[,valid_colnames])

rownames(median_sum_feature_c18) <- as.character(seq(1,nrow(median_sum_feature_c18)))
dim(median_sum_feature_c18)

# 787 1445



## Get chemical names and confidence of matching
 
median_sum_feature_c18$Annotation.confidence.score <- rep(NA_real_, nrow(median_sum_feature_c18))
median_sum_feature_c18$chem_name <- rep(NA_character_, nrow(median_sum_feature_c18))

for(i in 1:nrow(median_sum_feature_c18)){
  ffg <- unique(strsplit(median_sum_feature_c18$KEGG[i],"/")[[1]])
  median_sum_feature_c18$Annotation.confidence.score[i] <- paste(unique(stage4_c18$Annotation.confidence.score[stage4_c18$KEGG.ID %in% ffg]), collapse = ";")
  median_sum_feature_c18$chem_name[i] <- paste(unique(stage4_c18$Name[stage4_c18$KEGG.ID %in% ffg]), collapse = ";")
  
}


# subset to unique KEGG IDs 

# u_chems <- unique(median_sum_feature_c18$KEGG)
# r_names <- NA_character_
# for(i in 1:length(u_chems)){
#   r_names <- c(r_names,names(which.max(apply(median_sum_feature_c18[median_sum_feature_c18$KEGG %in% u_chems[i],valid_colnames], 1, 
#                                              function(x){sd(x, na.rm = T)/mean(x, na.rm = T)}))))
# }
# r_names <- r_names[-1]
# dim(median_sum_feature_c18[r_names,])
# median_sum_feature_c18 <- median_sum_feature_c18[r_names,]

data_c18 <- median_sum_feature_c18[is.na(median_sum_feature_c18$KEGG)==F,c(valid_colnames,"mz","time","KEGG","seq","Annotation.confidence.score","chem_name")]
dim(data_c18) 
#  787 1241

# head(data_c18)
View(data_c18)


### length of unique KEGG

# length(unique(median_sum_feature_c18$KEGG)) #109
# unique(median_sum_feature_c18$KEGG) 


#list of adduct
unique(stage4_c18$Adduct)

# "M-H"  "M+Cl"

toy <- median_sum_feature_c18
joy <- toy[1,]
for(i in 1: nrow(toy)){

  for(j in 1:length(strsplit(toy[i,c("KEGG")],"/")[[1]])){
    joy <- rbind(joy, toy[i,])
  }
}
joy <- joy[-1,]
dim(joy)
# 1449 1447

jc <- NA_character_
for(i in 1:nrow(toy)){
  jc <- c(jc, c(strsplit(toy[i,c("KEGG")],"/")[[1]]))
}
jc <- jc[-1]
joy$chemicals <- joy$KEGG
joy$KEGG <- jc


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

DF <- assign_hierarchy(count_data = joy, keep_unknowns = TRUE, identifier = "KEGG")
DF$Class <- as.character(DF$Class)
DF$Class <- ifelse(is.na(DF$Class),"Other", DF$Class)

DF$Subclass_1 <- as.character(DF$Subclass_1)
DF$Subclass_1 <- ifelse(is.na(DF$Subclass_1),"none", DF$Subclass_1)

DF$Subclass_2 <- as.character(DF$Subclass_2)
DF$Subclass_2 <- ifelse(is.na(DF$Subclass_2),"none", DF$Subclass_2)


toy$all_class <- rep(NA_character_, nrow(toy))
toy$all_Subclass_1 <- rep(NA_character_, nrow(toy))
toy$all_Subclass_2 <- rep(NA_character_, nrow(toy))

toy$most_freq_class <- rep(NA_character_, nrow(toy))
toy$most_freq_Subclass_1 <- rep(NA_character_, nrow(toy))
toy$most_freq_Subclass_2 <- rep(NA_character_, nrow(toy))


for(i in 1:length(toy$KEGG)){
  toy$all_class[i] <- knitr::combine_words(unique(DF$Class[which((DF$chemicals %in% toy$KEGG[i]))]), and = "", sep = "/")
  
  if(sum(DF$Class[which((DF$chemicals %in% toy$KEGG[i]))] %in% "Other") == 0){
    toy$most_freq_class[i] <- getmode(DF$Class[which((DF$chemicals %in% toy$KEGG[i]))])
  }
  else 
    toy$most_freq_class[i] <- getmode(DF$Class[which((DF$chemicals %in% toy$KEGG[i]))][DF$Class[which((DF$chemicals %in% toy$KEGG[i]))]!= "Other"])
}
toy$most_freq_class <- ifelse(is.na(toy$most_freq_class), "Other",toy$most_freq_class)


for(i in 1:length(toy$KEGG)){
  toy$all_Subclass_1[i] <- knitr::combine_words(unique(DF$Subclass_1[which((DF$chemicals %in% toy$KEGG[i]))]), and = "", sep = "/")
  
  if(sum(DF$Subclass_1[which((DF$chemicals %in% toy$KEGG[i]))] %in% "none") == 0){
    toy$most_freq_Subclass_1[i] <- getmode(DF$Subclass_1[which((DF$chemicals %in% toy$KEGG[i]))])
  }
  else 
    toy$most_freq_Subclass_1[i] <- getmode(DF$Subclass_1[which((DF$chemicals %in% toy$KEGG[i]))][DF$Subclass_1[which((DF$chemicals %in% toy$KEGG[i]))]!= "none"])
}
toy$most_freq_Subclass_1 <- ifelse(is.na(toy$most_freq_Subclass_1), "none",toy$most_freq_Subclass_1)



for(i in 1:length(toy$KEGG)){
  toy$all_Subclass_2[i] <- knitr::combine_words(unique(DF$Subclass_2[which((DF$chemicals %in% toy$KEGG[i]))]), and = "", sep = "/")
  
  if(sum(DF$Subclass_2[which((DF$chemicals %in% toy$KEGG[i]))] %in% "none") == 0){
    toy$most_freq_Subclass_2[i] <- getmode(DF$Subclass_2[which((DF$chemicals %in% toy$KEGG[i]))])
  }
  else 
    toy$most_freq_Subclass_2[i] <- getmode(DF$Subclass_2[which((DF$chemicals %in% toy$KEGG[i]))][DF$Subclass_2[which((DF$chemicals %in% toy$KEGG[i]))]!= "none"])
}
toy$most_freq_Subclass_2 <- ifelse(is.na(toy$most_freq_Subclass_2), "none",toy$most_freq_Subclass_2)



# toy[,c("KEGG","all_class","most_freq_class")]
data_c18$Class <- toy$most_freq_class
data_c18$Subclass_1 <- toy$most_freq_Subclass_1
data_c18$Subclass_2 <- toy$most_freq_Subclass_2

data_c18[rownames(data_c18)[data_c18$Subclass_1 == "Amino acids" & data_c18$Class == "Organic acids"],"Subclass_1"] = "none"
data_c18[rownames(data_c18)[data_c18$Subclass_1 == "Vitamins" & data_c18$Class == "Organic acids"],"Subclass_1"] = "none"
data_c18$Met_id <- paste0("Met",seq(1:nrow(data_c18)))


# data_c18[rownames(data_c18)[data_c18$Subclass_1 == "Fatty acyls" & data_c18$Class == "Peptides"],"Subclass_1"] = "none"
# data_c18[rownames(data_c18)[data_c18$Subclass_2 == "Amino sugars" & data_c18$Subclass_1 == "Carboxylic acids"],"Subclass_2"] = "none"
# data_c18[rownames(data_c18)[data_c18$Subclass_2 == "Fatty Acids and Conjugates" & data_c18$Subclass_1 == "none"],"Subclass_2"] = "none"
# data_c18[rownames(data_c18)[data_c18$Subclass_1 == "none"],"Subclass_2"] = "none"
# table(data_c18$Class)
# table(data_c18$Subclass_1)
# table(data_c18$Subclass_2)


write.csv(data_c18, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/data_c18.csv", row.names = F)

###################################################################################################################

data_c18 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/data_c18.csv", check.names = F)

# data_c18$Class[order(data_c18$Class,data_c18$Met_id)]
# data_c18$Met_id[order(data_c18$Class,data_c18$Met_id)]

t_data_c18 <- data.table::transpose(data_c18)
# get row and colnames in order
colnames(t_data_c18) <- data_c18$Met_id
t_data_c18$id_Year <- colnames(data_c18)
dim(t_data_c18)


# View(t_data_c18)
# t_data_c18[,c("chem_1","chem_2","id_Year")]

t_data_c18$id <- rep(NA_character_,nrow(t_data_c18))
t_data_c18$Year <- rep(NA_character_,nrow(t_data_c18))

for(i in 1:nrow(t_data_c18)){
  t_data_c18$id[i] <- strsplit(t_data_c18$id_Year[i],"_")[[1]][1]
  t_data_c18$Year[i] <- as.numeric(knitr::combine_words(strsplit(strsplit(t_data_c18$id_Year[i],"_")[[1]][2],"y")[[1]][1],and = "", sep =""))
}


t_data_c18 <- t_data_c18[,c(data_c18$Met_id[order(data_c18$Class,data_c18$Met_id)], "id","Year")]
# median_sum_feature_c18$seq # chemical sequence - very important
merged_omics <- merge(d_long,t_data_c18, by = c("id","Year"))
merged_omics$id <- as.character(merged_omics$id)
# dim(merged_omics)
# head(merged_omics)

## Subset for only 125 ids

merged_omics <- merged_omics[which(merged_omics$id %in% as.vector(t_data_c18[t_data_c18$Year == "7","id"])),]
dim(merged_omics)
# head(merged_omics)
merged_omics$Year <- as.character(merged_omics$Year)
merged_omics$Year <- factor(merged_omics$Year, levels  = c("7","14","22","28"))

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


merged_omics <- merge(merged_omics, d_subset[,c("pfoa_0","pfos_0","pfhxs_0","pfna_0","pfda_0","id",
                                                "cpfoa0","cpfos0","cpfhxs0","cpfna0","cpfda0",
                                                "cpfoa7","cpfos7","cpfhxs7","cpfna7","cpfda7",
                                                "cpfoa14","cpfos14","cpfhxs14","cpfna14","cpfda14",
                                                "cpfoa22","cpfos22","cpfhxs22","cpfna22","cpfda22",
                                                "cpfoa28","cpfos28","cpfhxs28","cpfna28","cpfda28")],by = "id")


merged_omics$cmbmi <- ifelse(merged_omics$mbmi >= 25, "High","Low")
merged_omics$cmage <- ifelse(merged_omics$mage >= as.numeric(quantile(merged_omics$mage)[4]), "High","Low")
merged_omics$cparity <- ifelse(merged_omics$parity >= 1, "1","0")
merged_omics$csmokepreg <- ifelse(merged_omics$smokepreg_2 == 1, "1","0")
merged_omics$cbreastfed_tot <- ifelse(merged_omics$breastfed_tot >= as.numeric(quantile(merged_omics$breastfed_tot)[2]), "High","Low")
merged_omics$cmatfishpreg <- as.numeric(ifelse(merged_omics$matfishpreg_cat2 >= 0.5, "1","0"))
merged_omics[,(colnames(merged_omics) %in% c(data_c18$Met_id[order(data_c18$Class,data_c18$Met_id)]))] <- (apply(merged_omics[,(colnames(merged_omics) %in% c(data_c18$Met_id[order(data_c18$Class,data_c18$Met_id)]))],2,as.numeric))


merged_omics$Carbohydrates <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Class == "Carbohydrates"])])
merged_omics$Hormones <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Class == "Hormones and transmitters"])])
merged_omics$Lipids <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Class == "Lipids"])])
merged_omics$Nucleic_acids <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Class == "Nucleic acids"])])
merged_omics$Organic_acids <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Class == "Organic acids"])])
merged_omics$Other <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Class == "Other"])])
merged_omics$Peptides <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Class == "Peptides"])])
merged_omics$Phytochemical <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Class == "Phytochemical compounds"])])
merged_omics$Vitamins <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Class == "Vitamins and Cofactors"])])

merged_omics$age7 <- merged_omics$age_7
merged_omics$age14 <- merged_omics$age_14

dim(merged_omics[merged_omics$Year == 7,(colnames(merged_omics) %in% c(data_c18$Met_id[order(data_c18$Class,data_c18$Met_id)]))])
dim(merged_omics[merged_omics$Year == 14,(colnames(merged_omics) %in% c(data_c18$Met_id[order(data_c18$Class,data_c18$Met_id)]))])
dim(merged_omics[merged_omics$Year == 22,(colnames(merged_omics) %in% c(data_c18$Met_id[order(data_c18$Class,data_c18$Met_id)]))])
dim(merged_omics[merged_omics$Year == 28,(colnames(merged_omics) %in% c(data_c18$Met_id[order(data_c18$Class,data_c18$Met_id)]))])
# 125 787

write.csv(merged_omics, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/merged_omics_c18.csv", row.names = F)



# Extra codes for fun

# c("Monosaccharides","Disaccharides","Aldoses","Amino_sugars",
#   "Purines","Pyrimidines","Ribonucleotides",
#   "Carboxylic_acids","Oxocarboxylic_acids",
#   "Alkaloids","Alkaloids_lysine","Alkaloids_nicotinic_acid",
#   "Fatty_Acids_Conjugates","Fatty_esters","Flavonoids","Isoprenoids","Sphingoid_bases","Steroids",
#   "Biogenic_amines","Amino_acids","Dipeptides")

# table(merged_omics$Year)
# table(merged_omics$id)

# merged_omics$Monosaccharides <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_1 == "Monosaccharides"])])
# merged_omics$Disaccharides <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Disaccharides"])])
# merged_omics$Aldoses <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Aldoses"])])
# merged_omics$Amino_sugars <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Amino sugars"])])
# 
# 
# merged_omics$Purines <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Purines"])])
# merged_omics$Pyrimidines <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Pyrimidines"])])
# merged_omics$Ribonucleotides <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Ribonucleotides"])])
# 
# 
# merged_omics$Alkaloids <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_1 == "Alkaloids"])])
# merged_omics$Alkaloids_lysine <- (merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Alkaloids derived from lysine"])])
# merged_omics$Alkaloids_nicotinic_acid <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Alkaloids derived from nicotinic acid"])])
# 
# merged_omics$Fatty_Acids_Conjugates <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Fatty Acids and Conjugates"])])
# merged_omics$Fatty_esters <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Fatty esters"])])
# merged_omics$Flavonoids <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Flavonoids"])])
# merged_omics$Isoprenoids <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Isoprenoids"])])
# merged_omics$Sphingoid_bases <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Sphingoid bases"])])
# merged_omics$Steroids <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Steroids"])])
# 
# merged_omics$Biogenic_amines  <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Biogenic amines"])])
# merged_omics$Amino_acids   <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_c18$Subclass_2 == "Amino acids"])])
# merged_omics$Dipeptides   <- rowMeans(merged_omics[,colnames(merged_omics) %in% c(data_c18$Met_id[data_hilic$Subclass_2 == "Dipeptides"])])















































###########################################################################################################################


get_beta_pval <- function(iterations, formula, data, exposure, y){
  model <- lm(formula ,data=data)
  s <- summary(model)
  ts <- s$coefficients[2,"Estimate"]
  pp <- NA_real_
  for(i in 1:iterations){
    
    data_permuted <- data
    
    data_permuted[,exposure] <- sample(data_permuted[,exposure])
    f_perm <- lm(formula ,data=data_permuted)
    s_f_perm <- summary(f_perm)
    pp <- c(pp,s_f_perm$coefficients[2,"Estimate"]) 
  }
  pp <- pp[-1]
  return(list(beta = round(s$coefficients[2,"Estimate"],3), model_pval = round(s$coefficients[2,"Pr(>|t|)"], 3), 
              simu_pval = round( length(pp[abs(pp) > abs(s$coefficients[2,"Estimate"])])/length(pp), 3)))
}


### FI

get_FI <- function(M, formula, data, exposure , y){
  
  model <- lm(formula ,data=data, weights = weights)
  s <- summary(model)
  ts <- s$coefficients[2,"Estimate"]
  
  func_val <- function(ts){
    tol = 0.15
    if(ts >0){
      return(seq(-max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
    }
    else 
      return(seq(-min(abs(c( -abs(ts) - tol, - abs(ts) + tol))), max(abs(c( -abs(ts) - tol, - abs(ts) + tol))), 0.005))
  }
  
  values <-   func_val(ts) 
  pval_FI <- rep(NA_real_,length(values))
  beta_FI <- rep(NA_real_,length(values))
  
  for(j in 1:length(values)){
    
    pp <- NA_real_
    data_permuted <- data
    data_permuted[,y][data_permuted[,exposure] == 1] <- data_permuted[,y][data_permuted[,exposure] == 1] + values[j]
    s <- summary(lm(formula,data=data_permuted, weights = weights))
    
    for(i in 1:M){
      
      
      data_permuted[,exposure] <- sample(data_permuted[,exposure])
      lmer_7 <-  suppressMessages(lm(formula ,data=data_permuted, weights = weights))
      s_lmer_7 <- summary(lmer_7)
      pp <- c(pp,s_lmer_7$coefficients[2,"Estimate"]) 
      
    }
    pp <- pp[-1]
    pval_FI[j] <- length(pp[abs(pp) > abs(s$coefficients[2,"Estimate"])])/length(pp)
    
  }
  
  return(paste0(round(min((values + ts)[pval_FI > 0.05]) + ts,4),", ",
               round( max((values + ts)[pval_FI > 0.05]) + ts, 4)))
  
}













# Plots


p_Carbohydrates <- ggplot(data = merged_omics[merged_omics$sex == 1,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids",
                                                                                                  "Other","Peptides","Phytochemical","Vitamins","id","Year")], aes(x = Year, y = Carbohydrates,group = id))
p_Carbohydrates <- (p_Carbohydrates + geom_line(alpha = 0.4)  + stat_smooth(aes(group = 1), method = "loess") 
                              + stat_summary(aes(group = 1), 
                            geom = "point", fun.y = mean, size = 3)+ theme_bw()) 
p_Carbohydrates <- p_Carbohydrates +  theme(plot.title=element_text(size=14,face="bold"),
                                            axis.text.x=element_text(face="bold"),
                                            axis.text.y = element_text(face = "bold"),
                                            axis.title=element_text(size=14,face="bold")) + ylab("Carbohydrates")
p_Carbohydrates




p_Hormones <- ggplot(data = merged_omics[merged_omics$sex == 1,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids","Other","Peptides","Phytochemical","Vitamins","id","Year")], aes(x = Year, y = Hormones,group = id))
p_Hormones <- (p_Hormones + geom_line(alpha = 0.4)  + stat_smooth(aes(group = 1), method = "loess") 
                 + stat_summary(aes(group = 1), geom = "point", fun.y = mean, size = 3)+ theme_bw()) 
p_Hormones <- p_Hormones +  theme(plot.title=element_text(size=14,face="bold"),
                                  axis.text.x=element_text(face="bold"),
                                  axis.text.y = element_text(face = "bold"),
                                  axis.title=element_text(size=14,face="bold")) + ylab("Hormones and transmitters")
p_Hormones





p_Lipids <- ggplot(data = merged_omics[merged_omics$sex == 1,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids","Other","Peptides","Phytochemical","Vitamins","id","Year")], aes(x = Year, y = Lipids,group = id))
p_Lipids <- (p_Lipids + geom_line(alpha = 0.4)  + stat_smooth(aes(group = 1), method = "loess") 
             + stat_summary(aes(group = 1), geom = "point", fun.y = mean, size = 3)+ theme_bw()) 
p_Lipids <- p_Lipids +  theme(plot.title=element_text(size=14,face="bold"),
                      axis.text.x=element_text(face="bold"),
                      axis.text.y = element_text(face = "bold"),
                      axis.title=element_text(size=14,face="bold")) +  ylab("Lipids") 
p_Lipids





p_Nucleic_acids <- ggplot(data = merged_omics[merged_omics$sex == 1,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids","Other","Peptides","Phytochemical","Vitamins","id","Year")], aes(x = Year, y = Nucleic_acids,group = id))
p_Nucleic_acids <- (p_Nucleic_acids + geom_line(alpha = 0.4)  + stat_smooth(aes(group = 1), method = "loess") 
                    + stat_summary(aes(group = 1), geom = "point", fun.y = mean, size = 3)+ theme_bw()) 
p_Nucleic_acids <- p_Nucleic_acids +  theme(plot.title=element_text(size=14,face="bold"),
                                            axis.text.x=element_text(face="bold"),
                                            axis.text.y = element_text(face = "bold"),
                                            axis.title=element_text(size=14,face="bold")) +  ylab("Nucleic Acids")
p_Nucleic_acids





p_Organic_acids <- ggplot(data = merged_omics[merged_omics$sex == 1,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids","Other","Peptides","Phytochemical","Vitamins","id","Year")], aes(x = Year, y = Organic_acids,group = id))
p_Organic_acids <- (p_Organic_acids + geom_line(alpha = 0.4)  + stat_smooth(aes(group = 1), method = "loess") 
                    + stat_summary(aes(group = 1), geom = "point", fun.y = mean, size = 3)+ theme_bw()) 
p_Organic_acids <- p_Organic_acids +  theme(plot.title=element_text(size=14,face="bold"),
                                            axis.text.x=element_text(face="bold"),
                                            axis.text.y = element_text(face = "bold"),
                                            axis.title=element_text(size=14,face="bold")) +  ylab("Organic Acids")
p_Organic_acids





p_Peptides <- ggplot(data = merged_omics[merged_omics$sex == 1,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids","Other","Peptides","Phytochemical","Vitamins","id","Year")], aes(x = Year, y = Peptides,group = id))
p_Peptides <- (p_Peptides + geom_line(alpha = 0.4)  + stat_smooth(aes(group = 1), method = "loess") 
               + stat_summary(aes(group = 1), geom = "point", fun.y = mean, size = 3)+ theme_bw()) 
p_Peptides <- p_Peptides +  theme(plot.title=element_text(size=14,face="bold"),
                                  axis.text.x=element_text(face="bold"),
                                  axis.text.y = element_text(face = "bold"),
                                  axis.title=element_text(size=14,face="bold")) + ylab("Peptides")
p_Peptides






p_Phytochemical <- ggplot(data = merged_omics[merged_omics$sex == 1,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids","Other","Peptides","Phytochemical","Vitamins","id","Year")], aes(x = Year, y = Phytochemical,group = id))
p_Phytochemical <- (p_Phytochemical + geom_line(alpha = 0.4)  + stat_smooth(aes(group = 1), method = "loess") 
                    + stat_summary(aes(group = 1), geom = "point", fun.y = mean, size = 3)+ theme_bw()) 
p_Phytochemical <- p_Phytochemical +  theme(plot.title=element_text(size=14,face="bold"),
                            axis.text.x=element_text(face="bold"),
                            axis.text.y = element_text(face = "bold"),
                            axis.title=element_text(size=14,face="bold")) +  ylab("Phytochemical compounds")
p_Phytochemical






p_Vitamins <- ggplot(data = merged_omics[merged_omics$sex == 1,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids","Other","Peptides","Phytochemical","Vitamins","id","Year")], aes(x = Year, y = Vitamins,group = id))
p_Vitamins <- (p_Vitamins + geom_line(alpha = 0.4)  + stat_smooth(aes(group = 1), method = "loess") 
               + stat_summary(aes(group = 1), geom = "point", fun.y = mean, size = 3)+ theme_bw()) 
p_Vitamins <- p_Vitamins  +  theme(plot.title=element_text(size=14,face="bold"),
                                   axis.text.x=element_text(face="bold"),
                                   axis.text.y = element_text(face = "bold"),
                                   axis.title=element_text(size=14,face="bold")) + ylab("Vitamins and Cofactors")
p_Vitamins






p_Other <- ggplot(data = merged_omics[merged_omics$sex == 0,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids","Other","Peptides","Phytochemical","Vitamins","id","Year")], aes(x = Year, y = Other,group = id))
p_Other <- (p_Other + geom_line(alpha = 0.4)  + stat_smooth(aes(group = 1), method = "loess") 
            + stat_summary(aes(group = 1), geom = "point", fun.y = mean, size = 3)+ theme_bw()) 
p_Other <-  p_Other +  theme(plot.title=element_text(size=14,face="bold"),
                 axis.text.x=element_text(face="bold"),
                 axis.text.y = element_text(face = "bold"),
                 axis.title=element_text(size=14,face="bold")) +  ylab("Others")
p_Other



jpeg("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/Figures and data for paper/hilic_pos_all_met.jpeg", units="in", width=12, height=12, res=600)
(LDL_plot <- ggpubr::ggarrange(p_Carbohydrates,p_Hormones, p_Lipids, 
                               p_Nucleic_acids, p_Organic_acids, p_Peptides,
                               p_Phytochemical, p_Vitamins, p_Other, ncol = 3,nrow = 3,common.legend = T, legend = "bottom") + theme_bw())

dev.off()


# # Heatmaps
# library("pheatmap")
# 
# D7 <- scale(merged_omics[merged_omics$Year == 7,!(colnames(merged_omics) %in% c("id_Year","id","Year"))])
# h1 <- pheatmap(D7, cutree_cols = 3, show_rownames =F, show_colnames = F)
# 
# D14 <- scale(as.data.frame(na.omit(merged_omics[merged_omics$Year == 14,!(colnames(merged_omics) %in% c("id_Year","id","Year","chem_946"))])))
# h2 <- pheatmap(D14, cutree_cols = 3, show_rownames =F, show_colnames = F)
# 
# D22 <- scale(merged_omics[merged_omics$Year == 22,!(colnames(merged_omics) %in% c("id_Year","id","Year"))])
# pheatmap(D22, cutree_cols = 3, show_rownames =F, show_colnames = F)
# 
# D28 <- scale(merged_omics[merged_omics$Year == 28,!(colnames(merged_omics) %in% c("id_Year","id","Year","chem_946"))])
# pheatmap(D28, cutree_cols = 3, show_rownames =F, show_colnames = F)


# Correlation plots

dt7 <- merged_omics[merged_omics$Year == 7,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids",
                                                                         "Peptides","Phytochemical","Vitamins","Other")]
M7 <- cor(dt7)
colnames(M7) <- c("Carbohydrates","Hormones","Lipids","Nucleic Acids","Organic Acids",
                  "Peptides","Phytochemical","Vitamins","Other")
rownames(M7) <- c("Carbohydrates","Hormones","Lipids","Nucleic Acids","Organic Acids",
                  "Peptides","Phytochemical","Vitamins","Other")

r7 <- corrplot::corrplot(M7, method = 'square', diag = FALSE, order = 'hclust', tl.col = 'black',
                        addrect = 3, rect.col = 'red',bg="gold2",  outline=TRUE,
                        cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))



dt14 <- merged_omics[merged_omics$Year == 14,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids",
                                                                         "Peptides","Phytochemical","Vitamins","Other")]
M14 <- cor(dt14)
colnames(M14) <- c("Carbohydrates","Hormones","Lipids","Nucleic Acids","Organic Acids",
                  "Peptides","Phytochemical","Vitamins","Other")
rownames(M14) <- c("Carbohydrates","Hormones","Lipids","Nucleic Acids","Organic Acids",
                  "Peptides","Phytochemical","Vitamins","Other")

corrplot::corrplot(M14, method = 'square', diag = FALSE, order = 'hclust', tl.col = 'black',
                   addrect = 3, rect.col = 'red',bg="gold2",  outline=TRUE,
                   cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))



dt22 <- merged_omics[merged_omics$Year == 22,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids",
                                                                           "Peptides","Phytochemical","Vitamins","Other")]
M22 <- cor(dt22)
colnames(M22) <- c("Carbohydrates","Hormones","Lipids","Nucleic Acids","Organic Acids",
                   "Peptides","Phytochemical","Vitamins","Other")
rownames(M22) <- c("Carbohydrates","Hormones","Lipids","Nucleic Acids","Organic Acids",
                   "Peptides","Phytochemical","Vitamins","Other")

corrplot::corrplot(M22, method = 'square', diag = FALSE, order = 'hclust', tl.col = 'black',
                   addrect = 3, rect.col = 'red',bg="gold2",  outline=TRUE,
                   cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))



dt28 <- merged_omics[merged_omics$Year == 28,colnames(merged_omics) %in% c("Carbohydrates","Hormones","Lipids","Nucleic_acids","Organic_acids",
                                                                           "Peptides","Phytochemical","Vitamins","Other")]
M28 <- cor(dt28)
colnames(M28) <- c("Carbohydrates","Hormones","Lipids","Nucleic Acids","Organic Acids",
                   "Peptides","Phytochemical","Vitamins","Other")
rownames(M28) <- c("Carbohydrates","Hormones","Lipids","Nucleic Acids","Organic Acids",
                   "Peptides","Phytochemical","Vitamins","Other")

corrplot::corrplot(M28, method = 'square', diag = FALSE, order = 'hclust', tl.col = 'black',
                   addrect = 3, rect.col = 'red',bg="gold2",  outline=TRUE,
                   cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))


































#####################################################################################################################

uid <- unique(t_data_hilic$id)
cv <- rep(NA_real_,length(median_sum_feature_hilic$seq))

for(j in 1:length(median_sum_feature_hilic$seq)){
  v  = NA_real_
    for(i in 1:length(uid)){
      v = c(v,sd(t_data_hilic[t_data_hilic$id == t_data_hilic$id[i],median_sum_feature_hilic$seq[j]], na.rm = T)/mean(t_data_hilic[t_data_hilic$id == t_data_hilic$id[i],median_sum_feature_hilic$seq[j]], na.rm = T))
    }
  cv[j] = mean(v[-1], na.rm=T)*100
}


hist(cv)
median_sum_feature_hilic$cv <- cv

# 01 Feb 2019
# c(04 May 2020,15 Dec 2020, 01 Aug 2020, 11 Sep 2020, 01 Aug 2020, 01 Mar 2020, 12 Aug 2020)
# 01 Sep 2021, 16 Dec 2021, 06 Aug 2021, 01 Sep 2021 , 01 Sep 2021 , 06 Aug 2021, 01 Aug 2021, 06 Mar 2021, 06 Dec 2021,
# 01 Sep 2021, 16 Apr 2021 , 22 Sep 2021, 07 Dec 2021, 01 Aug 2021 , 13 Jun 2021



# median_sum_feature_hilic$seq # chemical sequence - very important
merged_omics <- merge(d_long,t_data_hilic, by = c("id","Year"))
merged_omics$id <- as.character(merged_omics$id)
# dim(merged_omics)
# head(merged_omics)

## Subset for only 125 ids

merged_omics <- merged_omics[which(merged_omics$id %in% as.vector(t_data_hilic[t_data_hilic$Year == "7","id"])),]
dim(merged_omics)
head(merged_omics)

merged_omics <- merge(merged_omics, d_subset[,c("pfoa_0","pfos_0","pfhxs_0","pfna_0","pfda_0","id")],by = "id")


###########################################################################################################################

# PFOA

## At age 7

Pfoa_7 <- data.frame(mz = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     rt = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     pval = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     stat = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     chem_seq = rep(NA_character_,nrow(median_sum_feature_hilic)))

for(i in 1:nrow(median_sum_feature_hilic)){
  s <- rlm(merged_omics[merged_omics$Year == "7", "pfoa_0"] ~ merged_omics[merged_omics$Year == "7",median_sum_feature_hilic$seq[i]] + 
            merged_omics[merged_omics$Year == "7", "mbmi"] + merged_omics[merged_omics$Year == "7", "parity"] + 
            merged_omics[merged_omics$Year == "7", "sex"] + merged_omics[merged_omics$Year == "7", "mage"] +
             merged_omics[merged_omics$Year == "7", "matfishpreg_cat2"]  + merged_omics[merged_omics$Year == "7", "age_7"] , maxit = 1e4, acc = 1e-4, data = merged_omics)
  s <- summary(s)
  st <- as.data.frame(s$coefficients)
  st$pvalue = 1 -pt(abs(st[,"t value"]), df = s$df[2]) + pt(-abs(st[,"t value"]), df = s$df[2])
  
  Pfoa_7$mz[i] = median_sum_feature_hilic$mz[i]
  Pfoa_7$rt[i] = median_sum_feature_hilic$time[i]
  Pfoa_7$stat[i] = st[2,"t value"]
  Pfoa_7$pval[i] = st[2,"pvalue"]
  Pfoa_7$chem_seq[i] = median_sum_feature_hilic$seq[i]
  
}
colnames(Pfoa_7)[1:4] <- c("m/z","retention time","p-value","statistic")
Pfoa_7$p_adjust <- p.adjust(Pfoa_7$`p-value`, "fdr")
Pfoa_7[Pfoa_7$p_adjust < 0.05,]


# s[,"p-value"] = 1 -pt(abs(s$coefficients[,"t value"]), df = s$df[2]) + pt(-abs(s$coefficients[,"t value"]), df = s$df[2])

write.table(Pfoa_7[,!(colnames(Pfoa_7) %in% c("chem_seq"))], file = "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/PFAS at age 0/PFOA/pfoa_at_age_0_metaomics_at_7.txt", sep = "\t",
            row.names = FALSE)


## At age 14

Pfoa_14 <- data.frame(mz = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     rt = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     pval = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     stat = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     chem_seq = rep(NA_character_,nrow(median_sum_feature_hilic)))

for(i in 1:nrow(median_sum_feature_hilic)){
  s <- rlm(merged_omics[merged_omics$Year == "14", "pfoa_0"] ~ merged_omics[merged_omics$Year == "14",median_sum_feature_hilic$seq[i]] + 
             merged_omics[merged_omics$Year == "14", "mbmi"] + merged_omics[merged_omics$Year == "14", "parity"] + 
             merged_omics[merged_omics$Year == "14", "sex"] + merged_omics[merged_omics$Year == "14", "mage"] +
             merged_omics[merged_omics$Year == "14", "matfishpreg_cat2"] + merged_omics[merged_omics$Year == "14", "age_14"] , maxit = 1e4, acc = 1e-4, data = merged_omics)
  s <- summary(s)
  st <- as.data.frame(s$coefficients)
  st$pvalue = 1 -pt(abs(st[,"t value"]), df = s$df[2]) + pt(-abs(st[,"t value"]), df = s$df[2])
  
  Pfoa_14$mz[i] = median_sum_feature_hilic$mz[i]
  Pfoa_14$rt[i] = median_sum_feature_hilic$time[i]
  Pfoa_14$stat[i] = st[2,"t value"]
  Pfoa_14$pval[i] = st[2,"pvalue"]
  Pfoa_14$chem_seq[i] = median_sum_feature_hilic$seq[i]
  
}
colnames(Pfoa_14)[1:4] <- c("m/z","retention time","p-value","statistic")
Pfoa_14$p_adjust <- p.adjust(Pfoa_14$`p-value`, "fdr")
Pfoa_14[Pfoa_14$p_adjust < 0.05,]


# s[,"p-value"] = 1 -pt(abs(s$coefficients[,"t value"]), df = s$df[2]) + pt(-abs(s$coefficients[,"t value"]), df = s$df[2])

write.table(Pfoa_14[,!(colnames(Pfoa_14) %in% c("chem_seq"))], file = "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/PFAS at age 0/PFOA/pfoa_at_age_0_metaomics_at_14.txt", sep = "\t",
            row.names = FALSE)



## At age 22

Pfoa_22 <- data.frame(mz = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      rt = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      pval = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      stat = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      chem_seq = rep(NA_character_,nrow(median_sum_feature_hilic)))

for(i in 1:nrow(median_sum_feature_hilic)){
  s <- rlm(merged_omics[merged_omics$Year == "22", "pfoa_0"] ~ merged_omics[merged_omics$Year == "22",median_sum_feature_hilic$seq[i]] + 
             merged_omics[merged_omics$Year == "22", "mbmi"] + merged_omics[merged_omics$Year == "22", "parity"] + 
             merged_omics[merged_omics$Year == "22", "sex"] + merged_omics[merged_omics$Year == "22", "mage"] +
             merged_omics[merged_omics$Year == "22", "matfishpreg_cat2"] + merged_omics[merged_omics$Year == "22", "age22"] , maxit = 1e4, acc = 1e-4, data = merged_omics)
  s <- summary(s)
  st <- as.data.frame(s$coefficients)
  st$pvalue = 1 -pt(abs(st[,"t value"]), df = s$df[2]) + pt(-abs(st[,"t value"]), df = s$df[2])
  
  Pfoa_22$mz[i] = median_sum_feature_hilic$mz[i]
  Pfoa_22$rt[i] = median_sum_feature_hilic$time[i]
  Pfoa_22$stat[i] = st[2,"t value"]
  Pfoa_22$pval[i] = st[2,"pvalue"]
  Pfoa_22$chem_seq[i] = median_sum_feature_hilic$seq[i]
  
}
colnames(Pfoa_22)[1:4] <- c("m/z","retention time","p-value","statistic")
Pfoa_22$p_adjust <- p.adjust(Pfoa_22$`p-value`, "fdr")
Pfoa_22[Pfoa_22$p_adjust < 0.05,]


# s[,"p-value"] = 1 -pt(abs(s$coefficients[,"t value"]), df = s$df[2]) + pt(-abs(s$coefficients[,"t value"]), df = s$df[2])

write.table(Pfoa_22[,!(colnames(Pfoa_22) %in% c("chem_seq"))], file = "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/PFAS at age 0/PFOA/pfoa_at_age_0_metaomics_at_22.txt", sep = "\t",
            row.names = FALSE)



## At age 28

Pfoa_28 <- data.frame(mz = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      rt = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      pval = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      stat = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      chem_seq = rep(NA_character_,nrow(median_sum_feature_hilic)))

for(i in 1:nrow(median_sum_feature_hilic)){
  s <- rlm(merged_omics[merged_omics$Year == "28", "pfoa_0"] ~ merged_omics[merged_omics$Year == "28",median_sum_feature_hilic$seq[i]] + 
             merged_omics[merged_omics$Year == "28", "mbmi"] + merged_omics[merged_omics$Year == "28", "parity"] + 
             merged_omics[merged_omics$Year == "28", "sex"] + merged_omics[merged_omics$Year == "28", "mage"] +
             merged_omics[merged_omics$Year == "28", "matfishpreg_cat2"]  + merged_omics[merged_omics$Year == "28", "age28"] , maxit = 1e4, acc = 1e-4, data = merged_omics)
  s <- summary(s)
  st <- as.data.frame(s$coefficients)
  st$pvalue = 1 -pt(abs(st[,"t value"]), df = s$df[2]) + pt(-abs(st[,"t value"]), df = s$df[2])
  
  Pfoa_28$mz[i] = median_sum_feature_hilic$mz[i]
  Pfoa_28$rt[i] = median_sum_feature_hilic$time[i]
  Pfoa_28$stat[i] = st[2,"t value"]
  Pfoa_28$pval[i] = st[2,"pvalue"]
  Pfoa_28$chem_seq[i] = median_sum_feature_hilic$seq[i]
  
}
colnames(Pfoa_28)[1:4] <- c("m/z","retention time","p-value","statistic")
Pfoa_28$p_adjust <- p.adjust(Pfoa_28$`p-value`, "fdr")
Pfoa_28[Pfoa_28$p_adjust < 0.05,]
Pfoa_28[Pfoa_28$`p-value` < 0.05,]

# s[,"p-value"] = 1 -pt(abs(s$coefficients[,"t value"]), df = s$df[2]) + pt(-abs(s$coefficients[,"t value"]), df = s$df[2])

write.table(Pfoa_28[,!(colnames(Pfoa_28) %in% c("chem_seq"))], file = "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/PFAS at age 0/PFOA/pfoa_at_age_0_metaomics_at_28.txt", sep = "\t",
            row.names = FALSE)



###########################################################################################################################

# pfna

## At age 7

pfna_7 <- data.frame(mz = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     rt = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     pval = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     stat = rep(NA_real_,nrow(median_sum_feature_hilic)),
                     chem_seq = rep(NA_character_,nrow(median_sum_feature_hilic)))

for(i in 1:nrow(median_sum_feature_hilic)){
  s <- rlm(merged_omics[merged_omics$Year == "7", "pfna_0"] ~ merged_omics[merged_omics$Year == "7",median_sum_feature_hilic$seq[i]] + 
             merged_omics[merged_omics$Year == "7", "mbmi"] + merged_omics[merged_omics$Year == "7", "parity"] + 
             merged_omics[merged_omics$Year == "7", "sex"] + merged_omics[merged_omics$Year == "7", "mage"] +
             merged_omics[merged_omics$Year == "7", "matfishpreg_cat2"]  + merged_omics[merged_omics$Year == "7", "age_7"] , maxit = 1e4, acc = 1e-4, data = merged_omics)
  s <- summary(s)
  st <- as.data.frame(s$coefficients)
  st$pvalue = 1 -pt(abs(st[,"t value"]), df = s$df[2]) + pt(-abs(st[,"t value"]), df = s$df[2])
  
  pfna_7$mz[i] = median_sum_feature_hilic$mz[i]
  pfna_7$rt[i] = median_sum_feature_hilic$time[i]
  pfna_7$stat[i] = st[2,"t value"]
  pfna_7$pval[i] = st[2,"pvalue"]
  pfna_7$chem_seq[i] = median_sum_feature_hilic$seq[i]
  
}
colnames(pfna_7)[1:4] <- c("m/z","retention time","p-value","statistic")
pfna_7$p_adjust <- p.adjust(pfna_7$`p-value`, "fdr")
pfna_7[pfna_7$p_adjust < 0.05,]


# s[,"p-value"] = 1 -pt(abs(s$coefficients[,"t value"]), df = s$df[2]) + pt(-abs(s$coefficients[,"t value"]), df = s$df[2])

write.table(pfna_7[,!(colnames(pfna_7) %in% c("chem_seq"))], file = "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/PFAS at age 0/pfna/pfna_at_age_0_metaomics_at_7.txt", sep = "\t",
            row.names = FALSE)


## At age 14

pfna_14 <- data.frame(mz = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      rt = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      pval = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      stat = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      chem_seq = rep(NA_character_,nrow(median_sum_feature_hilic)))

for(i in 1:nrow(median_sum_feature_hilic)){
  s <- rlm(merged_omics[merged_omics$Year == "14", "pfna_0"] ~ merged_omics[merged_omics$Year == "14",median_sum_feature_hilic$seq[i]] + 
             merged_omics[merged_omics$Year == "14", "mbmi"] + merged_omics[merged_omics$Year == "14", "parity"] + 
             merged_omics[merged_omics$Year == "14", "sex"] + merged_omics[merged_omics$Year == "14", "mage"] +
             merged_omics[merged_omics$Year == "14", "matfishpreg_cat2"] + merged_omics[merged_omics$Year == "14", "age_14"] , maxit = 1e4, acc = 1e-4, data = merged_omics)
  s <- summary(s)
  st <- as.data.frame(s$coefficients)
  st$pvalue = 1 -pt(abs(st[,"t value"]), df = s$df[2]) + pt(-abs(st[,"t value"]), df = s$df[2])
  
  pfna_14$mz[i] = median_sum_feature_hilic$mz[i]
  pfna_14$rt[i] = median_sum_feature_hilic$time[i]
  pfna_14$stat[i] = st[2,"t value"]
  pfna_14$pval[i] = st[2,"pvalue"]
  pfna_14$chem_seq[i] = median_sum_feature_hilic$seq[i]
  
}
colnames(pfna_14)[1:4] <- c("m/z","retention time","p-value","statistic")
pfna_14$p_adjust <- p.adjust(pfna_14$`p-value`, "fdr")
pfna_14[pfna_14$p_adjust < 0.05,]


# s[,"p-value"] = 1 -pt(abs(s$coefficients[,"t value"]), df = s$df[2]) + pt(-abs(s$coefficients[,"t value"]), df = s$df[2])

write.table(pfna_14[,!(colnames(pfna_14) %in% c("chem_seq"))], file = "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/PFAS at age 0/pfna/pfna_at_age_0_metaomics_at_14.txt", sep = "\t",
            row.names = FALSE)



## At age 22

pfna_22 <- data.frame(mz = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      rt = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      pval = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      stat = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      chem_seq = rep(NA_character_,nrow(median_sum_feature_hilic)))

for(i in 1:nrow(median_sum_feature_hilic)){
  s <- rlm(merged_omics[merged_omics$Year == "22", "pfna_0"] ~ merged_omics[merged_omics$Year == "22",median_sum_feature_hilic$seq[i]] + 
             merged_omics[merged_omics$Year == "22", "mbmi"] + merged_omics[merged_omics$Year == "22", "parity"] + 
             merged_omics[merged_omics$Year == "22", "sex"] + merged_omics[merged_omics$Year == "22", "mage"] +
             merged_omics[merged_omics$Year == "22", "matfishpreg_cat2"] + merged_omics[merged_omics$Year == "22", "age22"] , maxit = 1e4, acc = 1e-4, data = merged_omics)
  s <- summary(s)
  st <- as.data.frame(s$coefficients)
  st$pvalue = 1 -pt(abs(st[,"t value"]), df = s$df[2]) + pt(-abs(st[,"t value"]), df = s$df[2])
  
  pfna_22$mz[i] = median_sum_feature_hilic$mz[i]
  pfna_22$rt[i] = median_sum_feature_hilic$time[i]
  pfna_22$stat[i] = st[2,"t value"]
  pfna_22$pval[i] = st[2,"pvalue"]
  pfna_22$chem_seq[i] = median_sum_feature_hilic$seq[i]
  
}
colnames(pfna_22)[1:4] <- c("m/z","retention time","p-value","statistic")
pfna_22$p_adjust <- p.adjust(pfna_22$`p-value`, "fdr")
pfna_22[pfna_22$p_adjust < 0.05,]


# s[,"p-value"] = 1 -pt(abs(s$coefficients[,"t value"]), df = s$df[2]) + pt(-abs(s$coefficients[,"t value"]), df = s$df[2])

write.table(pfna_22[,!(colnames(pfna_22) %in% c("chem_seq"))], file = "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/PFAS at age 0/pfna/pfna_at_age_0_metaomics_at_22.txt", sep = "\t",
            row.names = FALSE)



## At age 28

pfna_28 <- data.frame(mz = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      rt = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      pval = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      stat = rep(NA_real_,nrow(median_sum_feature_hilic)),
                      chem_seq = rep(NA_character_,nrow(median_sum_feature_hilic)))

for(i in 1:nrow(median_sum_feature_hilic)){
  s <- rlm(merged_omics[merged_omics$Year == "28", "pfna_0"] ~ merged_omics[merged_omics$Year == "28",median_sum_feature_hilic$seq[i]] + 
             merged_omics[merged_omics$Year == "28", "mbmi"] + merged_omics[merged_omics$Year == "28", "parity"] + 
             merged_omics[merged_omics$Year == "28", "sex"] + merged_omics[merged_omics$Year == "28", "mage"] +
             merged_omics[merged_omics$Year == "28", "matfishpreg_cat2"]  + merged_omics[merged_omics$Year == "28", "age28"] , maxit = 1e4, acc = 1e-4, data = merged_omics)
  s <- summary(s)
  st <- as.data.frame(s$coefficients)
  st$pvalue = 1 -pt(abs(st[,"t value"]), df = s$df[2]) + pt(-abs(st[,"t value"]), df = s$df[2])
  
  pfna_28$mz[i] = median_sum_feature_hilic$mz[i]
  pfna_28$rt[i] = median_sum_feature_hilic$time[i]
  pfna_28$stat[i] = st[2,"t value"]
  pfna_28$pval[i] = st[2,"pvalue"]
  pfna_28$chem_seq[i] = median_sum_feature_hilic$seq[i]
  
}
colnames(pfna_28)[1:4] <- c("m/z","retention time","p-value","statistic")
pfna_28$p_adjust <- p.adjust(pfna_28$`p-value`, "fdr")
pfna_28[pfna_28$p_adjust < 0.05,]
pfna_28[pfna_28$`p-value` < 0.05,]

# s[,"p-value"] = 1 -pt(abs(s$coefficients[,"t value"]), df = s$df[2]) + pt(-abs(s$coefficients[,"t value"]), df = s$df[2])

write.table(pfna_28[,!(colnames(pfna_28) %in% c("chem_seq"))], file = "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/PFAS at age 0/pfna/pfna_at_age_0_metaomics_at_28.txt", sep = "\t",
            row.names = FALSE)

set.seed(1234)
dfg <- data.frame(x1 = rnorm(100,0,1), x2 = rnorm(100,0,1), x3 = rnorm(100,0,1), x4 = rnorm(100,0,1))
err <- rnorm(100,0,1)
dfg$y <- ifelse(exp(0.5*dfg$x1 + 0.002*dfg$x2 + dfg$x3 +err)/(1+exp(0.5*dfg$x1 + 0.002*dfg$x2 + dfg$x3 + err))>0.5, 1,0)
summary(glm(y~ x1 + x2 + x3 + x4, data= dfg, family = "binomial"))


exposures_wqp <- c("x1","x2","x3","x4")
wqs_all <- gwqs(y ~ wqs , mix_name = exposures_wqp, data = dfg,
                 q = 10, validation = 0, b = 200, b1_pos = T, b1_constr = T,
                 family = "binomial")

summary(wqs_all)
wqs_all$final_weights

dfg$y
dfg$int <- ifelse(rownames(dfg) %in% which(dfg$y ==1)[1:45], 0.9, 0)

exposures_wqp <- c("x1","x2","x3","x4","int")
wqs_all <- gwqs(y ~ wqs , mix_name = exposures_wqp, data = dfg,
                q = 10, validation = 0, b = 200, b1_pos = T, b1_constr = T,
                family = "binomial")

summary(wqs_all)
wqs_all$final_weights



