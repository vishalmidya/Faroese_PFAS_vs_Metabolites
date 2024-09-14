library(car)
library(readr)
library(dplyr)
library(knitr)
library(xtable)
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
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(Hmisc)
library(gtsummary)
library(data.table)



pfas_met_tab_all_125 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/pfas_met_tab_all.csv", check.names = F)
pfas_met_tab_all_500 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/500/pfas_met_tab_all_500.csv", check.names = F)
names(pfas_met_tab_all_500)[4]<- "beta_500"
names(pfas_met_tab_all_500)[5]<- "model_pval_500"
names(pfas_met_tab_all_500)[6]<- "simu_pval_500"
names(pfas_met_tab_all_500)[7]<- "rand_adj_pval_500"
names(pfas_met_tab_all_500)[8]<- "fdr_500"

pfas_met_tab_all_125_sig<- pfas_met_tab_all_125 %>% 
                           filter(rand_adj_pval < 0.3)

pfas_met_tab_all_500_sig<- pfas_met_tab_all_500 %>% 
                           filter(rand_adj_pval_500 < 0.3)


## check metabolites significant in 125 samples
std125_sig125<- fread("~/Projects/Faroese_Minerva/check/std125_sig125.txt")
std500_sig125<- fread("~/Projects/Faroese_Minerva/check/std500_sig125.txt")
### EFFECT SIZES
pfas_met_tab_all_sig125<- pfas_met_tab_all_500 %>% 
                          select(Met_id, PFAS, Age, Mode, PFAS_age, beta_500, model_pval_500, simu_pval_500, rand_adj_pval_500, fdr_500) %>% 
                          right_join(pfas_met_tab_all_125_sig, by = c("Met_id", "PFAS", "Age", "Mode", "PFAS_age")) %>% 
                          inner_join(std125_sig125, by = c("Met_id", "PFAS", "Age", "Mode", "PFAS_age")) %>% 
                          inner_join(std500_sig125, by = c("Met_id", "PFAS", "Age", "Mode", "PFAS_age")) %>% 
                          mutate(es_diff = beta - beta_500,
                                 std_diff = std.err.125 - std.err.500) %>% 
                          select(PFAS, PFAS_age, Mode, Met_id, Age, Metabolite, Metabolite_multiple, mz, time, beta, beta_500, 
                                 std.err.125, std.err.500, model_pval, model_pval_500, simu_pval, simu_pval_500, 
                                 rand_adj_pval, rand_adj_pval_500, fdr, fdr_500, es_diff, std_diff)









