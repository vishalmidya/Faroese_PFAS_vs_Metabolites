library(foreach)
library(doParallel)
library(car)
library(MASS)
library(tidyverse)
library(data.table)

cl <- makeCluster(5) 
registerDoParallel(cl)

start.time <- Sys.time()


#!!!!!!!!!!!!!
keep_metabolites_hilic<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/keep_metabolites_hilic.csv", check.names = F)
keep_metabolites_c18<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/keep_metabolites_c18.csv", check.names = F)

hilic_metid<- keep_metabolites_hilic$Met_id
c18_metid<- keep_metabolites_c18$Met_id

hilic_resultsid<- gsub("Met", "result.", keep_metabolites_hilic$Met_id)
c18_resultsid<- gsub("Met", "result.", keep_metabolites_c18$Met_id)

d4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_hypothetical_test_stat_hilic.txt")

d8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_hypothetical_test_stat_c18.txt") 

test_stat_table <- as.data.frame(cbind(d4, d8))


p4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_hilic.txt")%>% filter(Met_id %in% hilic_metid)

p8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_c18.txt") %>% filter(Met_id %in% c18_metid)

fisher_p_val <- c(p4$simu_pval, p8$simu_pval)


hyp_p_val <- foreach(r = 1:ncol(test_stat_table), .combine = 'cbind') %dopar% {
  1 - rank(abs(test_stat_table[,r]))/nrow(test_stat_table) + 1/nrow(test_stat_table)
}
hyp_p_val <- as.data.frame(hyp_p_val)


min_p_nrep <- foreach(p = 1:dim(hyp_p_val)[1], .combine = 'c') %dopar% {
  min(hyp_p_val[p,], na.rm = T)
}

write.table(hyp_p_val,"/sc/arion/projects/Faroese/pfas_met/combined_pvalues/randomized_adj_pvalues/pfoa/pfoa_22/pfoa_22_hypothetical_pval.txt", row.names = FALSE)
# calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
adj_pval <- foreach(i = 1:length(fisher_p_val), .combine = 'c') %dopar% {
  mean(min_p_nrep <= fisher_p_val[i])
}


#!!!!!!!!!!!!!!!!!!!!!!!
p4$rand_adj_pval <- adj_pval[1:nrow(p4)]

p8$rand_adj_pval <- adj_pval[(nrow(p4)+1):(nrow(p4)+nrow(p8))]

# fdr !!!!!!!!!!!!!!!!!!!!!!!!!!!
p4$fdr <- p.adjust(fisher_p_val, "fdr")[1:nrow(p4)]

p8$fdr <- p.adjust(fisher_p_val, "fdr")[(nrow(p4)+1):(nrow(p4)+nrow(p8))]


write.table(p4,"/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_hilic.txt", row.names = FALSE)

write.table(p8,"/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_c18.txt", row.names = FALSE)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)



