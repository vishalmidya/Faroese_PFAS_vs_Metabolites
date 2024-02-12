library(foreach)
library(doParallel)
library(car)
library(MASS)
library(ggplot2)
library(data.table)

cl <- makeCluster(5) 
registerDoParallel(cl)

start.time <- Sys.time()

d4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_hypothetical_test_stat_hilic.txt")

d8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_hypothetical_test_stat_c18.txt")

test_stat_table <- as.data.frame(cbind(d4, d8))


p4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_hilic.txt")

p8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_c18.txt")

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

p4$rand_adj_pval <- adj_pval[1:1107]

p8$rand_adj_pval <- adj_pval[1108: 1759]


write.table(p4,"/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_hilic.txt", row.names = FALSE)

write.table(p8,"/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_c18.txt", row.names = FALSE)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)



