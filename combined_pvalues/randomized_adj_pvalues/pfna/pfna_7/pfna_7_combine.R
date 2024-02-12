library(foreach)
library(doParallel)
library(car)
library(MASS)
library(ggplot2)
library(data.table)

cl <- makeCluster(15) 
registerDoParallel(cl)

start.time <- Sys.time()

d2 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfna/pfna_7/minerva_data_pfna_7_metabolites_14/pfna_7_met_14_hypothetical_test_stat_hilic.txt")
d3 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfna/pfna_7/minerva_data_pfna_7_metabolites_22/pfna_7_met_22_hypothetical_test_stat_hilic.txt")
d4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/pfna_7_met_28_hypothetical_test_stat_hilic.txt")

d6 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_14/pfna_7_met_14_hypothetical_test_stat_c18.txt")
d7 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_22/pfna_7_met_22_hypothetical_test_stat_c18.txt")
d8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/pfna_7_met_28_hypothetical_test_stat_c18.txt")

test_stat_table <- as.data.frame(cbind(d2,d3,d4, d6,d7,d8))


p2 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfna/pfna_7/minerva_data_pfna_7_metabolites_14/pfna_7_met_14_beta_fisher_hilic.txt")
p3 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfna/pfna_7/minerva_data_pfna_7_metabolites_22/pfna_7_met_22_beta_fisher_hilic.txt")
p4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/pfna_7_met_28_beta_fisher_hilic.txt")

p6 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_14/pfna_7_met_14_beta_fisher_c18.txt")
p7 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_22/pfna_7_met_22_beta_fisher_c18.txt")
p8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/pfna_7_met_28_beta_fisher_c18.txt")

fisher_p_val <- c(p2$simu_pval,p3$simu_pval,p4$simu_pval,
                  
                 p6$simu_pval,p7$simu_pval,p8$simu_pval)


hyp_p_val <- foreach(r = 1:ncol(test_stat_table), .combine = 'cbind') %dopar% {
  1 - rank(abs(test_stat_table[,r]))/nrow(test_stat_table) + 1/nrow(test_stat_table)
}
hyp_p_val <- as.data.frame(hyp_p_val)


min_p_nrep <- foreach(p = 1:dim(hyp_p_val)[1], .combine = 'c') %dopar% {
  min(hyp_p_val[p,], na.rm = T)
}

write.table(hyp_p_val,"/sc/arion/projects/Faroese/pfas_met/combined_pvalues/randomized_adj_pvalues/pfna/pfna_7/pfna_7_hypothetical_pval.txt", row.names = FALSE)

# calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
adj_pval <- foreach(i = 1:length(fisher_p_val), .combine = 'c') %dopar% {
  mean(min_p_nrep <= fisher_p_val[i])
}

p2$rand_adj_pval <- adj_pval[1:1991]
p3$rand_adj_pval <- adj_pval[1992: 3982]
p4$rand_adj_pval <- adj_pval[3983: 5973]

p6$rand_adj_pval <- adj_pval[5974: 6760]
p7$rand_adj_pval <- adj_pval[6761: 7547]
p8$rand_adj_pval <- adj_pval[7548: 8334]


write.table(p2,"/sc/arion/projects/Faroese/pfas_met/hilic/pfna/pfna_7/minerva_data_pfna_7_metabolites_14/pfna_7_met_14_beta_fisher_hilic.txt", row.names = FALSE)
write.table(p3,"/sc/arion/projects/Faroese/pfas_met/hilic/pfna/pfna_7/minerva_data_pfna_7_metabolites_22/pfna_7_met_22_beta_fisher_hilic.txt", row.names = FALSE)
write.table(p4,"/sc/arion/projects/Faroese/pfas_met/hilic/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/pfna_7_met_28_beta_fisher_hilic.txt", row.names = FALSE)

write.table(p6,"/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_14/pfna_7_met_14_beta_fisher_c18.txt", row.names = FALSE)
write.table(p7,"/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_22/pfna_7_met_22_beta_fisher_c18.txt", row.names = FALSE)
write.table(p8,"/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/pfna_7_met_28_beta_fisher_c18.txt", row.names = FALSE)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)



