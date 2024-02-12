library(foreach)
library(doParallel)
library(car)
library(MASS)
library(ggplot2)
library(data.table)

cl <- makeCluster(15) 
registerDoParallel(cl)

start.time <- Sys.time()

d1 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_7/pfhxs_0_met_7_hypothetical_test_stat_hilic.txt")
d2 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_14/pfhxs_0_met_14_hypothetical_test_stat_hilic.txt")
d3 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_22/pfhxs_0_met_22_hypothetical_test_stat_hilic.txt")
d4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_28/pfhxs_0_met_28_hypothetical_test_stat_hilic.txt")

d5 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_7/pfhxs_0_met_7_hypothetical_test_stat_c18.txt")
d6 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_14/pfhxs_0_met_14_hypothetical_test_stat_c18.txt")
d7 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_22/pfhxs_0_met_22_hypothetical_test_stat_c18.txt")
d8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_28/pfhxs_0_met_28_hypothetical_test_stat_c18.txt")

test_stat_table <- as.data.frame(cbind(d1,d2,d3,d4,d5,d6,d7,d8))


p1 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_7/pfhxs_0_met_7_beta_fisher_hilic.txt")
p2 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_14/pfhxs_0_met_14_beta_fisher_hilic.txt")
p3 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_22/pfhxs_0_met_22_beta_fisher_hilic.txt")
p4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_28/pfhxs_0_met_28_beta_fisher_hilic.txt")

p5 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_7/pfhxs_0_met_7_beta_fisher_c18.txt")
p6 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_14/pfhxs_0_met_14_beta_fisher_c18.txt")
p7 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_22/pfhxs_0_met_22_beta_fisher_c18.txt")
p8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_28/pfhxs_0_met_28_beta_fisher_c18.txt")

fisher_p_val <- c(p1$simu_pval,p2$simu_pval,p3$simu_pval,p4$simu_pval,
              
              p5$simu_pval,p6$simu_pval,p7$simu_pval,p8$simu_pval)


hyp_p_val <- foreach(r = 1:ncol(test_stat_table), .combine = 'cbind') %dopar% {
  1 - rank(abs(test_stat_table[,r]))/nrow(test_stat_table) + 1/nrow(test_stat_table)
}
hyp_p_val <- as.data.frame(hyp_p_val)


min_p_nrep <- foreach(p = 1:dim(hyp_p_val)[1], .combine = 'c') %dopar% {
  min(hyp_p_val[p,], na.rm = T)
}

write.table(hyp_p_val,"/sc/arion/projects/Faroese/pfas_met/combined_pvalues/randomized_adj_pvalues/pfhxs/pfhxs_0/pfhxs_0_hypothetical_pval.txt", row.names = FALSE)

# calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
adj_pval <- foreach(i = 1:length(fisher_p_val), .combine = 'c') %dopar% {
  mean(min_p_nrep <= fisher_p_val[i])
}

p1$rand_adj_pval <- adj_pval[1:1107]
p2$rand_adj_pval <- adj_pval[1108: 2214]
p3$rand_adj_pval <- adj_pval[2215: 3321]
p4$rand_adj_pval <- adj_pval[3322: 4428]

p5$rand_adj_pval <- adj_pval[4429: 5080]
p6$rand_adj_pval <- adj_pval[5081: 5732]
p7$rand_adj_pval <- adj_pval[5733: 6384]
p8$rand_adj_pval <- adj_pval[6385: 7036]


write.table(p1,"/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_7/pfhxs_0_met_7_beta_fisher_hilic.txt", row.names = FALSE)
write.table(p2,"/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_14/pfhxs_0_met_14_beta_fisher_hilic.txt", row.names = FALSE)
write.table(p3,"/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_22/pfhxs_0_met_22_beta_fisher_hilic.txt", row.names = FALSE)
write.table(p4,"/sc/arion/projects/Faroese/pfas_met/hilic/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_28/pfhxs_0_met_28_beta_fisher_hilic.txt", row.names = FALSE)

write.table(p5,"/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_7/pfhxs_0_met_7_beta_fisher_c18.txt", row.names = FALSE)
write.table(p6,"/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_14/pfhxs_0_met_14_beta_fisher_c18.txt", row.names = FALSE)
write.table(p7,"/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_22/pfhxs_0_met_22_beta_fisher_c18.txt", row.names = FALSE)
write.table(p8,"/sc/arion/projects/Faroese/pfas_met/c18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_28/pfhxs_0_met_28_beta_fisher_c18.txt", row.names = FALSE)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)



