library(foreach)
library(doParallel)
library(car)
library(MASS)
library(ggplot2)
library(data.table)

cl <- makeCluster(10) 
registerDoParallel(cl)

start.time <- Sys.time()


#!!!!!!!!!!!!!
keep_metabolites_hilic<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/keep_metabolites_hilic.csv", check.names = F)
keep_metabolites_c18<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/keep_metabolites_c18.csv", check.names = F)

hilic_metid<- keep_metabolites_hilic$Met_id
c18_metid<- keep_metabolites_c18$Met_id

hilic_resultsid<- gsub("Met", "result.", keep_metabolites_hilic$Met_id)
c18_resultsid<- gsub("Met", "result.", keep_metabolites_c18$Met_id)

d3 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/pfoa_14_met_22_hypothetical_test_stat_hilic.txt")%>% select(all_of(hilic_resultsid))
d4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/pfoa_14_met_28_hypothetical_test_stat_hilic.txt")%>% select(all_of(hilic_resultsid))

d7 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/pfoa_14_met_22_hypothetical_test_stat_c18.txt") %>% select(all_of(c18_resultsid))
d8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/pfoa_14_met_28_hypothetical_test_stat_c18.txt") %>% select(all_of(c18_resultsid))

test_stat_table <- as.data.frame(cbind(d3,d4, d7,d8))


p3 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/pfoa_14_met_22_beta_fisher_hilic.txt")%>% filter(Met_id %in% hilic_metid)
p4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/pfoa_14_met_28_beta_fisher_hilic.txt")%>% filter(Met_id %in% hilic_metid)

p7 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/pfoa_14_met_22_beta_fisher_c18.txt") %>% filter(Met_id %in% c18_metid)
p8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/pfoa_14_met_28_beta_fisher_c18.txt") %>% filter(Met_id %in% c18_metid)

fisher_p_val <- c(p3$simu_pval,p4$simu_pval,
                  
                  p7$simu_pval,p8$simu_pval)


hyp_p_val <- foreach(r = 1:ncol(test_stat_table), .combine = 'cbind') %dopar% {
  1 - rank(abs(test_stat_table[,r]))/nrow(test_stat_table) + 1/nrow(test_stat_table)
}
hyp_p_val <- as.data.frame(hyp_p_val)


min_p_nrep <- foreach(p = 1:dim(hyp_p_val)[1], .combine = 'c') %dopar% {
  min(hyp_p_val[p,], na.rm = T)
}

write.table(hyp_p_val,"/sc/arion/projects/Faroese/pfas_met/combined_pvalues/randomized_adj_pvalues/pfoa/pfoa_14/pfoa_14_hypothetical_pval.txt", row.names = FALSE)


# calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
adj_pval <- foreach(i = 1:length(fisher_p_val), .combine = 'c') %dopar% {
  mean(min_p_nrep <= fisher_p_val[i])
}

#!!!!!!!!!!!!!!!!!!!!!!!
p3$rand_adj_pval <- adj_pval[1:nrow(p1)]
p4$rand_adj_pval <- adj_pval[(nrow(p1)+1):(nrow(p1)*2)]

p7$rand_adj_pval <- adj_pval[(nrow(p1)*2+1):(nrow(p1)*2+nrow(p5))]
p8$rand_adj_pval <- adj_pval[(nrow(p1)*2+nrow(p5)+1):(nrow(p1)*2+nrow(p5)*2)]

# fdr !!!!!!!!!!!!!!!!!!!!!!!!!!!
p3$fdr<- p.adjust(p3$simu_pval, "fdr")
p4$fdr<- p.adjust(p4$simu_pval, "fdr")

p7$fdr<- p.adjust(p7$simu_pval, "fdr")
p8$fdr<- p.adjust(p8$simu_pval, "fdr")


write.table(p3,"/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/pfoa_14_met_22_beta_fisher_hilic.txt", row.names = FALSE)
write.table(p4,"/sc/arion/projects/Faroese/pfas_met/hilic/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/pfoa_14_met_28_beta_fisher_hilic.txt", row.names = FALSE)

write.table(p7,"/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/pfoa_14_met_22_beta_fisher_c18.txt", row.names = FALSE)
write.table(p8,"/sc/arion/projects/Faroese/pfas_met/c18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/pfoa_14_met_28_beta_fisher_c18.txt", row.names = FALSE)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)



