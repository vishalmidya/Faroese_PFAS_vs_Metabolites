library(foreach)
library(doParallel)
library(car)
library(MASS)
library(tidyverse)
library(data.table)

cl <- makeCluster(15) 
registerDoParallel(cl)

start.time <- Sys.time()


#!!!!!!!!!!!!!
keep_metabolites_hilic<- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/keep_metabolites_hilic.csv", check.names = F)
keep_metabolites_c18<- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/keep_metabolites_c18.csv", check.names = F)

hilic_metid<- keep_metabolites_hilic$Met_id
c18_metid<- keep_metabolites_c18$Met_id

hilic_resultsid<- gsub("Met", "result.", keep_metabolites_hilic$Met_id)
c18_resultsid<- gsub("Met", "result.", keep_metabolites_c18$Met_id)

# !!!!!!!!!!!!!!!
d1 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/pfos_0_met_7_hypothetical_test_stat_hilic.txt")
d2 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/pfos_0_met_14_hypothetical_test_stat_hilic.txt")
d3 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/pfos_0_met_22_hypothetical_test_stat_hilic.txt")
d4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/pfos_0_met_28_hypothetical_test_stat_hilic.txt")

d5 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/pfos_0_met_7_hypothetical_test_stat_c18.txt") 
d6 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/pfos_0_met_14_hypothetical_test_stat_c18.txt") 
d7 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/pfos_0_met_22_hypothetical_test_stat_c18.txt") 
d8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/pfos_0_met_28_hypothetical_test_stat_c18.txt") 

test_stat_table <- as.data.frame(cbind(d1,d2,d3,d4,d5,d6,d7,d8))


p1 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/pfos_0_met_7_beta_fisher_hilic.txt")%>% filter(Met_id %in% hilic_metid)
p2 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/pfos_0_met_14_beta_fisher_hilic.txt")%>% filter(Met_id %in% hilic_metid)
p3 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/pfos_0_met_22_beta_fisher_hilic.txt")%>% filter(Met_id %in% hilic_metid)
p4 <- fread("/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/pfos_0_met_28_beta_fisher_hilic.txt")%>% filter(Met_id %in% hilic_metid)

p5 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/pfos_0_met_7_beta_fisher_c18.txt") %>% filter(Met_id %in% c18_metid)
p6 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/pfos_0_met_14_beta_fisher_c18.txt") %>% filter(Met_id %in% c18_metid)
p7 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/pfos_0_met_22_beta_fisher_c18.txt") %>% filter(Met_id %in% c18_metid)
p8 <- fread("/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/pfos_0_met_28_beta_fisher_c18.txt") %>% filter(Met_id %in% c18_metid)

fisher_p_val <- c(p1$simu_pval,p2$simu_pval,p3$simu_pval,p4$simu_pval,
              
              p5$simu_pval,p6$simu_pval,p7$simu_pval,p8$simu_pval)


hyp_p_val <- foreach(r = 1:ncol(test_stat_table), .combine = 'cbind') %dopar% {
  1 - rank(abs(test_stat_table[,r]))/nrow(test_stat_table) + 1/nrow(test_stat_table)
}
hyp_p_val <- as.data.frame(hyp_p_val)


min_p_nrep <- foreach(p = 1:dim(hyp_p_val)[1], .combine = 'c') %dopar% {
  min(hyp_p_val[p,], na.rm = T)
}

write.table(hyp_p_val,"/sc/arion/projects/Faroese/pfas_met/combined_pvalues/randomized_adj_pvalues/pfos/pfos_0/pfos_0_hypothetical_pval.txt", row.names = FALSE)

# calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
adj_pval <- foreach(i = 1:length(fisher_p_val), .combine = 'c') %dopar% {
  mean(min_p_nrep <= fisher_p_val[i])
}

#!!!!!!!!!!!!!!!!!!!!!!!
p1$rand_adj_pval <- adj_pval[1:nrow(p1)]
p2$rand_adj_pval <- adj_pval[(nrow(p1)+1):(nrow(p1)*2)]
p3$rand_adj_pval <- adj_pval[(nrow(p1)*2+1):(nrow(p1)*3)]
p4$rand_adj_pval <- adj_pval[(nrow(p1)*3+1):(nrow(p1)*4)]

p5$rand_adj_pval <- adj_pval[(nrow(p1)*4+1):(nrow(p1)*4+nrow(p5))]
p6$rand_adj_pval <- adj_pval[(nrow(p1)*4+nrow(p5)+1):(nrow(p1)*4+nrow(p5)*2)]
p7$rand_adj_pval <- adj_pval[(nrow(p1)*4+nrow(p5)*2+1):(nrow(p1)*4+nrow(p5)*3)]
p8$rand_adj_pval <- adj_pval[(nrow(p1)*4+nrow(p5)*3+1):(nrow(p1)*4+nrow(p5)*4)]


# fdr !!!!!!!!!!!!!!!!!!!!!!!!!!!
p1$fdr <- p.adjust(fisher_p_val, "fdr")[1:nrow(p1)]
p2$fdr <- p.adjust(fisher_p_val, "fdr")[(nrow(p1)+1):(nrow(p1)*2)]
p3$fdr <- p.adjust(fisher_p_val, "fdr")[(nrow(p1)*2+1):(nrow(p1)*3)]
p4$fdr <- p.adjust(fisher_p_val, "fdr")[(nrow(p1)*3+1):(nrow(p1)*4)]


p5$fdr <- p.adjust(fisher_p_val, "fdr")[(nrow(p1)*4+1):(nrow(p1)*4+nrow(p5))]
p6$fdr <- p.adjust(fisher_p_val, "fdr")[(nrow(p1)*4+nrow(p5)+1):(nrow(p1)*4+nrow(p5)*2)]
p7$fdr <- p.adjust(fisher_p_val, "fdr")[(nrow(p1)*4+nrow(p5)*2+1):(nrow(p1)*4+nrow(p5)*3)]
p8$fdr <- p.adjust(fisher_p_val, "fdr")[(nrow(p1)*4+nrow(p5)*3+1):(nrow(p1)*4+nrow(p5)*4)]




write.table(p1,"/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/pfos_0_met_7_beta_fisher_hilic.txt", row.names = FALSE)
write.table(p2,"/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/pfos_0_met_14_beta_fisher_hilic.txt", row.names = FALSE)
write.table(p3,"/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/pfos_0_met_22_beta_fisher_hilic.txt", row.names = FALSE)
write.table(p4,"/sc/arion/projects/Faroese/pfas_met/hilic/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/pfos_0_met_28_beta_fisher_hilic.txt", row.names = FALSE)

write.table(p5,"/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/pfos_0_met_7_beta_fisher_c18.txt", row.names = FALSE)
write.table(p6,"/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/pfos_0_met_14_beta_fisher_c18.txt", row.names = FALSE)
write.table(p7,"/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/pfos_0_met_22_beta_fisher_c18.txt", row.names = FALSE)
write.table(p8,"/sc/arion/projects/Faroese/pfas_met/c18/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/pfos_0_met_28_beta_fisher_c18.txt", row.names = FALSE)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)


pfas_met_tab_all <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/pfas_met_tab_all.csv", check.names = F)

pfhxs_0<- pfas_met_tab_all %>% 
          filter(PFAS == "PFHxS" & PFAS_age == 0)
pfhxs_0$FDR<- p.adjust(pfhxs_0$simu_pval, "fdr")
