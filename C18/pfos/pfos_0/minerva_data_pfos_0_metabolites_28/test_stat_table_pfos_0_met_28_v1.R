library(foreach)
library(doParallel)
library(car)
library(MASS)
library(ggplot2)
library(data.table)


cores=detectCores()
cl <- makeCluster(15) 
registerDoParallel(cl)

start.time <- Sys.time()

data_hilic <- read.csv("/sc/arion/work/midyav01/faroese/data_hilic.csv", check.names = F)
m.out1.pfos0_age28.matched <- read.csv("/sc/arion/work/midyav01/faroese/pfos/age0/matched_data_pfos_at_0_met_at_28.csv")
# m.out1.pfos0_age28.matched <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 0/minerva_data_pfos_0_metabolite_28/matched_data_pfos_at_0_met_at_28.csv")

data = m.out1.pfos0_age28.matched[,c(paste0("Met",seq(1:nrow(data_hilic))), 'cpfos0', 'sex',
                                    'mage',  'mbmi', 'smokepreg_2', 'cmatfishpreg', 'cparity', 'age28' )]
iterations = 50000
data.pfos.0.met_at_28 <- cbind(data_hilic[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Met_id")])
data.pfos.0.met_at_28$beta <- rep(NA_real_, nrow(data.pfos.0.met_at_28))
data.pfos.0.met_at_28$model_pval <- rep(NA_real_, nrow(data.pfos.0.met_at_28))
data.pfos.0.met_at_28$simu_pval <- rep(NA_real_, nrow(data.pfos.0.met_at_28))
observed_test_stat <- NA_real_
test_stat <- rep(NA_real_, iterations)


set.seed(runif(1,0,1e4))
data_permuted <- data
system.time(test_stat <- foreach(i=1:nrow(data.pfos.0.met_at_28), .combine='rbind', .multicombine=TRUE) %:% 
              foreach(j = 1:iterations, .combine = 'c') %dopar% {
                f_perm <- lm(sample(data_permuted[,i]) ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted)
                s_f_perm <- summary(f_perm)
                s_f_perm$coefficients[2,"Estimate"]
                
              }
)

# transpose

test_stat <- as.data.frame(test_stat)
colnames(test_stat) <- paste0("iter",seq(1,iterations))
test_stat_table <- data.table::transpose(test_stat)
colnames(test_stat_table) <- rownames(test_stat)
rownames(test_stat_table) <- colnames(test_stat)

write.table(test_stat_table,"/sc/arion/work/midyav01/faroese/pfos/age0/test_stat_table_pfos_0_met_28_v1.txt", row.names = FALSE)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)
