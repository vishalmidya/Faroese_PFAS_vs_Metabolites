library(foreach)
library(doParallel)
library(car)
library(MASS)
library(ggplot2)
library(data.table)

cl <- makeCluster(10) 
registerDoParallel(cl)

start.time <- Sys.time()

data_hilic <- read.csv("/sc/arion/work/midyav01/faroese/data_hilic.csv", check.names = F)
m.out1.pfos0_age7.matched <- read.csv("/sc/arion/work/midyav01/faroese/pfos/age0/matched_data_pfos_at_0_met_at_7.csv")
# m.out1.pfos0_age7.matched <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 0/minerva_data_pfos_0_metabolite_7/matched_data_pfos_at_0_met_at_7.csv")

data = m.out1.pfos0_age7.matched[,c(paste0("Met",seq(1:nrow(data_hilic))), 'cpfos0', 'sex',
                                    'mage',  'mbmi', 'smokepreg_2', 'cmatfishpreg', 'cparity', 'age7' )]

data.pfos.0.met_at_7 <- cbind(data_hilic[,c("mz","time","KEGG","Annotation.confidence.score","chem_name","Met_id")])
data.pfos.0.met_at_7$beta <- rep(NA_real_, nrow(data.pfos.0.met_at_7))
data.pfos.0.met_at_7$model_pval <- rep(NA_real_, nrow(data.pfos.0.met_at_7))
data.pfos.0.met_at_7$simu_pval <- rep(NA_real_, nrow(data.pfos.0.met_at_7))

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

oper <- foreach(i=1:nrow(data.pfos.0.met_at_7), .combine='comb', .multicombine=TRUE,
                .init=list(list(), list())) %dopar% {
                  exposure <- paste0("Met",i)
                  model <- (lm(data[,i] ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age7, data = data))
                  s <- summary(model)
                  list(s$coefficients[2,"Estimate"], s$coefficients[2,"Pr(>|t|)"])
                }

data.pfos.0.met_at_7$beta <- unlist(oper[[1]])
data.pfos.0.met_at_7$model_pval <- unlist(oper[[2]])
observed_test_stat <- unlist(oper[[1]])


d1 <- fread("/sc/arion/work/midyav01/faroese/pfos/age0/test_stat_table_pfos_0_met_7_v1.txt")
d2 <- fread("/sc/arion/work/midyav01/faroese/pfos/age0/test_stat_table_pfos_0_met_7_v2.txt")
d3 <- fread("/sc/arion/work/midyav01/faroese/pfos/age0/test_stat_table_pfos_0_met_7_v3.txt")
d4 <- fread("/sc/arion/work/midyav01/faroese/pfos/age0/test_stat_table_pfos_0_met_7_v4.txt")

l = list(d1,d2, d3, d4)
test_stat_table <- as.data.frame(data.table::rbindlist(l, use.names=TRUE))
# write.table(test_stat_table,"/sc/arion/work/midyav01/faroese/pfos/age0/test_stat_pfos_0_met_7.txt", row.names = FALSE)


observed_test_stat <- data.pfos.0.met_at_7$beta
p_values <- foreach(p = 1:dim(data.pfos.0.met_at_7)[1], .combine = 'c') %dopar% {
  mean(abs(test_stat_table[,p]) >= abs(observed_test_stat[p]))
}
data.pfos.0.met_at_7$simu_pval <- p_values


hyp_p_value <- foreach(r = 1:ncol(test_stat_table), .combine = 'cbind') %dopar% {
  1 - rank(abs(test_stat_table[,r]))/nrow(test_stat_table) + 1/nrow(test_stat_table)
}
hyp_p_value <- as.data.frame(hyp_p_value)


min_p_nrep <- foreach(p = 1:dim(hyp_p_value)[1], .combine = 'c') %dopar% {
  min(hyp_p_value[p,], na.rm = T)
}

wrt <- data.frame(q = c(0.05, 0.1, 0.2, 0.25, 1/3, 0.4, 0.5), val = quantile(min_p_nrep, c(0.05, 0.1, 0.2, 0.25, 1/3, 0.4, 0.5)))
write.table(wrt, "/sc/arion/work/midyav01/faroese/pfos/age0/minp_pfos_0_met_7.txt", row.names = FALSE)


# calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
data.pfos.0.met_at_7$adj_pval <- foreach(i = 1:length(p_values), .combine = 'c') %dopar% {
  mean(min_p_nrep <= p_values[i])
}


write.table(data.pfos.0.met_at_7,"/sc/arion/work/midyav01/faroese/pfos/age0/2e5_mwas_pfos_0_met_7.txt", row.names = FALSE)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)


###############################################


# d <- read.table("/sc/arion/work/midyav01/faroese/pfos/age0/15e4_mwas_pfos_0_met_7.txt", header = T)



