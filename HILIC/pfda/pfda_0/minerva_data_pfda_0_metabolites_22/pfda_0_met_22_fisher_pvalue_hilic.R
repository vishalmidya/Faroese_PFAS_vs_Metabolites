library(foreach)
library(doParallel)
library(car)
library(MASS)
library(ggplot2)
library(data.table)
library(iterators)
library(parallel)
library(carData)


cores=detectCores()
cl <- makeCluster(10) 
registerDoParallel(cl)

start.time <- Sys.time()

data_hilic <- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/data_hilic.csv", check.names = F)
m.out1.pfda_0_age22.matched <- read.csv("/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_0/minerva_data_pfda_0_metabolites_22/matched_data_pfda_at_0_met_at_22.csv")

met_name<- data_hilic$Met_id

data = m.out1.pfda_0_age22.matched[,c(met_name, 'cpfda0', 'sex',
                                    'mage',  'mbmi', 'smokepreg_2', 'cmatfishpreg', 'cparity', 'age22' )]

data.pfda_0.met_at_22 <- cbind(data_hilic[,c("mz","time","Met_id")])
data.pfda_0.met_at_22$beta <- rep(NA_real_, nrow(data.pfda_0.met_at_22))
data.pfda_0.met_at_22$model_pval <- rep(NA_real_, nrow(data.pfda_0.met_at_22))
data.pfda_0.met_at_22$simu_pval <- rep(NA_real_, nrow(data.pfda_0.met_at_22))

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# Model test statistic and raw p-value

oper <- foreach(i=1:nrow(data.pfda_0.met_at_22), .combine='comb', .multicombine=TRUE,
                .init=list(list(), list())) %dopar% {
                  exposure <- met_name[i]
                  model <- (lm(data[,i] ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data = data))
                  s <- summary(model)
                  list(s$coefficients[2,"Estimate"], s$coefficients[2,"Pr(>|t|)"])
                }

data.pfda_0.met_at_22$beta <- unlist(oper[[1]])
data.pfda_0.met_at_22$model_pval <- unlist(oper[[2]])
observed_test_stat <- unlist(oper[[1]])


# Fisher p-value

iterations = 50000
test_stat <- rep(NA_real_, iterations)

set.seed(runif(1,0,1e4))
data_permuted <- data
system.time(test_stat <- foreach(i=1:nrow(data.pfda_0.met_at_22), .combine='rbind', .multicombine=TRUE) %:% 
              foreach(j = 1:iterations, .combine = 'c') %dopar% {
                f_perm <- lm(sample(data_permuted[,i]) ~ cpfda0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22, data=data_permuted)
                s_f_perm <- summary(f_perm)
                s_f_perm$coefficients[2,"Estimate"]
                
              }
)

test_stat <- as.data.frame(test_stat)
colnames(test_stat) <- paste0("iter",seq(1,iterations))
test_stat_table <- data.table::transpose(test_stat)
colnames(test_stat_table) <- rownames(test_stat)
rownames(test_stat_table) <- colnames(test_stat)

write.table(test_stat_table,"/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_0/minerva_data_pfda_0_metabolites_22/pfda_0_met_22_hypothetical_test_stat_hilic.txt", row.names = FALSE)


p_values <- foreach(p = 1:dim(data.pfda_0.met_at_22)[1], .combine = 'c') %dopar% {
  mean(abs(test_stat_table[,p]) >= abs(observed_test_stat[p]))
}
data.pfda_0.met_at_22$simu_pval <- p_values

write.table(data.pfda_0.met_at_22,"/sc/arion/projects/Faroese/pfas_met/hilic/pfda/pfda_0/minerva_data_pfda_0_metabolites_22/pfda_0_met_22_beta_fisher_hilic.txt", row.names = FALSE)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)
