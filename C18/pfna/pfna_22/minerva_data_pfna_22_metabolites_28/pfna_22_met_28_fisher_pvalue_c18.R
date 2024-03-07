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

data_c18 <- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/data_c18.csv", check.names = F)
m.out1.pfna_22_age28.matched <- read.csv("/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_22/minerva_data_pfna_22_metabolites_28/matched_data_pfna_at_22_met_at_28.csv")

met_name<- data_c18$Met_id

data = m.out1.pfna_22_age28.matched[,c(met_name, 'cpfna22', 'sex',
                                    'mage',  'mbmi', 'smokepreg_2', 'cmatfishpreg', 'cparity', 'age28' )]

data.pfna_22.met_at_28 <- cbind(data_c18[,c("mz","time","Met_id")])
data.pfna_22.met_at_28$beta <- rep(NA_real_, nrow(data.pfna_22.met_at_28))
data.pfna_22.met_at_28$model_pval <- rep(NA_real_, nrow(data.pfna_22.met_at_28))
data.pfna_22.met_at_28$simu_pval <- rep(NA_real_, nrow(data.pfna_22.met_at_28))

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# Model test statistic and raw p-value

oper <- foreach(i=1:nrow(data.pfna_22.met_at_28), .combine='comb', .multicombine=TRUE,
                .init=list(list(), list())) %dopar% {
                  exposure <- paste0("Met",i)
                  model <- (lm(data[,i] ~ cpfna22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data = data))
                  s <- summary(model)
                  list(s$coefficients[2,"Estimate"], s$coefficients[2,"Pr(>|t|)"])
                }

data.pfna_22.met_at_28$beta <- unlist(oper[[1]])
data.pfna_22.met_at_28$model_pval <- unlist(oper[[2]])
observed_test_stat <- unlist(oper[[1]])


# Fisher p-value

iterations = 50000
test_stat <- rep(NA_real_, iterations)

set.seed(runif(1,0,1e4))
data_permuted <- data
system.time(test_stat <- foreach(i=1:nrow(data.pfna_22.met_at_28), .combine='rbind', .multicombine=TRUE) %:% 
              foreach(j = 1:iterations, .combine = 'c') %dopar% {
                f_perm <- lm(sample(data_permuted[,i]) ~ cpfna22 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28, data=data_permuted)
                s_f_perm <- summary(f_perm)
                s_f_perm$coefficients[2,"Estimate"]
                
              }
)

test_stat <- as.data.frame(test_stat)
colnames(test_stat) <- paste0("iter",seq(1,iterations))
test_stat_table <- data.table::transpose(test_stat)
colnames(test_stat_table) <- rownames(test_stat)
rownames(test_stat_table) <- colnames(test_stat)

write.table(test_stat_table,"/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_22/minerva_data_pfna_22_metabolites_28/pfna_22_met_28_hypothetical_test_stat_c18.txt", row.names = FALSE)


p_values <- foreach(p = 1:dim(data.pfna_22.met_at_28)[1], .combine = 'c') %dopar% {
  mean(abs(test_stat_table[,p]) >= abs(observed_test_stat[p]))
}
data.pfna_22.met_at_28$simu_pval <- p_values

write.table(data.pfna_22.met_at_28,"/sc/arion/projects/Faroese/pfas_met/c18/pfna/pfna_22/minerva_data_pfna_22_metabolites_28/pfna_22_met_28_beta_fisher_c18.txt", row.names = FALSE)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)
