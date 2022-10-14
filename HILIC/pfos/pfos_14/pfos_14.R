# merged_omics_pfos0 <- merged_omics[,c(colnames(merged_omics) %in% c(data_hilic$Met_id), )]
# matfishpreg_cat2
# cmatfishpreg

library(ggplot2)
library(MatchIt)
library(ggrepel)


pfos AT AGE 14 

# pfos14

# sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14


m.out1.pfos14_age14 <- matchit(cpfos14 ~    sex + mage + mbmi   + cmatfishpreg  + cparity + age14,
                              data = merged_omics[merged_omics$Year == 14,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.1)

summary(m.out1.pfos14_age14)
f <- love.plot(m.out1.pfos14_age14)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


m.out1.pfos14_age22 <- matchit(cpfos14 ~  sex  + mbmi  + cmatfishpreg  + cparity + age22, 
                              data = merged_omics[merged_omics$Year == 22,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.5)

summary(m.out1.pfos14_age22)
f <- love.plot(m.out1.pfos14_age22)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


m.out1.pfos14_age28 <- matchit(cpfos14 ~  sex + mage   + smokepreg_2 + cmatfishpreg  + cparity + age28, 
                              data = merged_omics[merged_omics$Year == 28,], discard = "both", method = "full", 
                              distance = "glm", caliper = 0.5)

summary(m.out1.pfos14_age28)
f <- love.plot(m.out1.pfos14_age28)
f + labs(title  = " ", x= "Standardized Mean Difference") + geom_vline(xintercept  = 0.1 , linetype="dotted",  color = "black", size=1.5) + geom_vline(xintercept  = -0.1 , linetype="dotted", color = "black", size=1.5)


m.out1.pfos14_age14.matched <- match.data(m.out1.pfos14_age14)
write.csv(m.out1.pfos14_age14.matched, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 14/minerva_data_pfos_14_metabolite_14/matched_data_pfos_at_14_met_at_14.csv",
          row.names = F)

m.out1.pfos14_age22.matched <- match.data(m.out1.pfos14_age22)
write.csv(m.out1.pfos14_age22.matched, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 14/minerva_data_pfos_14_metabolite_22/matched_data_pfos_at_14_met_at_22.csv",
          row.names = F)


m.out1.pfos14_age28.matched <- match.data(m.out1.pfos14_age28)
write.csv(m.out1.pfos14_age28.matched, "C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 14/minerva_data_pfos_14_metabolite_28/matched_data_pfos_at_14_met_at_28.csv",
          row.names = F)



########################################################################################################################################################################################

# pfos at 0

data_hilic <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/Figures and data for paper/data_hilic.csv", check.names = F)


########################

## Metabolites at 14

data.pfos.0.met_at_14 <- read.table("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 0/minerva_data_pfos_0_metabolite_14/mwas_pfos_0_met_14.txt", header=TRUE)
data.pfos.0.met_at_14$class <- data_hilic$Class
data.pfos.0.met_at_14$Subclass_1 <- data_hilic$Subclass_1
data.pfos.0.met_at_14$Subclass_2 <- data_hilic$Subclass_2

cutoff <- 0.001

data.pfos.0.met_at_14$simu_pval[which(data.pfos.0.met_at_14$simu_pval == 0)] <- cutoff

d_lm_pfas_meta_pos <- data.pfos.0.met_at_14
d_lm_pfas_plot <- d_lm_pfas_meta_pos
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$beta > 0 & d_lm_pfas_plot$simu_pval <= 0.05] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$beta < 0 & d_lm_pfas_plot$simu_pval <= 0.05] <- "Negative"

cut_label <- 0.01
d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$simu_pval < cut_label ] <- d_lm_pfas_plot$chem_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$simu_pval < cut_label]
for(i in 1:nrow(d_lm_pfas_plot)) d_lm_pfas_plot$delabel[i] <- strsplit(d_lm_pfas_plot$delabel,";")[[i]][1]


dft <- table(d_lm_pfas_plot$delabel[d_lm_pfas_plot$simu_pval < cut_label])
r <- names(dft)[as.numeric(dft) >1]
lst <- NA_character_
if(length(r)!= 0){
        for(i in 1:length(r)){
                m <-  which.min(d_lm_pfas_plot$simu_pval[d_lm_pfas_plot$delabel %in% r[i]])
                ml <- length(d_lm_pfas_plot$simu_pval[d_lm_pfas_plot$delabel %in% r[i]])
                lst <- c(lst, d_lm_pfas_plot$Met_id[d_lm_pfas_plot$delabel %in% r[i]][setdiff(seq(1,ml),m)])
        }
}
lst <- lst[-1]

vol <- (ggplot(d_lm_pfas_plot[!(d_lm_pfas_plot$Met_id %in% lst),], aes(x=beta, y=-log10(simu_pval), col=Association, label=delabel)) +# Show all points
                geom_point(size=2) +
                geom_hline(yintercept= -log(0.05, base = 10), color = "black", size = 1, ) + 
                geom_hline(yintercept= -log(cutoff, 10),
                           color = "green", size = 1, ) + 
                labs(x = "Beta estimates", y = "-log10(nomial p value)", title = "pfos at age 0 vs. Metabolite Features at age 14") + theme_bw()
        +  geom_label_repel(size = 3.5,
                            box.padding = unit(0.5, "lines"),
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                            force = 2, force_pull = 2, show.legend = FALSE)
        + scale_color_manual(values=c("blue", "black", "red"))  + ylim(c(0,4)))

# + annotate("text", x = 600, y = -log(max(d_lm_pfas_meta_pos$simu_pval[d_lm_pfas_meta_pos$qvalue <= 0.8] - 0.002), base = 10), 
#            label = "min(q-value) <= 0.8", colour = "black")

data.pfos.0.met_at_14_plot <- vol + theme(legend.position = "none",
                                         axis.text.x= element_text(size = 12, face = "bold"),
                                         axis.text.y = element_text(size = 12, face = "bold"),
                                         axis.title=element_text(size=12,face="bold"))
# data.pfos.0.met_at_14_plot


########################

## Metabolites at 22

data.pfos.0.met_at_22 <- read.table("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 0/minerva_data_pfos_0_metabolite_22/mwas_pfos_0_met_22.txt", header=TRUE)
data.pfos.0.met_at_22$class <- data_hilic$Class
data.pfos.0.met_at_22$Subclass_1 <- data_hilic$Subclass_1
data.pfos.0.met_at_22$Subclass_2 <- data_hilic$Subclass_2

cutoff <- 0.001

data.pfos.0.met_at_22$simu_pval[which(data.pfos.0.met_at_22$simu_pval == 0)] <- cutoff

d_lm_pfas_meta_pos <- data.pfos.0.met_at_22
d_lm_pfas_plot <- d_lm_pfas_meta_pos
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$beta > 0 & d_lm_pfas_plot$simu_pval <= 0.05] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$beta < 0 & d_lm_pfas_plot$simu_pval <= 0.05] <- "Negative"

cut_label <- 0.01
d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$simu_pval < cut_label ] <- d_lm_pfas_plot$chem_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$simu_pval < cut_label]
for(i in 1:nrow(d_lm_pfas_plot)) d_lm_pfas_plot$delabel[i] <- strsplit(d_lm_pfas_plot$delabel,";")[[i]][1]


dft <- table(d_lm_pfas_plot$delabel[d_lm_pfas_plot$simu_pval < cut_label])
r <- names(dft)[as.numeric(dft) >1]
lst <- NA_character_
if(length(r)!= 0){
        for(i in 1:length(r)){
                m <-  which.min(d_lm_pfas_plot$simu_pval[d_lm_pfas_plot$delabel %in% r[i]])
                ml <- length(d_lm_pfas_plot$simu_pval[d_lm_pfas_plot$delabel %in% r[i]])
                lst <- c(lst, d_lm_pfas_plot$Met_id[d_lm_pfas_plot$delabel %in% r[i]][setdiff(seq(1,ml),m)])
        }
}
lst <- lst[-1]

vol <- (ggplot(d_lm_pfas_plot[!(d_lm_pfas_plot$Met_id %in% lst),], aes(x=beta, y=-log10(simu_pval), col=Association, label=delabel)) +# Show all points
                geom_point(size=2) +
                geom_hline(yintercept= -log(0.05, base = 10), color = "black", size = 1, ) + 
                geom_hline(yintercept= -log(cutoff, 10),
                           color = "green", size = 1, ) + 
                labs(x = "Beta estimates", y = "-log10(nomial p value)", title = "pfos at age 0 vs. Metabolite Features at age 22") + theme_bw()
        +  geom_label_repel(size = 3.5,
                            box.padding = unit(0.5, "lines"),
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                            force = 2, force_pull = 2, show.legend = FALSE)
        + scale_color_manual(values=c("blue", "black", "red"))  + ylim(c(0,4)))

# + annotate("text", x = 600, y = -log(max(d_lm_pfas_meta_pos$simu_pval[d_lm_pfas_meta_pos$qvalue <= 0.8] - 0.002), base = 10), 
#            label = "min(q-value) <= 0.8", colour = "black")

data.pfos.0.met_at_22_plot <- vol + theme(legend.position = "none",
                                          axis.text.x= element_text(size = 12, face = "bold"),
                                          axis.text.y = element_text(size = 12, face = "bold"),
                                          axis.title=element_text(size=12,face="bold"))
# data.pfos.0.met_at_22_plot



########################

## Metabolites at 28

data.pfos.0.met_at_28 <- read.table("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 0/minerva_data_pfos_0_metabolite_28/mwas_pfos_0_met_28.txt", header=TRUE)
data.pfos.0.met_at_28$class <- data_hilic$Class
data.pfos.0.met_at_28$Subclass_1 <- data_hilic$Subclass_1
data.pfos.0.met_at_28$Subclass_2 <- data_hilic$Subclass_2

cutoff <- 0.001

data.pfos.0.met_at_28$simu_pval[which(data.pfos.0.met_at_28$simu_pval == 0)] <- cutoff

d_lm_pfas_meta_pos <- data.pfos.0.met_at_28
d_lm_pfas_plot <- d_lm_pfas_meta_pos
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$beta > 0 & d_lm_pfas_plot$simu_pval <= 0.05] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$beta < 0 & d_lm_pfas_plot$simu_pval <= 0.05] <- "Negative"

cut_label <- 0.01
d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$simu_pval < cut_label ] <- d_lm_pfas_plot$chem_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$simu_pval < cut_label]
for(i in 1:nrow(d_lm_pfas_plot)) d_lm_pfas_plot$delabel[i] <- strsplit(d_lm_pfas_plot$delabel,";")[[i]][1]


dft <- table(d_lm_pfas_plot$delabel[d_lm_pfas_plot$simu_pval < cut_label])
r <- names(dft)[as.numeric(dft) >1]
lst <- NA_character_
if(length(r)!= 0){
        for(i in 1:length(r)){
                m <-  which.min(d_lm_pfas_plot$simu_pval[d_lm_pfas_plot$delabel %in% r[i]])
                ml <- length(d_lm_pfas_plot$simu_pval[d_lm_pfas_plot$delabel %in% r[i]])
                lst <- c(lst, d_lm_pfas_plot$Met_id[d_lm_pfas_plot$delabel %in% r[i]][setdiff(seq(1,ml),m)])
        }
}
lst <- lst[-1]

vol <- (ggplot(d_lm_pfas_plot[!(d_lm_pfas_plot$Met_id %in% lst),], aes(x=beta, y=-log10(simu_pval), col=Association, label=delabel)) +# Show all points
                geom_point(size=2) +
                geom_hline(yintercept= -log(0.05, base = 10), color = "black", size = 1, ) + 
                geom_hline(yintercept= -log(cutoff, 10),
                           color = "green", size = 1, ) + 
                labs(x = "Beta estimates", y = "-log10(nomial p value)", title = "pfos at age 0 vs. Metabolite Features at age 28") + theme_bw()
        +  geom_label_repel(size = 3.5,
                            box.padding = unit(0.5, "lines"),
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                            force = 2, force_pull = 2, show.legend = FALSE)
        + scale_color_manual(values=c("blue", "black", "red"))  + ylim(c(0,4)))

# + annotate("text", x = 600, y = -log(max(d_lm_pfas_meta_pos$simu_pval[d_lm_pfas_meta_pos$qvalue <= 0.8] - 0.002), base = 10), 
#            label = "min(q-value) <= 0.8", colour = "black")

data.pfos.0.met_at_28_plot <- vol + theme(legend.position = "none",
                                          axis.text.x= element_text(size = 12, face = "bold"),
                                          axis.text.y = element_text(size = 12, face = "bold"),
                                          axis.title=element_text(size=12,face="bold"))
# data.pfos.0.met_at_28_plot



jpeg("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 0/volcano.jpeg",
     units="in", width=20, height=10, res=800)

ggpubr::ggarrange(data.pfos.0.met_at_14_plot,data.pfos.0.met_at_14_plot,
                  data.pfos.0.met_at_22_plot,data.pfos.0.met_at_28_plot,
                  ncol = 2,nrow = 2,common.legend = F)

dev.off()


pfos_met <- rbind(met_tab_pfos_28,met_tab_pfos_22,met_tab_pfos_14,met_tab_pfos_14)


################################################################################################

























































################################################################################################################################################

# Fit Cubic Splines

library(mgcv)

mod_gam14 <- gam(Aldoses ~ s(pfos_0, bs="cr") + s(age14, bs="cr") + s(mbmi, bs = 'cr') + s(mage, bs = 'cr') 
                 + parity + cmatfishpreg  + sex, data=m.out1.pfos0_age14.matched)
summary(mod_gam14)
# plot(mod_gam14)

mod_gam14 <- gam(Aldoses ~ s(pfos_0, bs="cr") + s(age14, bs="cr") + s(mbmi, bs = 'cr') + s(mage, bs = 'cr') 
                + parity + cmatfishpreg  + sex, data=m.out1.pfos0_age14.matched)
summary(mod_gam14)
# plot(mod_gam14)

mod_gam22 <- gam(Aldoses ~ s(pfos_0, bs="cr") + s(age22, bs="cr") + s(mbmi, bs = 'cr') + s(mage, bs = 'cr') 
                 + parity + cmatfishpreg  + sex, data=m.out1.pfos0_age22.matched)
summary(mod_gam22)
plot(mod_gam22)

mod_gam28 <- gam(Aldoses ~ s(pfos_0, bs="cr") + s(age28, bs="cr") + s(mbmi, bs = 'cr') + s(mage, bs = 'cr') 
                 + parity + cmatfishpreg  + sex, data=m.out1.pfos0_age28.matched)
summary(mod_gam28)
# plot(mod_gam22)


####################################################################################################################################################################

# carbohydrate

## Age 14

### Aldoses

start.time <- Sys.time()
get_beta_pval(M = 1e3, formula = as.formula(Aldoses ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Aldoses")
end.time <- Sys.time()
(time.taken <- end.time - start.time)

# $beta
# [1] 0.014
# 
# $model_pval
# [1] 0.1486
# 
# $simu_pval
# [1] 0.1494


start.time <- Sys.time()
get_FI(M = 5e2, formula = as.formula(Aldoses ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Aldoses")
end.time <- Sys.time()
(time.taken <- end.time - start.time)

### Amino_sugars

get_beta_pval(M = 1e4, formula = as.formula(Amino_sugars ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Amino_sugars")

# $beta
# [1] 0.084
# 
# $model_pval
# [1] 0.002
# 
# $simu_pval
# [1] 0.002

get_FI(M = 5e2, formula = as.formula(Amino_sugars ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Amino_sugars")


### Disaccharides

get_beta_pval(M = 1e3, formula = as.formula(Disaccharides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Disaccharides")

# $beta
# [1] 0.025
# 
# $model_pval
# [1] 0.619
# 
# $simu_pval
# [1] 0.6143


get_FI(M = 5e2, formula = as.formula(Disaccharides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Disaccharides")

#################################################################################

## Age 14

### Aldoses

start.time <- Sys.time()
get_beta_pval(M = 1e3, formula = as.formula(Aldoses ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Aldoses")
end.time <- Sys.time()
(time.taken <- end.time - start.time)

# $beta
# [1] 0.016
# 
# $model_pval
# [1] 0.1416
# 
# $simu_pval
# [1] 0.1438


start.time <- Sys.time()
get_FI(M = 5e2, formula = as.formula(Aldoses ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Aldoses")
end.time <- Sys.time()
(time.taken <- end.time - start.time)

### Amino_sugars

get_beta_pval(M = 1e3, formula = as.formula(Amino_sugars ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Amino_sugars")

# $beta
# [1] 0.04
# 
# $model_pval
# [1] 0.106
# 
# $simu_pval
# [1] 0.138

get_FI(M = 5e2, formula = as.formula(Amino_sugars ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Amino_sugars")


### Disaccharides

get_beta_pval(M = 1e3, formula = as.formula(Disaccharides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Disaccharides")

# $beta
# [1] 0.029
# 
# $model_pval
# [1] 0.48
# 
# $simu_pval
# [1] 0.5514


get_FI(M = 5e2, formula = as.formula(Disaccharides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Disaccharides")


#################################################################################

## Age 22

### Aldoses

start.time <- Sys.time()
get_beta_pval(M = 1e3, formula = as.formula(Aldoses ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + I(age22)),
              data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Aldoses")
end.time <- Sys.time()
(time.taken <- end.time - start.time)

# $beta
# [1] -0.01414
# 
# $model_pval
# [1] 0.259
# 
# $simu_pval
# [1] 0.266


start.time <- Sys.time()
get_FI(M = 5e2, formula = as.formula(Aldoses ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
       data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Aldoses")
end.time <- Sys.time()
(time.taken <- end.time - start.time)

### Amino_sugars

get_beta_pval(M = 1e4, formula = as.formula(Amino_sugars ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
              data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Amino_sugars")

# $beta
# [1] -0.04
# 
# $model_pval
# [1] 0.101
# 
# $simu_pval
# [1] 0.1214

get_FI(M = 5e2, formula = as.formula(Amino_sugars ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
       data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Amino_sugars")

# [1] "-0.3159, -0.0609"


### Disaccharides

get_beta_pval(M = 1e3, formula = as.formula(Disaccharides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
              data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Disaccharides")

# $beta
# [1] -0.046
# 
# $model_pval
# [1] 0.282
# 
# $simu_pval
# [1] 0.333


get_FI(M = 5e2, formula = as.formula(Disaccharides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
       data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Disaccharides")


#################################################################################

## Age 28

### Aldoses

start.time <- Sys.time()
get_beta_pval(M = 1e3, formula = as.formula(Aldoses ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
              data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Aldoses")
end.time <- Sys.time()
(time.taken <- end.time - start.time)

# $beta
# [1] 0.01
# 
# $model_pval
# [1] 0.843
# 
# $simu_pval
# [1] 0.86

start.time <- Sys.time()
get_FI(M = 5e2, formula = as.formula(Aldoses ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
       data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Aldoses")
end.time <- Sys.time()
(time.taken <- end.time - start.time)

### Amino_sugars

get_beta_pval(M = 1e3, formula = as.formula(Amino_sugars ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
              data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Amino_sugars")

# $beta
# [1] -0.01
# 
# $model_pval
# [1] 0.1412
# 
# $simu_pval
# [1] 0.1435

get_FI(M = 5e2, formula = as.formula(Amino_sugars ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
       data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Amino_sugars")


### Disaccharides

get_beta_pval(M = 1e3, formula = as.formula(Disaccharides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
              data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Disaccharides")

# $beta
# [1] -0.041
# 
# $model_pval
# [1] 0.24
# 
# $simu_pval
# [1] 0.289


get_FI(M = 5e2, formula = as.formula(Disaccharides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
       data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Disaccharides")




####################################################################################################################################################################

# Hormones

## Age 14

get_beta_pval(M = 1e3, formula = as.formula(Hormones ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Hormones")

# $beta
# [1] 0.0025
# 
# $model_pval
# [1] 0.6492
# 
# $simu_pval
# [1] 0.68

get_FI(M = 5e2, formula = as.formula(Hormones ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Hormones")

## Age 14

get_beta_pval(M = 1e3, formula = as.formula(Hormones ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Hormones")


# $beta
# [1] 0.0022
# 
# $model_pval
# [1] 0.6852
# 
# $simu_pval
# [1] 0.14214

get_FI(M = 5e2, formula = as.formula(Hormones ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Hormones")

## Age 22

get_beta_pval(M = 1e3, formula = as.formula(Hormones ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Hormones")


# $beta
# [1] 0.0023
# 
# $model_pval
# [1] 0.6806
# 
# $simu_pval
# [1] 0.1426

get_FI(M = 5e2, formula = as.formula(Hormones ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Hormones")

## Age 28

get_beta_pval(M = 1e3, formula = as.formula(Hormones ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Hormones")


# $beta
# [1] 0.0024
# 
# $model_pval
# [1] 0.65144
# 
# $simu_pval
# [1] 0.681



get_FI(M = 5e2, formula = as.formula(Hormones ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Hormones")



####################################################################################################################################################################

# Nucleic acids

## Age 14

### Purines

get_beta_pval(M = 1e3, formula = as.formula(Purines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Purines")

# $beta
# [1] -0.033
# 
# $model_pval
# [1] 0.3114
# 
# $simu_pval
# [1] 0.38

get_FI(M = 5e2, formula = as.formula(Purines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Purines")

## Age 14

get_beta_pval(M = 1e3, formula = as.formula(Purines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Purines")


# $beta
# [1] 0.022
# 
# $model_pval
# [1] 0.5314
# 
# $simu_pval
# [1] 0.583


get_FI(M = 5e2, formula = as.formula(Purines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Purines")

## Age 22

get_beta_pval(M = 1e3, formula = as.formula(Purines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
              data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Purines")

# $beta
# [1] -0.022
# 
# $model_pval
# [1] 0.5
# 
# $simu_pval
# [1] 0.523

get_FI(M = 5e2, formula = as.formula(Purines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
       data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Purines")

## Age 28

get_beta_pval(M = 1e3, formula = as.formula(Purines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
              data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Purines")


# $beta
# [1] 0.03
# 
# $model_pval
# [1] 0.333
# 
# $simu_pval
# [1] 0.358

get_FI(M = 5e2, formula = as.formula(Purines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
       data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Purines")


### Pyrimidines

get_beta_pval(M = 1e3, formula = as.formula(Pyrimidines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Pyrimidines")

# $beta
# [1] 0.114
# 
# $model_pval
# [1] 0
# 
# $simu_pval
# [1] 0


get_FI(M = 5e2, formula = as.formula(Pyrimidines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Pyrimidines")

## Age 14

get_beta_pval(M = 1e3, formula = as.formula(Pyrimidines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Pyrimidines")

# $beta
# [1] 0.022
# 
# $model_pval
# [1] 0.5314
# 
# $simu_pval
# [1] 0.594

get_FI(M = 5e2, formula = as.formula(Pyrimidines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Pyrimidines")

## Age 22

get_beta_pval(M = 1e3, formula = as.formula(Pyrimidines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
              data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Pyrimidines")


# $beta
# [1] 0.061
# 
# $model_pval
# [1] 0.01414
# 
# $simu_pval
# [1] 0.092

get_FI(M = 5e2, formula = as.formula(Pyrimidines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
       data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Pyrimidines")

## Age 28

get_beta_pval(M = 1e3, formula = as.formula(Pyrimidines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
              data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Pyrimidines")

# $beta
# [1] 0.0114
# 
# $model_pval
# [1] 0.621
# 
# $simu_pval
# [1] 0.639

get_FI(M = 5e2, formula = as.formula(Pyrimidines ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
       data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Pyrimidines")



### Ribonucleotides

get_beta_pval(M = 1e3, formula = as.formula(Ribonucleotides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Ribonucleotides")

# $beta
# [1] -0.01
# 
# $model_pval
# [1] 0.822
# 
# $simu_pval
# [1] 0.86


get_FI(M = 5e2, formula = as.formula(Ribonucleotides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Ribonucleotides")

## Age 14

get_beta_pval(M = 1e3, formula = as.formula(Ribonucleotides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Ribonucleotides")


# $beta
# [1] -0.044
# 
# $model_pval
# [1] 0.33
# 
# $simu_pval
# [1] 0.3314


get_FI(M = 5e2, formula = as.formula(Ribonucleotides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Ribonucleotides")

## Age 22

get_beta_pval(M = 1e3, formula = as.formula(Ribonucleotides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
              data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Ribonucleotides")


# $beta
# [1] -0.014
# 
# $model_pval
# [1] 0.1422
# 
# $simu_pval
# [1] 0.1431


get_FI(M = 5e2, formula = as.formula(Ribonucleotides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
       data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Ribonucleotides")

## Age 28

get_beta_pval(M = 1e3, formula = as.formula(Ribonucleotides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
              data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Ribonucleotides")


# $beta
# [1] 0.033
# 
# $model_pval
# [1] 0.4814
# 
# $simu_pval
# [1] 0.511
# 


get_FI(M = 5e2, formula = as.formula(Ribonucleotides ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
       data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Ribonucleotides")




####################################################################################################################################################################

# Vitamins

## Age 14

get_beta_pval(M = 1e3, formula = as.formula(Vitamins ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Vitamins")

# $beta
# [1] 0.044
# 
# $model_pval
# [1] 0.054
# 
# $simu_pval
# [1] 0.08


get_FI(M = 5e2, formula = as.formula(Vitamins ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Vitamins")

## Age 14

get_beta_pval(M = 1e3, formula = as.formula(Vitamins ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
              data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Vitamins")


# $beta
# [1] 0.0314
# 
# $model_pval
# [1] 0.113
# 
# $simu_pval
# [1] 0.136


get_FI(M = 5e2, formula = as.formula(Vitamins ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age_14),
       data = m.out1.pfos0_age14.matched , exposure = "cpfos0", y = "Vitamins")

## Age 22

get_beta_pval(M = 1e3, formula = as.formula(Vitamins ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
              data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Vitamins")

# $beta
# [1] 0.036
# 
# $model_pval
# [1] 0.13
# 
# $simu_pval
# [1] 0.1145


get_FI(M = 5e2, formula = as.formula(Vitamins ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age22),
       data = m.out1.pfos0_age22.matched , exposure = "cpfos0", y = "Vitamins")

## Age 28

get_beta_pval(M = 1e3, formula = as.formula(Vitamins ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
              data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Vitamins")

# $beta
# [1] 0.028
# 
# $model_pval
# [1] 0.269
# 
# $simu_pval
# [1] 0.319

get_FI(M = 5e2, formula = as.formula(Vitamins ~ cpfos0 + sex + mage + mbmi  + smokepreg_2 + cmatfishpreg  + cparity + age28),
       data = m.out1.pfos0_age28.matched , exposure = "cpfos0", y = "Vitamins")

