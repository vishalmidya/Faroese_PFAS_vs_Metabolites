##################################################################################################################################################################################################

# pfos at 0

data_hilic <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/Figures and data for paper/data_hilic.csv", check.names = F)

## Metabolites at 7

data.pfos.0.met_at_7 <- read.table("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 0/minerva_data_pfos_0_metabolite_7/mwas_pfos_0_met_7.txt", header=TRUE)
data.pfos.0.met_at_7$class <- data_hilic$Class
data.pfos.0.met_at_7$Subclass_1 <- data_hilic$Subclass_1
data.pfos.0.met_at_7$Subclass_2 <- data_hilic$Subclass_2

cutoff <- 0.001

data.pfos.0.met_at_7$simu_pval[which(data.pfos.0.met_at_7$simu_pval == 0)] <- cutoff

d_lm_pfas_meta_pos <- data.pfos.0.met_at_7
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
          labs(x = "Beta estimates", y = "-log10(nomial p value)", title = "pfos at age 0 vs. Metabolite Features at age 7") + theme_bw()
        +  geom_label_repel(size = 3.5,
                            box.padding = unit(0.5, "lines"),
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                            force = 2, force_pull = 2, show.legend = FALSE)
        + scale_color_manual(values=c("blue", "black", "red")) + ylim(c(0,4)) )

# + annotate("text", x = 600, y = -log(max(d_lm_pfas_meta_pos$simu_pval[d_lm_pfas_meta_pos$qvalue <= 0.8] - 0.002), base = 10), 
#            label = "min(q-value) <= 0.8", colour = "black")

data.pfos.0.met_at_7_plot <- vol + theme(legend.position = "none",
                                         axis.text.x= element_text(size = 12, face = "bold"),
                                         axis.text.y = element_text(size = 12, face = "bold"),
                                         axis.title=element_text(size=12,face="bold"))
# data.pfos.0.met_at_7_plot

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
data.pfos.0.met_at_14_plot



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
data.pfos.0.met_at_22_plot


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
data.pfos.0.met_at_28_plot



jpeg("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/pfos/pfos at age 0/volcano.jpeg",
     units="in", width=20, height=10, res=800)

ggpubr::ggarrange(data.pfos.0.met_at_7_plot,data.pfos.0.met_at_14_plot,
                  data.pfos.0.met_at_22_plot,data.pfos.0.met_at_28_plot,
                  ncol = 2,nrow = 2,common.legend = F)

dev.off()

