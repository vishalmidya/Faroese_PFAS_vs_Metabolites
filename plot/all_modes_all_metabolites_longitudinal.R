library(car)
library(readr)
library(lattice)
library(nlme)
library(ggplot2)
library(GGally)
library(nnet)
library(foreign)
library(biotools)
library(glmmML)
library(MASS)
library(lme4)
library(multcomp)
library(dplyr)
library(knitr)
library(xtable)
library(kableExtra)
library(DT)
library(glmnet)
library(corrplot)
library(ggpubr)
library("EnvStats")
library(lmerTest)
library("merTools")
library(reshape2)
library(ggplot2)
library(GGally)
library(mgcv)
library(gplots)
library(tidyr)
library(cluster)
library(factoextra)
library(psych)
library(xgboost)
library(SHAPforxgboost)
library(bkmr)
library(mi)
library(stargazer)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) 
library(spatstat)
library(Hmisc)
library(gtsummary)
library(blme)
library("labelled")
library(lavaan)
library(igraph)
library(semPlot)
library(grpreg)
library(robustHD)
library(gWQS)
library(gridExtra)
library(ggcorrplot)
library(BWQS)
library(qwraps2)
library(MatchIt)
library(data.table)
library(ggrepel)
library(mice)
library(ggrepel)

data_hilic <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/data_hilic.csv", check.names = F)
data_c18 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/data_c18.csv", check.names = F)

pfas_met_tab_all <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/pfas_met_tab_all.csv", check.names = F)

## AGE 0
pfas_met_tab <- pfas_met_tab_all [(pfas_met_tab_all$Met_id == "Met867" & pfas_met_tab_all$PFAS == "PFOA") & pfas_met_tab_all$Mode == "HILIC"| 
                                    (pfas_met_tab_all$Met_id == "Met258" & pfas_met_tab_all$PFAS == "PFOS" & pfas_met_tab_all$Mode == "HILIC") | 
                                    (pfas_met_tab_all$Met_id == "Met300" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode == "C18")| 
                                    (pfas_met_tab_all$Met_id == "Met314" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode == "C18")| 
                                    (pfas_met_tab_all$Met_id == "Met580" & pfas_met_tab_all$PFAS == "PFNA" & pfas_met_tab_all$Mode == "C18"),]
pfas_met_tab<- pfas_met_tab[pfas_met_tab$PFAS_age == "0", ]

dim(pfas_met_tab)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("7","14","22","28"))

table(pfas_met_tab$Age)

pfas0 <- pfas_met_tab

## AGE 7

pfas_met_tab <- pfas_met_tab_all [(pfas_met_tab_all$Met_id == "Met192" & pfas_met_tab_all$PFAS == "PFOS" & pfas_met_tab_all$Mode == "HILIC") | 
                                    (pfas_met_tab_all$Met_id == "Met389" & pfas_met_tab_all$PFAS == "PFOA" & pfas_met_tab_all$Mode == "HILIC") ,]
pfas_met_tab<- pfas_met_tab[pfas_met_tab$PFAS_age == "7", ]

dim(pfas_met_tab)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("14","22","28"))

table(pfas_met_tab$Age)


pfas7 <- pfas_met_tab

## AGE 14

pfas_met_tab <- pfas_met_tab_all [(pfas_met_tab_all$Met_id == "Met519" & pfas_met_tab_all$PFAS == "PFHxS" & pfas_met_tab_all$Mode =="HILIC") | 
                                    (pfas_met_tab_all$Met_id == "Met793" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC") | 
                                    (pfas_met_tab_all$Met_id == "Met438" & pfas_met_tab_all$PFAS == "PFNA" & pfas_met_tab_all$Mode =="C18"),]
pfas_met_tab<- pfas_met_tab[pfas_met_tab$PFAS_age == "14", ]
dim(pfas_met_tab)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("22","28"))

table(pfas_met_tab$Age)

pfas14 <- pfas_met_tab

## AGE 22

pfas_met_tab <- pfas_met_tab_all[(pfas_met_tab_all$Met_id == "Met120" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC") | 
                                   (pfas_met_tab_all$Met_id == "Met190" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC") | 
                                   (pfas_met_tab_all$Met_id == "Met206" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC")|
                                   (pfas_met_tab_all$Met_id == "Met267" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC") |
                                   (pfas_met_tab_all$Met_id == "Met348" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC") | 
                                   (pfas_met_tab_all$Met_id == "Met583" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC") | 
                                   (pfas_met_tab_all$Met_id == "Met676" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC")|
                                   (pfas_met_tab_all$Met_id == "Met677" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC")|
                                   (pfas_met_tab_all$Met_id == "Met728" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC") | 
                                   (pfas_met_tab_all$Met_id == "Met780" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC") | 
                                   (pfas_met_tab_all$Met_id == "Met610" & pfas_met_tab_all$PFAS == "PFOS" & pfas_met_tab_all$Mode =="HILIC")|
                                   (pfas_met_tab_all$Met_id == "Met742" & pfas_met_tab_all$PFAS == "PFOS" & pfas_met_tab_all$Mode =="HILIC")|
                                   (pfas_met_tab_all$Met_id == "Met79" & pfas_met_tab_all$PFAS == "PFOA" & pfas_met_tab_all$Mode =="HILIC")|
                                   (pfas_met_tab_all$Met_id == "Met86" & pfas_met_tab_all$PFAS == "PFOA" & pfas_met_tab_all$Mode =="HILIC") | 
                                   (pfas_met_tab_all$Met_id == "Met886" & pfas_met_tab_all$PFAS == "PFOA" & pfas_met_tab_all$Mode =="HILIC") | 
                                   (pfas_met_tab_all$Met_id == "Met977" & pfas_met_tab_all$PFAS == "PFOA" & pfas_met_tab_all$Mode =="HILIC")|
                                   (pfas_met_tab_all$Met_id == "Met452" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="C18"),]
pfas_met_tab<- pfas_met_tab[pfas_met_tab$PFAS_age == "22", ]

dim(pfas_met_tab)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("28"))

table(pfas_met_tab$Age)

pfas22 <- pfas_met_tab

all_sig_hits <- rbind(pfas0,pfas7,pfas14,pfas22)


# write.csv(all_sig_hits, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/sig_metabolites_longitudinal.csv", row.names = F )


#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

##### Figure

all_sig_hits <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/sig_metabolites_longitudinal.csv")
all_sig_hits$Age <- factor(all_sig_hits$Age,
                           levels = c(7, 14, 22, 28))

all_sig_hits$PFAS<-factor(all_sig_hits$PFAS,
                          levels=c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))

all_sig_hits$sig<-factor(ifelse(all_sig_hits$rand_adj_pval<0.2, 1, 0),
                         levels = c(0,1))

## AGE 0
vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 0,], aes(x=Age, y=beta, color=PFAS, label=Metabolite, shape = Mode, fill = sig)) +# Show all points
          geom_point(size = 2)   
        +  scale_shape_manual(values = c(15,17))
        + scale_fill_manual(values = c("white","yellow"))
        + scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "darkgreen", 
                                         "black",   "mediumpurple1"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))
        + geom_hline(yintercept= 0, color = "black", size = 1, ) + 
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients", 
               title = "Prenatal PFAS exposure and associated metabolites") + theme_bw()
        +  geom_label_repel(size = 7, family = 'serif',
                            fontface = 'italic',
                            box.padding = unit(0.5, "lines"),
                            max.iter = 2e4,
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                            force = 2, force_pull = 2, show.legend = F))

vol <-  (vol + theme(plot.title=element_text(size=16,face="bold"),
                     axis.title=element_text(size=14,face="bold"),
                     plot.tag = element_text(size = 14,face = "bold"),
                     axis.text.y = element_text(size=14,face="bold"),
                     axis.text.x = element_text(size=14,face="bold"),
                     strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
                     legend.title = element_text(face = "bold"), legend.position = 'bottom',
                     legend.text = element_text(size = 12),
                     legend.background = element_rect(fill="white", 
                                                      size=0.5, linetype="solid",  colour ="darkblue"),
                     plot.margin=unit(c(0.5,0.5,1,0.5), "cm")) 
         + ylim(c(-1,1))
         + annotate("rect", xmin = c(0.75, 1.75, 2.75, 3.75), 
                    xmax = c(1.25, 2.25, 3.25, 4.25), 
                    ymin = c(-1, -1, -1, -1), ymax =rep(1,4),
                    alpha = 0.05)
         + guides(color = guide_legend(override.aes = list(size = 4)),
                  shape = guide_legend(override.aes = list(size = 4))))
vol_0 <- vol
vol_0

#################################################################################
## AGE 7

vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 7,], aes(x=Age, y=beta, color=PFAS, label=Metabolite, shape = Mode, fill = sig)) +# Show all points
          geom_point(size = 2)   +  scale_shape_manual(values = c(15,17))
        + scale_fill_manual(values = c("white","yellow"))
        + scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "darkgreen", 
                                         "black",   "mediumpurple1"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))
        + geom_hline(yintercept= 0, color = "black", size = 1, ) + 
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients", 
               title = "PFAS exposure at age 7 and associated metabolites") + theme_bw()
        +  geom_label_repel(size = 7, family = 'serif',
                            fontface = 'italic',
                            box.padding = unit(0.5, "lines"),
                            max.iter = 1e4,
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                            force = 2, force_pull = 2, show.legend = F))

vol <-  (vol + theme(plot.title=element_text(size=16,face="bold"),
                     axis.title=element_text(size=14,face="bold"),
                     plot.tag = element_text(size = 14,face = "bold"),
                     axis.text.y = element_text(size=14,face="bold"),
                     axis.text.x = element_text(size=14,face="bold"),
                     strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
                     legend.title = element_text(face = "bold"), legend.position = 'bottom',
                     legend.text = element_text(size = 12),
                     legend.background = element_rect(fill="white", 
                                                      size=0.5, linetype="solid",  colour ="darkblue"),
                     plot.margin=unit(c(0.5,0.5,1,0.5), "cm")) 
         + ylim(c(-1,1.3))
         + annotate("rect", xmin = c(0.75, 1.75, 2.75), 
                    xmax = c(1.25, 2.25, 3.25), 
                    ymin = c(-1, -1, -1), 
                    ymax =rep(1.2,3),
                    alpha = 0.05)
         + guides(color = guide_legend(override.aes = list(size = 4)),
                  shape = guide_legend(override.aes = list(size = 4))))
vol_7 <- vol
vol_7

#################################################################################

## AGE 14

vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 14,], aes(x=Age, y=beta, color=PFAS, label=Metabolite, shape = Mode, fill = sig)) +# Show all points
          geom_point(size = 2)   +  scale_shape_manual(values = c(15,17))
        + scale_fill_manual(values = c("white","yellow"))
        + scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "darkgreen", 
                                         "black",   "mediumpurple1"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))
        + geom_hline(yintercept= 0, color = "black", size = 1, ) + 
          
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients", 
               title = "PFAS exposure at age 14 and associated metabolites") + theme_bw()
        +  geom_label_repel(size = 7, family = 'serif',
                            fontface = 'italic',
                            box.padding = unit(0.5, "lines"),
                            max.iter = 2e4,
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                            force = 2, force_pull = 2, show.legend = F))

vol <-  (vol + theme(plot.title=element_text(size=16,face="bold"),
                     axis.title=element_text(size=14,face="bold"),
                     plot.tag = element_text(size = 14,face = "bold"),
                     axis.text.y = element_text(size=14,face="bold"),
                     axis.text.x = element_text(size=14,face="bold"),
                     strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
                     legend.title = element_text(face = "bold"), legend.position = 'bottom',
                     legend.text = element_text(size = 12),
                     legend.background = element_rect(fill="white", 
                                                      size=0.5, linetype="solid",  colour ="darkblue"),
                     plot.margin=unit(c(0.5,0.5,1,0.5), "cm")) 
         + ylim(c(-1,1))
         + annotate("rect", xmin = c(0.75, 1.75), 
                    xmax = c(1.25, 2.25), 
                    ymin = c(-1, -1), 
                    ymax =rep(1,2),
                    alpha = 0.05)
         + guides(color = guide_legend(override.aes = list(size = 4)),
                  shape = guide_legend(override.aes = list(size = 4))))
vol_14 <- vol
vol_14

#################################################################################

## AGE 22

vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 22,], aes(x=Age, y=beta, color=PFAS, label=Metabolite, shape = Mode, fill = sig)) +# Show all points
          geom_point(size = 2)   +  scale_shape_manual(values = c(15,17))
        + scale_fill_manual(values = c("yellow"))
        + scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "darkgreen", 
                                         "black",   "mediumpurple1"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))
        + geom_hline(yintercept= 0, color = "black", size = 1, ) + 
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients", 
               title = "PFAS exposure at age 22 and associated metabolites") + theme_bw()
        +  geom_label_repel(size = 7, family = 'serif',
                            fontface = 'italic',
                            box.padding = unit(0.5, "lines"),
                            max.iter = 2e4,
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                            force = 2, force_pull = 2, show.legend = F))

vol <-  (vol + theme(plot.title=element_text(size=16,face="bold"),
                     axis.title=element_text(size=14,face="bold"),
                     plot.tag = element_text(size = 14,face = "bold"),
                     axis.text.y = element_text(size=14,face="bold"),
                     axis.text.x = element_text(size=14,face="bold"),
                     strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
                     legend.title = element_text(face = "bold"), legend.position = 'bottom',
                     legend.text = element_text(size = 12),
                     legend.background = element_rect(fill="white", 
                                                      size=0.5, linetype="solid",  colour ="darkblue"),
                     plot.margin=unit(c(0.5,0.5,1,0.5), "cm")) 
         + ylim(c(-1,1))
         + annotate("rect", xmin = c(0.75), 
                    xmax = c(1.25), 
                    ymin = c(-1), 
                    ymax =rep(1,1),
                    alpha = 0.05)
         + guides(color = guide_legend(override.aes = list(size = 4)),
                  shape = guide_legend(override.aes = list(size = 4))))
vol_22 <- vol
vol_22



# write.csv(df, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/New_faroese/sig_metabolites.csv", row.names = F )


jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/all_modes_all_metabolites_longitudinal_order.jpeg",
     units="in", width=24, height=28, res=600)

layout<- c(
  area(t=1, l=1, b = 1, r = 5),
  area(t=2, l=2, b = 2, r = 5),
  area(t=3, l=3, b = 3, r = 5),
  area(t=4, l=4, b = 4, r = 5)
)

plot1<- vol_0 + vol_7 + vol_14 + vol_22 + 
  plot_layout(design=layout,
              guides='collect',
  )+
  plot_annotation(title="Prospective associations between individual PFAS exposure and metabolic alterations throughout life-course",
                  theme = list(plot.title=element_text(size=20,face="bold"))) &
  theme(legend.position='bottom')
plot1


dev.off()


#############################################################

# plot<- ggpubr::ggarrange(vol_0, vol_7, vol_14, vol_22, nrow = 4, ncol = 1, common.legend = T, 
#                          legend = "bottom")
# annotate_figure(plot, top = text_grob("Prospective associations between individual PFAS exposure and metabolic alterations throughout life-course", 
#                                       color = "black", face = "bold", size = 18))
# 



# dft <- table(pfas_met_tab_28$Metabolite)
# r <- names(dft)[as.numeric(dft) >1]
# lst <- NA_character_
# if(length(r)!= 0){
#   for(i in 1:length(r)){
#     m <-  which.min(pfas_met_tab_28$simu_pval[pfas_met_tab_28$Metabolite %in% r[i]])
#     ml <- length(pfas_met_tab_28$simu_pval[pfas_met_tab_28$Metabolite %in% r[i]])
#     lst <- c(lst, pfas_met_tab_28$Met_id[pfas_met_tab_28$Metabolite %in% r[i]][setdiff(seq(1,ml),m)])
#   }
# }
# lst <- lst[-1]
# pfas_met_tab_28 <- pfas_met_tab_28[!(pfas_met_tab_28$Met_id %in% lst),]

