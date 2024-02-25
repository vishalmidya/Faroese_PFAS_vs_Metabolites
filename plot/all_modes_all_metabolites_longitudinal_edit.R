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
library(doParallel)
library(patchwork)


##### Figure

all_sig_hits <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/sig_metabolites_longitudinal.csv")
FI<- read.table("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/FI.txt", header = TRUE)
update<- read.table("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/update.txt", header = TRUE)


FI[FI$Met_id=="Met452",]<- update[1,]

all_sig_hits<- all_sig_hits %>% 
               inner_join(FI, by=c("PFAS_age", "Age", "PFAS", "Met_id", "Mode")) 

all_sig_hits$Age <- factor(all_sig_hits$Age,
                           levels = c(7, 14, 22, 28))

all_sig_hits$PFAS<-factor(all_sig_hits$PFAS,
                          levels=c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))


## AGE 0
vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 0,], aes(x=Age, y=beta, color=PFAS, label=Metabolite, shape = Mode)) +# Show all points
          geom_point(size = 1.5, position = position_dodge(width = 0.2)) +  
          geom_errorbar(aes(ymin = lower_FI, ymax = higher_FI), width = 0.03,size=0.8,position = position_dodge(width = 0.2))+
          scale_shape_manual(values = c(15,17))+ 
          scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "#1C8041", 
                                         "grey27",   "#E55709"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))+
          geom_hline(yintercept= 0, color = "black", size = 1, linetype="dashed") + 
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients w/ 95% Fiducial Interval", 
               title = "Prenatal PFAS exposure") + 
          theme_bw() +  
          geom_text_repel(size = 6.5, 
                           position = position_dodge(width = 0.2),
                           family = 'serif',
                           fontface = 'italic',
                           hjust = 1,
                           direction = "y",
                           box.padding = unit(0.5, "lines"),
                           max.iter = 2e4,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                           force = 2, force_pull = 2, show.legend = F)) 


vol <-  (vol + theme(plot.title=element_text(size=18,face="bold"),
                     axis.title=element_text(size=18,face="bold"),
                     plot.tag = element_text(size = 18,face = "bold"),
                     axis.text.y = element_text(size=16,face="bold"),
                     axis.text.x = element_text(size=16,face="bold"),
                     strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
                     legend.title = element_text(face = "bold"), legend.position = 'bottom',
                     legend.text = element_text(size = 18),
                     legend.background = element_rect(fill="white", 
                                                      size=0.5, linetype="solid",  colour ="darkblue"),
                     plot.margin=unit(c(0.5,0.5,1,0.5), "cm")) 
         + ylim(c(-1.2,1.2))
         + annotate("rect", xmin = c(0.75, 1.75, 2.75, 3.75), 
                    xmax = c(1.25, 2.25, 3.25, 4.25), 
                    ymin = rep(-1.2,4), ymax =rep(1.2,4),
                    alpha = 0.05)
         + guides(color = "none",
                  shape = "none"))
vol_0 <- vol
vol_0

#################################################################################
## AGE 7

vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 7,], aes(x=Age, y=beta, color=PFAS, label=Metabolite, shape = Mode, alpha=Met_id)) +# Show all points
          geom_point(size = 1.5, position = position_dodge(width = 0.2)) +  
          geom_errorbar(aes(ymin = lower_FI, ymax = higher_FI), width = 0.03,size=0.8,position = position_dodge(width = 0.2))+
          scale_alpha_manual(values = c(1,1,1))+
          scale_shape_manual(values = c(15,17))+ 
          scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "#1C8041", 
                                         "grey27",   "#E55709"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))+
          geom_hline(yintercept= 0, color = "black", size = 1, linetype="dashed") + 
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients w/ 95% Fiducial Interval", 
               title = "PFAS exposure at age 7") + 
          theme_bw() +  
          geom_text_repel(size = 6.5, 
                           position = position_dodge(width = 0.2),
                           family = 'serif',
                           fontface = 'italic',
                           hjust = 1,
                           direction = "y",
                           box.padding = unit(0.5, "lines"),
                           max.iter = 2e4,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                           force = 2, force_pull = 2, show.legend = F))


vol <-  (vol + theme(plot.title=element_text(size=18,face="bold"),
                     axis.title=element_text(size=18,face="bold"),
                     plot.tag = element_text(size = 18,face = "bold"),
                     axis.text.y = element_text(size=16,face="bold"),
                     axis.text.x = element_text(size=16,face="bold"),
                     strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
                     legend.title = element_text(face = "bold"), legend.position = 'bottom',
                     legend.text = element_text(size = 18),
                     legend.background = element_rect(fill="white", 
                                                      size=0.5, linetype="solid",  colour ="darkblue"),
                     plot.margin=unit(c(0.5,0.5,1,0.5), "cm")) 
         + ylim(c(-1,1.5))
         + annotate("rect", xmin = c(0.75, 1.75, 2.75), 
                    xmax = c(1.25, 2.25, 3.25), 
                    ymin = rep(-1,3), ymax =rep(1.5,3),
                    alpha = 0.05)
         + guides(color = "none",
                  shape = "none",
                  alpha = "none"))

vol_7 <- vol
vol_7

#################################################################################

## AGE 14

vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 14,], aes(x=Age, y=beta, color=PFAS, label=Metabolite, shape = Mode, alpha=Met_id)) +# Show all points
          geom_point(size = 1.5, position = position_dodge(width = 0.2)) +  
          geom_errorbar(aes(ymin = lower_FI, ymax = higher_FI), width = 0.03,size=0.8,position = position_dodge(width = 0.2))+
          scale_alpha_manual(values = rep(1,6))+
          scale_shape_manual(values = c(15,17))+ 
          scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "#1C8041", 
                                         "grey27",   "#E55709"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))+
          geom_hline(yintercept= 0, color = "black", size = 1, linetype="dashed") + 
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients w/ 95% Fiducial Interval", 
               title = "PFAS exposure at age 14") + 
          theme_bw() +  
          geom_text_repel(size =6.5, 
                          position = position_dodge(width = 0.2),
                          family = 'serif',
                          fontface = 'italic',
                          hjust = 1.2,
                          direction = "y",
                          box.padding = unit(0.5, "lines"),
                          max.iter = 2e4,
                          max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
                          force = 2, force_pull = 2, show.legend = F) )


vol <-  (vol + theme(plot.title=element_text(size=18,face="bold"),
                     axis.title=element_text(size=18,face="bold"),
                     plot.tag = element_text(size = 18,face = "bold"),
                     axis.text.y = element_text(size=16,face="bold"),
                     axis.text.x = element_text(size=16,face="bold"),
                     strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
                     legend.title = element_text(face = "bold"), legend.position = 'bottom',
                     legend.text = element_text(size = 18),
                     legend.background = element_rect(fill="white", 
                                                      size=0.5, linetype="solid",  colour ="darkblue"),
                     plot.margin=unit(c(0.5,0.5,1,0.5), "cm")) 
         + ylim(c(-1.2,0.5))
         + annotate("rect", xmin = c(0.75, 1.75), 
                    xmax = c(1.25, 2.25), 
                    ymin = rep(-1.2,2), ymax =rep(0.5,2),
                    alpha = 0.05)
         + guides(color = guide_legend(override.aes = list(size = 4)),
                  shape = guide_legend(override.aes = list(size = 4)),
                  alpha = "none"))

vol_14 <- vol
vol_14

#################################################################################

## AGE 22

vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 22,], aes(x=Age, y=beta, color=PFAS, label=Metabolite, shape = Mode, alpha=Met_id)) +# Show all points
          geom_point(size = 1.5, position = position_dodge(width = 0.2)) +  
          geom_errorbar(aes(ymin = lower_FI, ymax = higher_FI), width = 0.03,size=0.8,position = position_dodge(width = 0.2))+
          scale_alpha_manual(values = rep(1,17))+
          scale_shape_manual(values = c(15,17))+ 
          scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "#1C8041", 
                                         "grey27",   "#E55709"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))+
          geom_hline(yintercept= 0, color = "black", size = 1, linetype="dashed") + 
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients w/ 95% Fiducial Interval", 
               title = "PFAS exposure at age 22") + 
          theme_bw() +  
          geom_text_repel(size = 6.5, 
                           family = 'serif',
                           fontface = 'italic',
                           box.padding = unit(1.5, "lines"),
                           max.iter = 2e4,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                           force = 2, force_pull = 2, show.legend = F))



vol <-  (vol + theme(plot.title=element_text(size=18,face="bold"),
                     axis.title=element_text(size=18,face="bold"),
                     plot.tag = element_text(size = 18,face = "bold"),
                     axis.text.y = element_text(size=16,face="bold"),
                     axis.text.x = element_text(size=16,face="bold"),
                     strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
                     legend.title = element_text(face = "bold"), legend.position = 'bottom',
                     legend.text = element_text(18),
                     legend.background = element_rect(fill="white", 
                                                      size=0.5, linetype="solid",  colour ="darkblue"),
                     plot.margin=unit(c(0.5,0.5,1,0.5), "cm")) 
         + ylim(c(-1.25,0.5))
         + annotate("rect", xmin = c(0.75), 
                    xmax = c(1.25), 
                    ymin = rep(-1.25,1), ymax =rep(0.5,1),
                    alpha = 0.05)
         + guides(color = "none",
                  shape = "none",
                  alpha = "none"))

vol_22 <- vol
vol_22


# write.csv(df, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/New_faroese/sig_metabolites.csv", row.names = F )


jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/all_modes_all_metabolites_longitudinal_edit.jpeg",
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

