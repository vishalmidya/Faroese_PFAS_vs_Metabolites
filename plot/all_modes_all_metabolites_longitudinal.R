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

#### Combined plot

###### AGE 0

pfas_met_tab_C18 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfas_0_metabolites_all_C18.txt", header = T)
pfas_met_tab_HILIC <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\HILIC\\pfas_0_metabolites_all_HILIC.txt", header = T)
pfas_met_tab_0 <- rbind(pfas_met_tab_HILIC, pfas_met_tab_C18)
pfas_met_tab_0$PFAS_age <- rep("0", nrow(pfas_met_tab_0))

###### AGE 7

pfas_met_tab_C18 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfas_7_metabolites_all_C18.txt", header = T)
pfas_met_tab_HILIC <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\HILIC\\pfas_7_metabolites_all_HILIC.txt", header = T)
pfas_met_tab_7 <- rbind(pfas_met_tab_HILIC, pfas_met_tab_C18)
pfas_met_tab_7$PFAS_age <- rep("7", nrow(pfas_met_tab_7))

###### AGE 14

pfas_met_tab_C18 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfas_14_metabolites_all_C18.txt", header = T)
pfas_met_tab_HILIC <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\HILIC\\pfas_14_metabolites_all_HILIC.txt", header = T)
pfas_met_tab_14 <- rbind(pfas_met_tab_HILIC, pfas_met_tab_C18)
pfas_met_tab_14$PFAS_age <- rep("14", nrow(pfas_met_tab_14))


###### AGE 22

pfas_met_tab_C18 <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\C18\\pfas_22_metabolites_all_C18.txt", header = T)
pfas_met_tab_HILIC <- read.table("C:\\Users\\yaom03\\OneDrive - The Mount Sinai Hospital\\New_faroese\\HILIC\\pfas_22_metabolites_all_HILIC.txt", header = T)
pfas_met_tab_22 <- rbind(pfas_met_tab_HILIC, pfas_met_tab_C18)
pfas_met_tab_22$PFAS_age <- rep("22", nrow(pfas_met_tab_22))

pfas_met_tab_all <- rbind(pfas_met_tab_0,pfas_met_tab_7,pfas_met_tab_14,pfas_met_tab_22)

for(i in 1:nrow(pfas_met_tab_all)){
  if(pfas_met_tab_all$PFAS[i] == "pfoa"){pfas_met_tab_all$PFAS[i] = "PFOA"}
  else if(pfas_met_tab_all$PFAS[i] == "pfos"){pfas_met_tab_all$PFAS[i] = "PFOS"}
  else if(pfas_met_tab_all$PFAS[i] == "pfna"){pfas_met_tab_all$PFAS[i] = "PFNA"}
  else if(pfas_met_tab_all$PFAS[i] == "pfda"){pfas_met_tab_all$PFAS[i] = "PFDA"}
  else if(pfas_met_tab_all$PFAS[i] == "pfhxs"){pfas_met_tab_all$PFAS[i] = "PFHxS"}
}

##################################################################################################################################################################

## AGE 0
pfas_met_tab <- pfas_met_tab_all [(pfas_met_tab_all$Met_id == "Met87" & pfas_met_tab_all$PFAS == "PFNA") & pfas_met_tab_all$Mode == "C18"| 
                                    (pfas_met_tab_all$Met_id == "Met1269" & pfas_met_tab_all$PFAS == "PFOS" & pfas_met_tab_all$Mode == "HILIC") | 
                                    (pfas_met_tab_all$Met_id == "Met253" & pfas_met_tab_all$PFAS == "PFOS" & pfas_met_tab_all$Mode == "C18"),]
pfas_met_tab<- pfas_met_tab[pfas_met_tab$PFAS_age == "0", ]

dim(pfas_met_tab)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("7","14","22","28"))

table(pfas_met_tab$Age)

pfas_met_tab$Metabolite <- rep(NA_character_,nrow(pfas_met_tab))
for(i in 1:nrow(pfas_met_tab)) pfas_met_tab$Metabolite[i] <- strsplit(pfas_met_tab$chem_name,";")[[i]][1]

pfas_met_tab_7 <- pfas_met_tab[pfas_met_tab$Age == "7",]
pfas_met_tab_14 <- pfas_met_tab[pfas_met_tab$Age == "14",]
pfas_met_tab_22 <- pfas_met_tab[pfas_met_tab$Age == "22",]
pfas_met_tab_28 <- pfas_met_tab[pfas_met_tab$Age == "28",]

pfas_met_tab <- rbind(pfas_met_tab_7,pfas_met_tab_14, pfas_met_tab_22, pfas_met_tab_28)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("7","14","22","28"))
table(pfas_met_tab$Age)

pfas0 <- pfas_met_tab

## AGE 7

pfas_met_tab <- pfas_met_tab_all [(pfas_met_tab_all$Met_id == "Met695" & pfas_met_tab_all$PFAS == "PFHxS" & pfas_met_tab_all$Mode == "C18") | 
                                    (pfas_met_tab_all$Met_id == "Met707" & pfas_met_tab_all$PFAS == "PFHxS" & pfas_met_tab_all$Mode == "C18") | 
                                    (pfas_met_tab_all$Met_id == "Met1055" & pfas_met_tab_all$PFAS == "PFOA" & pfas_met_tab_all$Mode == "HILIC"),]
pfas_met_tab<- pfas_met_tab[pfas_met_tab$PFAS_age == "7", ]

dim(pfas_met_tab)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("14","22","28"))

table(pfas_met_tab$Age)

pfas_met_tab$Metabolite <- rep(NA_character_,nrow(pfas_met_tab))
for(i in 1:nrow(pfas_met_tab)) pfas_met_tab$Metabolite[i] <- strsplit(pfas_met_tab$chem_name,";")[[i]][1]

pfas_met_tab_14 <- pfas_met_tab[pfas_met_tab$Age == "14",]
pfas_met_tab_22 <- pfas_met_tab[pfas_met_tab$Age == "22",]
pfas_met_tab_28 <- pfas_met_tab[pfas_met_tab$Age == "28",]

pfas_met_tab <- rbind(pfas_met_tab_14, pfas_met_tab_22, pfas_met_tab_28)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("14","22","28"))
table(pfas_met_tab$Age)

pfas7 <- pfas_met_tab

## AGE 14

pfas_met_tab <- pfas_met_tab_all [(pfas_met_tab_all$Met_id == "Met1571" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC") | 
                                    (pfas_met_tab_all$Met_id == "Met363" & pfas_met_tab_all$PFAS == "PFHxS" & pfas_met_tab_all$Mode =="HILIC") | 
                                    (pfas_met_tab_all$Met_id == "Met459" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC")|
                                    (pfas_met_tab_all$Met_id == "Met1936" & pfas_met_tab_all$PFAS == "PFOA" & pfas_met_tab_all$Mode =="HILIC") | 
                                    (pfas_met_tab_all$Met_id == "Met407" & pfas_met_tab_all$PFAS == "PFHxS" & pfas_met_tab_all$Mode =="HILIC") | 
                                    (pfas_met_tab_all$Met_id == "Met676" & pfas_met_tab_all$PFAS == "PFNA" & pfas_met_tab_all$Mode =="C18"),]
pfas_met_tab<- pfas_met_tab[pfas_met_tab$PFAS_age == "14", ]
dim(pfas_met_tab)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("22","28"))

table(pfas_met_tab$Age)

pfas_met_tab$Metabolite <- rep(NA_character_,nrow(pfas_met_tab))
for(i in 1:nrow(pfas_met_tab)) pfas_met_tab$Metabolite[i] <- strsplit(pfas_met_tab$chem_name,";")[[i]][1]

pfas_met_tab_22 <- pfas_met_tab[pfas_met_tab$Age == "22",]
pfas_met_tab_28 <- pfas_met_tab[pfas_met_tab$Age == "28",]

pfas_met_tab <- rbind(pfas_met_tab_22, pfas_met_tab_28)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("22","28"))
table(pfas_met_tab$Age)

pfas14 <- pfas_met_tab

## AGE 22

pfas_met_tab <- pfas_met_tab_all[(pfas_met_tab_all$Met_id == "Met1771" & pfas_met_tab_all$PFAS == "PFHxS" & pfas_met_tab_all$Mode =="HILIC") | 
                                   (pfas_met_tab_all$Met_id == "Met389" & pfas_met_tab_all$PFAS == "PFDA" & pfas_met_tab_all$Mode =="HILIC") | 
                                   (pfas_met_tab_all$Met_id == "Met1473" & pfas_met_tab_all$PFAS == "PFNA" & pfas_met_tab_all$Mode =="HILIC")|
                                   (pfas_met_tab_all$Met_id == "Met398" & pfas_met_tab_all$PFAS == "PFOA" & pfas_met_tab_all$Mode =="HILIC"),]
pfas_met_tab<- pfas_met_tab[pfas_met_tab$PFAS_age == "22", ]

dim(pfas_met_tab)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("28"))

table(pfas_met_tab$Age)

pfas_met_tab$Metabolite <- rep(NA_character_,nrow(pfas_met_tab))
for(i in 1:nrow(pfas_met_tab)) pfas_met_tab$Metabolite[i] <- strsplit(pfas_met_tab$chem_name,";")[[i]][1]

pfas_met_tab_28 <- pfas_met_tab[pfas_met_tab$Age == "28",]

pfas_met_tab <- rbind(pfas_met_tab_28)
pfas_met_tab$Age <- factor(pfas_met_tab$Age, levels = c("28"))
table(pfas_met_tab$Age)

pfas22 <- pfas_met_tab

all_sig_hits <- rbind(pfas0,pfas7,pfas14,pfas22)
all_sig_hits <- all_sig_hits[,c("mz","time","beta","PFAS","Age","PFAS_age","Mode","Met_id","simu_pval","rand_adj_pval")]


##################################################################################################################################################################
# Annotations

conf_met <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/confirmed_metabolites.csv",fileEncoding = "Latin1")
conf_met_HILIC <- conf_met[conf_met$Method.RT == "HILIC+",]
conf_met_C18 <- conf_met[conf_met$Method.RT == "C18-",]

stage4_hilic <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/Stage4.csv")
stage4_c18 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/Stage4.csv")

all_sig_hits$chemical_ID <- rep(NA_character_, dim(all_sig_hits)[1])
all_sig_hits$Annotation.confidence.score <- rep(NA_character_, dim(all_sig_hits)[1])
all_sig_hits$Name <- rep(NA_character_, dim(all_sig_hits)[1])
all_sig_hits$Adduct <- rep(NA_character_, dim(all_sig_hits)[1])

for(i in 1:nrow(all_sig_hits)){
  
  if(all_sig_hits$Mode[i] == "HILIC"){
    
    rrt <- which(abs(all_sig_hits$mz[i] - stage4_hilic$mz) < 5 & abs(all_sig_hits$time[i] - stage4_hilic$time) < 30)
    dft <- stage4_hilic[unique(rrt[!is.na(rrt)]),]
    b <- dft[dft$Formula %in% intersect(conf_met_HILIC$Elemental.Composition, unique(dft$Formula)),]
    
    if(nrow(b)!= 0){
      
      score <- max(stage4_hilic$Annotation.confidence.score[stage4_hilic$chemical_ID %in% b$chemical_ID])
      dft2 <- stage4_hilic[stage4_hilic$chemical_ID %in% b$chemical_ID,c("chemical_ID","Annotation.confidence.score","Name","Adduct")]
      all_sig_hits$chemical_ID[i] = knitr::combine_words(unique(dft2$chemical_ID[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      all_sig_hits$Annotation.confidence.score[i] = knitr::combine_words(unique(dft2$Annotation.confidence.score[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      all_sig_hits$Name[i] = knitr::combine_words(unique(dft2$Name[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      
      names <- unlist(strsplit(all_sig_hits$Name[i],"/"))
      tp <- NA_character_
      for(j in 1:length(names)){
        
        tp <- c(tp,unique(dft2$Adduct[dft2$Name == names[j]]))
        
      }
      tp <- tp[-1]
      all_sig_hits$Adduct[i] = knitr::combine_words(tp, and = "", sep = "/")
      
    }
  }
  
  else if(all_sig_hits$Mode[i] == "C18"){
    
    rrt <- which(abs(all_sig_hits$mz[i] - stage4_c18$mz) < 5 & abs(all_sig_hits$time[i] - stage4_c18$time) < 30)
    dft <- stage4_c18[unique(rrt[!is.na(rrt)]),]
    b <- dft[dft$Formula %in% intersect(conf_met_C18$Elemental.Composition, unique(dft$Formula)),]
    
    if(nrow(b)!= 0){
      
      score <- max(stage4_c18$Annotation.confidence.score[stage4_c18$chemical_ID %in% b$chemical_ID])
      dft2 <- stage4_c18[stage4_c18$chemical_ID %in% b$chemical_ID,c("chemical_ID","Annotation.confidence.score","Name","Adduct")]
      all_sig_hits$chemical_ID[i] = knitr::combine_words(unique(dft2$chemical_ID[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      all_sig_hits$Annotation.confidence.score[i] = knitr::combine_words(unique(dft2$Annotation.confidence.score[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      all_sig_hits$Name[i] = knitr::combine_words(unique(dft2$Name[dft2$Annotation.confidence.score == score]), and = "", sep = "/")
      
      names <- unlist(strsplit(all_sig_hits$Name[i],"/"))
      tp <- NA_character_
      for(j in 1:length(names)){
        
        tp <- c(tp,unique(dft2$Adduct[dft2$Name == names[j]]))
        
      }
      tp <- tp[-1]
      all_sig_hits$Adduct[i] = knitr::combine_words(tp, and = "", sep = "/")
      
      
    }
  }
  
}


# write.csv(all_sig_hits, "C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/sig_metabolites_longitudinal.csv", row.names = F )


#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

##### Figure

all_sig_hits <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/sig_metabolites_longitudinal_ryan_updated.csv")
all_sig_hits$Age <- as.factor(all_sig_hits$Age)
all_sig_hits$Name[which(all_sig_hits$Name == "N-Acetylneuraminic acid/N-Acetyl-a-neuraminic acid")] <- "N-Acetylneuraminic acid/\nN-Acetyl-a-neuraminic acid"
all_sig_hits$Name[which(all_sig_hits$Name == "Asymmetric dimethylarginine/Symmetric dimethylarginine (Dimethyl-Arginine)")] <- "Asymmetric dimethylarginine/\nSymmetric dimethylarginine\n(Dimethyl-Arginine)"
all_sig_hits$Name[which(all_sig_hits$Name == "81114-Eicosatrienoic acid (linolenic acid) omega-6")] <- "81114-Eicosatrienoic acid\n(linolenic acid) omega-6"
all_sig_hits$Name[which(all_sig_hits$Name == "Betaine/L-Valine/Vaporole/N-Methyl-a-aminoisobutyric acid/5-Aminopentanoic acid/Norvaline/Amyl Nitrite")] <- "Betaine/L-Valine/Vaporole/\nN-Methyl-a-aminoisobutyric acid/5-Aminopentanoic acid/\nNorvaline/Amyl Nitrite"
all_sig_hits$Name[which(all_sig_hits$Name == "L-Aspartic acid/D-Aspartic acid (Aspartic Acid)")] <- "L-Aspartic acid/\nD-Aspartic acid (Aspartic Acid)"

all_sig_hits$PFAS<-factor(all_sig_hits$PFAS,
                          levels=c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))

all_sig_hits$sig<-ifelse(all_sig_hits$rand_adj_pval<0.2, "*", NA)

## AGE 0
vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 0,], aes(x=Age, y=beta, color=PFAS, label=Name, shape = Mode)) +# Show all points
          geom_point(size = 2)   +  scale_shape_manual(values = c(15,17))
        + scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "green", 
                                         "grey27",   "mediumpurple1"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))
        + geom_hline(yintercept= 0, color = "black", size = 1, ) + 
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients", 
               title = "Prenatal PFAS exposure and associated metabolites") + theme_bw()
        +  geom_label_repel(size = 4, family = 'serif',
                            fontface = 'italic',
                            box.padding = unit(0.5, "lines"),
                            max.iter = 2e4,
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                            force = 2, force_pull = 2, show.legend = F)+
          geom_text(aes(label = sig), nudge_y = 0.01, nudge_x=0.01))

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

vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 7,], aes(x=Age, y=beta, color=PFAS, label=Name, shape = Mode)) +# Show all points
          geom_point(size = 2)   +  scale_shape_manual(values = c(15,17))
        + scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "green", 
                                         "grey27",   "mediumpurple1"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))
        + geom_hline(yintercept= 0, color = "black", size = 1, ) + 
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients", 
               title = "PFAS exposure at age 7 and associated metabolites") + theme_bw()
        +  geom_label_repel(size = 4, family = 'serif',
                            fontface = 'italic',
                            box.padding = unit(0.5, "lines"),
                            max.iter = 1e4,
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                            force = 2, force_pull = 2, show.legend = F)+
          geom_text(aes(label = sig), nudge_y = 0.01, nudge_x=0.01))

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

vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 14,], aes(x=Age, y=beta, color=PFAS, label=Name, shape = Mode)) +# Show all points
          geom_point(size = 2)   +  scale_shape_manual(values = c(15,17))
        + scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "green", 
                                         "grey27",   "mediumpurple1"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))
        + geom_hline(yintercept= 0, color = "black", size = 1, ) + 
          
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients", 
               title = "PFAS exposure at age 14 and associated metabolites") + theme_bw()
        +  geom_label_repel(size = 4, family = 'serif',
                            fontface = 'italic',
                            box.padding = unit(0.5, "lines"),
                            max.iter = 2e4,
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                            force = 2, force_pull = 2, show.legend = F)+
          geom_text(aes(label = sig), nudge_y = 0.01, nudge_x=0.01))

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

vol <- (ggplot(all_sig_hits[all_sig_hits$PFAS_age == 22,], aes(x=Age, y=beta, color=PFAS, label=Name, shape = Mode)) +# Show all points
          geom_point(size = 2)   +  scale_shape_manual(values = c(17))
        + scale_color_manual(drop = FALSE,
                             values = c( "red",  "blue",  "green", 
                                         "grey27",   "mediumpurple1"),
                             labels = c("PFOS" ,"PFOA" ,"PFNA" ,"PFDA" ,"PFHxS"))
        + geom_hline(yintercept= 0, color = "black", size = 1, ) + 
          labs(x = "Age at metabolomic assessment", y = "Beta Coefficients", 
               title = "PFAS exposure at age 22 and associated metabolites") + theme_bw()
        +  geom_label_repel(size = 4, family = 'serif',
                            fontface = 'italic',
                            box.padding = unit(0.5, "lines"),
                            max.iter = 2e4,
                            max.overlaps = getOption("ggrepel.max.overlaps", default = 40),
                            force = 2, force_pull = 2, show.legend = F)+
          geom_text(aes(label = sig), nudge_y = 0.01, nudge_x=0.01))

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


jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/all_modes_all_metabolites_longitudinal.jpeg",
     units="in", width=15, height=20, res=600)


plot<- ggpubr::ggarrange(vol_0, vol_7, vol_14, vol_22, nrow = 4, ncol = 1, common.legend = T, 
                  legend = "bottom")
annotate_figure(plot, top = text_grob("Prospective associations between individual PFAS exposure and metabolic alterations throughout life-course", 
                                      color = "black", face = "bold", size = 18))


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

