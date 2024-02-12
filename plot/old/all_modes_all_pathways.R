library(ggplot2)
library(MatchIt)
library(ggrepel)
library(grid)
library(pBrackets)

# HILIC

###########################################################################################################################################################

## PFAS at age 0
### pfoa

em_pathways_pfoa0_meta7  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_7\\mummichog_pathway_enrichment.csv")
p7 <- em_pathways_pfoa0_meta7[em_pathways_pfoa0_meta7$Hits.sig >3 ,]
p7$age <- rep(7, nrow(p7))

em_pathways_pfoa0_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfoa0_meta14[em_pathways_pfoa0_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfoa0_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfoa0_meta22[em_pathways_pfoa0_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfoa0_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfoa0_meta28[em_pathways_pfoa0_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfoa_pathway <- rbind(p7,p14,p22,p28)
pfoa_pathway$age <- factor(pfoa_pathway$age, levels = c(7,14,22,28)) 



### pfos

em_pathways_pfos0_meta7  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_7\\mummichog_pathway_enrichment.csv")
p7 <- em_pathways_pfos0_meta7[em_pathways_pfos0_meta7$Hits.sig >3 ,]
p7$age <- rep(7, nrow(p7))

em_pathways_pfos0_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfos0_meta14[em_pathways_pfos0_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfos0_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfos0_meta22[em_pathways_pfos0_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfos0_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfos0_meta28[em_pathways_pfos0_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfos_pathway <- rbind(p7,p14,p22,p28)
pfos_pathway$age <- factor(pfos_pathway$age, levels = c(7,14,22,28)) 



### pfna

em_pathways_pfna0_meta7  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_7\\mummichog_pathway_enrichment.csv")
p7 <- em_pathways_pfna0_meta7[em_pathways_pfna0_meta7$Hits.sig >3 ,]
p7$age <- rep(7, nrow(p7))

em_pathways_pfna0_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfna0_meta14[em_pathways_pfna0_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfna0_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfna0_meta22[em_pathways_pfna0_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfna0_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfna0_meta28[em_pathways_pfna0_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfna_pathway <- rbind(p7,p14,p22,p28)
pfna_pathway$age <- factor(pfna_pathway$age, levels = c(7,14,22,28)) 


### pfda

em_pathways_pfda0_meta7  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_7\\mummichog_pathway_enrichment.csv")
p7 <- em_pathways_pfda0_meta7[em_pathways_pfda0_meta7$Hits.sig >3 ,]
p7$age <- rep(7, nrow(p7))

em_pathways_pfda0_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfda0_meta14[em_pathways_pfda0_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfda0_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfda0_meta22[em_pathways_pfda0_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfda0_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfda0_meta28[em_pathways_pfda0_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfda_pathway <- rbind(p7,p14,p22,p28)
pfda_pathway$age <- factor(pfda_pathway$age, levels = c(7,14,22,28)) 



### pfhxs

em_pathways_pfhxs0_meta7  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_7\\mummichog_pathway_enrichment.csv")
p7 <- em_pathways_pfhxs0_meta7[em_pathways_pfhxs0_meta7$Hits.sig >3 ,]
p7$age <- rep(7, nrow(p7))

em_pathways_pfhxs0_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfhxs0_meta14[em_pathways_pfhxs0_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfhxs0_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfhxs0_meta22[em_pathways_pfhxs0_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfhxs0_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfhxs0_meta28[em_pathways_pfhxs0_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfhxs_pathway <- rbind(p7,p14,p22,p28)
pfhxs_pathway$age <- factor(pfhxs_pathway$age, levels = c(7,14,22,28)) 


pfas_pathway <- rbind(pfoa_pathway,pfos_pathway,pfna_pathway,pfda_pathway,pfhxs_pathway)
pfas_pathway$PFAS <- c(rep("PFOA",nrow(pfoa_pathway)),rep("PFOS",nrow(pfos_pathway)),rep("PFNA",nrow(pfna_pathway)),
                       rep("PFDA",nrow(pfda_pathway)),rep("PFHxS",nrow(pfhxs_pathway)))
colnames(pfas_pathway)[1] <- "Pathways"

pfas_pathway_ex0 <- pfas_pathway
pfas_pathway_ex0$age <- paste0("0_",pfas_pathway_ex0$age) 


###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

## PFAS at age 7
### pfoa

em_pathways_pfoa7_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfoa\\pfoa_7\\minerva_data_pfoa_7_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfoa7_meta14[em_pathways_pfoa7_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfoa7_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfoa\\pfoa_7\\minerva_data_pfoa_7_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfoa7_meta22[em_pathways_pfoa7_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfoa7_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfoa\\pfoa_7\\minerva_data_pfoa_7_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfoa7_meta28[em_pathways_pfoa7_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfoa_pathway <- rbind(p14,p22,p28)
pfoa_pathway$age <- factor(pfoa_pathway$age, levels = c(14,22,28)) 


### pfos

em_pathways_pfos7_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfos\\pfos_7\\minerva_data_pfos_7_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfos7_meta14[em_pathways_pfos7_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfos7_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfos\\pfos_7\\minerva_data_pfos_7_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfos7_meta22[em_pathways_pfos7_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfos7_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfos\\pfos_7\\minerva_data_pfos_7_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfos7_meta28[em_pathways_pfos7_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfos_pathway <- rbind(p14,p22,p28)
pfos_pathway$age <- factor(pfos_pathway$age, levels = c(14,22,28)) 



### pfna

em_pathways_pfna7_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfna\\pfna_7\\minerva_data_pfna_7_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfna7_meta14[em_pathways_pfna7_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfna7_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfna\\pfna_7\\minerva_data_pfna_7_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfna7_meta22[em_pathways_pfna7_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfna7_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfna\\pfna_7\\minerva_data_pfna_7_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfna7_meta28[em_pathways_pfna7_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfna_pathway <- rbind(p14,p22,p28)
pfna_pathway$age <- factor(pfna_pathway$age, levels = c(14,22,28)) 



### pfda

em_pathways_pfda7_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfda\\pfda_7\\minerva_data_pfda_7_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfda7_meta14[em_pathways_pfda7_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfda7_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfda\\pfda_7\\minerva_data_pfda_7_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfda7_meta22[em_pathways_pfda7_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfda7_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfda\\pfda_7\\minerva_data_pfda_7_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfda7_meta28[em_pathways_pfda7_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfda_pathway <- rbind(p14,p22,p28)
pfda_pathway$age <- factor(pfda_pathway$age, levels = c(14,22,28)) 


### pfhxs

em_pathways_pfhxs7_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\HILIC\\pfhxs\\pfhxs_7\\minerva_data_pfhxs_7_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfhxs7_meta14[em_pathways_pfhxs7_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfhxs7_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfhxs\\pfhxs_7\\minerva_data_pfhxs_7_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfhxs7_meta22[em_pathways_pfhxs7_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfhxs7_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfhxs\\pfhxs_7\\minerva_data_pfhxs_7_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfhxs7_meta28[em_pathways_pfhxs7_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfhxs_pathway <- rbind(p14,p22,p28)
pfhxs_pathway$age <- factor(pfhxs_pathway$age, levels = c(14,22,28)) 


pfas_pathway <- rbind(pfoa_pathway,pfos_pathway,pfna_pathway,pfda_pathway,pfhxs_pathway)
pfas_pathway$PFAS <- c(rep("PFOA",nrow(pfoa_pathway)),rep("PFOS",nrow(pfos_pathway)),rep("PFNA",nrow(pfna_pathway)),
                       rep("PFDA",nrow(pfda_pathway)),rep("PFHxS",nrow(pfhxs_pathway)))
colnames(pfas_pathway)[1] <- "Pathways"

pfas_pathway_ex7 <- pfas_pathway
pfas_pathway_ex7$age <- paste0("7_",pfas_pathway_ex7$age) 


###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

## PFAS at age 14
### pfoa

em_pathways_pfoa14_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfoa\\pfoa_14\\minerva_data_pfoa_14_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfoa14_meta22[em_pathways_pfoa14_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfoa14_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfoa\\pfoa_14\\minerva_data_pfoa_14_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfoa14_meta28[em_pathways_pfoa14_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfoa_pathway <- rbind(p22,p28)
pfoa_pathway$age <- factor(pfoa_pathway$age, levels = c(22,28)) 


### pfos

em_pathways_pfos14_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfos\\pfos_14\\minerva_data_pfos_14_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfos14_meta22[em_pathways_pfos14_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfos14_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfos\\pfos_14\\minerva_data_pfos_14_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfos14_meta28[em_pathways_pfos14_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfos_pathway <- rbind(p22,p28)
pfos_pathway$age <- factor(pfos_pathway$age, levels = c(22,28)) 



### pfna

em_pathways_pfna14_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfna\\pfna_14\\minerva_data_pfna_14_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfna14_meta22[em_pathways_pfna14_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfna14_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfna\\pfna_14\\minerva_data_pfna_14_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfna14_meta28[em_pathways_pfna14_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfna_pathway <- rbind(p22,p28)
pfna_pathway$age <- factor(pfna_pathway$age, levels = c(22,28)) 



### pfda

em_pathways_pfda14_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfda\\pfda_14\\minerva_data_pfda_14_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfda14_meta22[em_pathways_pfda14_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfda14_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfda\\pfda_14\\minerva_data_pfda_14_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfda14_meta28[em_pathways_pfda14_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfda_pathway <- rbind(p22,p28)
pfda_pathway$age <- factor(pfda_pathway$age, levels = c(22,28)) 



### pfhxs

em_pathways_pfhxs14_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfhxs\\pfhxs_14\\minerva_data_pfhxs_14_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfhxs14_meta22[em_pathways_pfhxs14_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfhxs14_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfhxs\\pfhxs_14\\minerva_data_pfhxs_14_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfhxs14_meta28[em_pathways_pfhxs14_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfhxs_pathway <- rbind(p22,p28)
pfhxs_pathway$age <- factor(pfhxs_pathway$age, levels = c(22,28)) 

pfas_pathway <- rbind(pfoa_pathway,pfos_pathway,pfna_pathway,pfda_pathway,pfhxs_pathway)
pfas_pathway$PFAS <- c(rep("PFOA",nrow(pfoa_pathway)),rep("PFOS",nrow(pfos_pathway)),rep("PFNA",nrow(pfna_pathway)),
                       rep("PFDA",nrow(pfda_pathway)),rep("PFHxS",nrow(pfhxs_pathway)))
colnames(pfas_pathway)[1] <- "Pathways"

pfas_pathway_ex14 <- pfas_pathway
pfas_pathway_ex14$age <- paste0("14_",pfas_pathway_ex14$age) 

###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

## PFAS at age 22
### pfoa

em_pathways_pfoa22_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfoa\\pfoa_22\\minerva_data_pfoa_22_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfoa22_meta28[em_pathways_pfoa22_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfoa_pathway <- rbind(p28)
pfoa_pathway$age <- factor(pfoa_pathway$age, levels = c(28)) 


### pfos

em_pathways_pfos22_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfos\\pfos_22\\minerva_data_pfos_22_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfos22_meta28[em_pathways_pfos22_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfos_pathway <- rbind(p28)
pfos_pathway$age <- factor(pfos_pathway$age, levels = c(28)) 


### pfna

em_pathways_pfna22_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfna\\pfna_22\\minerva_data_pfna_22_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfna22_meta28[em_pathways_pfna22_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfna_pathway <- rbind(p28)
pfna_pathway$age <- factor(pfna_pathway$age, levels = c(28)) 

### pfda

em_pathways_pfda22_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfda\\pfda_22\\minerva_data_pfda_22_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfda22_meta28[em_pathways_pfda22_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfda_pathway <- rbind(p28)
pfda_pathway$age <- factor(pfda_pathway$age, levels = c(28)) 

### pfhxs

em_pathways_pfhxs22_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\HILIC\\pfhxs\\pfhxs_22\\minerva_data_pfhxs_22_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfhxs22_meta28[em_pathways_pfhxs22_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfhxs_pathway <- rbind(p28)
pfhxs_pathway$age <- factor(pfhxs_pathway$age, levels = c(28)) 

pfas_pathway <- rbind(pfoa_pathway,pfos_pathway,pfna_pathway,pfda_pathway,pfhxs_pathway)
pfas_pathway$PFAS <- c(rep("PFOA",nrow(pfoa_pathway)),rep("PFOS",nrow(pfos_pathway)),rep("PFNA",nrow(pfna_pathway)),
                       rep("PFDA",nrow(pfda_pathway)),rep("PFHxS",nrow(pfhxs_pathway)))
colnames(pfas_pathway)[1] <- "Pathways"

pfas_pathway_ex22 <- pfas_pathway
pfas_pathway_ex22$age <- paste0("22_",pfas_pathway_ex22$age) 


pfas_pathway <- rbind(pfas_pathway_ex0,pfas_pathway_ex7,pfas_pathway_ex14,pfas_pathway_ex22)
pfas_pathway$age <- factor(pfas_pathway$age, levels = c("0_7","0_14","0_22","0_28",
                                                        "7_7","7_14","7_22","7_28",
                                                        "14_14","14_22","14_28",
                                                        "22_22","22_28", "28_28"))

pfas_pathway$Mode <- rep("HILIC",nrow(pfas_pathway))
pfas_pathway_HILIC <- pfas_pathway


################################################################################################################################################################################

# C18
## PFAS at age 0
### pfoa

em_pathways_pfoa0_meta7  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_7\\mummichog_pathway_enrichment.csv")
p7 <- em_pathways_pfoa0_meta7[em_pathways_pfoa0_meta7$Hits.sig >3 ,]
p7$age <- rep(7, nrow(p7))

em_pathways_pfoa0_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfoa0_meta14[em_pathways_pfoa0_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfoa0_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfoa0_meta22[em_pathways_pfoa0_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfoa0_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfoa\\pfoa_0\\minerva_data_pfoa_0_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfoa0_meta28[em_pathways_pfoa0_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfoa_pathway <- rbind(p7,p14,p22,p28)
pfoa_pathway$age <- factor(pfoa_pathway$age, levels = c(7,14,22,28)) 



### pfos

em_pathways_pfos0_meta7  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_7\\mummichog_pathway_enrichment.csv")
p7 <- em_pathways_pfos0_meta7[em_pathways_pfos0_meta7$Hits.sig >3 ,]
p7$age <- rep(7, nrow(p7))

em_pathways_pfos0_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfos0_meta14[em_pathways_pfos0_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfos0_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfos0_meta22[em_pathways_pfos0_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfos0_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfos\\pfos_0\\minerva_data_pfos_0_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfos0_meta28[em_pathways_pfos0_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfos_pathway <- rbind(p7,p14,p22,p28)
pfos_pathway$age <- factor(pfos_pathway$age, levels = c(7,14,22,28)) 



### pfna

em_pathways_pfna0_meta7  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_7\\mummichog_pathway_enrichment.csv")
p7 <- em_pathways_pfna0_meta7[em_pathways_pfna0_meta7$Hits.sig >3 ,]
p7$age <- rep(7, nrow(p7))

em_pathways_pfna0_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfna0_meta14[em_pathways_pfna0_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfna0_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfna0_meta22[em_pathways_pfna0_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfna0_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfna\\pfna_0\\minerva_data_pfna_0_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfna0_meta28[em_pathways_pfna0_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfna_pathway <- rbind(p7,p14,p22,p28)
pfna_pathway$age <- factor(pfna_pathway$age, levels = c(7,14,22,28)) 


### pfda

em_pathways_pfda0_meta7  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_7\\mummichog_pathway_enrichment.csv")
p7 <- em_pathways_pfda0_meta7[em_pathways_pfda0_meta7$Hits.sig >3 ,]
p7$age <- rep(7, nrow(p7))

em_pathways_pfda0_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfda0_meta14[em_pathways_pfda0_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfda0_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfda0_meta22[em_pathways_pfda0_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfda0_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfda\\pfda_0\\minerva_data_pfda_0_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfda0_meta28[em_pathways_pfda0_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfda_pathway <- rbind(p7,p14,p22,p28)
pfda_pathway$age <- factor(pfda_pathway$age, levels = c(7,14,22,28)) 



### pfhxs

em_pathways_pfhxs0_meta7  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_7\\mummichog_pathway_enrichment.csv")
p7 <- em_pathways_pfhxs0_meta7[em_pathways_pfhxs0_meta7$Hits.sig >3 ,]
p7$age <- rep(7, nrow(p7))

em_pathways_pfhxs0_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfhxs0_meta14[em_pathways_pfhxs0_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfhxs0_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfhxs0_meta22[em_pathways_pfhxs0_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfhxs0_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfhxs\\pfhxs_0\\minerva_data_pfhxs_0_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfhxs0_meta28[em_pathways_pfhxs0_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfhxs_pathway <- rbind(p7,p14,p22,p28)
pfhxs_pathway$age <- factor(pfhxs_pathway$age, levels = c(7,14,22,28)) 


pfas_pathway <- rbind(pfoa_pathway,pfos_pathway,pfna_pathway,pfda_pathway,pfhxs_pathway)
pfas_pathway$PFAS <- c(rep("PFOA",nrow(pfoa_pathway)),rep("PFOS",nrow(pfos_pathway)),rep("PFNA",nrow(pfna_pathway)),
                       rep("PFDA",nrow(pfda_pathway)),rep("PFHxS",nrow(pfhxs_pathway)))
colnames(pfas_pathway)[1] <- "Pathways"

pfas_pathway_ex0 <- pfas_pathway
pfas_pathway_ex0$age <- paste0("0_",pfas_pathway_ex0$age) 


###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

## PFAS at age 7
### pfoa

em_pathways_pfoa7_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfoa\\pfoa_7\\minerva_data_pfoa_7_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfoa7_meta14[em_pathways_pfoa7_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfoa7_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfoa\\pfoa_7\\minerva_data_pfoa_7_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfoa7_meta22[em_pathways_pfoa7_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfoa7_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfoa\\pfoa_7\\minerva_data_pfoa_7_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfoa7_meta28[em_pathways_pfoa7_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfoa_pathway <- rbind(p14,p22,p28)
pfoa_pathway$age <- factor(pfoa_pathway$age, levels = c(14,22,28)) 


### pfos

em_pathways_pfos7_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfos\\pfos_7\\minerva_data_pfos_7_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfos7_meta14[em_pathways_pfos7_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfos7_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfos\\pfos_7\\minerva_data_pfos_7_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfos7_meta22[em_pathways_pfos7_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfos7_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfos\\pfos_7\\minerva_data_pfos_7_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfos7_meta28[em_pathways_pfos7_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfos_pathway <- rbind(p14,p22,p28)
pfos_pathway$age <- factor(pfos_pathway$age, levels = c(14,22,28)) 



### pfna

em_pathways_pfna7_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfna\\pfna_7\\minerva_data_pfna_7_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfna7_meta14[em_pathways_pfna7_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfna7_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfna\\pfna_7\\minerva_data_pfna_7_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfna7_meta22[em_pathways_pfna7_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfna7_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfna\\pfna_7\\minerva_data_pfna_7_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfna7_meta28[em_pathways_pfna7_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfna_pathway <- rbind(p14,p22,p28)
pfna_pathway$age <- factor(pfna_pathway$age, levels = c(14,22,28)) 



### pfda

em_pathways_pfda7_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfda\\pfda_7\\minerva_data_pfda_7_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfda7_meta14[em_pathways_pfda7_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfda7_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfda\\pfda_7\\minerva_data_pfda_7_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfda7_meta22[em_pathways_pfda7_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfda7_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfda\\pfda_7\\minerva_data_pfda_7_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfda7_meta28[em_pathways_pfda7_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfda_pathway <- rbind(p14,p22,p28)
pfda_pathway$age <- factor(pfda_pathway$age, levels = c(14,22,28)) 


### pfhxs

em_pathways_pfhxs7_meta14  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\New_faroese\\C18\\pfhxs\\pfhxs_7\\minerva_data_pfhxs_7_metabolites_14\\mummichog_pathway_enrichment.csv")
p14 <- em_pathways_pfhxs7_meta14[em_pathways_pfhxs7_meta14$Hits.sig >3 ,]
p14$age <- rep(14, nrow(p14))

em_pathways_pfhxs7_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfhxs\\pfhxs_7\\minerva_data_pfhxs_7_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfhxs7_meta22[em_pathways_pfhxs7_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfhxs7_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfhxs\\pfhxs_7\\minerva_data_pfhxs_7_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfhxs7_meta28[em_pathways_pfhxs7_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfhxs_pathway <- rbind(p14,p22,p28)
pfhxs_pathway$age <- factor(pfhxs_pathway$age, levels = c(14,22,28)) 


pfas_pathway <- rbind(pfoa_pathway,pfos_pathway,pfna_pathway,pfda_pathway,pfhxs_pathway)
pfas_pathway$PFAS <- c(rep("PFOA",nrow(pfoa_pathway)),rep("PFOS",nrow(pfos_pathway)),rep("PFNA",nrow(pfna_pathway)),
                       rep("PFDA",nrow(pfda_pathway)),rep("PFHxS",nrow(pfhxs_pathway)))
colnames(pfas_pathway)[1] <- "Pathways"

pfas_pathway_ex7 <- pfas_pathway
pfas_pathway_ex7$age <- paste0("7_",pfas_pathway_ex7$age) 


###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

## PFAS at age 14
### pfoa

em_pathways_pfoa14_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfoa\\pfoa_14\\minerva_data_pfoa_14_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfoa14_meta22[em_pathways_pfoa14_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfoa14_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfoa\\pfoa_14\\minerva_data_pfoa_14_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfoa14_meta28[em_pathways_pfoa14_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfoa_pathway <- rbind(p22,p28)
pfoa_pathway$age <- factor(pfoa_pathway$age, levels = c(22,28)) 


### pfos

em_pathways_pfos14_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfos\\pfos_14\\minerva_data_pfos_14_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfos14_meta22[em_pathways_pfos14_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfos14_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfos\\pfos_14\\minerva_data_pfos_14_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfos14_meta28[em_pathways_pfos14_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfos_pathway <- rbind(p22,p28)
pfos_pathway$age <- factor(pfos_pathway$age, levels = c(22,28)) 



### pfna

em_pathways_pfna14_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfna\\pfna_14\\minerva_data_pfna_14_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfna14_meta22[em_pathways_pfna14_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfna14_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfna\\pfna_14\\minerva_data_pfna_14_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfna14_meta28[em_pathways_pfna14_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfna_pathway <- rbind(p22,p28)
pfna_pathway$age <- factor(pfna_pathway$age, levels = c(22,28)) 



### pfda

em_pathways_pfda14_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfda\\pfda_14\\minerva_data_pfda_14_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfda14_meta22[em_pathways_pfda14_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfda14_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfda\\pfda_14\\minerva_data_pfda_14_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfda14_meta28[em_pathways_pfda14_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfda_pathway <- rbind(p22,p28)
pfda_pathway$age <- factor(pfda_pathway$age, levels = c(22,28)) 



### pfhxs

em_pathways_pfhxs14_meta22  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfhxs\\pfhxs_14\\minerva_data_pfhxs_14_metabolites_22\\mummichog_pathway_enrichment.csv")
p22 <- em_pathways_pfhxs14_meta22[em_pathways_pfhxs14_meta22$Hits.sig >3 ,]
p22$age <- rep(22, nrow(p22))

em_pathways_pfhxs14_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfhxs\\pfhxs_14\\minerva_data_pfhxs_14_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfhxs14_meta28[em_pathways_pfhxs14_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfhxs_pathway <- rbind(p22,p28)
pfhxs_pathway$age <- factor(pfhxs_pathway$age, levels = c(22,28)) 

pfas_pathway <- rbind(pfoa_pathway,pfos_pathway,pfna_pathway,pfda_pathway,pfhxs_pathway)
pfas_pathway$PFAS <- c(rep("PFOA",nrow(pfoa_pathway)),rep("PFOS",nrow(pfos_pathway)),rep("PFNA",nrow(pfna_pathway)),
                       rep("PFDA",nrow(pfda_pathway)),rep("PFHxS",nrow(pfhxs_pathway)))
colnames(pfas_pathway)[1] <- "Pathways"

pfas_pathway_ex14 <- pfas_pathway
pfas_pathway_ex14$age <- paste0("14_",pfas_pathway_ex14$age) 

###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

## PFAS at age 22
### pfoa

em_pathways_pfoa22_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfoa\\pfoa_22\\minerva_data_pfoa_22_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfoa22_meta28[em_pathways_pfoa22_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfoa_pathway <- rbind(p28)
pfoa_pathway$age <- factor(pfoa_pathway$age, levels = c(28)) 


### pfos

em_pathways_pfos22_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfos\\pfos_22\\minerva_data_pfos_22_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfos22_meta28[em_pathways_pfos22_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfos_pathway <- rbind(p28)
pfos_pathway$age <- factor(pfos_pathway$age, levels = c(28)) 


### pfna

em_pathways_pfna22_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfna\\pfna_22\\minerva_data_pfna_22_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfna22_meta28[em_pathways_pfna22_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfna_pathway <- rbind(p28)
pfna_pathway$age <- factor(pfna_pathway$age, levels = c(28)) 

### pfda

em_pathways_pfda22_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfda\\pfda_22\\minerva_data_pfda_22_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfda22_meta28[em_pathways_pfda22_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfda_pathway <- rbind(p28)
pfda_pathway$age <- factor(pfda_pathway$age, levels = c(28)) 

### pfhxs

em_pathways_pfhxs22_meta28  <- read.csv("C:\\Users\\midyav01\\OneDrive - The Mount Sinai Hospital\\MSSM Projects\\METABOLOMICS\\\\New_faroese\\C18\\pfhxs\\pfhxs_22\\minerva_data_pfhxs_22_metabolites_28\\mummichog_pathway_enrichment.csv")
p28 <- em_pathways_pfhxs22_meta28[em_pathways_pfhxs22_meta28$Hits.sig >3 ,]
p28$age <- rep(28, nrow(p28))

pfhxs_pathway <- rbind(p28)
pfhxs_pathway$age <- factor(pfhxs_pathway$age, levels = c(28)) 

pfas_pathway <- rbind(pfoa_pathway,pfos_pathway,pfna_pathway,pfda_pathway,pfhxs_pathway)
pfas_pathway$PFAS <- c(rep("PFOA",nrow(pfoa_pathway)),rep("PFOS",nrow(pfos_pathway)),rep("PFNA",nrow(pfna_pathway)),
                       rep("PFDA",nrow(pfda_pathway)),rep("PFHxS",nrow(pfhxs_pathway)))
colnames(pfas_pathway)[1] <- "Pathways"

pfas_pathway_ex22 <- pfas_pathway
pfas_pathway_ex22$age <- paste0("22_",pfas_pathway_ex22$age) 


pfas_pathway <- rbind(pfas_pathway_ex0,pfas_pathway_ex7,pfas_pathway_ex14,pfas_pathway_ex22)
pfas_pathway$age <- factor(pfas_pathway$age, levels = c("0_7","0_14","0_22","0_28",
                                                        "7_7","7_14","7_22","7_28",
                                                        "14_14","14_22","14_28",
                                                        "22_22","22_28", "28_28"))

pfas_pathway$Mode <- rep("C18",nrow(pfas_pathway))
pfas_pathway_C18 <- pfas_pathway


#####################################################################################

pfas_pathway <- rbind(pfas_pathway_HILIC, pfas_pathway_C18)
pfas_pathway$age <- factor(pfas_pathway$age, levels = c("0_7","0_14","0_22","0_28",
                                                        "7_7","7_14","7_22","7_28",
                                                        "14_14","14_22","14_28",
                                                        "22_22","22_28", "28_28"))

pfas_pathway$fdr <- p.adjust(pfas_pathway$FET, "fdr")
# pfas_pathway_subset <- pfas_pathway[pfas_pathway$fdr < 0.2 & pfas_pathway$FET < 0.05,]
pfas_pathway_subset <- pfas_pathway[pfas_pathway$fdr < 0.05,]


all_pathways <- read.csv("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/Faroese data/all pfas/all_pathways.csv")
pfas_pathway_subset <- merge(pfas_pathway_subset, all_pathways, by.x = "Pathways")
pfas_pathway_subset$Group <- factor(pfas_pathway_subset$Group, levels = c("Carbohydrate metabolism",
                                                                          "Lipid metabolism", "Metabolism of amino acids and derivatives",
                                                                          "Nucleotide metabolism" ,
                                                                          "Metabolism of vitamins and cofactors", 
                                                                          "Steroid hormone biosynthesis", 
                                                                          "Energy metabolism" , "Xenobiotics biodegradation and metabolism",
                                                                          "Metabolism of terpenoids and polyketides"))

pfas_pathway_subset <- pfas_pathway_subset[order(pfas_pathway_subset$Group),]
uni_paths <- pfas_pathway_subset$Pathways
uni_paths <- uni_paths[!duplicated(uni_paths)]
pfas_pathway_subset$Pathways <- factor(pfas_pathway_subset$Pathways, levels = uni_paths)

ret <- c(length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[1]])),
         length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[2]])),
         length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[3]])),
         length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[4]])),
         length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[5]])),
         length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[6]])),
         length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[7]])),
         length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[8]])),
         length(unique(pfas_pathway_subset$Pathways[pfas_pathway_subset$Group ==  unique(pfas_pathway_subset$Group)[9]]))
         )

gp <- (ggplot(pfas_pathway_subset, aes(x = age,  y = Pathways, color = PFAS)) + 
         geom_point(aes(size = -log(FET , base = 10), shape = Mode, color = PFAS), 
                    alpha  = 5, 
                    position = position_dodge2(width = 0.5, preserve  = "total", padding = 0.5, reverse = T)) + 
         theme_bw()   + scale_shape_manual(values = c(15,16)) +
         ylab(NULL) + xlab("Age at Metabolomic Assessment") +
         ggtitle("Enriched pathways associated with serum PFAS concentrations over lifecourse\n (Both HILIC +ve and C18 -ve Mode with false discovery rate <20%)") + labs(size = "-log(raw pvalue, 10)") 
       + scale_x_discrete(labels= c("7","14","22","28","14","22","28",
                                    "22","28","28"))
       + theme(plot.title=element_text(size=14,face="bold"), axis.title=element_text(size=14,face="bold"),
               plot.tag = element_text(size = 13,face = "bold"),
               axis.text.y = element_text(size=11,face="bold"),
               axis.text.x = element_text(size=12,face="bold"),
               strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
               legend.title = element_text(face = "bold"), legend.position = 'right',
               legend.text = element_text(size = 12),
               legend.background = element_rect(fill="white", 
                                                size=0.5, linetype="solid",  colour ="darkblue"),
               plot.margin=unit(c(0.5,0.5,1,0.5), "cm"))
       +  annotate("rect", xmin = seq(0.75,9.75), 
                   xmax = seq(1.25,10.25), 
                   ymin = rep(0,10), ymax =rep(36.5,10),
                   alpha = .1)  
       + guides(color = guide_legend(override.aes = list(size = 4)),
                shape = guide_legend(override.aes = list(size = 4))))

gp <- (gp +  geom_hline(yintercept = c(cumsum(ret))+0.5)
       + geom_hline(yintercept = 36.6, size = 2)
       + geom_vline(xintercept = c(4.5, 7.5, 9.5, 11.1), size = 1.5)
       + annotate("text", x = rep(12.9, length(unique(pfas_pathway_subset$Group))), 
                  y = c(5.7,13,26,33.5,37,40, 41.5,43.4,45), 
                  label = c("Carbohydrate\nmetabolism",
                            "Lipid metabolism", "Metabolism of amino\nacids and derivatives",
                            "Nucleotide metabolism" ,
                            "Metabolism of vitamins\nand cofactors", 
                            "Steroid hormone biosynthesis", 
                            "Energy metabolism" , "Xenobiotics biodegradation and metabolism",
                            "Metabolism of terpenoids and polyketides"), color = "red",
                  size = c(5,5,5,4,4,4,4,4,4), fontface  = 'bold')
       + annotate("text", label = c("PFAS exposure\nPrenatal", "PFAS exposure\nMid-childhood: Age 7", 
                                    "PFAS exposure\nPuberty: Age 14", "PFAS exposure\nAdulthood: Age 22"),
                  x = c(2.5, 6, 8.5, 10.3), y = c(rep(40.4, 4)), fontface = "bold", size = 4)
       + coord_cartesian(xlim = c(1.2, 14), ylim = c(0.5,41), clip = "off") )


jpeg("C:/Users/midyav01/OneDrive - The Mount Sinai Hospital/MSSM Projects/METABOLOMICS/New_faroese/all_modes_all_pathways.jpeg",
     units="in", width=21, height=12, res=800)

gp

dev.off()


