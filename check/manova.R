library(readr)
library(car)
library(MASS)
library(ggplot2)
library(data.table)
library(iterators)
library(tidyverse)
library(ggcorrplot)
library(factoextra)
library(ggpubr)
library(vegan)

############################## HILIC
## import data
merged_omics <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/merged_omics_hilic.csv", check.names = F)

all_met<- merged_omics %>% 
  select(Year, starts_with("Met"), sex, mage, mbmi, smokepreg_2, cmatfishpreg, cparity, age7, age14, age22, age28)

### age 7
all_met_age7<- all_met %>% 
  filter(Year == 7) %>% 
  select(starts_with("Met"), sex, mage, mbmi, smokepreg_2, cmatfishpreg, cparity, age7)


### age 14
all_met_age14<- all_met %>% 
  filter(Year == 14) %>% 
  select(starts_with("Met"), sex, mage, mbmi, smokepreg_2, cmatfishpreg, cparity, age14)


### age 22
all_met_age22<- all_met %>% 
  filter(Year == 22) %>% 
  select(starts_with("Met"), sex, mage, mbmi, smokepreg_2, cmatfishpreg, cparity, age22)


### age 28
all_met_age28<- all_met %>% 
  filter(Year == 28) %>% 
  select(starts_with("Met"), sex, mage, mbmi, smokepreg_2, cmatfishpreg, cparity, age28)


## manova

### age 7
met_age7 <- as.matrix(all_met_age7[,1:1991])

manova_age7_model <- adonis2(met_age7 ~ sex+ mage+ mbmi+ smokepreg_2+ 
                             cmatfishpreg+ cparity+ age7, data = all_met_age7, 
                             by='margin', permutations = 1e4, method="euclidean")

### age 14

met_age14 <- as.matrix(all_met_age14[,1:1991])
# library(compositions)ilr(clo(met_age7))
manova_age14_model <- adonis2(met_age14 ~ sex+ mage+ mbmi+ smokepreg_2+ 
                              cmatfishpreg+ cparity+ age14, data = all_met_age14, 
                              by='margin', permutations = 1e4, method="euclidean")



### age 22

met_age22 <- as.matrix(all_met_age22[,1:1991])
# library(compositions)ilr(clo(met_age7))
manova_age22_model <- adonis2(met_age22 ~ sex+ mage+ mbmi+ smokepreg_2+ 
                              cmatfishpreg+ cparity+ age22, data = all_met_age22, 
                              by='margin', permutations = 1e4, method="euclidean")

### age 28

met_age28 <- as.matrix(all_met_age28[,1:1991])
# library(compositions)ilr(clo(met_age7))
manova_age28_model <- adonis2(met_age28 ~ sex+ mage+ mbmi+ smokepreg_2+ 
                                cmatfishpreg+ cparity+ age28, data = all_met_age28, 
                              by='margin', permutations = 1e4, method="euclidean")

manova_age7_model
manova_age14_model
manova_age22_model
manova_age28_model



############################## C18
## import data
merged_omics <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/merged_omics_c18.csv", check.names = F)

all_met<- merged_omics %>% 
  select(Year, starts_with("Met"), sex, mage, mbmi, smokepreg_2, cmatfishpreg, cparity, age7, age14, age22, age28)

### age 7
all_met_age7<- all_met %>% 
  filter(Year == 7) %>% 
  select(starts_with("Met"), sex, mage, mbmi, smokepreg_2, cmatfishpreg, cparity, age7)


### age 14
all_met_age14<- all_met %>% 
  filter(Year == 14) %>% 
  select(starts_with("Met"), sex, mage, mbmi, smokepreg_2, cmatfishpreg, cparity, age14)


### age 22
all_met_age22<- all_met %>% 
  filter(Year == 22) %>% 
  select(starts_with("Met"), sex, mage, mbmi, smokepreg_2, cmatfishpreg, cparity, age22)


### age 28
all_met_age28<- all_met %>% 
  filter(Year == 28) %>% 
  select(starts_with("Met"), sex, mage, mbmi, smokepreg_2, cmatfishpreg, cparity, age28)


## manova

### age 7
met_age7 <- as.matrix(all_met_age7[,1:787])

manova_age7_model <- adonis2(met_age7 ~ sex+ mage+ mbmi+ smokepreg_2+ 
                               cmatfishpreg+ cparity+ age7, data = all_met_age7, 
                             by='margin', permutations = 1e4, method="euclidean")

### age 14

met_age14 <- as.matrix(all_met_age14[,1:787])
# library(compositions)ilr(clo(met_age7))
manova_age14_model <- adonis2(met_age14 ~ sex+ mage+ mbmi+ smokepreg_2+ 
                                cmatfishpreg+ cparity+ age14, data = all_met_age14, 
                              by='margin', permutations = 1e4, method="euclidean")



### age 22

met_age22 <- as.matrix(all_met_age22[,1:787])
# library(compositions)ilr(clo(met_age7))
manova_age22_model <- adonis2(met_age22 ~ sex+ mage+ mbmi+ smokepreg_2+ 
                                cmatfishpreg+ cparity+ age22, data = all_met_age22, 
                              by='margin', permutations = 1e4, method="euclidean")

### age 28

met_age28 <- as.matrix(all_met_age28[,1:787])
# library(compositions)ilr(clo(met_age7))
manova_age28_model <- adonis2(met_age28 ~ sex+ mage+ mbmi+ smokepreg_2+ 
                                cmatfishpreg+ cparity+ age28, data = all_met_age28, 
                              by='margin', permutations = 1e4, method="euclidean")


manova_age7_model
manova_age14_model
manova_age22_model
manova_age28_model












