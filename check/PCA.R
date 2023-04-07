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
library(ggrepel)

# HILIC
## import data
merged_omics <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/merged_omics_hilic.csv", check.names = F)


## focus on metabolites only
all_met<- merged_omics %>% 
          select(Year, starts_with("Met"))

### age 7
all_met_age7<- all_met %>% 
               filter(Year == 7) %>% 
               select(starts_with("Met"))


### age 14
all_met_age14<- all_met %>% 
                filter(Year == 14) %>% 
                select(starts_with("Met"))


### age 22
all_met_age22<- all_met %>% 
                filter(Year == 22) %>% 
                select(starts_with("Met"))
 

### age 28
all_met_age28<- all_met %>% 
                filter(Year == 28) %>% 
                select(starts_with("Met"))

## check missing values
colSums(is.na(all_met_age7))
colSums(is.na(all_met_age14))
colSums(is.na(all_met_age22))
colSums(is.na(all_met_age28))



## check correlation & PCA

### age7
corr_met_age7<- cor(all_met_age7)
# ggcorrplot(corr_met_age7)

pca_met_age7<- princomp(corr_met_age7)
pca_met_age7_summary<- summary(pca_met_age7)

pca_met_age7_loading<- pca_met_age7$loadings[, 1:2]
pca_met_age7_screenplot<- fviz_eig(pca_met_age7, addlabels = TRUE)

pca_met_age7_varplot<- fviz_pca_var(pca_met_age7, col.var = "black", 
                                    repel = TRUE, geom = c("point", "text"),
                                    select.var = list(name = c('Met1269',
                                                               'Met1055',
                                                               'Met1571',
                                                               'Met363',
                                                               'Met459',
                                                               'Met1936',
                                                               'Met407',
                                                               'Met1771',
                                                               'Met389',
                                                               'Met1473',
                                                               'Met398'
                                                               
                                    )))  +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  labs(title = "Met at age 7")+
  theme(panel.border = element_rect(color = "grey",
                                    fill = NA,
                                    size = 5))



### age14
corr_met_age14<- cor(all_met_age14)
# ggcorrplot(corr_met_age14)

pca_met_age14<- princomp(corr_met_age14)
pca_met_age14_summary<- summary(pca_met_age14)

pca_met_age14_loading<- pca_met_age14$loadings[, 1:2]
pca_met_age14_screenplot<- fviz_eig(pca_met_age14, addlabels = TRUE)

pca_met_age14_varplot<- fviz_pca_var(pca_met_age14, col.var = "black", 
                                     repel = TRUE, geom = c("point", "text"),
                                     select.var = list(name = c('Met1269',
                                                                'Met1055',
                                                                'Met1571',
                                                                'Met363',
                                                                'Met459',
                                                                'Met1936',
                                                                'Met407',
                                                                'Met1771',
                                                                'Met389',
                                                                'Met1473',
                                                                'Met398'
                                     )))  +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  labs(title = "Met at age 14")+
  theme(panel.border = element_rect(color = "grey",
                                    fill = NA,
                                    size = 5))



### age22
corr_met_age22<- cor(all_met_age22)
# ggcorrplot(corr_met_age22)

pca_met_age22<- princomp(corr_met_age22)
pca_met_age22_summary<- summary(pca_met_age22)

pca_met_age22_loading<- pca_met_age22$loadings[, 1:2]
pca_met_age22_screenplot<- fviz_eig(pca_met_age22, addlabels = TRUE)

pca_met_age22_varplot<- fviz_pca_var(pca_met_age22, col.var = "black", 
                                     repel = TRUE, geom = c("point", "text"),
                                     select.var = list(name = c('Met1269',
                                                                'Met1055',
                                                                'Met1571',
                                                                'Met363',
                                                                'Met459',
                                                                'Met1936',
                                                                'Met407',
                                                                'Met1771',
                                                                'Met389',
                                                                'Met1473',
                                                                'Met398'
                                     )))  +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  labs(title = "Met at age 22")+
  theme(panel.border = element_rect(color = "grey",
                                    fill = NA,
                                    size = 5))


### age28
corr_met_age28<- cor(all_met_age28)
# ggcorrplot(corr_met_age28)

pca_met_age28<- princomp(corr_met_age28)
pca_met_age28_summary<- summary(pca_met_age28)

pca_met_age28_loading<- pca_met_age28$loadings[, 1:2]
pca_met_age28_screenplot<- fviz_eig(pca_met_age28, addlabels = TRUE)

pca_met_age28_varplot<- fviz_pca_var(pca_met_age28, col.var = "black", 
                                    geom = c("point", "text"),
                                    select.var = list(name = c('Met1269',
                                                               'Met1055',
                                                               'Met1571',
                                                               'Met363',
                                                               'Met459',
                                                               'Met1936',
                                                               'Met407',
                                                               'Met1771',
                                                               'Met389',
                                                               'Met1473',
                                                               'Met398'
                                    )))  +
                      xlim(-0.2, 0.2) +
                      ylim(-0.2, 0.2) +
                      labs(title = "Met at age 28")+
                      theme(panel.border = element_rect(color = "grey",
                                                         fill = NA,
                                                         size = 5))






### combine
jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/pca_met_hilic.jpeg",
     units="in", width=24, height=24, res=600)

ggarrange(pca_met_age7_varplot, pca_met_age14_varplot,
          pca_met_age22_varplot,pca_met_age28_varplot,
          ncol=2, nrow=2)


dev.off()


# C18
## import data
merged_omics <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/merged_omics_c18.csv", check.names = F)


## focus on metabolites only
all_met<- merged_omics %>% 
  select(Year, starts_with("Met"))

### age 7
all_met_age7<- all_met %>% 
  filter(Year == 7) %>% 
  select(starts_with("Met"))


### age 14
all_met_age14<- all_met %>% 
  filter(Year == 14) %>% 
  select(starts_with("Met"))


### age 22
all_met_age22<- all_met %>% 
  filter(Year == 22) %>% 
  select(starts_with("Met"))


### age 28
all_met_age28<- all_met %>% 
  filter(Year == 28) %>% 
  select(starts_with("Met"))

## check missing values
colSums(is.na(all_met_age7))
colSums(is.na(all_met_age14))
colSums(is.na(all_met_age22))
colSums(is.na(all_met_age28))



## check correlation & PCA

### age7
corr_met_age7<- cor(all_met_age7)
# ggcorrplot(corr_met_age7)

pca_met_age7<- princomp(corr_met_age7)
pca_met_age7_summary<- summary(pca_met_age7)

pca_met_age7_loading<- pca_met_age7$loadings[, 1:2]
pca_met_age7_screenplot<- fviz_eig(pca_met_age7, addlabels = TRUE)

pca_met_age7_varplot<- fviz_pca_var(pca_met_age7, col.var = "black", 
                                    repel = TRUE, geom = c("point", "text"),
                                    select.var = list(name = c('Met87',
                                                               'Met253',
                                                               'Met695',
                                                               'Met707',
                                                               'Met676'
                                                               
                                    )))  +
  xlim(-0.3, 0.3) +
  ylim(-0.3, 0.3) +
  labs(title = "Met at age 7")+
  theme(panel.border = element_rect(color = "grey",
                                    fill = NA,
                                    size = 5))



### age14
corr_met_age14<- cor(all_met_age14)
# ggcorrplot(corr_met_age14)

pca_met_age14<- princomp(corr_met_age14)
pca_met_age14_summary<- summary(pca_met_age14)

pca_met_age14_loading<- pca_met_age14$loadings[, 1:2]
pca_met_age14_screenplot<- fviz_eig(pca_met_age14, addlabels = TRUE)

pca_met_age14_varplot<- fviz_pca_var(pca_met_age14, col.var = "black", 
                                     repel = TRUE, geom = c("point", "text"),
                                     select.var = list(name = c('Met87',
                                                                'Met253',
                                                                'Met695',
                                                                'Met707',
                                                                'Met676'
                                     )))  +
  xlim(-0.3, 0.3) +
  ylim(-0.3, 0.3) +
  labs(title = "Met at age 14")+
  theme(panel.border = element_rect(color = "grey",
                                    fill = NA,
                                    size = 5))



### age22
corr_met_age22<- cor(all_met_age22)
# ggcorrplot(corr_met_age22)

pca_met_age22<- princomp(corr_met_age22)
pca_met_age22_summary<- summary(pca_met_age22)

pca_met_age22_loading<- pca_met_age22$loadings[, 1:2]
pca_met_age22_screenplot<- fviz_eig(pca_met_age22, addlabels = TRUE)

pca_met_age22_varplot<- fviz_pca_var(pca_met_age22, col.var = "black", 
                                     repel = TRUE, geom = c("point", "text"),
                                     select.var = list(name = c('Met87',
                                                                'Met253',
                                                                'Met695',
                                                                'Met707',
                                                                'Met676'
                                     )))  +
  xlim(-0.3, 0.3) +
  ylim(-0.3, 0.3) +
  labs(title = "Met at age 22")+
  theme(panel.border = element_rect(color = "grey",
                                    fill = NA,
                                    size = 5))


### age28
corr_met_age28<- cor(all_met_age28)
# ggcorrplot(corr_met_age28)

pca_met_age28<- princomp(corr_met_age28)
pca_met_age28_summary<- summary(pca_met_age28)

pca_met_age28_loading<- pca_met_age28$loadings[, 1:2]
pca_met_age28_screenplot<- fviz_eig(pca_met_age28, addlabels = TRUE)

pca_met_age28_varplot<- fviz_pca_var(pca_met_age28, col.var = "black", 
                                     repel = TRUE, geom = c("point", "text"),
                                     select.var = list(name = c('Met87',
                                                                'Met253',
                                                                'Met695',
                                                                'Met707',
                                                                'Met676'
                                     )))  +
  xlim(-0.3, 0.3) +
  ylim(-0.3, 0.3) +
  labs(title = "Met at age 28")+
  theme(panel.border = element_rect(color = "grey",
                                    fill = NA,
                                    size = 5))



### combine
jpeg("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/pca_met_c18.jpeg",
     units="in", width=24, height=24, res=600)

ggarrange(pca_met_age7_varplot, pca_met_age14_varplot,
          pca_met_age22_varplot,pca_met_age28_varplot,
          ncol=2, nrow=2)


dev.off()


