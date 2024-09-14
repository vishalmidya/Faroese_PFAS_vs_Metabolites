

library(tidyverse)
library(ggpubr)
library(data.table)
library(readr)
library(readxl)
library(formattable)
library(ComplexHeatmap)
library(OlinkAnalyze)
library(gtsummary)
library(flextable)
library(ggplotify)
library(pheatmap)
library(plotly)
library(gapminder)
library(ggfortify)
library(naniar)
library(reactable)
library(sandwich)
library(ggrepel)
library(mice)
library(corrplot)



merged_omics_500 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/merged_omics_hilic.csv", check.names = F)
merged_omics_125 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/merged_omics_hilic_125.csv", check.names = F)

merged_omics_500$Year<- factor(merged_omics_500$Year, 
                               levels = c(7, 14, 22, 28))


merged_omics_125$Year<- factor(merged_omics_125$Year, 
                               levels = c(7, 14, 22, 28))

### PFAS
merged_omics_500_long<- merged_omics_500 %>% 
                        pivot_longer(
                          cols = cpfoa0:cpfda28,
                          names_to = c("PFAS_name", "PFAS_age"),
                          names_pattern = "c([a-zA-Z]+)([0-9]+)",
                          values_to = "PFAS_value"
                        )
merged_omics_500_long$PFAS_value<- factor(merged_omics_500_long$PFAS_value,
                                          levels = c(0, 1),
                                          labels = c("Lower", "Higher"))

PFAS_plot<- ggplot(merged_omics_500_long, aes(fill=PFAS_value, x=PFAS_age)) +
            geom_bar(position="stack", stat = "count") +
            facet_grid(~PFAS_name, scale="free")+ 
            theme_bw() +
            labs(title = "") +
            xlab("PFAS Age") + 
            ylab("Count") + 
            theme(plot.title = element_text(hjust = 0.5, size = 18),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.text = element_text(size = 12, face="bold"),
                  axis.title = element_text(size = 14, face="bold", vjust=0.5),
                  axis.title.x = element_text(size = 14, face="bold", vjust=-0.8),
                  legend.title = element_blank(),
                  panel.grid = element_blank()) 

PFAS_plot

### metabolites
ggplot(merged_omics_500) +
  geom_density(aes(x = Met196)) 

ggplot(merged_omics_125) +
  geom_density(aes(x = Met196)) 

Met196_500<- ggplot(merged_omics_500) +
             geom_density(aes(x = Met196, color=Year)) 

Met196_125<- ggplot(merged_omics_125) +
             geom_density(aes(x = Met196, color=Year)) 

Met196_500
Met196_125


ggplot(merged_omics_500) +
  geom_density(aes(x = Met1)) 

ggplot(merged_omics_125) +
  geom_density(aes(x = Met1)) 

ggplot(merged_omics_500) +
  geom_density(aes(x = Met1))
Met1_500<- ggplot(merged_omics_500) +
           geom_density(aes(x = Met1, color=Year)) 

Met1_125<- ggplot(merged_omics_125) +
           geom_density(aes(x = Met1, color=Year)) 

Met1_500
Met1_125

Met3_500<- ggplot(merged_omics_500) +
           geom_density(aes(x = Met3, color=Year)) 

Met3_125<- ggplot(merged_omics_125) +
           geom_density(aes(x = Met3, color=Year)) 

Met3_500
Met3_125





