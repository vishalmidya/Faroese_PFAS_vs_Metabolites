library(tidyverse)
library(ggpubr)
library(data.table)
library(dplyr)



## only keep four columns used in MetaboAnalyst
format<- function(input_path, output_path) {
  input<- fread(input_path, header = TRUE)
  output<- input %>%
    select(mz, time, beta, simu_pval) %>%
    dplyr::rename(m.z=mz, r.t=time, t.score=beta, p.value=simu_pval) %>%
    arrange(p.value)
  write.table(output, output_path, row.names = FALSE)
  
}

##----------------------------------------- pfos
### HILIC
#### at 0
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/pfos_0_met_7_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/pfos_0_met_7_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/pfos_0_met_14_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/pfos_0_met_14_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/pfos_0_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/pfos_0_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/pfos_0_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/pfos_0_met_28_beta_fisher_hilic_mummichog.txt")

#### at 7
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_7/minerva_data_pfos_7_metabolites_14/pfos_7_met_14_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_7/minerva_data_pfos_7_metabolites_14/pfos_7_met_14_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_7/minerva_data_pfos_7_metabolites_22/pfos_7_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_7/minerva_data_pfos_7_metabolites_22/pfos_7_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_7/minerva_data_pfos_7_metabolites_28/pfos_7_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_7/minerva_data_pfos_7_metabolites_28/pfos_7_met_28_beta_fisher_hilic_mummichog.txt")

#### at 14
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_14/minerva_data_pfos_14_metabolites_22/pfos_14_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_14/minerva_data_pfos_14_metabolites_22/pfos_14_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_14/minerva_data_pfos_14_metabolites_28/pfos_14_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_14/minerva_data_pfos_14_metabolites_28/pfos_14_met_28_beta_fisher_hilic_mummichog.txt")

#### at 22
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_22/minerva_data_pfos_22_metabolites_28/pfos_22_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfos/pfos_22/minerva_data_pfos_22_metabolites_28/pfos_22_met_28_beta_fisher_hilic_mummichog.txt")

### C18
#### at 0
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/pfos_0_met_7_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_0/minerva_data_pfos_0_metabolites_7/pfos_0_met_7_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/pfos_0_met_14_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_0/minerva_data_pfos_0_metabolites_14/pfos_0_met_14_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/pfos_0_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_0/minerva_data_pfos_0_metabolites_22/pfos_0_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/pfos_0_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_0/minerva_data_pfos_0_metabolites_28/pfos_0_met_28_beta_fisher_c18_mummichog.txt")

#### at 7
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_7/minerva_data_pfos_7_metabolites_14/pfos_7_met_14_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_7/minerva_data_pfos_7_metabolites_14/pfos_7_met_14_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_7/minerva_data_pfos_7_metabolites_22/pfos_7_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_7/minerva_data_pfos_7_metabolites_22/pfos_7_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_7/minerva_data_pfos_7_metabolites_28/pfos_7_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_7/minerva_data_pfos_7_metabolites_28/pfos_7_met_28_beta_fisher_c18_mummichog.txt")

#### at 14
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_14/minerva_data_pfos_14_metabolites_22/pfos_14_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_14/minerva_data_pfos_14_metabolites_22/pfos_14_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_14/minerva_data_pfos_14_metabolites_28/pfos_14_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_14/minerva_data_pfos_14_metabolites_28/pfos_14_met_28_beta_fisher_c18_mummichog.txt")

#### at 22
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_22/minerva_data_pfos_22_metabolites_28/pfos_22_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfos/pfos_22/minerva_data_pfos_22_metabolites_28/pfos_22_met_28_beta_fisher_c18_mummichog.txt")



##----------------------------------------- pfoa
### HILIC
#### at 0
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_7/pfoa_0_met_7_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_7/pfoa_0_met_7_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_14/pfoa_0_met_14_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_14/pfoa_0_met_14_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_22/pfoa_0_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_22/pfoa_0_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_28/pfoa_0_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_28/pfoa_0_met_28_beta_fisher_hilic_mummichog.txt")

#### at 7
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_14/pfoa_7_met_14_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_14/pfoa_7_met_14_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_22/pfoa_7_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_22/pfoa_7_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_28/pfoa_7_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_28/pfoa_7_met_28_beta_fisher_hilic_mummichog.txt")

#### at 14
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/pfoa_14_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/pfoa_14_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/pfoa_14_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/pfoa_14_met_28_beta_fisher_hilic_mummichog.txt")

#### at 22
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_hilic_mummichog.txt")

### C18
#### at 0
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_7/pfoa_0_met_7_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_7/pfoa_0_met_7_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_14/pfoa_0_met_14_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_14/pfoa_0_met_14_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_22/pfoa_0_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_22/pfoa_0_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_28/pfoa_0_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_0/minerva_data_pfoa_0_metabolites_28/pfoa_0_met_28_beta_fisher_c18_mummichog.txt")

#### at 7
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_14/pfoa_7_met_14_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_14/pfoa_7_met_14_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_22/pfoa_7_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_22/pfoa_7_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_28/pfoa_7_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_7/minerva_data_pfoa_7_metabolites_28/pfoa_7_met_28_beta_fisher_c18_mummichog.txt")

#### at 14
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/pfoa_14_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_22/pfoa_14_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/pfoa_14_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_14/minerva_data_pfoa_14_metabolites_28/pfoa_14_met_28_beta_fisher_c18_mummichog.txt")

#### at 22
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfoa/pfoa_22/minerva_data_pfoa_22_metabolites_28/pfoa_22_met_28_beta_fisher_c18_mummichog.txt")


##----------------------------------------- pfna
### HILIC
#### at 0
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_0/minerva_data_pfna_0_metabolites_7/pfna_0_met_7_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_0/minerva_data_pfna_0_metabolites_7/pfna_0_met_7_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_0/minerva_data_pfna_0_metabolites_14/pfna_0_met_14_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_0/minerva_data_pfna_0_metabolites_14/pfna_0_met_14_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_0/minerva_data_pfna_0_metabolites_22/pfna_0_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_0/minerva_data_pfna_0_metabolites_22/pfna_0_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_0/minerva_data_pfna_0_metabolites_28/pfna_0_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_0/minerva_data_pfna_0_metabolites_28/pfna_0_met_28_beta_fisher_hilic_mummichog.txt")

#### at 7
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_7/minerva_data_pfna_7_metabolites_14/pfna_7_met_14_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_7/minerva_data_pfna_7_metabolites_14/pfna_7_met_14_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_7/minerva_data_pfna_7_metabolites_22/pfna_7_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_7/minerva_data_pfna_7_metabolites_22/pfna_7_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/pfna_7_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/pfna_7_met_28_beta_fisher_hilic_mummichog.txt")

#### at 14
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_14/minerva_data_pfna_14_metabolites_22/pfna_14_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_14/minerva_data_pfna_14_metabolites_22/pfna_14_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_14/minerva_data_pfna_14_metabolites_28/pfna_14_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_14/minerva_data_pfna_14_metabolites_28/pfna_14_met_28_beta_fisher_hilic_mummichog.txt")

#### at 22
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_22/minerva_data_pfna_22_metabolites_28/pfna_22_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfna/pfna_22/minerva_data_pfna_22_metabolites_28/pfna_22_met_28_beta_fisher_hilic_mummichog.txt")

### C18
#### at 0
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_0/minerva_data_pfna_0_metabolites_7/pfna_0_met_7_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_0/minerva_data_pfna_0_metabolites_7/pfna_0_met_7_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_0/minerva_data_pfna_0_metabolites_14/pfna_0_met_14_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_0/minerva_data_pfna_0_metabolites_14/pfna_0_met_14_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_0/minerva_data_pfna_0_metabolites_22/pfna_0_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_0/minerva_data_pfna_0_metabolites_22/pfna_0_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_0/minerva_data_pfna_0_metabolites_28/pfna_0_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_0/minerva_data_pfna_0_metabolites_28/pfna_0_met_28_beta_fisher_c18_mummichog.txt")

#### at 7
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_7/minerva_data_pfna_7_metabolites_14/pfna_7_met_14_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_7/minerva_data_pfna_7_metabolites_14/pfna_7_met_14_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_7/minerva_data_pfna_7_metabolites_22/pfna_7_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_7/minerva_data_pfna_7_metabolites_22/pfna_7_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/pfna_7_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_7/minerva_data_pfna_7_metabolites_28/pfna_7_met_28_beta_fisher_c18_mummichog.txt")

#### at 14
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_14/minerva_data_pfna_14_metabolites_22/pfna_14_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_14/minerva_data_pfna_14_metabolites_22/pfna_14_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_14/minerva_data_pfna_14_metabolites_28/pfna_14_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_14/minerva_data_pfna_14_metabolites_28/pfna_14_met_28_beta_fisher_c18_mummichog.txt")

#### at 22
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_22/minerva_data_pfna_22_metabolites_28/pfna_22_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfna/pfna_22/minerva_data_pfna_22_metabolites_28/pfna_22_met_28_beta_fisher_c18_mummichog.txt")


##----------------------------------------- pfda
### HILIC
#### at 0
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_0/minerva_data_pfda_0_metabolites_7/pfda_0_met_7_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_0/minerva_data_pfda_0_metabolites_7/pfda_0_met_7_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_0/minerva_data_pfda_0_metabolites_14/pfda_0_met_14_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_0/minerva_data_pfda_0_metabolites_14/pfda_0_met_14_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_0/minerva_data_pfda_0_metabolites_22/pfda_0_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_0/minerva_data_pfda_0_metabolites_22/pfda_0_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_0/minerva_data_pfda_0_metabolites_28/pfda_0_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_0/minerva_data_pfda_0_metabolites_28/pfda_0_met_28_beta_fisher_hilic_mummichog.txt")

#### at 7
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_7/minerva_data_pfda_7_metabolites_14/pfda_7_met_14_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_7/minerva_data_pfda_7_metabolites_14/pfda_7_met_14_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_7/minerva_data_pfda_7_metabolites_22/pfda_7_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_7/minerva_data_pfda_7_metabolites_22/pfda_7_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_7/minerva_data_pfda_7_metabolites_28/pfda_7_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_7/minerva_data_pfda_7_metabolites_28/pfda_7_met_28_beta_fisher_hilic_mummichog.txt")

#### at 14
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_14/minerva_data_pfda_14_metabolites_22/pfda_14_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_14/minerva_data_pfda_14_metabolites_22/pfda_14_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_14/minerva_data_pfda_14_metabolites_28/pfda_14_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_14/minerva_data_pfda_14_metabolites_28/pfda_14_met_28_beta_fisher_hilic_mummichog.txt")

#### at 22
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/pfda_22_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/pfda_22_met_28_beta_fisher_hilic_mummichog.txt")

### C18
#### at 0
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_0/minerva_data_pfda_0_metabolites_7/pfda_0_met_7_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_0/minerva_data_pfda_0_metabolites_7/pfda_0_met_7_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_0/minerva_data_pfda_0_metabolites_14/pfda_0_met_14_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_0/minerva_data_pfda_0_metabolites_14/pfda_0_met_14_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_0/minerva_data_pfda_0_metabolites_22/pfda_0_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_0/minerva_data_pfda_0_metabolites_22/pfda_0_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_0/minerva_data_pfda_0_metabolites_28/pfda_0_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_0/minerva_data_pfda_0_metabolites_28/pfda_0_met_28_beta_fisher_c18_mummichog.txt")

#### at 7
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_7/minerva_data_pfda_7_metabolites_14/pfda_7_met_14_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_7/minerva_data_pfda_7_metabolites_14/pfda_7_met_14_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_7/minerva_data_pfda_7_metabolites_22/pfda_7_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_7/minerva_data_pfda_7_metabolites_22/pfda_7_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_7/minerva_data_pfda_7_metabolites_28/pfda_7_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_7/minerva_data_pfda_7_metabolites_28/pfda_7_met_28_beta_fisher_c18_mummichog.txt")

#### at 14
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_14/minerva_data_pfda_14_metabolites_22/pfda_14_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_14/minerva_data_pfda_14_metabolites_22/pfda_14_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_14/minerva_data_pfda_14_metabolites_28/pfda_14_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_14/minerva_data_pfda_14_metabolites_28/pfda_14_met_28_beta_fisher_c18_mummichog.txt")

#### at 22
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/pfda_22_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfda/pfda_22/minerva_data_pfda_22_metabolites_28/pfda_22_met_28_beta_fisher_c18_mummichog.txt")


##----------------------------------------- pfhxs
### HILIC
#### at 0
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_7/pfhxs_0_met_7_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_7/pfhxs_0_met_7_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_14/pfhxs_0_met_14_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_14/pfhxs_0_met_14_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_22/pfhxs_0_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_22/pfhxs_0_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_28/pfhxs_0_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_28/pfhxs_0_met_28_beta_fisher_hilic_mummichog.txt")

#### at 7
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_14/pfhxs_7_met_14_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_14/pfhxs_7_met_14_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_22/pfhxs_7_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_22/pfhxs_7_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_28/pfhxs_7_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_28/pfhxs_7_met_28_beta_fisher_hilic_mummichog.txt")

#### at 14
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_22/pfhxs_14_met_22_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_22/pfhxs_14_met_22_beta_fisher_hilic_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_28/pfhxs_14_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_28/pfhxs_14_met_28_beta_fisher_hilic_mummichog.txt")

#### at 22
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_22/minerva_data_pfhxs_22_metabolites_28/pfhxs_22_met_28_beta_fisher_hilic.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/HILIC/pfhxs/pfhxs_22/minerva_data_pfhxs_22_metabolites_28/pfhxs_22_met_28_beta_fisher_hilic_mummichog.txt")

### C18
#### at 0
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_7/pfhxs_0_met_7_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_7/pfhxs_0_met_7_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_14/pfhxs_0_met_14_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_14/pfhxs_0_met_14_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_22/pfhxs_0_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_22/pfhxs_0_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_28/pfhxs_0_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_0/minerva_data_pfhxs_0_metabolites_28/pfhxs_0_met_28_beta_fisher_c18_mummichog.txt")

#### at 7
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_14/pfhxs_7_met_14_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_14/pfhxs_7_met_14_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_22/pfhxs_7_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_22/pfhxs_7_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_28/pfhxs_7_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_7/minerva_data_pfhxs_7_metabolites_28/pfhxs_7_met_28_beta_fisher_c18_mummichog.txt")

#### at 14
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_22/pfhxs_14_met_22_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_22/pfhxs_14_met_22_beta_fisher_c18_mummichog.txt")

format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_28/pfhxs_14_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_14/minerva_data_pfhxs_14_metabolites_28/pfhxs_14_met_28_beta_fisher_c18_mummichog.txt")

#### at 22
format(input_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_22/minerva_data_pfhxs_22_metabolites_28/pfhxs_22_met_28_beta_fisher_c18.txt", 
       output_path="C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/New_faroese/C18/pfhxs/pfhxs_22/minerva_data_pfhxs_22_metabolites_28/pfhxs_22_met_28_beta_fisher_c18_mummichog.txt")

