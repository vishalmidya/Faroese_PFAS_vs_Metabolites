minerva.hpc.mssm.edu

##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################

## pfos at 0 

### met_7

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_7_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_7_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_7_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_7_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_7_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_7_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_7_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_7_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J comb_pfos_0_met_7 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_0_met_7.R


################################################################################

## pfos at 0 

### met_14

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_14_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_14_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_14_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_14_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_14_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_14_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_14_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_14_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J comb_pfos_0_met_14 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_0_met_14.R


################################################################################

## pfos at 0 

### met_22

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_22_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_22_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_22_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_22_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_22_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_22_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_22_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_22_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J comb_pfos_0_met_22 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_0_met_22.R


################################################################################

## pfos at 0 

### met_28

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_28_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_28_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_28_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash


module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_28_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_28_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash


module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_28_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J teststat_pfos_0_met_28_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_0_met_28_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age0
bsub -J comb_pfos_0_met_28 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_0_met_28.R


##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################

## pfos at 7

### met_7

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_7_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_7_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_7_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_7_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_7_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_7_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_7_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_7_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J comb_pfos_7_met_7 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_7_met_7.R


################################################################################

### met_14

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_14_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_14_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_14_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_14_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_14_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_14_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_14_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_14_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J comb_pfos_7_met_14 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_7_met_14.R


################################################################################

### met_22

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_22_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_22_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_22_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_22_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_22_v3 -P acc_damasv01a  -q premium  -n 5 -R rusage[mem=17G] -R span[ptile=5] -W 9:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_22_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_22_v4 -P acc_damasv01a  -q premium  -n 5 -R rusage[mem=16G] -R span[ptile=5] -W 9:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_22_v4.R


#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J comb_pfos_7_met_22 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_7_met_22.R


################################################################################


### met_28

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_28_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_28_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_28_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_28_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_28_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_28_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J teststat_pfos_7_met_28_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_7_met_28_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age7
bsub -J comb_pfos_7_met_28 -P acc_damasv01a  -q premium -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_7_met_28.R




##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################

## pfos at 14

### met_14

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_14_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_14_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_14_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_14_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_14_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_14_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_14_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_14_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J comb_pfos_14_met_14 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_14_met_14.R


################################################################################

### met_22

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_22_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_22_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_22_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_22_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_22_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_22_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_22_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_22_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J comb_pfos_14_met_22 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_14_met_22.R


################################################################################


### met_28

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_28_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_28_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_28_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_28_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_28_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_28_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J teststat_pfos_14_met_28_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_14_met_28_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age14
bsub -J comb_pfos_14_met_28 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_14_met_28.R




##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################

## pfos at 22

### met_22

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age22
bsub -J teststat_pfos_22_met_22_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_22_met_22_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age22
bsub -J teststat_pfos_22_met_22_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_22_met_22_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age22
bsub -J teststat_pfos_22_met_22_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_22_met_22_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age22
bsub -J teststat_pfos_22_met_22_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_22_met_22_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age22
bsub -J comb_pfos_22_met_22 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_22_met_22.R


################################################################################


### met_28

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age22
bsub -J teststat_pfos_22_met_28_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_22_met_28_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age22
bsub -J teststat_pfos_22_met_28_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_22_met_28_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age22
bsub -J teststat_pfos_22_met_28_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_22_met_28_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age22
bsub -J teststat_pfos_22_met_28_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_22_met_28_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age22
bsub -J comb_pfos_22_met_28 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_22_met_28.R




##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################
##############################################################################################################################################################################

## pfos at 28

### met_28

#### v1

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age28
bsub -J teststat_pfos_28_met_28_v1 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_28_met_28_v1.R

#### v2

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age28
bsub -J teststat_pfos_28_met_28_v2 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_28_met_28_v2.R

#### v3

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age28
bsub -J teststat_pfos_28_met_28_v3 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_28_met_28_v3.R

#### v4

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age28
bsub -J teststat_pfos_28_met_28_v4 -P acc_damasv01a  -q premium  -n 15 -R rusage[mem=6G] -R span[ptile=15] -W 14:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript test_stat_table_pfos_28_met_28_v4.R

#### combine

cd /
  cd sc/arion/work/midyav01/faroese/pfos/age28
bsub -J comb_pfos_28_met_28 -P acc_damasv01a  -q premium  -n 10 -R rusage[mem=9G] -R span[ptile=10] -W 01:00 -o %J.stdout -eo %J.stderr -L /bin/bash

module load gcc/8.3.0 
module load R/4.1.0
Rscript comb_pfos_28_met_28.R



