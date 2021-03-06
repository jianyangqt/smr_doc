
## MeCS

#### MeCS: Meta-analysis of cis-eQTL in Correlated Samples

#### \# Overview

MeCS is a method that only requires summary-level cis-eQTL data to
perform a meta-analysis of cis-eQTLs from multiple cohorts (or tissues)
with sample overlaps. It estimates the proportion of sample overlap from
null SNPs in the cis-regions and meta-analysed the eQTL effects using a
generalized least squares approach. The method can be applied to data
from genetic studies of molecular phenotypes (e.g. DNA methylation and
histone modification).

Bug reports or questions to Jian Yang (<jian.yang@uq.edu.au>) at
Institute for Molecular Bioscience, The University of Queensland.

#### \# Tutorial

>Example

```
smr --besd-flist my_file.list --mecs --thread-num 5 --out mecs_result 
```

**\--mecs** implements the MeCS analysis.

Specify a p-value threshold to exclude the significant SNPs from calculating the cohort correlation matrix.
```
osca --besd-flist my_file.list --mecs --pmecs 0.01 --out mymecs
```
**\--pmecs** reads a p-value threshold to exclude the significant SNPs from calculating the cohort correlation matrix. The default value is 0.01. 

Specify a minimum number of null SNPs for calculating the cohort correlation matrix.
```
osca --besd-flist my_file.list --mecs --nmecs 100 --out mymecs
```
**\--nmecs** reads a minimum number of null SNPs  that are required to calculate the cohort correlation matrix. The default value is 100. 


>Example

```
smr --besd-flist my_file.list --meta --thread-num 5 --out meta_result 
```
**\--meta** implements the conventional inverse-variance-weighted meta-analysis assuming all the cohorts are independent.


#### Citation
Qi T, Wu Y, Zeng J, Zhang F, Xue A, Jiang L, Zhu Z, Kemper K, Yengo L, Zheng Z, eQTLGen Consortium, Marioni RE, Montgomery GW, Deary IJ, Wray NR, Visscher PM, McRae AF & Yang J (2018) Identifying gene targets for brain-related traits using transcriptomic and methylomic data from blood. [Nature Communications, 9: 2282.](https://www.nature.com/articles/s41467-018-04558-1)
