
## MeCS

### MeCS: Meta-analysis of cis-eQTL in Correlated Samples

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

#### \# Citation

Qi et al. Identifying gene targets for brain-related traits using
transcriptomic and methylomic data from blood. *Submitted*.

#### \# Tutorial

>Example

``` {.r}
smr --besd-flist my_file.list --mecs --thread-num 5 --out mecs_result 
```

**\--mecs** implements the MeCS analysis.

**\--meta** implements the conventional inverse-variance-weighted
meta-analysis assuming all the cohorts are independent.

>Example

``` {.r}
smr --besd-flist my_file.list --meta --thread-num 5 --out meta_result 
```
