
## SMR & HEIDI analysis {: .expand}

### SMR

#### \# run SMR and HEIDI test

```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --out mysmr --thread-num 10 
```
**\--bfile** reads individual-level SNP genotype data (in PLINK
binary format) from a reference sample for LD estimation, i.e. .bed,
.bim, and .fam files.

**\--gwas-summary** reads summary-level data from GWAS. The input
format follows that for GCTA-COJO analysis (
<http://cnsgenomics.com/software/gcta/#COJO>).

***mygwas.ma***
```
SNP	A1	A2	freq	b	se	p	n
rs1001	A	G	0.8493 	0.0024 	0.0055 	0.6653	129850
rs1002	C	G	0.03606	0.0034	0.0115	0.7659	129799
rs1003	A	C	0.5128	0.045	0.038	0.2319	129830
......
```
Columns are SNP, the coded allele (also called the effect allele or the reference allele), the other allele (also called the alternative allele), frequency of the effect allele, effect size, standard error, p-value and sample size. The headers are not keywords and will be omitted by the program. <font color='red'>Important: “A1” needs to be the effect allele with “A2” being the other allele and “freq” needs to be the frequency of “A1”</font>. **NOTE:1) For a case-control study, the effect size should be log(odds ratio) with its corresponding standard error. 2) We use the GCTA-COJO format here for compatibility considerations. It is of note that the columns "freq" and "n" will not be used in either SMR or HEIDI analysis and thus can be replaced by "NA" if they are not available.**

**\--beqtl-summary** reads summary-level data from a eQTL study in
binary format. We store eQTL summary data in three separate files
.esi (SNP information, in the same format as the PLINK .bim file),
.epi (probe information) and .besd (eQTL summary statistics in
binary format). See [Data Management](#DataManagement) for more
information. We have prepared the data from the Westra study (Westra
et al. 2013 Nat Genet) in this format, which is available for
download at [Download](#Download).

**\--out** saves the results from the SMR analysis in .smr file
(text format).

***mysmr.smr***
```
ProbeID	Probe_Chr	Gene	Probe_bp	SNP	SNP_Chr	SNP_bp	A1	A2	Freq	b_GWAS	se_GWAS	p_GWAS	b_eQTL	se_eQTL	p_eQTL	b_SMR	se_SMR	p_SMR	p_HEIDI	nsnp_HEIDI
prb01	1	Gene1	1001	rs01	1	1011	C	T	0.95	-0.024	0.0063	1.4e-04	0.36	0.048	6.4e-14	-0.0668	0.0197	6.8e-04	NA	NA
prb02	1	Gene2	2001	rs02	1	2011	G	C	0.0747	0.0034	0.0062	5.8e-01	0.62	0.0396	2e-55	0.0055	0.01	5.8e-01	4.17e-01	28
......
```
Columns are probe ID, probe chromosome, gene name, probe position, SNP name,SNP chromosome, SNP position, the effect (coded) allele, the other allele, frequency of the effect allele (estimated from the reference samples), effect size from GWAS, SE from GWAS, p-value from GWAS, effect size from eQTL study, SE from eQTL study, p-value from eQTL study, effect size from SMR, SE from SMR, p-value from SMR, p-value from HEIDI (HEterogeneity In Depedent Instruments) test, and number of SNPs used in the HEIDI test.

**Missing Value** is represented by "<font color='red'>NA</font>".

**\--thread-num** specifies the number of OpenMP threads for
parallel computing. The default value is 1.

#### \# Specify a method for HEIDI test

```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --heidi-mtd 0 --out mysmr 
```

**\--heidi-mtd** specifies a method for HEIDI test. 0 for the
original HEIDI test approach as in Zhu et al. ([2016 Nature
Genetics](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3538.html)),
and 1 for a new HEIDI test ( beta version for testing). The default
value is 1. The new approach uses up to the top 20 SNPs in the
cis-eQTL region (including the top cis-eQTL) for heterogeneity test
because our latest simulation shows that the power of HEIDI test
increases initially but then decreases with increasing number of
SNPs (m) with a peak at m = \~20.

#### \# Filter SNPs by MAF (in the reference sample)

```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --maf 0.01 --out mysmr  
```

**\--maf** removes SNPs based on a minor allele frequency (MAF)
     threshold in the reference sample.


#### \# Include or exclude a subset of individuals

```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --keep myindi.list --out mysmr 
```

**\--keep** includes a subset of individuals in the reference sample
for analysis.

**\--remove** excludes a subset of individuals in the reference
sample from the analysis.

***myindi.list***

```
 F001 S001
 F002 S002
 F003 S001
 ... 
```

#### \# Include or exclude a subset of eQTL summary data

```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --extract-snp mysnp.list --extract-probe myprobe.list --out mysmr 
```

**\--extract-snp** extracts a subset of SNPs for analysis.

**\--exclude-snp** excludes a subset of SNPs from analysis.

***mysnp.list***
```
rs1001
rs1002
rs1003
...
```

**\--extract-probe** extracts a subset of probes for analysis.

**\--exclude-probe** excludes a subset of probes from analysis.

***myprobe.list***
```
probe1001
probe1002
probe1003
...
```

#### \# Other parameters

```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --peqtl-smr 5e-8 --ld-upper-limit 0.9 --ld-lower-limit 0.05 --peqtl-heidi 1.57e-3 --heidi-min-m 3 --heidi-max-m 20 --cis-wind 2000 --thread-num 5 --out mysmr   
```

**\--peqtl-smr** p-value threshold to select the top associated eQTL
for the SMR test. The default value is 5.0e-8. By default, we only
run the SMR analysis in the cis regions. Please see below for the
SMR analysis in trans regions.

**\--peqtl-heidi** threshold of eQTL p-value to select eQTLs for the
HEIDI test. The default value is 1.57e-3, which is equivalent to a
chi-squared value (df=1) of 10.

**\--ld-upper-limit** LD r-squared threshold used to prune SNPs (eQTLs) in the
HEIDI test, i.e. removing SNPs in very strong LD with the top associated eQTL.
The default value is 0.9.

**\--ld-lower-limit** LD r-squared threshold used to prune SNPs (eQTLs) in the
HEIDI test, i.e. removing SNPs in low LD or not in LD with the top associated eQTL.
The default value is 0.05.

**\--heidi-min-m** minimum requirement of the number of cis-SNPs used in
the HEIDI test. We will skip the HEIDI test if the number of SNPs is
smaller than the threshold. This is because if the number of SNPs is
too small, HEIDI test has little power to detect heterogeneity and
possibly generates misleading result. The default value is 3.

**\--heidi-max-m** maximum number of eQTLs used in
the HEIDI test. If the number of cis-SNPs included in the HEIDI test after LD pruning is larger than m, then only the top m SNPs (ranked by eQTL p-value) will be used in the test. The default value is 20.

**\--cis-wind** defines a window centred around the probe to select
cis-eQTLs (passing a p-value threshold) for the SMR analysis. The
default value is 2000Kb.

#### \# Specify a target SNP for the SMR and HEIDI tests

By default, we use the top cis-eQTL as a target in the SMR analysis,
i.e. using the top cis-eQTL in the SMR test and then using the top
cis-eQTL to test against the other cis-eQTLs in the region for
heterogeneity in the HEIDI test. You can also specific the target by the
following option. Note that this option will ignore p-value specified by
the \--peqtl-smr option (\--peqtl-heidi still applies).

```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --target-snp rs12345 --out mysmr 
```

**\--target-snp** specifies a SNP as the target for the SMR and
HEIDI tests as described above.

#### \# Turn off the HEIDI test

```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --heidi-off --out mysmr 
```
**\--heidi-off** turns off the HEIDI test.


#### \# Run SMR_SO (SMR Sample Overlap)
```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --sample-overlap --pmecs 0.01 --mmecs 2 --out mysmr 
```
**\--sample-overlap** turns on SMR_SO.

**\--pmecs** specifies a p-value threshold to select insignificant SNPs to calculate the correlation between the outcome and the exposure. The default value is 0.01.

**\--mmecs** specifies a minimum SNP number to calculate the correlation. The default value is 2.


### SMR and HEIDI tests in trans regions
The trans-eQTLs are defined as the eQTLs that are more than 5Mb away
from the probe.

```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --out mytrans --trans --trans-wind 1000 
```

**\--trans** turns on SMR and HEIDI tests in trans regions.

**\--trans-wind** defines the size of a window on either side of the top associated
trans-eQTL to select SNPs (passing a p-value threshold) for the SMR
and HEIDI test. The default value is 1000 Kb (i.e. a whole region of 2000 Kb).

***mytrans.smr***

```
ProbeID	Probe_Chr	Gene	Probe_bp	trans_chr	trans_leftBound	trans_rightBound	SNP	SNP_Chr	SNP_bp	A1	A2	Freq	b_GWAS	se_GWAS	p_GWAS	b_eQTL	se_eQTL	p_eQTL	b_SMR	se_SMR	p_SMR	p_HEIDI	nsnp_HEIDI
prb01	1	Gene1	1001	16	5349752	7350902	rs01	16	6349942	C	T	0.131	0.0021	0.0152	8.8e-01	-0.214	0.038	3.26e-08	-0.0098	0.071	8.9e-01	1.73e-1	19
prb01	1	Gene1	1001	21	6443018	8459725	rs02	21	7460164	G	C	0.0747	0.0034	0.0062	5.8e-01	0.62	0.0396	2e-55	0.0055	0.01	5.8e-01	4.17e-01	8
......					
```

Columns are probe ID, probe chromosome, gene name, probe position,
tans-eQTL chromosome, left boundary of the trans-region, right boundary
of the trans-region, SNP name, SNP chromosome, SNP position, the effect
(coded) allele, the other allele, frequency of the effect allele
(estimated from the reference samples), effect size from GWAS, SE from
GWAS, p-value from GWAS, effect size from eQTL study, SE from eQTL
study, p-value from eQTL study, effect size from SMR, SE from SMR,
p-value from SMR, p-value from HEIDI test, and number of SNPs used in
the HEIDI test.

### Multi-SNP-based SMR test

Below shows an option to combine the information from all the SNPs in a
region that pass a p-value threshold (the default value is 5.0e-8 which
can be modified by the flag **\--peqtl-smr**) to conduct a multi-SNP-based SMR
analysis ([Wu et al. 2018 Nature Communications](https://www.nature.com/articles/s41467-018-03371-0)).

The SNPs are pruned for LD using a weighted vertex coverage algorithm
with a LD r2 threshold (the default value is 0.9 which can be modified
by the flag **--ld-pruning**) and eQTL p-value as the weight.

```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --out mymulti --smr-multi 
```

**\--smr-multi** turns on set-based SMR test in the cis-region.

```
smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --out mymulti --smr-multi --set-wind 500 
```

**\--set-wind** defines a window width (Kb) centred around the top
associated cis-eQTL to select SNPs in the cis-region. The default
value is -9 resulting in selecting SNPs in the whole cis-region if
this option is not specified.

### SMR analysis of two molecular traits
Here we provide an option to test the pleotropic association between two
molecular traits using summary data. Take the analysis of DNA
methylation and gene expression data as an example. In this case, we
will need mQTL and eQTL summary data in BESD format. The current version of the program focuses only on the analysis 
in the cis-region, i.e. only testing for associations between genes and DNA methylation sites that are in less than 2 Mb distance.

```
smr --bfile mydata --beqtl-summary myexposure --beqtl-summary myoutcome  --out myomics

```

**\--beqtl-summary** the first one reads mQTL summary data as the
exposure. The second one reads eQTL summary data from as the
outcome.

```
smr --bfile mydata --beqtl-summary myexposure --beqtl-summary myoutcome --extract-exposure-probe myeprobein.list --out myomics

```

**\--extract-exposure-probe** extracts a subset of exposure probes
for analysis.

**\--extract-outcome-probe** extracts a subset of outcome probes for
analysis.

**\--exclude-exposure-probe** excludes a subset of exposure probes
from analysis.

**\--exclude-outcome-probe** excludes a subset of outcome probes
from analysis.

```
smr --bfile mydata --beqtl-summary myexposure --beqtl-summary myoutcome --extract-single-exposure-probe eprobe1 --extract-single-outcome-probe oprobe1 --out myomics

```

**\--extract-single-exposure-probe** extracts a single exposure
probe for analysis.

**\--extract-single-outcome-probe** extracts a single outcome probe
for analysis.

```
smr --bfile mydata --beqtl-summary myexposure --beqtl-summary myoutcome --exclude-single-exposure-probe eprobe1 --exclude-single-outcome-probe oprobe1 --out myomics

```

**\--exclude-single-outcome-probe** excludes a single outcome probe
from analysis.

**\--exclude-single-exposure-probe** excludes a single exposure
probe from analysis.

