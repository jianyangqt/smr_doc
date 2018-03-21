
## Data Management {: .expand}

eQTL summary data are usually generated from association tools such as
PLINK and stored in separate files in text format (usually one file for
each probe) with a very large file size in total. Here we provide an
efficient way to store the eQTL summary data in binary format (BESD),
with flexible options to query the data for any subset of SNPs and/or
probes (see [Query eQTL Results](#QueryeQTLResults)). We further provide
a sparse version of the BESD format (used as default), which is
extremely storage-efficient without losing too much information. The
basic idea is that we store summary data of all SNPs within 2Mb of a
probe in either direction, all SNPs within 1Mb of any trans-eQTL in
either direction, and all SNPs with p &lt; 1e-5 in the rest of the genome
(note that all these parameters can be specified by users). We also
provide several options to import data in various of formats (e.g.
PLINK, GEMMA, BOLT-LMM and other text formats).


### BESD format

#### \# BESD format: an efficient format to store eQTL summary data

We store eQTL summary data in three separate files .esi (SNP
information, similar as a PLINK .bim file), .epi (probe information) and
.besd (a binary file to store the summary statistics).

***myeqtl.esi***
```
1	rs1001	0	744055	A	G	0.23
1	rs1002	0	765522 	C	G	0.06
1	rs1003	0	995669	T	C	0.11
......
```
Columns are chromosome, SNP, genetic distance (can be any arbitary value
since it will not be used in the SMR analysis), basepair position, the
effect (coded) allele, the other allele and frequency of the effect
allele.

***myeqtl.epi***
```
1	probe1001	0	924243	Gene01	+
1	probe1002	0	939564	Gene02	-
1	probe1003	0	1130681	Gene03	-
......
```
Columns are chromosome, probe ID(can be the ID of an exon or a
transcript for RNA-seq data), genetic distance (can be any arbitary
value), physical position, gene ID and gene orientation (this
information will only be used for graphic presentation, please see [[SMR plot]](#SMRlocusplot)).

***myeqtl.besd***

eQTL summary-level statistics (effect size and SE) in binary format.
<font color='red'>Please do not try to open this file with a text editor</font>.

Given the large numbers of SNPs and probes, the size of a .besd file
will still be very large. Since the eQTL association signals are highly
enriched in the cis-region and often there are not many trans-eQTLs, we
could reduce the size of the .besd file by orders of magnitude if we
only store the data for SNPs within 2Mb of a probe in either direction,
SNPs within 1Mb of any trans-eQTL in either direction, and SNPs with p &lt; 1e-5 
in the rest of the genome (see below for options to change these
paramters). We call this the sparse BESD format.

We only store effect size (b) and SE in the BESD file, and re-calculate
p-value for analysis when necessary, assuming b / SE follows a standard
normal distribution, i.e. N(0, 1). Strictly speaking, b / SE follows a
t-distribution which is approximately N(0, 1) if sample size is large.
For data sets with small sample sizes (e.g. GTEx), this might lead to a bias
in p-value. In this scenario, we recommend users to compute z\* based on the original p-value from
a standard normal distribution, and adjust the standard error as SE = b
/ z\*. This adjustment guarantees that the re-computed p-value from b
and SE being exact the same as the original p-value.

See [below](#MakeaBESDfile) for options to make a BESD file from data in several different
formats.

### Make a BESD file

We provide eight different ways of converting cis-eQTL data in other formats to BESD format.

#### 1. Make a BESD file from eQTL summary data in ESD format

To compile data in sparse BESD format

```
smr --eqtl-flist my.flist --make-besd --out mybesd 
```

**\--eqtl-flist** reads a file to get probe information and file
paths of the eQTL summary data.

**\--make-besd** saves summary data in BESD format. By default, the
data will be stored in sparse BESD format (See below for the option
**\--make-besd-dense** to store the data in dense BESD format). By
default, the data will be stored in sparse BESD format if the
sparsity given the parameters (by default, ±2Mb of the cis-region,
±1Mb of any trans-eQTL and all SNP at p &lt; 1e-5) is lower than 0.4. It
will also output a text file (.summary) to summarise the genomic
regions stored in the .besd file (sparse format) for each probe.

***my.flist***
```
Chr	ProbeID	GeneticDistance	ProbeBp	Gene	Orientation	PathOfEsd
9	cg00000658	0	139997924	MAN1B1	-	path/my01.esd
20	cg26036652	0	33735834 	NA	NA	path/my02.esd
1	cg00489772	0	3775078	NA	NA	path/my03.esd
......
```
This is a text file **with headers**. The first 6 columns are the same as in .epi. The last column is the full path of an eQTL summary data file (.esd file, see below for the format of a .esd file). 

***my01.esd***
```
Chr	SNP	Bp	A1	A2	Freq	Beta	se	p
9	rs12349815	150048	T	A	0.968	0.019	0.016	0.2434
20	rs141129176	62955484	G	A	0.89	0.012	0.009	0.2156
......
```
This is a text file **with headers**. Columns are chromosome, SNP, the
effect (coded) allele, the other allele, basepair position, frequency of
the effect allele, effect size, standard error and p-value.

**HINT** : if the SNPs in all of the .esd files are identical, the
efficiency of the analysis can be largely improved by adding the
\--geno-uni option. This option call be used in all the commands of this
section.

```
smr --eqtl-flist my.flist --make-besd --geno-uni --out mybesd 
```

**\--geno-uni** indicates all the input .esd files are identical.

To compile eQTL summary data in sparse BESD format with user-specified
parameters

```
smr --eqtl-flist my.flist --make-besd  --cis-wind 2000 --trans-wind 1000 --peqtl-trans 5.0e-8 --peqtl-other 1.0e-5 --out mybesd 
```

**\--cis-wind** specifies a window (in Kb unit) to store all the
SNPs within the window of the probe in either direction. The default
value is 2000Kb.

**\--trans-wind** specifies a window (in Kb unit) to store all the
SNPs in a trans-region. If there is a trans-eQTL with p-value &lt; the
specified threshold (\--peqtl-trans), it will store all the SNPs
within the window of the top associated trans-eQTL in either
direction. The default value is 1000Kb.

**\--peqtl-trans** p-value threshold for trans-eQTLs. The default
value is 5.0e-8.

**\--peqtl-other** Apart from the cis and trans regions, it will
also store all SNPs with eQTL p-values &lt; this threshold. The
default value is 1.0e-5 .

To compile the eQTL summary data in dense BESD format

```
smr --eqtl-flist my.flist --make-besd-dense --out mybesd 
```
**\--make-besd-dense** saves summary data of all SNPs for all probes.

**WARNING** : This will generate a huge file.

**NOTE** : the **\--make-besd-dense** option can be used in all the commands above and below.

#### 2. Make a BESD file from Matrix eQTL output 

```
smr --eqtl-summary mateQTL.txt --matrix-eqtl-format --make-besd --out mybesd 
```
**\--eqtl-summary** reads eQTL summary statistics in text format or compressed text format (e.g. *.tar.gz file).

**\--matrix-eqtl-format** indicates eQTL summary data in Matrix eQTL output format.

***mateQTL.txt***
```
SNP	gene	beta	t-stat	p-value	FDR
rs13258200	ENSG00000071894.10	-1.00783189089702	-16.641554315712	2.3556801409867e-24	1.12905157909337e-18
rs6599528	ENSG00000071894.10	-1.06253739134798	-15.8412867110849	2.73027622294589e-23	5.51886367106636e-18
rs2272666	ENSG00000071894.10	-1.04810713295621	-15.6736668186937	4.6058755123246e-23	5.51886367106636e-18
rs4313195	ENSG00000071894.10	-1.04810713295621	-15.6736668186937	4.6058755123246e-23	5.51886367106636e-18
rs2280836	ENSG00000071894.10	-1.00773805984667	-15.2332537951202	1.84980830591084e-22	1.77318554626341e-17
...
```
This file has headers. The six columns are SNP, gene, beta (i.e. SNP effect on gene expression), t-statistic, p-value and q-value
([[http://www.bios.unc.edu/research/genomic\_software/Matrix\_eQTL/]](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/)).

**NOTE** : 1) The program is able to read \*.tar.gz file. 2) The SNP and probe information in the SMR eQTL output files (.esi and .epi) converted from Matrix eQTL output are not complete and need to be updated using the options in [Update a BESD file](#UpdateaBESDfile).

#### 3. Make a BESD file from FastQTL output

```
smr --eqtl-summary fastqtlnomi.txt --fastqtl-nominal-format --make-besd --out mybesd 
```
**\--fastqtl-nominal-format** indicates eQTL summary data in FastQTL "nominal pass" output format.

***fastqtlnomi.txt***
```
ENSG00000237438.1 indel:4D_22_16518157 -999303 0.542909 -0.0510761
ENSG00000237438.1 snp_22_16519289 -998171 0.432764 0.124424
ENSG00000237438.1 snp_22_16520141 -997319 0.0945196 -0.251906
ENSG00000237438.1 snp_22_16520948 -996512 0.102846 -0.274157
ENSG00000237438.1 snp_22_16523696 -993764 0.0676318 -0.324492
ENSG00000237438.1 snp_22_16523730 -993730 0.0674215 -0.206578
...
```
This file has **no header**. the five columns are gene, SNP, distance in bp between the SNP and the gene, p-value and beta (i.e. SNP effect on gene expression) 
([[http://fastqtl.sourceforge.net/pages/cis_nominal.html]](http://fastqtl.sourceforge.net/pages/cis_nominal.html)).

**NOTE** : 1) The FastQTL output file should be generated by a FastQTL version higher than v2.184. 2) The program is able to read \*.tar.gz file. 3) The SNP and probe information in the SMR output files (.esi and .epi) converted from FastQTL output are not complete and need to be updated using the options in [Update a BESD file](#UpdateaBESDfile).

#### 4. Make a BESD file from eQTL summary data in PLINK-qassoc output format

The output file from a PLINK \--assoc analysis does not contain allele
information. We therefore need to read the alleles from a PLINK .bim
file. The file path of the PLINK .bim file needs to be added as the last
column of the .flist file (see the example below).

```
smr --eqtl-flist my.flist  --plink-qassoc-format --make-besd --out mybesd 
```

**\--plink-qassoc-format** reads eQTL summary data in PLINK-qassoc output format (output file from a PLINK \--assoc analysis for a quantitative trait).


 ***my.flist***
```
Chr	ProbeID	GeneticDistance	ProbeBp	Gene	Orientation	PathOfEsd	PathOfBim
9	cg00000658	0	139997924	MAN1B1	-	path_assoc/my01.qassoc	path_genotype/chr9
20	cg26036652	0	33735834 	NA	NA	path_assoc/my02.qassoc	path_genotype/chr20
1	cg00489772	0	3775078	NA	NA	path_assoc/my03.qassoc	path_genotype/chr19
......	
```
**NOTE** : The program is able to read \*.tar.gz file, e.g.
path\_assoc/my03.qassoc.tar.gz

#### 5. Make a BESD file from eQTL summary data in GEMMA output format

```
smr --eqtl-flist my.flist --gemma-format --make-besd --out mybesd 
```

**\--gemma-format** reads eQTL summary data in GEMMA association output format

 
```
chr	rs	ps	n_miss	allel1	allel0	af	beta	se	l_remle	p_wald
1	rs3683945	3197400	0 	A	G	0.443	-7.788665e-02	6.193502e-02	4.317993e+00	2.087616e-01
1	rs3707673	3407393	0	G	A	0.443	-6.654282e-02	6.210234e-02	4.316144e+00	2.841271e-01
1	rs6269442	3492195	0	A	G	0.365	-5.344241e-02	5.377464e-02	4.323611e+00	3.204804e-01
......	
```
The 11 columns are: chromosome, SNP ID, basepair position, number of
missing values for a given SNP, the effect (coded) allele, the other
allele, frequency of the effect allele, effect size, standard error,
lambda and p-value
([[http://www.xzlab.org/software.html]](http://www.xzlab.org/software.html)).

#### 6. Make a BESD file from eQTL summary data in BOLT-LMM output format

```
smr --eqtl-flist my.flist  --bolt-assoc-format --make-besd --out mybesd 
```

**\--bolt-assoc-format** reads eQTL summary data in BOLT\_LMM output
format
```
SNP	CHR	BP	GENPOS	ALLELE1	ALLELE0	A1FREQ	F_MISS	BETA	SE	P_BOLT_LMM_INF	P_BOLT_LMM
rs58108140	1	10583	0.000000 	A	G	0.109810	0.011935	0.074942	0.045043	9.6E-02	9.7E-02
rs180734498	1	13302	0.000000	T	C	0.061042	0.007595	0.084552	0.058078	1.5E-01	1.4E-01
rs151118460	1	91581	0.000000	A	G	0.399377	0.013382	0.024344	0.034394	4.8E-01	4.8E-01
......
```
The 12 columns are: SNP ID, chromosome, basepair position, genetic
position, the effect (coded) allele, the other allele, frequency of the
effect allele, fraction of individuals with missing genotype at the SNP,
effect size, standard error, infinitesimal model (mixture model)
association test p-value, and non-infinitesimal model association test
p-value
([[https://data.broadinstitute.org/alkesgroup/BOLT-LMM/\#x1-440008.1]](https://data.broadinstitute.org/alkesgroup/BOLT-LMM/#x1-440008.1)).

#### 7. Make a BESD file from SMR query output

```
smr --qfile myquery.txt --make-besd --out mybesd 
```

**\--qfile** reads eQTL summary data in SMR query output format (see
[Query eQTL Results for the format of a query output
file](#QueryeQTLResults)).

#### 8. Make a BESD file from BESD file(s)

To make a sparse BESD file from a single dense BESD file

```
smr --beqtl-summary my_beqtl --make-besd --out my_sparse 
```

```
smr --beqtl-summary my_beqtl --cis-wind 2000 --trans-wind 1000 --peqtl-trans 5.0e-8 --peqtl-other 1.0e-5 --make-besd --out my_sparse
```

To make a sparse BESD file from multiple sparse or dense BESD files (can
be a mixture of both types)

```
smr --besd-flist my_file.list --make-besd --out my_sparse  
```

**\--besd-flist** reads a file to get the full paths of the BESD
files.

***my_file.list***
```
path1/my_besd1
path2/my_besd2
path3/my_besd3
...
```
**NOTE** : this command can be used to merge multiple BESD files.

**HINT** : if the SNPs in all the .esi files are identical, you can speed up
the analysis using the **\--geno-uni** option.

### Update a BESD file

Some of the information such as SNP chromosome, SNP position, the effect allele, the other allele, probe chromosome and probe position might be missing in the output files generated by some of the eQTL analysis tools (e.g.  Matrix eQTL and FastQTL). These information are necessary for the SMR analysis and some of the data management operations in the SMR tool, and thus need to be updated in the .esi or .epi file.

Users can update the .esi file and .epi file manually, but please be aware of not changing the order of SNPs in the .esi file or probes in the .epi file because they are associated with the information in the .besd file. We also provide the options below to update the information in the .esi  and .epi files.

#### \# Update the .esi or .epi files

```
smr --beqtl-summary my_beqtl --update-esi mybigpool.esi 
```
**\--update-esi** reads a .esi file for updating and backup the original .esi file.

```
smr --beqtl-summary my_beqtl --update-epi mybigpool.epi
```
**\--update-epi** reads a .epi file for updating and backup the original .epi file.

#### \# Add or update the frequencies of the effect alleles

```
smr --beqtl-summary myeqtl --update-freq mysnp.freq 
```

**\--update-freq** reads an input file with allele frequency
information and adds a new column (i.e. frequency the effect allele)
to the .esi file.

***mysnp.freq***
```
rs12349815	T	A	0.968
rs141129176	G	A	0.89
......
```
The input is a text file **without headers**. Columns are SNP, the
effect allele, the other allele and frequency of the effect allele.

**NOTE** : the SMR program is compatible with .esi files with or without
frequency information.

### Mange the sample size

#### \# Add or update the sample size in the BESD file

```
smr --beqtl-summary myeqtl --add-n 1000 --make-besd --out mybesd

```
**\--add-n** reads the sample size.

**NOTE** : The flag is valid for all the options to make a BESD file. For example: 

```
smr --qfile myquery.txt --add-n 100 --make-besd --out mybesd 
```

#### \# Show the sample size

```
smr --beqtl-summary myeqtl --show-n 
```
**\--show-n** shows the sample size on the screen or in the log output.


### Extract/remove a subset of data

#### \# Extract a subset of SNPs and/or probes

```
smr --beqtl-summary myeqtl --extract-snp mysnp.list --extract-probe myprobe.list  --make-besd --out mybesd 
```

#### \# Remove a subset of SNPs and/or probes

```
smr --beqtl-summary myeqtl --exclude-snp mysnp.list --exclude-probe myprobe.list  --make-besd --out mybesd 
```

#### \# Extract a subset of SNPs with p < a threshold

```
smr --beqtl-summary myeqtl --extract-snp-p 1e-5  --make-besd --out mybesd 
```
**\--extract-snp-p** reads a p-value threshold to extract SNPs

#### \# Remove a subset of SNPs with p < a threshold

```
smr --beqtl-summary myeqtl --exclude-snp-p 1e-5  --make-besd --out mybesd 
```
**\--exclude-snp-p** reads a p-value threshold to exclude SNPs

#### \# Extract cis-regions of eQTL summary data

```
smr --beqtl-summary myeqtl --extract-cis --make-besd --out mybesd
```
**\--extract-cis** extracts the cis-eQTL summary data.

### Remove technical eQTLs

#### \# Remove technical eQTLs

Filtering out eQTLs for which there is a significant cis-eQTL in the
hybridization region of the probe. This option will remove all the SNPs
in the cis-region of the probe and save the removed data in a file in
SMR Query format (see [Query eQTL Results for the format of a query
output file](#QueryeQTLResults)). The default p-value threshold is 5e-8,
which can be changed by the **\--p-technical** (see below).

```
smr --beqtl-summary myeqtl --rm-technical probe_hybrid.txt --make-besd --out mybesd

```

**\--rm-technical** specifies the probe hybridization region and
excludes the technical eQTLs.

***probe_hybrid.txt***
```
19	probe0	50310094	50310143
19	probe1	406496	406545
10	probe2	119293020	119293069
......	
```
This is a text file **without headers**. Columns are chromosome, probe
ID, start of the hybridization region and end of the hybridization
region.

```

smr --beqtl-summary myeqtl --rm-technical probe_hybrid.txt --p-technical 5e-8 --make-besd --out mybesd

```

**\--p-technical** reads a p-value threshold to select technical
eQTLs. The default value is 5e-8.
