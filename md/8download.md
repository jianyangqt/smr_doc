
## Download {: .expand}
### Executable Files (version 1.0) 

[smr\_Linux.zip](download/smr_Linux.zip)

[smr\_Mac.zip](download/smr_Mac.zip)

[smr\_Win.zip](download/smr_Win.zip)


The executable files (binary code) are released under MIT license.

### Source code Files (version 1.0) 
[smr\_src.zip](download/smr_src.tar.gz)

The source code are released under GPL v3.

### R script for SMR locus plot 

R script and sample file for SMR locus plot:

[plot.zip](download/plot.zip)


### Update log 
**31.**   Version 1.0 (15 January, 2019): 1) Formally released SMR version 1.0; 2) Added flags to query descriptive summary of  the cis-eQTL and trans-eQTL data.

**30.**   Version 0.712 (6 September, 2018): Added an option to make and query LD information in binary format (i.e. .bld, and .esi files) and an interface to read LD matrix (in BLD format) as the reference for the HEIDI test.

**29.**   Version 0.711 (30 August, 2018): Added flags --qtltoos-nominal-format and --qtltoos-permu-format to transform eQTL summary statistics in QTLtools output format to SMR BESD format.

**28.**   Version 0.710 (21 June, 2018): improved the allele frequency checking step of SMR (more flexible than the previous version). Note that we also updated the frequencies of the effect alleles in the McRae et al. mQTL data, GTEx eQTL data and Brain-mMeta mQTL data.

**27.**   Version 0.709 (14 June, 2018): added a QC step to check the differences in allele frequency among eQTL, GWAS and LD reference data for the SMR analysis.

**26.**   Version 0.708 (31 May, 2018): added two flags --extract-target-snp-probe and --extract-snp-probe to extract specific SNP-probe pairs for the SMR analysis.

**25.**   Version 0.707 (11 May, 2018): changed the default LD pruning r2 threshold for SMR-multi from 0.9 to 0.1.

**24.**   Version 0.706 (1 April, 2018): fixed a problem of linking a library in Linux and a bug with --target-snp. 

**23.**   Version 0.705 (16 Feburary, 2018): Updated MeCS and removed some confusing log information when making BESD files. 

**22.**   Version 0.704 (01 Feburary, 2018): Fixed a bug in making BESD files. 

**21.**   Version 0.703 (19 January, 2018): 1) Added flags --add-n to add sample size in a BESD file and --show-n to display the sample size on screen or in a log output. 2) Added flag --matrix-eqtl-format to transform the eQTL summary statistics in Matrix eQTL output format to SMR BESD format. 3) Added flag --fastqtl-nominal-format to transform the eQTL summary statistics in FastQTL outformat to SMR BESD format. 4) Added flags --update-epi and --update-esi to update or complete the information in .epi file and .esi file respectively.

**20.**   Version 0.702 (07 January, 2018): Added flags --extract-snp-p and --exclude-snp-p to make a subset of BESD with a p-value threshold. 

**19.**   (02 January, 2018): Released the GTEx eQTL summary data in SMR binary (BESD) format.

**18.**   Version 0.701 (22 December, 2017): 1) Updated some parameters used in the HEIDI test. A lower limit of LD r-squared threshold (the default value is 0.05) has been added to remove SNPs that are not in LD or in low LD with the top eQTL. 2) Added a flag --heidi-max-m to specify the maximum number of SNPs used in the HEIDI test. 

**17.**   Version 0.69 (7 October, 2017): added features to run
multi-SNP based SMR and SMR analysis of two molecular traits. Also add a
feature to remove technical eQTLs.

**16.**   (12 September, 2017): Luke R. Lloyd-Jones et al. released CAGE
eQTL summary statistics for SMR analysis.

**15.**   Version 0.68 (11 August, 2017): updated the SMR and HEIDI
tests in the trans regions (the previous version focuses only on the top
trans-eQTL locus and the new version will run the tests for all the
trans-eQTL loci one at a time).

**14.**   Version 0.67 (22 June, 2017): updated the functions to make
BESD file by the following strategy: 1) Z\* from N(0, 1) given the
p-value. 2) SE\* = b / Z\* . 3) store b and SE\* in BESD. This
adjustment guarantees that the re-computed p-value from b and SE being
exact the same as the original p-value, useful for data with small
sample size.

**13.**   Version 0.66 (10 January, 2017): updated the function to
generate the file for locus plot. The new version is able to read a gene
list with/without strand information.

**12.**   Version 0.65 (12 December, 2016): added a flag (\--heidi-mtd)
for users to choose the original approach or a new approach for HEIDI
test.

**11.**   Version 0.64 (8 August, 2016): updated the .esi file format;
updated the HEIDI test (a new method that improves the power of the
HEIDI test); updated the SMR query output format; improved the analysis
to combine multiple BESD files.

**10.**   Version 0.632 (28 June, 2016): added a feature to make a BESD
file from BOLT-LMM output format.

**9.**   Version 0.631 (23 June, 2016): more options to make BESD files
and more memory-efficient when making binary besd files.

**8.**   Version 0.630 (23 May, 2016): updated features to make binary
besd file from plain text file(s).

**7.**   Version 0.628 (11 May, 2016): added a feature to visualize SMR
results.

**6.**   Version 0.620 (12 Aril, 2016): added a feature to deal with
duplicate IDs.

**5.**   Version 0.619 (4 Aril, 2016): updated sparse besd format;
updated features to make sparse verison of BESD; added features to query
eQTL summary results; added features to combine BESD files.

**4.**   Version 0.6 (10 Nov, 2015): added features of SMR and HEIDI
test for the trans regions.

**3.**  12 Oct, 2015: Eigen library and OpenMP were used.

**2.**  17 Sept, 2015: updated the format of sparse besd file; added a
function to make sparse besd file by extracting information from full
dense besd file; added a function to check quickly how many probes are
associated with a SNP at p &lt; a threshold(e.g. 5e-8).

**1.**  24 Aug, 2015: first release.



