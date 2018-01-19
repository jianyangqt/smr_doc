
## Download {: .expand}
### Executable Files (version 0.703) 

[smr\_Linux.zip](download/smr_Linux.zip)

[smr\_Mac.zip](download/smr_Mac.zip)

[smr\_Win.zip](download/smr_Win.zip)

The executable files (binary code) are released under MIT lincense.

### eQTL summary data 

#### \# Westra eQTL summary data

Westra eQTL summary data [(Westra et al. 2013 Nat Genet)](http://www.ncbi.nlm.nih.gov/pubmed/24013639) in SMR
binary (BESD) format:

[westra\_eqtl\_data\_hg18.zip (hg18)](https://nextcloud.qriscloud.org.au/index.php/s/U6oh8SKSWfyfhT4) (10.3 MB)

[westra\_eqtl\_data\_hg19.zip (hg19)](https://nextcloud.qriscloud.org.au/index.php/s/b4uHQzoKAbGbSFc) (10.3 MB)


#### \# CAGE eQTL summary data

CAGE eQTL summary data [(Luke R. Lloyd-Jones et al. 2017 AJHG)](http://www.cell.com/ajhg/abstract/S0002-9297(16)30532-8) in SMR binary (BESD) format:

[cage\_eqtl\_data\_hg19.tgz (hg19)](https://nextcloud.qriscloud.org.au/index.php/s/52V62aEFNxrf8h7) (3.8 GB)

[cage\_eqtl\_data\_lite\_hg19.tgz (hg19)](https://nextcloud.qriscloud.org.au/index.php/s/nvGmXefCxiHPAz4) Lite version of the CAGE data (only SNPs with P < 1e-5 are included; 86.1 MB)

The CAGE eQTL results have finer coverage than the Westra et al. 2013
results with comparable power. Please note that the EGCUT cohort is
common to both the Westra et al. 2013 and CAGE data sets. Please see the
above link to the CAGE paper that outlines how these eQTL results were
generated.

Please see the associated [Shiny App](http://cnsgenomics.com/shiny/CAGE/) for further interactive
interrogation of the results generated in the CAGE analysis.

#### \# McRae et al. mQTL summary data

McRae et al. mQTL summary data [(Allan McRae et al. 2017 bioRxi)](https://www.biorxiv.org/content/early/2017/07/21/166710) in SMR binary (BESD) format:

[LBC\_BSGS\_meta.tar.gz (hg19)](https://nextcloud.qriscloud.org.au/index.php/s/v5WoHlCHGt2mm3r) (7.4 GB)

This is a set of mQTL summary data from a meta-analysis of the BSGS and LBC data (McRae et al. 2017 bioRxiv). Only the DNA methylation probes with at least an cis-mQTL at P < 5e-8 and only SNPs within 2Mb distance from each probe are available. 

#### \# V7 release of the GTEx eQTL summary data

V7 release of the GTEx eQTL summary data [(GTEx Consortium 2017 Nature)](https://www.nature.com/articles/nature24277) in SMR binary (BESD) format:

[GTEx\_V7\_cis_eqtl\_summary.tar.gz (hg19)](https://nextcloud.qriscloud.org.au/index.php/s/ppnjr9vIpcCeMsR) (53.6 GB)

This is a set of cis-eQTL summary data across 48 human tissues from the GTEx project. Only SNPs within 1Mb of the transcription start site are available. The standard errors in the BESD files were re-computed from the observed effect sizes and p-values based on a chi-squared distribution with 1 degree of freedom. The forth column of the *.epi file is the middle position of the probe sequence rather than the transcription start site. See [GTEx Portal](http://www.gtexportal.org/) for details about the eQTL analysis.

### R script for SMR locus plot 

R script and sample file for SMR locus plot:

[plot.zip](download/plot.zip)


### Update log 


**1.**  24 Aug, 2015: first release.

**2.**  17 Sept, 2015: updated the format of sparse besd file; added a
function to make sparse besd file by extracting information from full
dense besd file; added a function to check quickly how many probes are
associated with a SNP at p &lt; a threshold(e.g. 5e-8).

**3.**  12 Oct, 2015: Eigen library and OpenMP were used.

**4.**   Version 0.6 (10 Nov, 2015): added features of SMR and HEIDI
test for the trans regions.

**5.**   Version 0.619 (4 Aril, 2016): updated sparse besd format;
updated features to make sparse verison of BESD; added features to query
eQTL summary results; added features to combine BESD files.

**6.**   Version 0.620 (12 Aril, 2016): added a feature to deal with
duplicate IDs.

**7.**   Version 0.628 (11 May, 2016): added a feature to visualize SMR
results.

**8.**   Version 0.630 (23 May, 2016): updated features to make binary
besd file from plain text file(s).

**9.**   Version 0.631 (23 June, 2016): more options to make BESD files
and more memory-efficient when making binary besd files.

**10.**   Version 0.632 (28 June, 2016): added a feature to make a BESD
file from BOLT-LMM output format.

**11.**   Version 0.64 (8 August, 2016): updated the .esi file format;
updated the HEIDI test (a new method that improves the power of the
HEIDI test); updated the SMR query output format; improved the analysis
to combine multiple BESD files.

**12.**   Version 0.65 (12 December, 2016): added a flag (\--heidi-mtd)
for users to choose the original approach or a new approach for HEIDI
test.

**13.**   Version 0.66 (10 January, 2017): updated the function to
generate the file for locus plot. The new version is able to read a gene
list with/without strand information.

**14.**   Version 0.67 (22 June, 2017): updated the functions to make
BESD file by the following strategy: 1) Z\* from N(0, 1) given the
p-value. 2) SE\* = b / Z\* . 3) store b and SE\* in BESD. This
adjustment guarantees that the re-computed p-value from b and SE being
exact the same as the original p-value, useful for data with small
sample size.

**15.**   Version 0.68 (11 August, 2017): updated the SMR and HEIDI
tests in the trans regions (the previous version focuses only on the top
trans-eQTL locus and the new version will run the tests for all the
trans-eQTL loci one at a time).

**16.**   (12 September, 2017): Luke R. Lloyd-Jones et al. released CAGE
eQTL summary statistics for SMR analysis.

**17.**   Version 0.69 (7 October, 2017): added features to run
multi-SNP based SMR and SMR analysis of two molecular traits. Also add a
feature to remove technical eQTLs.

**18.**   Version 0.701 (22 December, 2017): 1) Updated some parameters used in the HEIDI test. A lower limit of LD r-squared threshold (the default value is 0.05) has been added to remove SNPs that are not in LD or in low LD with the top eQTL. 2) Added a flag --heidi-max-m to specify the maximum number of SNPs used in the HEIDI test. 

**19.**   (02 January, 2018): Released the GTEx eQTL summary data in SMR binary (BESD) format.

**20.**   Version 0.702 (07 January, 2018): added flags --extract-snp-p and --exclude-snp-p to make a subset of BESD with a p-value threshold. 

**21.**   Version 0.703 (13 January, 2018): 1) added flags --add-n to add the sample size to BESD file and --show-n to display the sample size on screen. 2) added flag --matrix-eqtl-format to transform the eQTL summary statistics in Matrix eQTL output format to SMR BESD format. 3) added flags --fastqtl-nominal-format and --fastqtl-permu-format to transform the eQTL summary statistics in FastQTL outformat to SMR BESD format. 4) added flags --update-epi and --update-esi to update or complete the information in .epi file and .esi file respectively.

