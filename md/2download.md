
## Download {: .expand}
### Executable Files (version 0.70) 

[smr\_Linux.zip](download/smr_Linux.zip)

[smr\_Mac.zip](download/smr_Mac.zip)

[smr\_Win.zip](download/smr_Win.zip)

The executable files (binary code) are released under MIT lincense.

### eQTL summary data 

#### \# Westra eQTL summary data

Westra eQTL summary data [([Westra et al. 2013 Nat Genet])](http://www.ncbi.nlm.nih.gov/pubmed/24013639) in SMR
binary (BESD) format:

[westra\_eqtl\_data\_hg18.zip (hg18)](download/westra_eqtl_hg18.zip)

[westra\_eqtl\_data\_hg19.zip (hg19)](download/westra_eqtl_hg19.zip)


#### \# CAGE eQTL summary data

CAGE eQTL summary data [([Luke R. Lloyd-Jones et al. 2017 AJHG])](http://www.cell.com/ajhg/abstract/S0002-9297(16)30532-8) in SMR binary (BESD) format:

[cage\_eqtl\_data\_hg19.tgz (hg19)](http://www.cnsgenomics.com/data/CAGE/cage_eqtl_data_hg19.tgz)

The CAGE eQTL results have finer coverage than the Westra et al. 2013
results with comparable power. Please note that the EGCUT cohort is
common to both the Westra et al. 2013 and CAGE data sets. Please see the
above link to the CAGE paper that outlines how these eQTL results were
generated.

Please see the associated [Shiny App](http://cnsgenomics.com/shiny/CAGE/) for further interactive
interrogation of the results generated in the CAGE analysis.

#### \# McRae et al. mQTL summary data

McRae et al. mQTL summary data [([Allan McRae et al. 2017 bioRxiv])](https://www.biorxiv.org/content/early/2017/07/21/166710) in SMR binary (BESD) format:

[LBC\_BSGS\_meta.tar.gz (hg19)](http://cnsgenomics.com/data/SMR/LBC_BSGS_meta.tar.gz)

This is a set of mQTL summary data from a meta-analysis of the BSGS and LBC data (McRae et al. 2017 bioRxiv). Only the DNA methylation probes with at least an cis-mQTL at P < 5e-8 and only SNPs within 2Mb distance from each probe are available. 

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

**18.**   Version 0.70 (8 December, 2017): updated some parameters in HEIDI test. 
A lower limit of LD R-square threshold (the default value is 0.05) has been implmented. The best number of SNPs which contain the top SNP for the HEIDI test is changed to 31.

