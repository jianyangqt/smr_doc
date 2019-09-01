
## Data Resource {: .expand}

### eQTL summary data 

#### \# Westra eQTL summary data
Westra eQTL summary data [(Westra et al. 2013 Nat Genet)](http://www.ncbi.nlm.nih.gov/pubmed/24013639) in SMR
binary (BESD) format:<br/>
[westra\_eqtl\_data\_hg18.zip (hg18)](../../data/SMR/westra_eqtl_hg18.zip) (10.3 MB)<br/>
[westra\_eqtl\_data\_hg19.zip (hg19)](../../data/SMR/westra_eqtl_hg19.zip) (10.3 MB)


#### \# CAGE eQTL summary data
CAGE eQTL summary data [(Luke R. Lloyd-Jones et al. 2017 AJHG)](http://www.cell.com/ajhg/abstract/S0002-9297(16)30532-8) in SMR binary (BESD) format:<br/>
[cage\_eqtl\_data\_hg19.tgz (hg19)](../../data/SMR/cage_eqtl_data_hg19.tgz) (3.8 GB)<br/>
[cage\_eqtl\_data\_lite\_hg19.tgz (hg19)](../../data/SMR/cage_eqtl_data_lite_hg19.tar.gz) Lite version of the CAGE data (only SNPs with P < 1e-5 are included; 86.1 MB)

The CAGE eQTL results have finer coverage than the Westra et al. 2013
results with comparable power. Please note that the EGCUT cohort is
common to both the Westra et al. 2013 and CAGE data sets. Please see the
above link to the CAGE paper that outlines how these eQTL results were
generated.

Please see the associated [Shiny App](http://cnsgenomics.com/shiny/CAGE/) for further interactive interrogation of the results generated in the CAGE analysis.

#### \# V7 release of the GTEx eQTL summary data
V7 release of the GTEx eQTL summary data [(GTEx Consortium 2017 Nature)](https://www.nature.com/articles/nature24277) in SMR binary (BESD) format:<br/>
[GTEx\_V7\_cis\_eqtl\_summary.tar.gz (hg19)](../../data/SMR/GTEx_V7_cis_eqtl_summary.tar.gz) (55 GB)<br/>
[GTEx\_V7\_cis\_eqtl\_summary\_lite.tar.gz (hg19)](../../data/SMR/GTEx_V7_cis_eqtl_summary_lite.tar.gz) Lite version of the GTEx V7 data (only SNPs with P < 1e-5 are included; 5.3 GB)

This is a set of cis-eQTL summary data across 48 human tissues from the GTEx project. Only SNPs within 1Mb of the transcription start site are available. The standard errors in the BESD files were re-computed from the observed effect sizes and p-values based on a chi-squared distribution with 1 degree of freedom. The forth column of the *.epi file is the middle position of the probe sequence rather than the transcription start site. See [GTEx Portal](http://www.gtexportal.org/) for details about the eQTL analysis.

#### \# Qi et al. brain eQTL summary data
**GTEx-brain eQTL data** (estimated effective n = 233)

GTEx-brain eQTL summary data ([Qi et al. 2018 Nat Commun](https://www.nature.com/articles/s41467-018-04558-1)) in SMR binary (BESD) format: [GTEx-brain.tar.gz](../../data/SMR/GTEx-brain.tar.gz) (1.1 GB)

This is a set of eQTL data from a meta-analysis of 10 brain regions in GTEx v6 ([GTEx Consortium 2017 Nature](https://www.nature.com/articles/nature24277)) correcting for sample overlap by the [MeCS](#MeCS) method. Only SNPs within 1Mb distance from each probe are available. 

**Brain-eMeta eQTL data** (estimated effective n = 1,194)

Brain-eMeta eQTL summary data ([Qi et al. 2018 Nat Commun](https://www.nature.com/articles/s41467-018-04558-1)) in SMR binary (BESD) format: [Brain-eMeta.tar.gz](../../data/SMR/Brain-eMeta.tar.gz) (1.1 GB)

This is a set of eQTL data from a meta-analysis of GTEx brain ([GTEx Consortium 2017 Nature](https://www.nature.com/articles/nature24277)), CMC ([Fromer et al. 2016 Nat Neurosci](https://www.nature.com/articles/nn.4399)), and ROSMAP ([Ng et al. 2017 Nat Neurosci](https://www.nature.com/articles/nn.4632)) by [MeCS](#MeCS). Only SNPs within 1Mb distance from each probe are available.


#### \# Geuvadis eQTL summary data
Geuvadis eQTL summary data [(Lappalainen et al. 2013 Nature)](https://www.nature.com/articles/nature12531) in SMR binary (BESD) format:<br/>
[geuvadis\_EUR\_rsid.tar.gz (hg19)](../../data/SMR/geuvadis_EUR_rsid.tar.gz) (650.5MB)

Unthresholded Geuvadis eQTL data for lymphoblastoid cell lines isolated from 373 EUR individuals were used, with YRI individuals excluded. The eQTL summary data were from [EUR373.gene.K10.noplim.cis\_assembled.txt.gz](http://jungle.unige.ch/~lappalainen/geuvadis/EUR373.gene.K10.noplim.cis_assembled.txt.gz). Since betas are not available for the unthresholded dataset, they were estimated from t-values, allele freqs and n=373 according to Formula (6) in the SMR paper [(Zhu et al. 2016 Nat Genet)](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3538.html). Geuvadis SNP ids were converted to rsids for compatibility with external GWAS data and plink 1KG files. The data are based on GRCh37 assembly, and gene IDs deprecated in the GRCh37 version of Ensembl were removed. We would like to acknowledge [Mikhail Spivakov](mailto:Mikhail.Spivakov@babraham.ac.uk) for his effort in transforming the data to SMR BESD format and writing the description above.

#### \# PsychENCODE eQTL summary data
Here are two sets of cis-eQTL summary data in the prefrontal cortex from the [PsychENCODE project](http://resource.psychencode.org/). Only the data of SNPs in a 1 Mb window around each gene are available.

a.    PsychENCODE eQTL summary data (correcting for 50 PEER factors) in SMR binary (BESD) format: 
[PsychENCODE\_cis\_eqtl\_PEER50\_summary.tar.gz (hg19)](http://cnsgenomics.com/data/SMR/PsychENCODE_cis_eqtl_PEER50_summary.tar.gz) (33 MB)

Only SNPs with FDR < 0.05 are available for each gene. The eQTL analyses were performed with 50 probabilistic estimation of expression residuals (PEER) factors included as covariates (see [Wang et al. 2018 Science](http://science.sciencemag.org/content/362/6420/eaat8464) for details about data generation and analysis). 

b.    PsychENCODE eQTL summary data (correcting for 100 HCP factors) in SMR binary (BESD) format:
[PsychENCODE\_cis\_eqtl\_HCP100\_summary.tar.gz (hg19)](http://cnsgenomics.com/data/SMR/PsychENCODE_cis_eqtl_HCP100_summary.tar.gz) (65 MB)

Similar to the data set above but the eQTL analyses were performed including100 hidden covariate (HCP) factors as covariates (see [Gandal et al. 2018 Science](http://science.sciencemag.org/content/362/6420/eaat8127) for details about data generation and analysis). 


### mQTL summary data 

#### \# McRae et al. mQTL summary data

McRae et al. mQTL summary data ([McRae et al. 2018 Sci Rep](https://www.nature.com/articles/s41598-018-35871-w); [Wu et al. 2018 Nat Commun](https://www.nature.com/articles/s41467-018-03371-0)) in SMR binary (BESD) format (n = 1,980):  
[LBC\_BSGS\_meta.tar.gz (hg19)](../../data/SMR/LBC_BSGS_meta.tar.gz) (7.5 GB)  
[LBC\_BSGS\_meta\_lite.tar.gz (hg19)](../../data/SMR/LBC_BSGS_meta_lite.tar.gz) Lite version of the McRae et al. mQTL data (only SNPs with P < 1e-5 are included; 241 MB)

The original mQTL data were generated in two cohorts BSGS (n = 614) and LBC (n = 1366) in peripheral blood ([McRae et al. 2018 Sci Rep](https://www.nature.com/articles/s41598-018-35871-w)). The methylation states of all the samples which are of European descent were measured based onIllumina HumanMethylation450 chips. The mQTL summary data available here were a meta-analysis of the BSGS and LBC data ([Wu et al. 2018 Nat Commun](https://www.nature.com/articles/s41467-018-03371-0)). Only the DNA methylation probes withat least a cis-mQTL at P < 5e-8 and only SNPs within 2Mb distance from each probe are available.

#### \# Qi et al. brain mQTL summary data

**Brain-mMeta mQTL data** (estimated effective n = 1,160)

Brain-mMeta mQTL summary data ([Qi et al. 2018 Nat Commun](https://www.nature.com/articles/s41467-018-04558-1)) in SMR binary (BESD) format: [Brain-mMeta.tar.gz](../../data/SMR/Brain-mMeta.tar.gz) (893 MB)

This is a set of mQTL data from a meta-analysis of ROSMAP ([Ng et al. 2017 Nat Neurosci](https://www.nature.com/articles/nn.4632)), Hannon et al. ([Hannon et al. 2016 Nat Neurosci](https://www.nature.com/articles/nn.4182)) and Jaffe et al. ([Jaffe et al. 2016 Nat Neurosci](https://www.nature.com/articles/nn.4181)). In the ROSMAP data, only SNPs within 5Kb of each DNA methylation probe are available. In the Hannon et al. data, only SNPs within 500Kb distance from each probe and with PmQTL < 1.0e-10 are available. In the Jaffe et al. data, only SNPs within 20Kb distance from each probe and with FDR < 0.1 are available. 

#### \# Hannon et al. mQTL summary data
Whole blood mQTL data set used in Hannon et al. ([2018 AJHG](https://www.sciencedirect.com/science/article/pii/S0002929718303185?via=ihub)).

[Hannon et al. WholeBlood dataset.zip](https://www.dropbox.com/s/os4cgkb4519wbvn/US_mQTLS_SMR_format.zip?dl=0) (sample size = 1,175; 121MB)


Three mQTL datasets used in Hannon et al. ([2017 AJHG](https://www.sciencedirect.com/science/article/pii/S0002929717301581?via%3Dihub)). All the files are in SMR BESD format.

Two blood mQTL datasets ([Hannon et al. 2016 Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1041-x)) <br/>
[Hannon et al. Blood dataset1.zip](../../data/SMR/Hannon_Blood_dataset1.zip) (sample size = 639; 42MB)<br/>
[Hannon et al. Blood dataset2.zip](../../data/SMR/Hannon_Blood_dataset2.zip) (sample size = 665; 25MB)

Fetal brain mQTL data ([Hannon et al. 2015 Nat Neurosci](https://www.nature.com/articles/nn.4182)) <br/>
[Hannon et al. FetalBrain.zip](../../data/SMR/Hannon_FetalBrain.zip) (sample size = 166; 4.8MB)

Note: 1) The SNPs are coded as chr:bp (based on the genome build hg19)
rather than with rsIDs. 2) SNPs with mQTL p-values > 1e-10 are not included. 3) Any question regarding to these datasets should be addressed to [Eilis Hannon](mailto:E.J.Hannon@exeter.ac.uk) or [Jonathan Mill](mailto:J.Mill@exeter.ac.uk). 

