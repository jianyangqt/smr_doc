
## Query eQTL Results

### Query eQTL Summary Results

Since the eQTL summary are stored in binary format for a large number of
probes and SNPs, we provide the options below to query the eQTL summary
results from command line options or from file-list options given a
specified eQTL p-value threshold.

#### \# Command line options for SNPs

To query the eQTL resutls for a single SNP, we could use this command

```
smr --beqtl-summary myeqtl --query 5.0e-8 --snp rs123 --out myquery 
```

**\--query** saves in text format a subset of the eQTL summary
dataset based on the specified eQTL p-value threshold. The default
value is 5.0e-8.

**\--snp** specifies a single SNP.

***myquery.txt***

```
SNP	Chr	BP	A1	A2	Freq	Probe	Probe_Chr	Probe_bp	Gene	Orientation	b	se	p
rs01	1	1001	A	G	0.23	cg01	1	1101	gene1	+	-0.033	0.006	3.8e-08
rs01	1	1001	A	G	0.06	cg02	1	1201	gene2	-	0.043	0.007	8.1e-10
......	
```
To query eQTL resutls for a range of SNPs in a genomic region

```
smr --beqtl-summary myeqtl --query 5.0e-8 --from-snp rs123 --to-snp rs456 --out myquery 
```

**\--from-snp** specifies the start SNP.

**\--to-snp** specifies the end SNP.

**NOTE** : All SNPs should be on the same chromosome.

To query eQTL results for all SNP on a chromosome

```
smr --beqtl-summary myeqtl --query 5.0e-8 --snp-chr 1 
```

**\--snp-chr** specifies a chromosome to select SNPs.

**NOTE** : The probes in the result could be on the other chromosomes if
there are *trans*-eQTLs.

To query SNPs based on physical positions

```
smr --beqtl-summary myeqtl --query 5.0e-8 --snp-chr 1 --from-snp-kb 100 --to-snp-kb 200 --out myquery 
```

**\--from-snp-kb** specifies the start physical position of the
region.

**\--to-snp-kb** specifies the end physical position of the region.

**NOTE** : You will need to specify a chromosome (using the '\--snp-chr'
option) when using this option.

To query based on a flanking region of a SNP

```
smr --beqtl-summary myeqtl --query 5.0e-8 --snp rs123 --snp-wind 50 --out myquery 
```

**\--snp-wind** defines a window (in Kb unit) centred on a specified SNP.

#### \# Command line options for probes

To query based on a single probe

```
smr --beqtl-summary myeqtl --query 5.0e-8 --probe cg123 --out myquery 
```

**\--probe** specifies a single probe.

To query based on a range of probes

```
smr --beqtl-summary myeqtl --query 5.0e-8 --from-probe cg123 --to-probe cg456 --out myquery 
```

**\--from-probe** specifies the start probe.

**\--to-probe** specifies the end probe.

NOTE : All probes should be on the same chromosome.

To query based on a chromosome

```
smr --beqtl-summary myeqtl --query 5.0e-8 --probe-chr 1 
```

**\--probe-chr** specifies a chromosome to select probes.

**NOTE** : The SNPs in the result could be on the other chromosomes if there
are *trans*-eQTLs.

To query based on physical positions of the probes

```
smr --beqtl-summary myeqtl --query 5.0e-8 --probe-chr 1 --from-probe-kb 1000 --to-probe-kb 2000 --out myquery 
```

**\--from-probe-kb** specifies the start physical position of the
probes.

**\--to-probe-kb** specifies the end physical position of the
probes.

**NOTE** : You will need to specify a chromosome (using the '\--probe-chr'
option) when using this option.

To query based on a flanking region of a probe

```
smr --beqtl-summary myeqtl --query 5.0e-8 --probe cg123 --probe-wind 1000 --out myquery 
```

**\--probe-wind** defines a window (in Kb unit) centred on a specified probe.

To query based on a gene

```
smr --beqtl-summary myeqtl --query 5.0e-8 --gene gene1 --out myquery 
```

**\--gene** specifies a single gene to select probes.

#### \# Command line option for cis-region

```
smr --beqtl-summary myeqtl --query 5.0e-8 --probe cg123 --cis-wind 2000 --out myquery 
```

#### \# File-list options

To query based on a list of SNPs

```
smr --beqtl-summary myeqtl --extract-snp snp.list --query 5.0e-8 --out myquery 
```

To query based on a list of probes

```
smr --beqtl-summary myeqtl --extract-probe probe.list --query 5.0e-8 --out myquery
```

To qurey based on a list of genes

```
smr --beqtl-summary myeqtl --genes gene.list --query 5.0e-8 --out myquery 
```

**\--genes** extracts a subset of probes which tag the genes in the
list.

***gene.list***
```
gene1
gene2
gene3
...
```
### Descriptive summary of the  cis-eQTL data

```
smr  --descriptive-cis --beqtl-summary myeqtl --out mydescriptive 
```
**\--descriptive-cis** outputs the descriptive summary of the cis-eQTL data (including the number of probes, the boundaries of the each cis-region, the top cis-eQTL in each cis-eQTL region and the summary statistics of the top cis-eQTL) at p-value threshold (the default value is 5.0e-8; can be changed by --peqtl-cis). The size of the cis-region can be changed by --cis-wind (the default value is 2000Kb).


```
smr  --descriptive-cis --beqtl-summary myeqtl --cis-wind 2000 --out mydescriptive 
```


### Descriptive summary of the  trans-eQTL data

```
smr --descriptive-trans --beqtl-summary myeqtl --out mydescriptive 
```
**\--descriptive-trans**  outputs the descriptive summary of the trans-eQTL data at p-value threshold (the default value is 5.0e-8; can be changed by --peqtl-trans). The size of a trans-region can be changed by --trans-wind (the default value is 1000Kb). Note that the trans-eQTLs are defined those that are more than 5Mb away from the probe (this distance threshold can be changed by --cis-wind).


```
smr --descriptive-trans --beqtl-summary myeqtl --trans-wind 1000 --out mydescriptive 
```

