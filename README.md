# ifpan-kukla-methylation

## Cytosine level analysis

### Methods:
Bismark genome-wide methylation reports were used to compute methylation ratio for each cytosine (n = 53698126) for each sample. Cytosines were annotated with closest gene within 2000 base pairs, using protein coding genes from Ensembl (version 104). Cytosines were filtered for coverage (minimum 10 in each sample, n = 1259063), for being annotated with gene symbol (n = 810573) and for variance (greter than zero, n = 795936). For statistical analysis we used t test within each brain structure (FCx, Hp) separately. For multiple testing correction False Discovery Rate was used.

### Results:
None of the analysed cytosines passed false discovery rate threshold FDR < 0.1.

## CpG islands level analysis

### Methods:
CpG islands were downloaded from UCSC Genome Browser using Table Browser. CpG islands were annotated with overlapping protein coding genes (extended by 2000 bases from each side) downloaded from Ensembl (version 104). In situation of multiple genes overlapping CpG island the gene which center was closest to CpG island was selected. CpG islands were filtered for being annotated with genes and having at least 10 cytosines with computed methylation ratio (n = 13378). For statistical analysis two way anova with additive model was used with treatment and cytosine site as factors. Each brain structure was analyzed separately.

### Results
Out of 13378 analyzed CpG islands 130 have altered methylation (FDR < 0.1) in FCx and 134 in Hp.

## Gene level analysis
### Methods
Protein coding genes were downloaded from Ensembl (version 104). Genomic ranges for genes were extended by 2000 bases (range for core promoter). Genomic range were filtered for having at least 10 cytosines with computed methylation ration (n = 19939). For statistical analysis two way anova with additive model was used with treatment and cytosine site as factors. Each brain structure was analyzed separately.

### results
Out of 19939 analyzed genes 362 have altered methylation (FDR < 0.1) in FCx and 779 in Hp.