## Cleans and Analyzes 16S data from Caribean Coral

### Required Packages
* [dada2](https://bioconductor.org/packages/release/bioc/html/dada2.html)
* [phyloseq](https://joey711.github.io/phyloseq/) 
* [vegan](https://cran.r-project.org/package=vegan)
* [CoDaSeq](https://github.com/ggloor/CoDaSeq)

### Description
Sequence the data (physical)
Pull the data from the Hypergator
DADA2 was used for filtering, error estimation, merging reads, depreplicaiton, removing chimeras, and selection of the amplicon squence variants (ASVs). the SILVA reference data set was used with DADA to perform taxonmy classification, exluding any data that could not be assingned as bacteria or archae. The taxonomy tables are imported into Phyloseq for community analysis based on Count, Taxonomy, and ASVs. The zCompisition package performed the count xzero multiplicative method  to transform zero counts so that  Aitchison distance metric could be calculated with CoDaSeq. Principle component anlysis of the distance was performed with the prcomp package and plotted with ggplot3. Vegan was used to compare the similarities between the distances and communities. Differential ASV abundance within genotype communities was calculated using DESeq2. Red and Green genotype communities were found to be relatively similar and therefor were compared against the yellow communities.

### WINDOWS USERS 
Code was developed using IOS. Make sure to change line 37 and set multithreading to false. You can also search for Windows users in your editor/IDE for easier editing. Results may vary when run on Windows. 

### ERRORS
Alternative solution if an error appears when computing the log geometric means
``` R
cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- DESeq(dds,test="Wald", fitType="parametric")
```
