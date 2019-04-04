## Cleans and Analyzes 16S data from Caribean Coral

### Required Packages
* [dada2](https://bioconductor.org/packages/release/bioc/html/dada2.html)
* [ShortRead](http://bioconductor.org/packages/release/bioc/html/ShortRead.html) 
* [ggplot2](https://cran.r-project.org/package=ggplot2)
* [phyloseq](https://joey711.github.io/phyloseq/) 
* [vegan](https://cran.r-project.org/package=vegan)
* [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html)
* [randomcoloR](https://cran.r-project.org/package=randomcoloR)
* [knitr](https://github.com/yihui/knitr#readme)
* [ALDEx2](https://bioconductor.org/packages/release/bioc/html/ALDEx2.html)
* [CoDaSeq](https://github.com/ggloor/CoDaSeq)
* [zCompositions](https://cran.r-project.org/package=zCompositions)
* [igraph](https://igraph.org/r/)
* [car](https://cran.r-project.org/package=car)
* [propr](https://cran.r-project.org/package=propr)
* [dendextend](https://cran.r-project.org/package=dendextend)
* [tibble](https://tibble.tidyverse.org/)
* [ggfortify](https://cran.r-project.org/package=ggfortify)
* [ggbiplot](https://github.com/vqv/ggbiplot)

Majority of these packages are part of CRAN

### Description

Sequence the data (physical)
Pull the data from the Hypergator
DADA2 was used for filtering, error estimation, merging reads, depreplicaiton, removing chimeras, and selection of the amplicon squence variants (ASVs). the SILVA reference data set was used with DADA to perform taxonmy classification, exluding any data that could not be assingned as bacteria or archae. The taxonomy tables are imported into Phyloseq for community analysis based on Count, Taxonomy, and ASVs. The zCompisition package performed the count xzero multiplicative method  to transform zero counts so that  Aitchison distance metric could be calculated with CoDaSeq. Principle component anlysis of the distance was performed with the prcomp package and plotted with ggplot3. Vegan was used to compare the similarities between the distances and communities. Differential ASV abundance within genotype communities was calculated using DESeq2. Red and Green genotype communities were found to be relatively similar and therefor were compared against the yellow communities.

## ERRORS
alternative if getting error cannot compute log geometric means
cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- DESeq(dds,test="Wald", fitType="parametric")
