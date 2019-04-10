library(phyloseq)
library(ggplot2)
library(cowplot)
createPsObject <- function (otuTable,taxaTable, metadata){
otu <- read.table(otuTable,sep="\t",header=TRUE, row.names=1)
taxon <- read.table(taxaTable,sep="\t",header=TRUE,row.names=1)
samples<-read.table(metadata,sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
class(OTU)
taxa_names(OTU)####IS THIS IMPERATIVE? ONLY PRINTS HALF OF SMAPLES BEFORE MAXING OUT
try(class(taxmat))
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
taxa_names(TAX) ####IS THIS IMPERATIVE? ONLY PRINTS HALF OF SMAPLES BEFORE MAXING OUT
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps <- prune_samples(sample_names(ps), ps)
print(ps)
return (ps)
}