source('16sDataCleaner-Nicole/createPS.R')
library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(vegan)

####put parsed, adaptors & primers removed, unjoined (R1 and R2 separate) fastq files
# into directory for DADA2 & make sure the full path is updated in the next line:

path <- "16sDataCleaner-Nicole/Data"
list.files(path)

##dada2 v1.6.0

# Samplename is everything before the first underscore
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

##Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[1:6])
#fixediterror#
^#plotQualityProfile(fnRs[1:6]) #Error from cutadapt (sequences with length 0) fix @ http://cutadapt.readthedocs.io/en/stable/guide.html#filtering-reads with parameter --minimum-length N

#Perform filtering and trimming
#Assign the filenames for the filtered fastq.gz files.
#Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter the forward and reverse reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#Learn the Error Rates, it TAKES TIME! do first forward and then reverse
# Forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
# Reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)

#visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)

#Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#Inspecting the dada-class object returned by dada:
dadaFs[[1]]

#Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers) ## The sequences being tabled vary in length.
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
#As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "dada_read_stats.txt",sep="\t",col.names=NA)

#####SAVE THIS FILE SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS, adjust name
saveRDS(seqtab.nochim, file="16sDataCleaner-Nicole/seqtab.nochim.rds")

# RELOAD THE SAVED INFO FROM HERE (if you have closed the project):
seqtab.nochim <- readRDS("16sDataCleaner-Nicole/seqtab.nochim.rds")

#Assign taxonomy
#Make sure the appropriate database is available in the DADA2 directory
taxa <- assignTaxonomy(seqtab.nochim, "16sDataCleaner-Nicole/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

#### FIX the NAs in the taxa table (I'm writing it out here, then reading back in for phyloseq)
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)
taxon$Phylum[is.na(taxon$Phylum)] <- taxon$Kingdom[is.na(taxon$Phylum)]
taxon$Class[is.na(taxon$Class)] <- taxon$Phylum[is.na(taxon$Class)]
taxon$Order[is.na(taxon$Order)] <- taxon$Class[is.na(taxon$Order)]
taxon$Family[is.na(taxon$Family)] <- taxon$Order[is.na(taxon$Family)]
taxon$Genus[is.na(taxon$Genus)] <- taxon$Family[is.na(taxon$Genus)]
write.table(taxon,"Acropora_silva_taxa_table.txt",sep="\t",col.names=NA)
write.table(seqtab.nochim, "Acropora_silva_otu_table.txt",sep="\t",col.names=NA)

#Removing Mitochondrial and Chloroplast Reads:
#####read in otu and taxonomy tables from dada2, sample data.
ps = createPsObject("Acropora_silva_otu_table.txt",
						   "Acropora_silva_taxa_table.txt",
						   "Acropora_metadata.txt")
#remove chloroplasts and mitochondria and Eukaryota
get_taxa_unique(ps, "Family") #559
get_taxa_unique(ps, "Order") #331
get_taxa_unique(ps, "Kingdom") #4
ps2 <- subset_taxa(ps, Family !="Mitochondria")
ps2 <- subset_taxa(ps2, Order !="Chloroplast")
ps2 <- subset_taxa(ps2, Kingdom !="Eukaryota")
ps2 <- subset_taxa(ps2, Kingdom !="NA")
get_taxa_unique(ps2, "Family") #555
get_taxa_unique(ps2, "Order") #327
get_taxa_unique(ps2, "Kingdom") #2
ps2
# filtered taxa with phyloseq, now export otu and taxa tables from phyloseq
otu = as(otu_table(ps2), "matrix")
taxon = as(tax_table(ps2), "matrix")
metadata = as(sample_data(ps2), "matrix")
write.table(otu,"Acropora_silva_nochloronomito_otu_table.txt",sep="\t",col.names=NA)
write.table(taxon,"Acropora_silva_nochloronomito_taxa_table.txt",sep="\t",col.names=NA)

# look at data and chose filtering method for very low abundance ASVs
ntaxa(ps) #11078
ps10<-filter_taxa(ps, function(x) mean(x) >10, TRUE)
ntaxa(ps10) #265
ps5<-filter_taxa(ps, function(x) mean(x) >5, TRUE)
ntaxa(ps5) #453
get_taxa_unique(ps, "Genus") # 1096
get_taxa_unique(ps5, "Genus") #139
get_taxa_unique(ps10, "Genus") #93

# filtered ASVs with very low abundance with phyloseq, now export otu and taxa tables from phyloseq for codaseq
otu = as(otu_table(ps5), "matrix")
taxon = as(tax_table(ps5), "matrix")
metadata = as(sample_data(ps5), "matrix")
write.table(otu,"Acropora_ps5_silva_nochloronomito_otu_table.txt",sep="\t",col.names=NA)
write.table(taxon,"Acropora_ps5_silva_nochloronomito_taxa_table.txt",sep="\t",col.names=NA)
write.table(metadata,"Acropora_ps5_silva_metadata.txt",sep="\t",col.names=NA)

#Creating Community Composition Bar Charts
##### bar charts using ps5 with 81 samples (373 taxa)
ps = createPsObject("Acropora_ps5_silva_nochloronomito_otu_table.txt",
						   "Acropora_ps5_silva_nochloronomito_taxa_table.txt",
						   "Acropora_ps5_silva_metadata.txt")
ps_ra<-transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
#figure out how many colors you need
get_taxa_unique(ps_ra, "Class") #27
get_taxa_unique(ps_ra, "Order") #65
get_taxa_unique(ps_ra, "Genus") #139

#you can make n any number of colors you want; with as much difference between the colors as possible (distinct colors)
n <- 64
palette <- distinctColorPalette(n)
#you can rerun the previous line to get a new selection of colors

ps_ra_red = subset_samples(ps_ra, Genotype == "R")
ps_ra_green = subset_samples(ps_ra, Genotype == "G")
ps_ra_yellow = subset_samples(ps_ra, Genotype == "Y")

p1=plot_bar(ps_ra_red, fill="Order")+
  geom_bar(aes(fill=Order), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+
  ggtitle("Genotype R")+
  facet_grid(.~Colony,scales="free",space="free")+
  theme(plot.title = element_text(face="italic"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

p2=plot_bar(ps_ra_green, fill="Order")+
  geom_bar(aes(fill=Order), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+
  ggtitle("Genotype G")+
  facet_grid(.~Colony,scales="free",space="free")+
  theme(plot.title = element_text(face="italic"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

p3=plot_bar(ps_ra_yellow, fill="Order")+
  geom_bar(aes(fill=Order), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+
  ggtitle("Genotype Y")+
  facet_grid(.~Colony,scales="free",space="free")+
  theme(plot.title = element_text(face="italic"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

####adjust width and height until it looks right for double columns
pdf("Acropora_BarCharts_Order_forthelegend.pdf",width=24, height=10)
plot_grid(p1,p2,p3,labels=c("A","B","C"), ncol=2, nrow=2)
dev.off()

#Bar Plots Individual:
# bar charts using ps5 with 62 samples (683 taxa)
ps = createPsObject("Acropora_ps5_silva_nochloronomito_otu_table.txt",
						    "Acropora_ps5_silva_nochloronomito_taxa_table.txt",
						    "Acropora_ps5_silva_metadata.txt")
ps_ra<-transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps_ra_G = subset_samples(ps_ra, Genotype == "G")
ps_ra_R = subset_samples(ps_ra, Genotype == "R")
ps_ra_Y = subset_samples(ps_ra, Genotype == "Y")
#figure out how many colors you need
get_taxa_unique(ps_ra, "Order") #65
get_taxa_unique(ps_ra, "Class") #27
#you can make n any number of colors you want; with as much difference between the colors as possible (distinct colors)
n <- 99
palette <- distinctColorPalette(n)
#you can rerun the previous line to get a new selection of colors
p1=plot_bar(ps_ra_G, fill="Order")
p1+geom_bar(aes(fill=Order), stat="identity",position="stack")+
ftheme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+theme(legend.position = "bottom")+
  facet_grid(.~Colony,scales="free",space="free")

p2=plot_bar(ps_ra_R, fill="Order")
p2+geom_bar(aes(fill=Order), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+theme(legend.position = "bottom")+
  facet_grid(.~Colony,scales="free",space="free")

p3=plot_bar(ps_ra_Y, fill="Order")
p3+geom_bar(aes(fill=Order), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+theme(legend.position = "bottom")+
  facet_grid(.~Colony,scales="free",space="free")

#PCA Analysis and ANOSIM
#Microbiome Datasets are Compositional: and this is not optional
#Gloor et al 2017 & other ggloor papers and codes on github (CoDA microbiome tutorial)
library(knitr)
library(ALDEx2)
library(CoDaSeq)
library(zCompositions)
library(igraph)
library(car)
library(grDevices)
library(propr)
library(vegan)
library(dendextend)
library(tibble)
library(ggfortify)
library(ggbiplot)
#sessionInfo()

#### READ IN ***FILTERED***, cleaned DADA2 OTU TABLE (chloroplasts and mitochondria removed) and taxonomy table
otu <- read.table("Acropora_ps5_silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("Acropora_ps5_silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("Acropora_ps5_silva_metadata.txt",sep="\t",header=T,row.names=1)
genus<-as.character(taxon$Genus)

# First, replace 0 values with an estimate (because normalization is taking log, can't have 0)
# Also transposing here, need samples as rows
# First, replace 0 values with an estimate (because normalization is taking log, can't have 0)
# Also transposing here, need samples as rows
d.czm <- cmultRepl(t(otu), method="CZM", label=0)
# Perform the center-log-ratio (CLR) transformation 
d.clr <- codaSeq.clr(d.czm)
# transpose matrix of CLR transformed data for ordination and dendrogram
E.clr <- t(d.clr)
# plot compositional PCA biplot (perform a singular value decomposition)
d.pcx <- prcomp(E.clr)
# calculate percent variance explained for the axis labels
pc1 <- round(d.pcx$sdev[1]^2/sum(d.pcx$sdev^2),2)
pc2 <- round(d.pcx$sdev[2]^2/sum(d.pcx$sdev^2),2)
xlab <- paste("PC1: ", pc1, sep="")
ylab <- paste("PC2: ", pc2, sep="")
biplot(d.pcx, cex=c(0.6,0.4), var.axes=F,scale=1, xlab=xlab, ylab=ylab, ylabs=genus)
summary(d.pcx)
str(d.pcx)
screeplot(d.pcx)

##### replot PCA with ggplot2 (showing samples only)
df_out <- as.data.frame(d.pcx$x)
theme_set(theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))
pdf("Acropora_PCA_Genotype_Branch.pdf",width=8.5)
p<-ggplot(df_out,aes(x=PC1,y=PC2,color=samples$Genotype,shape=samples$Branch))
p<-p+geom_point(size=3)+theme(axis.title = element_text(size=14))+theme(axis.text=element_text(size=12))+
  theme(legend.title = element_text(size=14))+theme(legend.text = element_text(size=12))+
  scale_color_manual(values=c("#009E73","#990000","#e79f00","#56B4E9"))
p + labs(x=xlab, y=ylab, color="Genotype", shape="Branch") + coord_fixed()
dev.off()

# set metadata as factors for anosim
conds<-as.character(samples$Genotype)
site<-as.character(samples$Branch)
coral<-as.character(samples$Colony)

# anosim between groups using Aitchison distance
dist.clr <- dist(E.clr)
ano <- anosim(dist.clr, conds, permutations=999)
plot(ano)

ano <- anosim(dist.clr, site, permutations=999)
png("Acropora_ANOSIM_Branch.pdf")
plot(ano)
dev.off()

ano <- anosim(dist.clr,conds, permutations=999)
pdf("Acropora_ANOSIM_Genotype.pdf",width=8.5)
plot(ano)
dev.off()

ano <- anosim(dist.clr, coral, permutations=999)
pdf("Acropora_ANOSIM_Colony.pdf",width=8.5)
plot(ano)
dev.off()


####DeSeq#####
###################### Differential abundance
### first get only samples below, adjust otu, taxa, and metadata tables.

############################################# Green V Red ############################################
ps = createPsObject(,"Acropora_ps5_silva_nochloronomito_taxa_table.txt",
							"Acropora_ps5_silva_nochloronomito_otu_table.txt",
							"Acropora_metadata_GR.txt")

#Subset data into TWO conditions (needed for Deseq)
ps_GR= subset_samples(ps, Genotype != "Y")
ps_GR # should be 54 samples
otu = as(otu_table(ps_GR), "matrix")
taxon = as(tax_table(ps_GR), "matrix")
write.table(otu,"Acropora_silva_nochloronomito_otu_table_Green_Red.txt",sep="\t",col.names=NA)
write.table(taxon,"Acropora_silva_nochloronomito_taxa_table_Green_Red.txt",sep="\t",col.names=NA)
#use already fixed metadata (new column)
#clear data and read back in
ps = createPsObject("Acropora_silva_nochloronomito_otu_table_Green_Red.txt",
						   "Acropora_silva_nochloronomito_taxa_table_Green_Red.txt",
						   "Acropora_metadata_Green_Red.txt")
#Define the order of the conditions for testing
#In this order, the positive fold change values are what increased in the cultures compared to roots
sample_data(ps)$Genotype<-factor(sample_data(ps)$Genotype,levels=c("G","R"))
head(sample_data(ps)$Genotype, 10)

#DESEQ2 analysis; no filtering of low abundance OTUs
library("DESeq2")
packageVersion("DESeq2")
dds = phyloseq_to_deseq2(ps, ~ Genotype)
dds
#filter rows with very few counts
dds <- dds[ rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds,test="Wald", fitType="parametric")

#######If error check read me

res = results(dds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

#save table of results
sig = as(sigtab, "matrix")
write.table(sig,"DESeq2_results_Green_Red.txt",sep="\t",col.names=NA)
##rplot from DESeq2 results
sigtab<-read.table("DESeq2_results_Green_Red.txt",sep="\t",header=TRUE,row.names=1)
#ggplot2 summary of the results
library(ggplot2)
theme_set(theme_bw())
my20colors<-c("#c26162","#d689c2","#d64142","#6db643","#9c58cb","#cdab39","#626fdf","#5cba78","#bd3fa0","#497535","#de76dd","#949345","#6e58a2","#cd6529","#46aed7","#d94681","#49b39c","#9e4e77","#bd814c","#798fd5")
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
pdf(file="DESeq_Acropora_Green_Red_Phylum.pdf",width=8.5)
ggplot(sigtab, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + ggtitle("Green Red") + geom_point(size=4) + coord_flip() +scale_color_manual(values=my20colors)
dev.off()
# Genus order
x2 = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x2) max(x2))
x2 = sort(x2, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x2))
pdf(file="DeSeq_Acropora_Green_Red_Genus.pdf",width=8.5)
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + ggtitle("Green Red") + geom_point(size=4) + coord_flip() +scale_color_manual(values=my20colors)
dev.off()

############################################# Red Yellow  ############################################
ps = createPsObject("Acropora_ps5_silva_nochloronomito_otu_table.txt",
						   "Acropora_ps5_silva_nochloronomito_taxa_table.txt",
						   "Acropora_metadata_Red_Yellow.txt")

#Subset data into TWO conditions (needed for Deseq)
ps_RY= subset_samples(ps, Genotype != "G")
ps_RY # should be 54 samples
otu = as(otu_table(ps_RY), "matrix")
taxon = as(tax_table(ps_RY), "matrix")
write.table(otu,"Acropora_silva_nochloronomito_otu_table_Red_Yellow.txt",sep="\t",col.names=NA)
write.table(taxon,"Acropora_silva_nochloronomito_taxa_table_Red_Yellow.txt",sep="\t",col.names=NA)
#use already fixed metadata (new column)
#clear data and read back in
ps = createPsObject("Acropora_silva_nochloronomito_otu_table_Red_Yellow.txt",
						   "Acropora_silva_nochloronomito_taxa_table_Red_Yellow.txt",
						   "Acropora_metadata_Red_Yellow.txt")

#Define the order of the conditions for testing
#In this order, the positive fold change values are what increased in the cultures compared to roots
sample_data(ps)$Genotype<-factor(sample_data(ps)$Genotype,levels=c("R","Y"))
head(sample_data(ps)$Genotype, 10)

#DESEQ2 analysis; no filtering of low abundance OTUs
library("DESeq2")
packageVersion("DESeq2")
dds = phyloseq_to_deseq2(ps, ~ Genotype)
dds
#filter rows with very few counts
dds <- dds[ rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds,test="Wald", fitType="parametric")

#######aif error check README
res = results(dds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

#save table of results
sig = as(sigtab, "matrix")
write.table(sig,"DESeq2_results_Red_Yellow.txt",sep="\t",col.names=NA)
##rplot from DESeq2 results
sigtab<-read.table("DESeq2_results_Red_Yellow.txt",sep="\t",header=TRUE,row.names=1)
#ggplot2 summary of the results
library(ggplot2)
theme_set(theme_bw())
my20colors<-c("#c26162","#d689c2","#d64142","#6db643","#9c58cb","#cdab39","#626fdf","#5cba78","#bd3fa0","#497535","#de76dd","#949345","#6e58a2","#cd6529","#46aed7","#d94681","#49b39c","#9e4e77","#bd814c","#798fd5")
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
pdf(file="DESeq_Acropora_Red_Yellow_Phylum.pdf",width=8.5)
ggplot(sigtab, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + ggtitle("Red Yellow") + geom_point(size=4) + coord_flip() +scale_color_manual(values=my20colors)
dev.off()
# Genus order
x2 = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x2) max(x2))
x2 = sort(x2, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x2))
pdf(file="DeSeq_Acropora_Red_Yellow_Genus.pdf",width=8.5)
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + ggtitle("Red Yellow") + geom_point(size=4) + coord_flip() +scale_color_manual(values=my20colors)
dev.off()

############################################# Green/Red Yellow ############################################
ps = createPsObject("Acropora_ps5_silva_nochloronomito_otu_table.txt",
						   "Acropora_ps5_silva_nochloronomito_taxa_table.txt",
						   "Acropora_metadata_GR.txt")
#Subset data into TWO conditions (needed for Deseq)
ps # should be 81 samples
otu = as(otu_table(ps), "matrix")
taxon = as(tax_table(ps), "matrix")
write.table(otu,"Acropora_silva_nochloronomito_otu_table_GreenRed_Yellow.txt",sep="\t",col.names=NA)
write.table(taxon,"Acropora_silva_nochloronomito_taxa_table_GreenRed_Yellow.txt",sep="\t",col.names=NA)
#use already fixed metadata (new column)
#clear data and read back in
ps = createPsObject("Acropora_silva_nochloronomito_otu_table_GreenRed_Yellow.txt",
						   "Acropora_silva_nochloronomito_taxa_table_GreenRed_Yellow.txt",
						   "Acropora_metadata_GR.txt")
#Define the order of the conditions for testing
#In this order, the positive fold change values are what increased in the cultures compared to roots
sample_data(ps)$Geno<-factor(sample_data(ps)$Geno,levels=c("GR","Y"))
head(sample_data(ps)$Geno, 10)

#DESEQ2 analysis; no filtering of low abundance OTUs
library("DESeq2")
packageVersion("DESeq2")
dds = phyloseq_to_deseq2(ps, ~ Geno)
dds
#filter rows with very few counts
dds <- dds[ rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds,test="Wald", fitType="parametric")

#######If there is an error check README

res = results(dds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

#save table of results
sig = as(sigtab, "matrix")
write.table(sig,"DESeq2_results_GreenRed_Yellow.txt",sep="\t",col.names=NA)
##rplot from DESeq2 results
sigtab<-read.table("DESeq2_results_GreenRed_Yellow.txt",sep="\t",header=TRUE,row.names=1)
#ggplot2 summary of the results
library(ggplot2)
theme_set(theme_bw())
my20colors<-c("#c26162","#d689c2","#d64142","#6db643","#9c58cb","#cdab39","#626fdf","#5cba78","#bd3fa0","#497535","#de76dd","#949345","#6e58a2","#cd6529","#46aed7","#d94681","#49b39c","#9e4e77","#bd814c","#798fd5")
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
pdf(file="DESeq_Acropora_GreenRed_Yellow_Phylum.pdf",width=8.5)
ggplot(sigtab, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + ggtitle("Green&Red Yellow") + geom_point(size=4) + coord_flip() +scale_color_manual(values=my20colors)
dev.off()
# Genus order
x2 = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x2) max(x2))
x2 = sort(x2, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x2))
pdf(file="DeSeq_Acropora_GreenRed_Yellow_Genus.pdf",width=8.5)
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + xlab("Genera") + geom_point(size=4) + coord_flip() +scale_color_manual(values= c("#33CCCC","mediumpurple1","deeppink4")) +geom_abline(m = 0, color = "red", slope = 0)
dev.off()

plotMA(res, ylim=c(-5,5))
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

#Distance to Centroid
#Microbiome Datasets are Compositional: and this is not optional
#Gloor et al 2017 & other ggloor papers and codes on github (CoDA microbiome tutorial)
library(knitr)
library(ALDEx2)
library(CoDaSeq)
library(zCompositions)
library(igraph)
library(car)
library(grDevices)
library(propr)
library(vegan)
library(ggplot2)
library(Rmisc)
library(reshape)
library(tibble)
library(ggplot2)
#### READ IN OTU data that has been filtered for very low abundance sequences
otu <- read.table("Acropora_ps5_silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("Acropora_ps5_silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("Acropora_ps5_silva_metadata.txt",sep="\t",header=T,row.names=1)
genus<-as.character(taxon$Genus)

# First, replace 0 values with an estimate (because normalization is taking log, can't have 0)
# Also transposing here, need samples as rows
d.czm <- cmultRepl(t(otu), method="CZM", label=0)
# Perform the center-log-ratio (CLR) transformation 
d.clr <- codaSeq.clr(d.czm)
# transpose matrix of CLR transformed data for ordination and dendrogram
E.clr <- t(d.clr)
# calculate Aitchison distance
dist.clr <- dist(E.clr)

# set metadata as factors for vegan
conds<-as.character(samples$Genotype)
site<-as.character(samples$Branch)
coral<-as.character(samples$Colony)

#calculate multivariate dispersions based on condition
mod <-betadisper(dist.clr, conds)
#one way anova
anova(mod)
#boxplots
plot(mod)
boxplot(mod)

## Compute mean distance to centroid per group
#this just prints values on the console
tapply(mod$distances, conds, mean)
## Same, but variance instead
tapply(mod$distances, conds, var)

#Get the distances to centroid from the model
mod$distances
dis <- mod$distances
#melt
dis.melt <- melt(dis)
#move rownames to columns so we can merge the dispersion values and metadata
dis.melt$Sample <- rownames(dis.melt)
samples$Sample <- rownames(samples)
#merge metadata and dispersion 
dis.treat <- merge(samples, dis.melt)
#rename column
colnames(dis.treat)[5] <- "distance"

#run linear model to test significance
distlm <-lm(distance~Genotype*Colony, data=dis.treat)
summary(distlm)
anova(distlm)

#### plot average dispersion by group, with all points shown
p1<-ggplot(dis.treat, aes(x=Colony,y=distance))+
  geom_boxplot()+
  theme_bw()+
  geom_jitter(position=position_jitter(width=.1, height=0),aes(color=Genotype),size=3)+
  scale_color_manual(values=c("#009E73","#990000","#e79f00","#56B4E9"))+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(face="italic"))+
  theme(text=element_text(size=16))+
  ylab("Distance to Centroid")
pdf("Figure4_Acropora_DistanceToCentroid.pdf",width=11,height=11)
plot(p1)
dev.off()

#Supplementary 
f <- function(x) {
  ans <- boxplot.stats(x)
  data.frame(ymin = ans$conf[1], ymax = ans$conf[2], y = ans$stats[3])
}

p2<-ggplot(dis.treat, aes(x=Colony,y=distance))+
  geom_boxplot()+
  theme_bw()+
  stat_summary(fun.data = f, geom = "crossbar", 
               colour = NA, fill = "skyblue", width = 0.8, alpha = 0.5)
geom_jitter(position=position_jitter(width=.1, height=0),aes(color=Genotype),size=3)+
  scale_color_manual(values=c("#009E73","#990000","#e79f00","#56B4E9"))+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(face="italic"))+
  theme(text=element_text(size=16))+
  ylab("Distance to Centroid")
pdf("Acropora_Notched_Colony_DistanceToCentroid.pdf",width=11,height=11)
plot(p2)
dev.off()

p3<-ggplot(dis.treat, aes(x=Genotype,y=distance))+
  geom_boxplot()+
  theme_bw()+
  stat_summary(fun.data = f, geom = "crossbar", 
               colour = NA, fill = "skyblue", width = 0.8, alpha = 0.5)
geom_jitter(position=position_jitter(width=.1, height=0),aes(color=Genotype),size=3)+
  scale_color_manual(values=c("#009E73","#990000","#e79f00","#56B4E9"))+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(face="italic"))+
  theme(text=element_text(size=16))+
  ylab("Distance to Centroid")
pdf("Acropora_Notched_Genotype_DistanceToCentroid.pdf",width=11,height=11)
plot(p3)
dev.off()
#END#