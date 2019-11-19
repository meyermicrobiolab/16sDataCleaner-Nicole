library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(vegan)
library(knitr)
library(ALDEx2)
library(CoDaSeq)
library(zCompositions)
library(igraph)
library(car)
library(grDevices)
library(propr)
library(cowplot)
library(randomcoloR)
library(DESeq2)
library(dplyr)
library(reshape2)
library(tibble)
library(exactRankTests)
library(nlme)
library(ggplot2)
library(data.table)
library(scales)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

###### Quality-filter reads and create Amplicon Sequence Variant tables

###adjust path name as appropriate
path <- "~/Documents/NurseryAcropora/cutadapt/"
list.files(path)

# Samplename is everything before the first underscore
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_*.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Perform filtering and trimming
# Assign the filenames for the filtered fastq.gz files.
# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter the forward and reverse reads
# WINDOWS USERS: set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn the Error Rates, it TAKES TIME! do first forward and then reverse
# Forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
# Reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)

# visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)

# Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Inspecting the dada-class object returned by dada:
dadaFs[[1]]

# Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "dada_read_stats.txt",sep="\t",col.names=NA)

# SAVE THIS FILE SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS, adjust name
saveRDS(seqtab.nochim, file="~/Documents/NurseryAcropora/seqtab.nochim.rds")
# RELOAD THE SAVED INFO FROM HERE (if you have closed the project):
# seqtab.nochim <- readRDS("~/Documents/NurseryAcropora/seqtab.nochim.rds")

######################ASSIGNING THE TAXONOMY###############################################################################
# Make sure the appropriate database is available in the DADA2 directory
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/NurseryAcropora/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

# FIX the NAs in the taxa table
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)
taxon$Phylum[is.na(taxon$Phylum)] <- taxon$Kingdom[is.na(taxon$Phylum)]
taxon$Class[is.na(taxon$Class)] <- taxon$Phylum[is.na(taxon$Class)]
taxon$Order[is.na(taxon$Order)] <- taxon$Class[is.na(taxon$Order)]
taxon$Family[is.na(taxon$Family)] <- taxon$Order[is.na(taxon$Family)]
taxon$Genus[is.na(taxon$Genus)] <- taxon$Family[is.na(taxon$Genus)]
write.table(taxon,"silva_taxa_table.txt",sep="\t",col.names=NA)
write.table(seqtab.nochim, "silva_otu_table.txt",sep="\t",col.names=NA)

#####################REMOVING MITOCHONDIRAL AND CHLOROPLAST READS############################################################
# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata.
otu <- read.table("silva_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps
# 11268 taxa and 81 samples
# remove chloroplasts and mitochondria and Eukaryota
get_taxa_unique(ps, "Family") #561
get_taxa_unique(ps, "Order") #331
get_taxa_unique(ps, "Kingdom") #4
ps <- subset_taxa(ps, Family !="Mitochondria")
ps <- subset_taxa(ps, Order !="Chloroplast")
ps <- subset_taxa(ps, Kingdom !="Eukaryota")
ps <- subset_taxa(ps, Kingdom !="NA")
get_taxa_unique(ps, "Family") #557
get_taxa_unique(ps, "Order") #327
get_taxa_unique(ps, "Kingdom") #2
ps # 10153 taxa and 81 samples

# filtered taxa with phyloseq, now export cleaned otu and taxa tables from phyloseq
otu = as(otu_table(ps), "matrix")
taxon = as(tax_table(ps), "matrix")
metadata = as(sample_data(ps), "matrix")
write.table(otu,"silva_nochloronomito_otu_table.txt",sep="\t",col.names=NA)
write.table(taxon,"silva_nochloronomito_taxa_table.txt",sep="\t",col.names=NA)

# look at data and chose filtering method for very low abundance ASVs
ntaxa(ps) #10153
ps5<-filter_taxa(ps, function(x) mean(x) >5, TRUE)
ntaxa(ps5) #372
ps10<-filter_taxa(ps, function(x) mean(x) >10, TRUE)
ntaxa(ps10) #234
get_taxa_unique(ps, "Genus") #1098
get_taxa_unique(ps5, "Genus") #136
get_taxa_unique(ps10, "Genus") #92

# filtered ASVs with very low abundance with phyloseq, now export otu and taxa tables from phyloseq for codaseq
otu = as(otu_table(ps5), "matrix")
taxon = as(tax_table(ps5), "matrix")
write.table(otu,"ps5_silva_nochloronomito_otu_table.txt",sep="\t",col.names=NA)
write.table(taxon,"ps5_silva_nochloronomito_taxa_table.txt",sep="\t",col.names=NA)

######### Perform center-log-ratio transformation on ASVs and calculate Aitchison Distance and principal components
# READ IN OTU data that has been filtered for very low abundance sequences; do not clear data here. Keep phyloseq object ps5 for anosim/permanova
otu <- read.table("ps5_silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("ps5_silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata.txt",sep="\t",header=T,row.names=1)

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
biplot(d.pcx, cex=c(0.6,0.4), var.axes=F,scale=1, xlab=xlab, ylab=ylab)
summary(d.pcx)
str(d.pcx)
screeplot(d.pcx)

######### FIGURE 2 PCA ######################################################################################################
# replot PCA with ggplot2 (showing samples only)
df_out <- as.data.frame(d.pcx$x)
theme_set(theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))
cols<-c("G"="#009E73","R"="#D55E00","Y"="#F0E442")
pdf("PCA.pdf",width=8.5)
p<-ggplot(df_out,aes(x=PC1,y=PC2,fill=samples$Genotype,shape=samples$Branch))
p<-p+geom_point(size=3)+theme(axis.title = element_text(size=14))+theme(axis.text=element_text(size=12))+
  theme(legend.title = element_text(size=14))+theme(legend.text = element_text(size=12))+
  scale_fill_manual(values=cols)+
  scale_shape_manual(values=c(21,22,24))+
  guides(fill = guide_legend(override.aes=list(shape=21)))
p + labs(x=xlab, y=ylab, fill="Genotype", shape="Branch") + coord_fixed()
dev.off()

####### Use phyloseq/vegan to perform ANOSIM/PERMANOVA
# set metadata as factors for anosim
genotype<-as.character(samples$Genotype)
branch<-as.character(samples$Branch)
colony<-as.character(samples$Colony)

# anosim between groups using Aitchison distance
dist.clr <- dist(E.clr)

ano <- anosim(dist.clr, genotype, permutations=999)
pdf("Acropora_ANOSIM_genotype.pdf")
plot(ano)
dev.off()

ano <- anosim(dist.clr,branch, permutations=999)
pdf("Acropora_ANOSIM_branch.pdf",width=8.5)
plot(ano)
dev.off()

ano <- anosim(dist.clr, colony, permutations=999)
pdf("Acropora_ANOSIM_colony.pdf",width=8.5)
plot(ano)
dev.off()

# permanova between groups using Aitchison distance
perm<-adonis(dist.clr~genotype*colony,as(sample_data(ps5),"data.frame"))
print(perm)

############ FIGURE 3 Stacked bar charts of bacterial orders ################################################################
ps_ra<-transform_sample_counts(ps5, function(OTU) OTU/sum(OTU))
#figure out how many colors you need
get_taxa_unique(ps_ra, "Class") #26
get_taxa_unique(ps_ra, "Order") #62
get_taxa_unique(ps_ra, "Genus") #136

#you can make n any number of colors you want; with as much difference between the colors as possible (distinct colors)
n <- 62
palette <- distinctColorPalette(n)
#you can rerun the previous line to get a new selection of colors
# keep list of colors used in palette that is most appealing
sink("palette5.txt")
print(palette)
sink()

ps_ra_green = subset_samples(ps_ra, Genotype == "G")
ps_ra_red = subset_samples(ps_ra, Genotype == "R")
ps_ra_yellow = subset_samples(ps_ra, Genotype == "Y")

p1=plot_bar(ps_ra_green, fill="Order")+
  geom_bar(aes(fill=Order), stat="identity",position="stack")+
  theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+
  scale_fill_manual(values=palette)+
  ggtitle("Genotype G")+
  facet_grid(.~Colony,scales="free",space="free")+
  theme(plot.title = element_text(face="italic"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")
p1
p2=plot_bar(ps_ra_red, fill="Order")+
  geom_bar(aes(fill=Order), stat="identity",position="stack")+
  theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+
  scale_fill_manual(values=palette)+
  ggtitle("Genotype R")+
  facet_grid(.~Colony,scales="free",space="free")+
  theme(plot.title = element_text(face="italic"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

p3=plot_bar(ps_ra_yellow, fill="Order")+
  geom_bar(aes(fill=Order), stat="identity",position="stack")+
  theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+
  scale_fill_manual(values=palette)+
  ggtitle("Genotype Y")+
  facet_grid(.~Colony,scales="free",space="free")+
  theme(plot.title = element_text(face="italic"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

####adjust width and height until it looks right for double columns
pdf("Acropora_BarCharts_Order5.pdf",width=24, height=10)
plot_grid(p1,p2,p3,labels=c("A","B","C"), ncol=2, nrow=2)
dev.off()
### additional edits performed in inkscape: designation of colonies on x-axis and manually adding legend in bottom right
# to get legend, plot p3 and change legend.position to right

pdf("for_the_legend5.pdf",width=24)
p3=plot_bar(ps_ra_yellow, fill="Order")+
  geom_bar(aes(fill=Order), stat="identity",position="stack")+
  theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+
  scale_fill_manual(values=palette)+
  ggtitle("Genotype Y")+
  facet_grid(.~Colony,scales="free",space="free")+
  theme(plot.title = element_text(face="italic"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "right")
p3
dev.off()

########################### FIGURE 4 - MD3-55 ASVs #########################################################################

rick<-subset_taxa(ps_ra, Genus=="MD3-55")
rick
otu.rick = as(otu_table(rick), "matrix")
taxon.rick = as(tax_table(rick), "matrix")
meta.rick = as(sample_data(rick), "matrix")
otu.rick<-as.data.frame(otu.rick)
otu.rick<-rownames_to_column(otu.rick,var="Sample")
#export ASV sequences for supplemental table
write.table(otu.rick,"MD3-55_ASVs.txt",sep="\t",col.names=NA)
#make short ASV names
names(otu.rick)[2]<-"ASV1"
names(otu.rick)[3]<-"ASV2"
names(otu.rick)[4]<-"ASV3"
names(otu.rick)[5]<-"ASV4"
names(otu.rick)[6]<-"ASV5"
meta.rick<-as.data.frame(meta.rick)
meta.rick<-rownames_to_column(meta.rick,var="Sample")
otu.rick.meta<-merge(meta.rick,otu.rick,"Sample")
otu_long<-melt(otu.rick.meta,id.vars=c("Sample","Genotype","Colony","Branch"),variable.name="ASV",value.name="Proportion")
cols2<-c("ASV1"="#56B4E9","ASV2"="#0072B2","ASV3"="#999999","ASV4"="#E69F00","ASV5"="#009E73")

pdf("MD3-55_bars.pdf",width=8.5)
p1<-ggplot(otu_long, aes(x=Sample,y=Proportion))+
  geom_bar(aes(fill=ASV), color="#333333", stat="identity",position="stack")+
  facet_grid(.~Genotype,scales="free",space="free")+
  theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90, size=8))+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=14))+
  scale_fill_manual(values=cols2)
p1
dev.off()


################################# ANCOM TEST OF DIFFERENTIALLY ABUNDANT FAMILIES ############################################
#ANCOM Function - compare across multiple treatments groups using a compositional appproach
#https://sites.google.com/site/siddharthamandal1985/research

###Need to run this first in order to run ANCOM on data
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}


#####ANCOM results
otu <- read.table("silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps
#10150 taxa and 81 samples
get_taxa_unique(ps, "Family") #557


#Need to run ANCOM function first
#ANCOM
#from phyloseq -use all of the data, not relative abundance, and  not filtered beyond mitochondria and choloroplasts 
#Make the otutable
#group the otus based on family
#this is probaby a little overkill, but it's how I've done it in the past - all you need is the otutable really, so probably that could have been done in one step, but i like this for future plot-making
dat <- tax_glom(ps, taxrank = "Family") #at the Family level

#melt the data, so it's like a dataframe
datm <- psmelt(dat)

#Cast the new datatable with columns that are of interest
datc <- data.table::dcast(datm, Sample + Genotype + Colony ~ Family, value.var = 'Abundance', fun.aggregate = sum)

dim(datc) #dimensions of the table

otud <- datc[,c(1,4:560)] #select the first column, and then all of the taxa columns  
colnames(otud)[1] <- "Sample.ID" #rename the first column to Sample.ID - this is to match ANCOM syntax

metadat <- sample_data(ps) #get the sample data
metadat <- as.data.frame(as.matrix(metadat)) #make into into a matrix
# at this point, my sample id numbers are the row names, not a separate column, move row names to column with dplylr

metadat <- tibble::rownames_to_column(metadat, "Sample.ID") #make sure the sample names column is called Sample.ID

names(otud) <- make.names(names(otud)) #get the names from the table for 
otu_test <- otud #rename otud to otu_test, for syntax in ANCOM

metadat <- select(metadat, c("Sample.ID","Genotype","Colony")) # use select to only use treatment columns of interest
map_test <- metadat #rename map_TEst
Vardat <- map_test #specify that this for Vardat - ANCOM syntax
#### ANCOM test - not adjusted, more than 2 levels = Kruskal Wallis
comparison_test_treat=ANCOM.main(OTUdat=otu_test, #calling the OTU table
                                 Vardat=map_test, #calling the metadata
                                 adjusted=FALSE, #true if covariates are to be included for adjustment
                                 repeated=FALSE, #repeated measure
                                 main.var="Genotype", #main variable or fator
                                 adj.formula= NULL, #other factors to include
                                 repeat.var=FALSE, #repeated measure
                                 long = FALSE, #longitudinal study
                                 multcorr=2,
                                 sig=0.05, #significance level
                                 prev.cut=0.90) #OTUs with proportion of zeroes greater than prev.cut are not included in the analysis

res <- comparison_test_treat$W.taxa #taxa that significantly vary across factor level of interest
write.table(res,"ANCOM_family_KruskallWallis_Genotype.txt",sep="\t",col.names=NA)
res2 <- res[which(res$detected_0.7==TRUE),] 

#### ANCOM test - Adjusted by Coral species, ANOVA
comparison_test_treat=ANCOM.main(OTUdat=otu_test, #calling the OTU table
                                 Vardat=map_test, #calling the metadata
                                 adjusted=TRUE, #true if covariates are to be included for adjustment
                                 repeated=FALSE, #repeated measure
                                 main.var="Genotype", #main variable or fator
                                 adj.formula= "Colony", #other factors to include
                                 repeat.var=FALSE, #repeated measure
                                 long = FALSE, #longitudinal study
                                 multcorr=2,
                                 sig=0.05, #significance level
                                 prev.cut=0.90) #OTUs with proportion of zeroes greater than prev.cut are not included in the analysis

res3 <- comparison_test_treat$W.taxa #taxa that significantly vary across factor level of interest
write.table(res3,"ANCOM_family_ANOVA_Genotype_adjColony.txt",sep="\t",col.names=NA)
res4 <- res[which(res3$detected_0.7==TRUE),] 

sig_sites <- glue::glue_collapse(droplevels(factor(res2$otu.names)), sep = ", ") #this is to get a list of the families that are different
print(sig_sites)
#Alphaproteobacteria, Francisellaceae, Desulfobacteraceae, JGI_0000069.P22, Pirellulaceae, Pseudomonadaceae, Halomonadaceae, Woesearchaeia, Alteromonadaceae, Clade_III, Marinimicrobia_.SAR406_clade., Puniceicoccaceae, S25.593, Fusobacteriaceae, Xenococcaceae, Phycisphaeraceae, PB19

#Calculate relative abundance
datc_relabund <-  sweep(datc[,4:560], 1, rowSums(datc[,4:560]), '/')
datc_relnames <- cbind(datc[,1:3],datc_relabund)

#only selet the significant families
sig_dis <- select(datc_relnames, Sample, Genotype, Colony, Alphaproteobacteria, Francisellaceae, Desulfobacteraceae, "JGI_0000069-P22", Pirellulaceae, Pseudomonadaceae, Halomonadaceae, Woesearchaeia, Clade_III, "Marinimicrobia_(SAR406_clade)", Puniceicoccaceae, "S25-593", Fusobacteriaceae, Xenococcaceae, Phycisphaeraceae, PB19)
sig_long <- melt(sig_dis, id.vars=c("Sample","Genotype","Colony"),variable.name="Family",value.name="Proportion")
sum_sig <- Rmisc::summarySE(sig_long, measurevar = "Proportion", groupvars = c("Genotype","Family"), na.rm=TRUE)

cols<-c("G"="#009E73","R"="#D55E00","Y"="#F0E442")
sum_sig$Genotype<-factor(sum_sig$Genotype, levels=c("G","R","Y"))

############ FIGURE 5 Differentially abundant families by coral genotype #####################################################
pdf("ANCOM_Families_Genotype.pdf",width=8.5)
fams <- ggplot(sum_sig, aes(x=Family, y=Proportion+0.001))+
  geom_point(size=4, aes(fill=Genotype,shape=Genotype))+
  scale_fill_manual(values=cols)+
  scale_shape_manual(values=c(21,22,24))+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.title.y=element_text(size=14))+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_errorbar(aes(ymin=Proportion+0.001-se, ymax=Proportion+0.001+se), width=.1)+
  theme(legend.title = element_blank())+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme(legend.text = element_text(size=12))
fams
dev.off()


