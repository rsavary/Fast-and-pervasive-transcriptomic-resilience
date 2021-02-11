###### Four plot from Plink for SNPS
library("ggplot2")
# Smicro CBASS:
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/40_New_SNPs_calling_and_plotting_all_samples_CBASS_RSS/Spis/")

library("adegenet")
library("ade4")
df<-read.PLINK(file = "plink.raw", chunkSize=1000)


#### not neeeded if missing genotype were removed before in plink with option --geno 0
toRemove <- is.na(glMean(df, alleleAsUnit = FALSE))
which(toRemove) # position of entirely non-typed loci
b <- df[, !toRemove]
pca1SpisCBASS<-glPca(b)
2

SNPsimportance<-pca1SpisCBASS$loadings
loadingplot(SNPsimportance)

SNPS<-as.data.frame(SNPsimportance)
SNPS2<-SNPS[order(-SNPS$Axis1),]

#### get percentage for each axis: https://groups.google.com/forum/#!topic/poppr/yUdJycNYrhc
var_fracpca1SpisCBASS<- pca1SpisCBASS$eig/sum(pca1SpisCBASS$eig)


pca1SpisCBASS<-as.data.frame(pca1SpisCBASS$scores)


### ggplot
library(ggplot2)
require(gridExtra)


names<-df$ind.names
names2<-gsub("29-5","295" , names)
names2<-gsub("34-5","345" , names2)
a<-strsplit(names2,"-")
df <- data.frame(matrix(NA,nrow=71,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][4]
}
gen<-df[,1]
gen<-gsub('.bam', '', gen)

gen<-gsub("G13A","22" , gen)
gen<-gsub("G14A","21" , gen)
gen<-gsub("G15A","24" , gen)
gen<-gsub("G8A","23" , gen)
gen<-gsub("G9A","25" , gen)

gen2<-as.numeric(as.character(gen))

g2<-ggplot(pca1SpisCBASS, aes(x=PC1, y=PC2)) +
  geom_point(col="black", shape=gen2, size=5) + ggtitle("PCA genomic var. Spis (76 samples / 43'194 SNPs)") +
  xlab(paste("PC1: ",round(var_fracpca1SpisCBASS[1]*100,1),"%",sep="")) + ylab(paste("PC2: ",round(var_fracpca1SpisCBASS[2]*100,1),"%",sep="")) + theme(
    plot.title = element_text(size=15, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12))


##### Smic SNPs plotting



a<-c(22,21,24,23,25)
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/40_New_SNPs_calling_and_plotting_all_samples_CBASS_RSS/Smic")

library("adegenet")
library("ade4")
df<-read.PLINK(file = "plink.raw", chunkSize=1000)


#### not neeeded if missing genotype were removed before in plink with option --geno 0
toRemove <- is.na(glMean(df, alleleAsUnit = FALSE))
which(toRemove) # position of entirely non-typed loci
b <- df[, !toRemove]
pca1SmicroCBASS<-glPca(b)
2

SNPsimportance<-pca1SmicroCBASS$loadings
loadingplot(SNPsimportance)

SNPS<-as.data.frame(SNPsimportance)
SNPS2<-SNPS[order(-SNPS$Axis1),]

#### get percentage for each axis: https://groups.google.com/forum/#!topic/poppr/yUdJycNYrhc
var_fracpca1SmicroCBASS<- pca1SmicroCBASS$eig/sum(pca1SmicroCBASS$eig)


pca1SmicroCBASS<-as.data.frame(pca1SmicroCBASS$scores)


### ggplot
library(ggplot2)
require(gridExtra)


names<-df$ind.names
names2<-gsub("29-5","295" , names)
names2<-gsub("34-5","345" , names2)
a<-strsplit(names2,"-")
df <- data.frame(matrix(NA,nrow=71,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][4]
}
gen<-df[,1]
gen<-gsub('.bam', '', gen)

gen<-gsub("G13A","22" , gen)
gen<-gsub("G14A","21" , gen)
gen<-gsub("G15A","24" , gen)
gen<-gsub("G8A","23" , gen)
gen<-gsub("G9A","25" , gen)

gen2<-as.numeric(as.character(gen))

g1<-ggplot(pca1SmicroCBASS, aes(x=PC1, y=PC2)) +
  geom_point(col="black", shape=gen2, size=5) + ggtitle("PCA genomic var. Smicro (71 samples / 19023 SNPs)") +
  xlab(paste("PC1: ",round(var_fracpca1SmicroCBASS[1]*100,1),"%",sep="")) + ylab(paste("PC2: ",round(var_fracpca1SmicroCBASS[2]*100,1),"%",sep="")) + theme(
    plot.title = element_text(size=15, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12))



grid.arrange(g2,g1,nrow=2)    

