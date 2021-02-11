#bacterial community PCA two plot Fig2 E and F CBASS and RSS and PERMANOVA
######


###### graph with all libraries but without field samples
#install.packages("ggplot2")
library("ggplot2")
#install.packages("gridExtra")
library(gridExtra)

#CBASS
bactdat<-read.table("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/02_CBASS_vs_RSS_T1-T2/01_Paper_Versions/26_Revisions_from_Anny_with_bacteria_ASV/PCA_ASVs_CBASS.txt",header=T)
#### number of sequence for size of symboles:
bactsize<-read.table("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/02_CBASS_vs_RSS_T1-T2/01_Paper_Versions/26_Revisions_from_Anny_with_bacteria_ASV/Bacteria_reads_number_ASV.txt", header=T)
bactsizeCBASS<-subset(bactsize,bactsize$experiment=="CBASS")

#### check if both dataset are correctly ordered:
paste(rownames(bactdat),bactsizeCBASS$Sample, sep=" = ")

#genotypes info
geno<-rownames(bactdat)
geno<-gsub("29-5", "295", geno)
geno<-gsub("34-5", "345", geno)
a<-strsplit(geno,"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df

a<-strsplit(geno,"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df


a<-strsplit(geno,"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][4]
}
geno<-df

genotime<-paste(geno$matrix.NA..nrow...40..ncol...1.,time$matrix.NA..nrow...40..ncol...1.,sep="-")

genotime<-gsub("G13A-T1",15,genotime)
genotime<-gsub("G14A-T1",16,genotime)
genotime<-gsub("G15A-T1",17,genotime)
genotime<-gsub("G8A-T1",23,genotime)
genotime<-gsub("G9A-T1",25,genotime)
genotime<-gsub("G13A-T2",0,genotime)
genotime<-gsub("G14A-T2",1,genotime)
genotime<-gsub("G15A-T2",2,genotime)
genotime<-gsub("G8A-T2",5,genotime)
genotime<-gsub("G9A-T2",6,genotime)

temp<-temp$matrix.NA..nrow...40..ncol...1.
temp<-gsub("27","blue",temp)
temp<-gsub("295","yellow",temp)
temp<-gsub("32","orange",temp)
temp<-gsub("345","red",temp)

g5<-ggplot(bactdat, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(size=2),fill=temp, col=temp, shape=as.numeric(genotime),size=log(as.numeric(gsub("'","",bactsizeCBASS$nonchimeric_sequences)))/2) + ggtitle("PCA Bacterial community short-term heat stress") +
  xlab("PC1: 20.1% variance") + ylab("PC2: 9.9% variance") + theme(
    plot.title = element_text(size=15, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12))

g5


##### RSS
bactdat<-read.table("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/02_CBASS_vs_RSS_T1-T2/01_Paper_Versions/26_Revisions_from_Anny_with_bacteria_ASV/PCA_ASVs_RSS.txt",header=T)
#### number of sequence for size of symboles:
bactsize<-read.table("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/02_CBASS_vs_RSS_T1-T2/01_Paper_Versions/26_Revisions_from_Anny_with_bacteria_ASV/Bacteria_reads_number_ASV.txt", header=T)
bactsizeRSS<-subset(bactsize,bactsize$experiment=="RSS")

#### check if both dataset are correctly ordered:
paste(rownames(bactdat),bactsizeRSS$Sample, sep=" = ")

#genotypes info
geno<-rownames(bactdat)
geno<-gsub("29-5", "295", geno)
geno<-gsub("34-5", "345", geno)
a<-strsplit(geno,"-")
df <- data.frame(matrix(NA,nrow=28,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df

a<-strsplit(geno,"-")
df <- data.frame(matrix(NA,nrow=28,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df


a<-strsplit(geno,"-")
df <- data.frame(matrix(NA,nrow=28,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][4]
}
geno<-df

genotime<-paste(geno$matrix.NA..nrow...28..ncol...1.,time$matrix.NA..nrow...28..ncol...1.,sep="-")

genotime<-gsub("G13A-T1",15,genotime)
genotime<-gsub("G14A-T1",16,genotime)
genotime<-gsub("G15A-T1",17,genotime)
genotime<-gsub("G8A-T1",23,genotime)
genotime<-gsub("G9A-T1",25,genotime)
genotime<-gsub("G13A-T2",0,genotime)
genotime<-gsub("G14A-T2",1,genotime)
genotime<-gsub("G15A-T2",2,genotime)
genotime<-gsub("G8A-T2",5,genotime)
genotime<-gsub("G9A-T2",6,genotime)

temp<-temp$matrix.NA..nrow...28..ncol...1.
temp<-gsub("27","blue",temp)
temp<-gsub("295","yellow",temp)
temp<-gsub("32","orange",temp)
temp<-gsub("345","red",temp)

g6<-ggplot(bactdat, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(size=2),fill=temp, col=temp, shape=as.numeric(genotime),size=log(as.numeric(gsub("'","",bactsizeRSS$nonchimeric_sequences)))/2) + ggtitle("PCA Bacterial community long-term heat stress") +
  xlab("PC1: 26.3% variance") + ylab("PC2: 10.6% variance") + theme(
    plot.title = element_text(size=15, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12))

g6


grid.arrange(g5,g6,nrow=1,ncol=2) 

### run script of RNAseq PCA plot and add to grid:



grid.arrange(g1,g2,g3,g4,g5,g6,nrow=3,ncol=2) 


