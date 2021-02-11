######
#Script of RNA-seq with Kallisto; Stylophora pistillata 
###Library needed
library("DESeq2")
library("limma")
library("edgeR")
library("ggplot2")
library("tximport")
library("adegenet")
library("pheatmap")

#install.packages("devtools") for post hoc test ADONIS 
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

########################################
########## Kallisto for Stylophora pistillata ####
########################################
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/02_Kallisto_abundance_Smicro-Spist/Abundance_Kallisto_RSS_vs_CBASS_without_field/")

filesToProcessSK <- dir(pattern = "*_abundance.tsv$")  #files to process.
names(filesToProcessSK)<-gsub('_abundance.tsv$','',filesToProcessSK,perl=TRUE)
samplesSK<-read.table('/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/02_Kallisto_abundance_Smicro-Spist/Abundance_Kallisto_RSS_vs_CBASS_without_field/CB-T1-27-G13A_abundance.tsv',h=T)

tx2geneSK<-cbind.data.frame(samplesSK$target_id,gsub('.[0-9].v6.1$','',samplesSK$target_id,perl=TRUE))
colnames(tx2geneSK)<-c('TXNAME','GENEID')
#add tximport function manualy, as the tximport is not made for this version. 
txiSK <- tximport(filesToProcessSK, type="kallisto", tx2gene=tx2geneSK, countsFromAbundance="lengthScaledTPM") # normalized library size and transcript length # countsFromAbundance="lengthScaledTPM" or countsFromAbundance="no"
summary(txiSK)

#######following protocol manual deseq2 page 8
#labels and design
countsSK<-txiSK[[2]]
header<-colnames(countsSK) # header<-gsub("\\d", "", colnames(counts) )


# get temp from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]

# get T1/T3 from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
temp<-df[,1]


#exp
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][1]
}
temp<-df[,1]


#genotype
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
temp<-df[,1]
temp<-gsub('[-]', '', temp)

###### summary matrix of pseudoaligned reads with Kallisto ####

#for Smicro
mat <- matrix(, nrow = dim(countsSK)[2], ncol = 8)
mat2<-as.data.frame(mat)
colnames(mat2)<-c("Nb_mapped_reads","mean_mapped_reads_per_transcript","max_mapped_reads_toatranscript","total_nb_transcript","nb_transcript>0reads","nb_transcript>10reads","nb_transcript>100reads", "nb_transcript>1000reads")
rownames(mat2)<-colnames(countsSK)
mat3<-mat2
for(i in 1:length(colnames(countsSK))){
  mat3[i,1]<-sum(countsSK[1:49109,i])
  mat3[i,2]<-mean(countsSK[1:49109,i])
  mat3[i,3]<-max(countsSK[1:49109,i])
  mat3[i,4]<-sum(countsSK[1:49109,i]>=0)
  mat3[i,5]<-sum(countsSK[1:49109,i]>0)
  mat3[i,6]<-sum(countsSK[1:49109,i]>10)
  mat3[i,7]<-sum(countsSK[1:49109,i]>100)
  mat3[i,8]<-sum(countsSK[1:49109,i]>1000)
  
}

Smicro_summary_Kallisto<-mat3

write.csv(as.data.frame(Smicro_summary_Kallisto),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/02_Kallisto_abundance_Smicro-Spist/Abundance_Kallisto_RSS_vs_CBASS_without_field/Smicroadriaticum_summary_Kallisto")
summary<-read.csv(file = "/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/02_Kallisto_abundance_Smicro-Spist/Summary_Kallisto_pseudoalignement.txt", header=F)
summary<-summary[6:85,]

#for Spist
mat <- matrix(, nrow = dim(countsSK)[2], ncol = 8)
mat2<-as.data.frame(mat)
colnames(mat2)<-c("Nb_mapped_reads","mean_mapped_reads_per_transcript","max_mapped_reads_toatranscript","total_nb_transcript","nb_transcript>0reads","nb_transcript>10reads","nb_transcript>100reads", "nb_transcript>1000reads")
rownames(mat2)<-colnames(countsSK)
mat3<-mat2
for(i in 1:length(colnames(countsSK))){
  mat3[i,1]<-sum(countsSK[49110:74878,i])
  mat3[i,2]<-mean(countsSK[49110:74878,i])
  mat3[i,3]<-max(countsSK[49110:74878,i])
  mat3[i,4]<-sum(countsSK[49110:74878,i]>=0)
  mat3[i,5]<-sum(countsSK[49110:74878,i]>0)
  mat3[i,6]<-sum(countsSK[49110:74878,i]>10)
  mat3[i,7]<-sum(countsSK[49110:74878,i]>100)
  mat3[i,8]<-sum(countsSK[49110:74878,i]>1000)
  
}

Spist_summary_Kallisto<-mat3


write.csv(as.data.frame(Spist_summary_Kallisto),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/02_Kallisto_abundance_Smicro-Spist/Abundance_Kallisto_RSS_vs_CBASS_without_field/Spistillata_summary_Kallisto")
summary<-read.csv(file = "/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/02_Kallisto_abundance_Smicro-Spist/Summary_Kallisto_pseudoalignement.txt", header=F)
summary<-summary[6:85,]



######## Fig S1
Totalsummmary<-cbind(summary,Spist_summary_Kallisto)
new<-Totalsummmary$V2-Spist_summary_Kallisto$Nb_mapped_reads
new2<-new-Smicro_summary_Kallisto$Nb_mapped_reads
new3<-t(cbind(Spist_summary_Kallisto$Nb_mapped_reads, Smicro_summary_Kallisto$Nb_mapped_reads, new2))

colnames(new3)<-Totalsummmary$V1

par(mfrow=c(1,1),mar=c(7,4,3,2))
b<-barplot(new3, las=2, col=c("blue", "green", "red"), cex.names = 0.6, cex.axis=0.5)

options(digits=1)
spistper<-new3[1,]/colSums(new3)*100
smicroper<-new3[2,]/colSums(new3)*100
unmap<-new3[3,]/colSums(new3)*100
text(b,new3[1,]/2, labels = round(spistper, digits=0),las=2, cex=0.6)
text(b,new3[1,]+((new3[2,])/2), labels = round(smicroper, digits=0),las=2, cex=0.6)
text(b,new3[1,]+new3[2,]+((new3[3,])/2), labels = round(unmap, digits=0),las=2, cex=0.6)
legend(1,42000000,legend = c("S. pistillata","S. microadriaticum (Clade A, A1)","non-pseudo-align"),fill = c("blue", "green", "red"),col = c("blue", "green", "red"), cex=0.6)


###########
#### 1 ####
###########

#######################################
#design only CBASS symbiodinium 
#######################################
#######################################
# get temp from header
# get temp from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]

# get T1/T3 from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]


#exp
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][1]
}
exp<-df[,1]


#genotype
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
gen<-df[,1]
gen<-gsub('[-]', '', gen)

temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime)
# test inculding genotype
temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # for CTL and AMF
#####get rid of library of the RSS
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[1:49109,-c(41:80)] #to remove lib
design2<-design[-c(41:80),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
#ddsSKtemprlog<-rlog(ddsSKtemp)

#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="gen")

vstCBASSsymb<-assay(ddsSKtempvst)
CBsymb<-plotPCA(ddsSKtempvst,intgroup="temp")
#plotPCA(ddsSKtempvst, ntop=100)

###########
#### 2 ####
###########

#######################################
#design only RSS symbiodinium 
#######################################
#######################################
# get temp from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]

# get T1/T3 from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]


#exp
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][1]
}
exp<-df[,1]


#genotype
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
gen<-df[,1]
gen<-gsub('[-]', '', gen)

temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime)
# test inculding genotype
temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # for CTL and AMF
#####get rid of library of the RSS
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[1:49109,-c(1:40)] #to remove lib
design2<-design[-c(1:40),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
#ddsSKtemprlog<-rlog(ddsSKtemp)
vstRSSsymb<-assay(ddsSKtempvst)
#write.csv(vstRSSsymb,"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/15_Kallisto-VST_normalized_counts/Smicro_RSS_vst_normalized_counts.csv")
#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst,intgroup="temp")
RSSsymb<-plotPCA(ddsSKtempvst, intgroup="temp")

###########
#### 3 ####
###########

#######################################
#design only CBASS spistillata
#######################################
#######################################
# get temp from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]

# get T1/T3 from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]


#exp
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][1]
}
exp<-df[,1]


#genotype
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
gen<-df[,1]
gen<-gsub('[-]', '', gen)

temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime)
# test inculding genotype
temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # for CTL and AMF
#####get rid of library of the RSS
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[49110:74878,-c(41:80)] #to remove lib
design2<-design[-c(41:80),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
#ddsSKtemprlog<-rlog(ddsSKtemp)
vstCBASScoral<-assay(ddsSKtempvst)
#write.csv(vstCBASScoral,"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis//15_Kallisto-VST_normalized_counts/Spist_CBASS_vst_normalized_counts.csv")

#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="temp")
CBpist<-plotPCA(ddsSKtempvst, intgroup="temp")


###########
#### 4 ####
###########

#######################################
#design only RSS spistillata
#######################################
#######################################
# get temp from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]

# get T1/T3 from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]


#exp
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][1]
}
exp<-df[,1]


#genotype
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
gen<-df[,1]
gen<-gsub('[-]', '', gen)

temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime)
# test inculding genotype
temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # for CTL and AMF
#####get rid of library of the RSS
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[49110:74878,-c(1:40)] #to remove lib
design2<-design[-c(1:40),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
#ddsSKtemprlog<-rlog(ddsSKtemp)
vstRSScoral<-assay(ddsSKtempvst)
#write.csv(vstRSScoral,"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis//15_Kallisto-VST_normalized_counts/Spist_RSS_vst_normalized_counts.csv")

#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="temp")
RSSpist<-plotPCA(ddsSKtempvst, intgroup="temp")





####### 4 plots togethere:

CBpistdat<-CBpist$data
RSSpistdat<-RSSpist$data
CBsymbdat<-CBsymb$data
RSSsymbdat<-RSSsymb$data
######################### col
#generate color
#colfunc<-colorRampPalette(c("blue", "red"))
col1<-c("blue","yellow","orange","red")
#############
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
col<-df[,1]
col<-col[1:40]
col<-gsub("27", col1[1], col)
col<-gsub("29", col1[2], col)
col<-gsub("32", col1[3], col)
col<-gsub("34", col1[4], col)


colCBASS<-col[1:40]
colRSSSpist<-col[1:40]
colRSSSmicro<-col[1:40]

rownames(CBpistdat)
#### add several graph in same window, list of color and shape to decide.
#### shape full for T1 empty for T3, square, rond triangle and other for genotypes
#genotype
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
temp<-df[,1]
temp<-gsub('[-]', '', temp)
gen<-temp
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
temp<-df[,1]


gen<-paste(gen,temp,sep="_")

# if we want to add genetic info for the supp graph, to do!!!!!
#gen[1:5]<-c("15","16","17","23","25") # full square, full circle, triangle up full, diamond full, triangle down
gen[1:20]<-c("15","16","17","23","25")
gen[21:40]<-c("0","1","2","5","6")
gen[41:60]<-c("15","16","17","23","25")
gen[61:80]<-c("0","1","2","5","6")

#gen[1:20]<-c("16","16","16","16","16")
#gen[21:40]<-c("1","1","1","1","1")
#gen[41:60]<-c("16","16","16","16","16")
#gen[61:80]<-c("1","1","1","1","1")

#CBgen<-gen[1:45]
#CBgen<-as.numeric(CBgen)
#RSSgen<-gen[c(1:5,46:85)]
# ggplot
library(ggplot2)
require(gridExtra)

g1<-ggplot(CBpistdat, aes(x=PC1, y=PC2)) +
  geom_point(fill=col, col=col, shape=as.numeric(gen[1:40]), size=log(Spist_summary_Kallisto$Nb_mapped_reads/5000)[1:40]) + ggtitle("PCA CBASS S. pistillata") +
  xlab(CBpist$labels[2]) + ylab(CBpist$labels[1]) + theme(
    plot.title = element_text(size=15, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12))

g2<-ggplot(RSSpistdat, aes(x=PC1, y=PC2)) +
  geom_point(fill=col, col=col, shape=as.numeric(gen[1:40]), size=log(Spist_summary_Kallisto$Nb_mapped_reads/5000)[c(41:80)]) + ggtitle("PCA RSS S. pistillata") +
  xlab(RSSpist$labels[2]) + ylab(RSSpist$labels[1]) + theme(
    plot.title = element_text(size=15, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12))

g3<-ggplot(CBsymbdat, aes(x=PC1, y=PC2)) +
  geom_point(aes(size=2),fill=col, col=col, shape=as.numeric(gen[1:40]),size=log(Smicro_summary_Kallisto$Nb_mapped_reads/5000)[1:40]) + ggtitle("PCA CBASS S. microadriaticum") +
  xlab(CBsymb$labels[2]) + ylab(CBsymb$labels[1]) + theme(
    plot.title = element_text(size=15, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12))

g4<-ggplot(RSSsymbdat, aes(x=PC1, y=PC2)) +
  geom_point(aes(size=2),fill=col, col=col, shape=as.numeric(gen[1:40]),size=log(Smicro_summary_Kallisto$Nb_mapped_reads/5000)[c(41:80)]) + ggtitle("PCA RSS S. microadriaticum") +
  xlab(RSSsymb$labels[2]) + ylab(RSSsymb$labels[1]) + theme(
    plot.title = element_text(size=15, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12))

grid.arrange(g1,g2,g3,g4,nrow=2)                                                                                                             





##### PERMANOVA 1:
a<-strsplit(as.character(CBpistdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]


# get T1/T3 from header
a<-strsplit(as.character(CBpistdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]

#exp
a<-strsplit(as.character(CBpistdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][1]
}
exp<-df[,1]


#genotype
a<-strsplit(as.character(CBpistdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
gen<-df[,1]
gen<-gsub('[-]', '', gen)


##### PERMANOVA 1:

library(vegan)
CBpistdat
distCBpistPCA<-dist(CBpistdat[,1:2], method = "euclidean")
#adonis(distCBpistPCA ~ gen+temp+paste(time,temp,"_"))
adonis(distCBpistPCA ~ gen+temp*time)


#### test post hoc test
Y2<-data.frame(gen,temp,time,temp_t=paste(temp,"_",time,sep=""))
pairwise.adonis2(CBpistdat[,1:2] ~ temp_t, data = Y2)



##### PERMANOVA 2:
a<-strsplit(as.character(RSSpistdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]

# get T1/T3 from header
a<-strsplit(as.character(RSSpistdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]

#exp
a<-strsplit(as.character(RSSpistdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][1]
}
exp<-df[,1]


#genotype
a<-strsplit(as.character(RSSpistdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
gen<-df[,1]
gen<-gsub('[-]', '', gen)


RSSpistdat
distRSSpistPCA<-dist(RSSpistdat[,1:2], method = "euclidean")

adonis(distRSSpistPCA ~ gen+temp*time)
#adonis(distRSSpistPCA ~ gen+temp+paste(time,temp,"_"))


### post hoc test:
Y2<-data.frame(gen,temp,time,temp_t=paste(temp,"_",time,sep=""))
pairwise.adonis2(RSSpistdat[,1:2] ~ temp_t, data = Y2)


##### PERMANOVA 3:
a<-strsplit(as.character(CBsymbdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]

# get T1/T3 from header
a<-strsplit(as.character(CBsymbdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]

#exp
a<-strsplit(as.character(CBsymbdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][1]
}
exp<-df[,1]


#genotype
a<-strsplit(as.character(CBsymbdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
gen<-df[,1]
gen<-gsub('[-]', '', gen)


CBsymbdat
distCBsymPCA<-dist(CBsymbdat[,1:2], method = "euclidean")

adonis(distCBsymPCA ~ gen+temp*time)
#adonis(distCBsymPCA ~ gen+temp+paste(time,temp,"_"))

#### post hoc test:
### post hoc test:
Y2<-data.frame(gen,temp,time,temp_t=paste(temp,"_",time,sep=""))
pairwise.adonis2(CBsymbdat[,1:2] ~ temp_t, data = Y2)


##### PERMANOVA 4:
a<-strsplit(as.character(RSSsymbdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]

# get T1/T3 from header
a<-strsplit(as.character(RSSsymbdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]

#exp
a<-strsplit(as.character(RSSsymbdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][1]
}
exp<-df[,1]


#genotype
a<-strsplit(as.character(RSSsymbdat$name),"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
gen<-df[,1]
gen<-gsub('[-]', '', gen)

distRSSsymPCA<-dist(RSSsymbdat[,1:2], method = "euclidean")

adonis(distRSSsymPCA ~ gen+temp*time)
#adonis(distRSSsymPCA ~ gen+temp+paste(time,temp,"_"))
