######
#Script of RNA-seq comparison Kallisto Stylophora Pistillata 
###Library needed
library("DESeq2")
library("limma")
library("edgeR")
library("ggplot2")
library("tximport")
library("adegenet")
library("pheatmap")



########################################
##########Kallisto for Stylophora pistillata ####
########################################
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/19_Symportal_ITS2_mapping_results/Abundance_Kallisto_RSS_vs_CBASS_all/")
filesToProcessSK <- dir(pattern = "*_abundance.tsv$")  #files to process.
names(filesToProcessSK)<-gsub('_abundance.tsv$','',filesToProcessSK,perl=TRUE)
samplesSK<-read.table('/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/19_Symportal_ITS2_mapping_results/Abundance_Kallisto_RSS_vs_CBASS_all/CB-RSS-field-23-G13A_abundance.tsv',h=T)

tx2geneSK<-cbind.data.frame(samplesSK$target_id,gsub('.[0-9].v6.1$','',samplesSK$target_id,perl=TRUE))
colnames(tx2geneSK)<-c('TXNAME','GENEID')
#add tximport function manualy, as the tximport is not made for this version. 
txiSK <- tximport(filesToProcessSK, type="kallisto", tx2gene=tx2geneSK, countsFromAbundance="lengthScaledTPM") # normalized library size and transcript length # countsFromAbundance="lengthScaledTPM" or countsFromAbundance="no"
summary(txiSK)
#######following protocol manual deseq2 page 8
#labels and design
countsSK<-txiSK[[2]]
header<-colnames(countsSK) # header<-gsub("\\d", "", colnames(counts) )
#temp=c(23,23,23,32,34,29,29,34,34,34,27,27,27,27,29,29,29,29,32,32,27,27,27,27,27,29,29,29,32,32,32,34,34) #temp
#temp=c(4,4,4,4,4,4,4,2,2,1,1,1,2,1,1,2,1,1,4,4,1,2,1,1,2,2,1,2,2,2,2,4,4) # lane
#temp=c("T0","T0","T0","T1","T1","T3","T3","T3","T3","T3","T1","T1","T1","T1","T1","T1","T1","T1","T1","T1","T3","T3","T3","T3","T3","T3","T3","T3","T3","T3","T3","T3","T3") # time
#temp<-c(1,1,1,1,34,1,1,34,1,1,1,1,1,1,1,1,1,1,1,1,34,34)
#temp<-c("Field","Field","Field","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS")

######## modification to the data after first analyze:
# 1) exchange the two field samples that were exchanged.
# 2) remove low counts libraries.



#for ITS database 
mat <- matrix(, nrow = dim(countsSK)[2], ncol = 8)
mat2<-as.data.frame(mat)
colnames(mat2)<-c("Nb_mapped_reads","mean_mapped_reads_per_transcript","max_mapped_reads_toatranscript","total_nb_transcript","nb_transcript>0reads","nb_transcript>10reads","nb_transcript>100reads", "nb_transcript>1000reads")
rownames(mat2)<-colnames(countsSK)
mat3<-mat2
for(i in 1:length(colnames(countsSK))){
  mat3[i,1]<-sum(countsSK[,i])
  mat3[i,2]<-mean(countsSK[,i])
  mat3[i,3]<-max(countsSK[,i])
  mat3[i,4]<-sum(countsSK[,i]>=0)
  mat3[i,5]<-sum(countsSK[,i]>0)
  mat3[i,6]<-sum(countsSK[,i]>10)
  mat3[i,7]<-sum(countsSK[,i]>100)
  mat3[i,8]<-sum(countsSK[,i]>1000)
  
}

SymITS_summary_Kallisto<-mat3

write.csv(as.data.frame(SymITS_summary_Kallisto),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/19_Symportal_ITS2_mapping_results/Symbiodinium_ITSBENHUME_summary_Kallisto")
summary<-read.csv(file = "/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/08_CBASS_vs_RSS_Symbiodinium_ITS2_BEN_HUME_Symportal/Symbiodinium_ITSBENHUME_summary_Kallisto", header=F)



a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=11,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
temp<-df[,1]
temp<-gsub('[-]', '', temp)


##### temp
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]


data<-txiSK$counts
data[is.nan(data)] <- 0
data<-data[rowSums(data) >1,]

q<-data[,-c(1:5)]
dim(q)
df <- data.frame(matrix(NA,nrow=dim(q)[1],ncol=dim(q)[2]))
colnames(df)<-colnames(q)
rownames(df)<-rownames(q)

for(i in 1:dim(q)[2]){
  for(j in 1:dim(q)[1]){
    df[j,i]<-q[j,i]/sum(q[,i])*100
  }
}

colfunc<-colorRampPalette(c("blue", "darkblue"))
colfunc2<-colorRampPalette(c("red", "darkred"))
col1<-colfunc(dim(q)[1])
col1[66]<-"pink"
col1[67]<-"orange"
col1[68]<-colfunc2(10)[3]
col1[69]<-colfunc2(10)[6]
col1[70]<-"firebrick4"
col1[71]<-colfunc2(10)[8]
col1[72:75]<-"red"
col1[c(2:3)]<-c("#FFFF00","#C3FF00")
col1[c(4:5)]<-c("lightblue","lightblue")
par(mfrow=c(1,1),mar=c(10,4,3,7))
par(xpd=TRUE)
b<-barplot(as.matrix(df), col=col1, las=2, cex.names = 0.7, cex.axis=0.7, bty="L", main="Barplot of RNAseq reads pseudoaligning to ITS2 database from Symportal 2019 ", ylab="% of reads aligned to the Symbiodinium ITS database (Symportal) for each type",cex.main=0.7)
legend(100,110,legend = rownames(df),fill =col1,col = col1, cex=0.2)
text(b,101, labels = round(colSums(q)),las=1, cex=0.4)











