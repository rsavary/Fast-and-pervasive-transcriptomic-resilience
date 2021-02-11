#
#Script of RNA-seq comparison Kallisto Stylophora Pistillata 
###Library needed
library("DESeq2")
library("limma")
library("edgeR")
library("ggplot2")
library("tximport")
library("adegenet")
library("pheatmap")
library("ape")




########################################
##########Kallisto for Stylophora pistillata ####
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

# get T1/T2 from header
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




#
#### 1  design only CBASS Smicro ####
#
# get temp from header

temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # for CTL and AMF
#####get rid of library CTL8,Q11, G4, G12 for CTL and AMF
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[1:49109,-c(41:80)] #to remove lib
design2<-design[-c(41:80),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene

ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
#ddsSKtemprlog<-rlog(ddsSKtemp)

#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="temp")


#Make a counts table that is scaled by the size factors
resSKtemp<-DESeq(ddsSKtemp,modelMatrixType="expanded", betaPrior=T )


temp = t(sizeFactors(resSKtemp))
sizematrix<-matrix(data=temp, nrow=nrow(counts(ddsSKtemp)), ncol=ncol(temp), byrow=TRUE)
scaledcounts = counts(ddsSKtemp)/sizematrix
head(scaledcounts)
head(sizematrix)

table<-cbind(rowMeans(scaledcounts[,1:5]),
      rowMeans(scaledcounts[,6:10]),
      rowMeans(scaledcounts[,11:15]),
      rowMeans(scaledcounts[,16:20]),
      rowMeans(scaledcounts[,21:25]),
      rowMeans(scaledcounts[,26:30]),
      rowMeans(scaledcounts[,31:35]),
      rowMeans(scaledcounts[,36:40]))
colnames(table)<-c("CB-T1-27","CB-T1-29","CB-T1-32","CB-T1-34","CB-T2-27","CB-T2-29","CB-T2-32","CB-T2-34")

new<-scaledcounts[,1:2]
for(i in 1:dim(table)[1]){
  new[i,1]<-(sum(table[i,] > 5)/8)*100   #### keep only site that have at least 5reads in 88%
}
new2<-subset(new,new[,1]>=80)
dim(new2)


resSKtemp2<-resSKtemp[rownames(counts(resSKtemp)) %in% rownames(new2), ]
########## log fold change shrinkage to remove gene and have a better visualization
resultsNames(ddsSKtemp)


library("gtools")
annotGOSm<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME/02_Symbiodinium/01_CLADE_A_Smicroadriaticum/smic_tabulated_annots-NEW.csv")
annotGOSm2<-transform(annotGOSm, Query = colsplit(Query, split = " ", names = c('Name','Namelong')))
colnames(annotGOSm2)[1]

#####results contrast "34T1","27T1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic12947", intgroup="temp",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)
### save results 0.05 sign from 34T1 vs 27T1

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_CBASS_34T1_vs_27T1_temp_and_gen_model.csv")

#####results contrast "32T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)
### save results 0.05 sign from 32T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_CBASS_32T1_vs_27T1_temp_and_gen_model.csv")


#####results contrast "29T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_CBASS_29T1_vs_27T1_temp_and_gen_model.csv")




#####results contrast "34T2","27T2" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic12947", intgroup="temp",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)
### save results 0.05 sign from 34T2 vs 27T2
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_CBASS_34T2_vs_27T2_temp_and_gen_model.csv")

#####results contrast "32T2","27T2" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)
### save results 0.05 sign from 32T2 vs 27T2
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_CBASS_32T2_vs_27T2_temp_and_gen_model.csv")


#####results contrast "29T2","27T2" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_CBASS_29T2_vs_27T2_temp_and_gen_model.csv")


#######
####### Stress vs recovery
#######
#####results contrast "34T2","34T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T2","34T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign 
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_CBASS_34T2_vs_34T1_temp_and_gen_model.csv")

#####results contrast "32T2","32T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T2","32T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T2 vs 32T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_CBASS_32T2_vs_32T1_temp_and_gen_model.csv")

#####results contrast "29T2","29T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T2","29T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 29T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_CBASS_29T2_vs_29T1_temp_and_gen_model.csv")

#####results contrast "27T2","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"27T2","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 27T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_CBASS_27T2_vs_27T1_temp_and_gen_model.csv")




#
#### 1  design only RSS Smicro ####
#
# get temp from header

a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]

a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]


temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # for CTL and AMF
#####get rid of library CTL8,Q11, G4, G12 for CTL and AMF
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[1:49109,-c(1:40)] #to remove lib
design2<-design[-c(1:40),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
#ddsSKtemprlog<-rlog(ddsSKtemp)

#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="temp")



#to do DEG which one with model and follow up on the simple model deseq2 and then more difficult model limma
#### reorder factor and model

resSKtemp<-DESeq(ddsSKtemp,modelMatrixType="expanded", betaPrior=T )

temp = t(sizeFactors(resSKtemp))
sizematrix<-matrix(data=temp, nrow=nrow(counts(ddsSKtemp)), ncol=ncol(temp), byrow=TRUE)
scaledcounts = counts(ddsSKtemp)/sizematrix
head(scaledcounts)
head(sizematrix)

table<-cbind(rowMeans(scaledcounts[,1:5]),
             rowMeans(scaledcounts[,6:10]),
             rowMeans(scaledcounts[,11:15]),
             rowMeans(scaledcounts[,16:20]),
             rowMeans(scaledcounts[,21:25]),
             rowMeans(scaledcounts[,26:30]),
             rowMeans(scaledcounts[,31:35]),
             rowMeans(scaledcounts[,36:40]))
colnames(table)<-c("CB-T1-27","CB-T1-29","CB-T1-32","CB-T1-34","CB-T2-27","CB-T2-29","CB-T2-32","CB-T2-34")

new<-scaledcounts[,1:2]
for(i in 1:dim(table)[1]){
  new[i,1]<-(sum(table[i,] > 5)/8)*100   #### keep only site that have at least 5reads in ~90%
}
new2<-subset(new,new[,1]>=80)
dim(new2)

#write.csv(table, "/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/05_New_filtering_following_Dan_advices_table/Smicro_RSS_table_rowMeans_per_treatments")

########## log fold change shrinkage to remove gene and have a better visualization
resultsNames(ddsSKtemp)
resSKtemp2<-resSKtemp[rownames(counts(resSKtemp)) %in% rownames(new2), ]
############


#####results contrast "32T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_RSS_32T1_vs_27T1_temp_and_gen_model.csv")


#####results contrast "29T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_RSS_29T1_vs_27T1_temp_and_gen_model.csv")


#####results contrast "32T2","27T2" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_RSS_32T2_vs_27T2_temp_and_gen_model.csv")


#####results contrast "29T2","27T2" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_RSS_29T2_vs_27T2_temp_and_gen_model.csv")

#######
####### stress vs recovery
#######
#####results contrast "32T2","32T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T2","32T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T2 vs 32T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_RSS_32T2_vs_32T1_temp_and_gen_model.csv")

#####results contrast "29T2","29T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T2","29T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 29T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_RSS_29T2_vs_29T1_temp_and_gen_model.csv")

#####results contrast "27T2","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"27T2","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 27T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_RSS_27T2_vs_27T1_temp_and_gen_model.csv")



###### for RSS Smicro, stricter filtering every comparison with 34T1 or 34T2, gene with only at least a mean of 5 reads in 34T1/34T2 are kept


temp = t(sizeFactors(resSKtemp))
sizematrix<-matrix(data=temp, nrow=nrow(counts(ddsSKtemp)), ncol=ncol(temp), byrow=TRUE)
scaledcounts = counts(ddsSKtemp)/sizematrix
head(scaledcounts)
head(sizematrix)

table<-cbind(rowMeans(scaledcounts[,1:5]),
             rowMeans(scaledcounts[,6:10]),
             rowMeans(scaledcounts[,11:15]),
             rowMeans(scaledcounts[,16:20]),
             rowMeans(scaledcounts[,21:25]),
             rowMeans(scaledcounts[,26:30]),
             rowMeans(scaledcounts[,31:35]),
             rowMeans(scaledcounts[,36:40]))
colnames(table)<-c("CB-T1-27","CB-T1-29","CB-T1-32","CB-T1-34","CB-T2-27","CB-T2-29","CB-T2-32","CB-T2-34")

new<-scaledcounts[,1:2]
for(i in 1:dim(table)[1]){
  new[i,1]<-(sum(table[i,] > 5)/8)*100   #### keep only site that have at least 5reads in 90%
}
new2<-subset(new,new[,1]>=100)
dim(new2)

resSKtemp3<-resSKtemp[rownames(counts(resSKtemp)) %in% rownames(new2), ]


#####results contrast "34T1","27T1" ########
resutlsSKtemp<-results(resSKtemp3,contrast = c('temp',"34T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 34T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_RSS_34T1_vs_27T1_temp_and_gen_model.csv")



#####results contrast "34T2","27T2" ########
resutlsSKtemp<-results(resSKtemp3,contrast = c('temp',"34T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 34T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_RSS_34T2_vs_27T2_temp_and_gen_model.csv")


#######
####### stress vs recovery
#######
#####results contrast "34T2","34T1" ####
resutlsSKtemp<-results(resSKtemp3,contrast = c('temp',"34T2","34T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Smicro_RSS_34T2_vs_34T1_temp_and_gen_model.csv")






#
#### 2  design only CBASS spistillata ####
#
# get temp from header

a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]

a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]

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

#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="temp")



#to do DEG which one with model and follow up on the simple model deseq2 and then more difficult model limma
#### reorder factor and model
#ddsCTLAMF$condition<-factor(ddsCTLAMF$condition, levels=c("CTL","AMF"))
resSKtemp<-DESeq(ddsSKtemp,modelMatrixType="expanded", betaPrior=T )

temp = t(sizeFactors(resSKtemp))
sizematrix<-matrix(data=temp, nrow=nrow(counts(ddsSKtemp)), ncol=ncol(temp), byrow=TRUE)
scaledcounts = counts(ddsSKtemp)/sizematrix
head(scaledcounts)
head(sizematrix)

table<-cbind(rowMeans(scaledcounts[,1:5]),
             rowMeans(scaledcounts[,6:10]),
             rowMeans(scaledcounts[,11:15]),
             rowMeans(scaledcounts[,16:20]),
             rowMeans(scaledcounts[,21:25]),
             rowMeans(scaledcounts[,26:30]),
             rowMeans(scaledcounts[,31:35]),
             rowMeans(scaledcounts[,36:40]))
colnames(table)<-c("CB-T1-27","CB-T1-29","CB-T1-32","CB-T1-34","CB-T2-27","CB-T2-29","CB-T2-32","CB-T2-34")

new<-scaledcounts[,1:2]
for(i in 1:dim(table)[1]){
  new[i,1]<-(sum(table[i,] > 1)/8)*100   #### keep only site that have at least 10reads in 90%
}
new2<-subset(new,new[,1]>=80)
dim(new2)

resSKtemp2<-resSKtemp[rownames(counts(resSKtemp)) %in% rownames(new2), ]
#write.csv(table, "/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/05_New_filtering_following_Dan_advices_table/Spist_CBASS_table_rowMeans_per_treatments")

########## log fold change shrinkage to remove gene and have a better visualization
resultsNames(ddsSKtemp)
counts(resSKtemp)
############


library("gtools")
annotGOSm<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME_reefgenomics/01_Coral/Stylophora_pistillata_genome_reefgenomics/spis_tabulated_annots_csv.csv")
annotGOSm2<-transform(annotGOSm, Query = colsplit(Query, split = " ", names = c('Name','Namelong')))
annotGOSm3<-as.data.frame(cbind(as.character(annotGOSm2$Query[,1]),as.character(annotGOSm2$GO.terms)))


par(mar=c(3,3,3,3))
#####results contrast "34T1","27T1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Spis4494", intgroup="temp",main="Spist CBASS Spis4494 heat shock protein Hsp-16.2")
###plot COPI 5 genes enriched in Spist CBASS-RSS common recovery, 1002 genes, 34t2-34t1
plotCounts(ddsSKtemp, gene="Spis3585", intgroup="temp",main="Spist CBASS Spis3585 Coatomer subunit alpha")
plotCounts(ddsSKtemp, gene="Spis14812.t1", intgroup="temp",main="Spist CBASS Spis14812.t1 Coatomer subunit zeta-1")
plotCounts(ddsSKtemp, gene="Spis20868", intgroup="temp",main="Spist CBASS Spis20868 Coatomer subunit beta'")
plotCounts(ddsSKtemp, gene="Spis5648", intgroup="temp",main="Spist CBASS Spis5648 Transmembrane emp24 domain-containing protein 3")
plotCounts(ddsSKtemp, gene="Spis20417", intgroup="temp",main="Spist CBASS Spis20417 Coatomer subunit epsilon")


which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 34T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_34T1_vs_27T1_temp_and_gen_model.csv")


#####results contrast "32T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_32T1_vs_27T1_temp_and_gen_model.csv")

#####results contrast "29T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_29T1_vs_27T1_temp_and_gen_model.csv")





#####results contrast "34T2","27T2" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Spis1000.t1", intgroup="temp",main="Spist CBASS Spis1000.t1 Common gene in recovery for CBASS and RSS for 29,32,34")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 34T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_34T2_vs_27T2_temp_and_gen_model.csv")


#####results contrast "32T2","27T2" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_32T2_vs_27T2_temp_and_gen_model.csv")

#####results contrast "29T2","27T2" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_29T2_vs_27T2_temp_and_gen_model.csv")




#######
####### Normal vs recovery
#######
#####results contrast "34T2","34T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T2","34T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_34T2_vs_34T1_temp_and_gen_model.csv")

#####results contrast "32T2","32T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T2","32T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Spis18245", intgroup="temp",main="Oxidative stress-responsive serine-rich protein 1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T2 vs 32T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_32T2_vs_32T1_temp_and_gen_model.csv")

#####results contrast "29T2","29T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T2","29T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 29T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_29T2_vs_29T1_temp_and_gen_model.csv")

#####results contrast "27T2","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"27T2","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 27T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_27T2_vs_27T1_temp_and_gen_model.csv")


#
#### 2  design only RSS spistillata ####
#
# get temp from header

a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]

a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]


temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # for CTL and AMF
#####get rid of library of the CBASS
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[49110:74878,-c(1:40)] #to remove lib
design2<-design[-c(1:40),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA

plotPCA(ddsSKtempvst, intgroup="temp")

#### reorder factor and model
resSKtemp<-DESeq(ddsSKtemp,modelMatrixType="expanded", betaPrior=T )

temp = t(sizeFactors(resSKtemp))
sizematrix<-matrix(data=temp, nrow=nrow(counts(ddsSKtemp)), ncol=ncol(temp), byrow=TRUE)
scaledcounts = counts(ddsSKtemp)/sizematrix
head(scaledcounts)
head(sizematrix)

table<-cbind(rowMeans(scaledcounts[,1:5]),
             rowMeans(scaledcounts[,6:10]),
             rowMeans(scaledcounts[,11:15]),
             rowMeans(scaledcounts[,16:20]),
             rowMeans(scaledcounts[,21:25]),
             rowMeans(scaledcounts[,26:30]),
             rowMeans(scaledcounts[,31:35]),
             rowMeans(scaledcounts[,36:40]))
colnames(table)<-c("CB-T1-27","CB-T1-29","CB-T1-32","CB-T1-34","CB-T2-27","CB-T2-29","CB-T2-32","CB-T2-34")

new<-scaledcounts[,1:2]
for(i in 1:dim(table)[1]){
  new[i,1]<-(sum(table[i,] > 1)/8)*100   #### keep only site that have at least a mean of 1reads
}
new2<-subset(new,new[,1]>=80) ###### in 100 of treatments
dim(new2)

#write.csv(table, "/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/05_New_filtering_following_Dan_advices_table/Spist_RSS_table_rowMeans_per_treatments")

resSKtemp2<-resSKtemp[rownames(counts(resSKtemp)) %in% rownames(new2), ]
########## log fold change shrinkage to remove gene and have a better visualization
resultsNames(ddsSKtemp)

#####results contrast "34T1","27T1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Spis1000.t1", intgroup="temp",main="Spist RSS Spis1000.t1 Common gene in recovery for CBASS and RSS for 29,32,34")

###plot COPI 5 genes enriched in Spist CBASS-RSS common recovery, 1002 genes, 34t2-34t1
plotCounts(ddsSKtemp, gene="Spis3585", intgroup="temp",main="Spist RSS Spis3585 Coatomer subunit alpha")
plotCounts(ddsSKtemp, gene="Spis14812.t1", intgroup="temp",main="Spist RSS Spis14812.t1 Coatomer subunit zeta-1")
plotCounts(ddsSKtemp, gene="Spis20868", intgroup="temp",main="Spist RSS Spis20868 Coatomer subunit beta'")
plotCounts(ddsSKtemp, gene="Spis5648", intgroup="temp",main="Spist RSS Spis5648 Transmembrane emp24 domain-containing protein 3")
plotCounts(ddsSKtemp, gene="Spis20417", intgroup="temp",main="Spist RSS Spis20417 Coatomer subunit epsilon")
### plot some gene of the UPR
plotCounts(ddsSKtemp, gene="Spis3351", intgroup="temp",main="Spist RSS Spis20417 Coatomer subunit epsilon")


which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 34T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_34T1_vs_27T1_temp_and_gen_model.csv")


#####results contrast "32T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_32T1_vs_27T1_temp_and_gen_model.csv")


#####results contrast "29T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_29T1_vs_27T1_temp_and_gen_model.csv")





#####results contrast "34T2","27T2" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Spis1000.t1", intgroup="temp",main="Spist RSS Spis1000.t1 Common gene in recovery for CBASS and RSS for 29,32,34")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 34T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_34T2_vs_27T2_temp_and_gen_model.csv")


#####results contrast "32T2","27T2" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_32T2_vs_27T2_temp_and_gen_model.csv")


#####results contrast "29T2","27T2" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T2","27T2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_29T2_vs_27T2_temp_and_gen_model.csv")

#######
####### Stress vs recovery
#######
#####results contrast "34T2","34T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T2","34T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_recovery_34T2_vs_34T1_temp_and_gen_model.csv")

#####results contrast "32T2","32T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T2","32T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 32T2
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_recovery_32T2_vs_32T1_temp_and_gen_model.csv")

#####results contrast "29T2","29T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T2","29T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 29T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_recovery_29T2_vs_29T1_temp_and_gen_model.csv")

#####results contrast "27T2","27T1" ####
resutlsSKtemp<-results(resSKtemp,contrast = c('temp',"27T2","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 27T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_recovery_27T2_vs_27T1_temp_and_gen_model.csv")


























