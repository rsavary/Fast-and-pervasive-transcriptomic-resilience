#### Spist Venne diagramm following mail Dan, venne diagramm for the same change in temp across all experiement CBASS, RSS t1 and t3


####### Venn diagram for GENE DEG between recovery treatments
library(VennDiagram)

####make venn plot with list of  GENE deg for 27 vs 29.5 in Spist across CBASS,t1,T2,RSS, t1,t2
#UP
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/")

SCt1<-read.csv("Spist_CBASS_29T1_vs_27T1_temp_and_gen_model.csv",header=T)
SCt2<-read.csv("Spist_CBASS_29T3_vs_27T3_temp_and_gen_model.csv",header=T)
SRt1<-read.csv("Spist_RSS_29T1_vs_27T1_temp_and_gen_model.csv",header=T)
SRt2<-read.csv("Spist_RSS_29T3_vs_27T3_temp_and_gen_model.csv",header=T)

SCt1<-SCt1[SCt1$log2FoldChange>0,]
SCt2<-SCt2[SCt2$log2FoldChange>0,]
SRt1<-SRt1[SRt1$log2FoldChange>0,]
SRt2<-SRt2[SRt2$log2FoldChange>0,]

SCt1<-SCt1$X
SCt2<-SCt2$X
SRt1<-SRt1$X
SRt2<-SRt2$X
par(mar=c(0,0,0,0),mfrow=c(1,1),bty="n")
plot.new()
venn.plot<-venn.diagram(list(SCt1 ,SRt2, SCt2,SRt1), height = 100, width = 300,NULL, fill=c("green", "blue", "yellow","red"), alpha=c(0.2, 0.2,0.2,0.2), cex = 2, cat.fontface=4, category.names=c("Spist CBASS T1", "Spist RSS T2","Spist CBASS T2","Spist RSS T1"), main="Gene DEG for 27 vs 29.5°C")
grid.draw(venn.plot)





###DOWN
SCt1<-read.csv("Spist_CBASS_29T1_vs_27T1_temp_and_gen_model.csv",header=T)
SCt2<-read.csv("Spist_CBASS_29T3_vs_27T3_temp_and_gen_model.csv",header=T)
SRt1<-read.csv("Spist_RSS_29T1_vs_27T1_temp_and_gen_model.csv",header=T)
SRt2<-read.csv("Spist_RSS_29T3_vs_27T3_temp_and_gen_model.csv",header=T)

SCt1<-SCt1[SCt1$log2FoldChange<0,]
SCt2<-SCt2[SCt2$log2FoldChange<0,]
SRt1<-SRt1[SRt1$log2FoldChange<0,]
SRt2<-SRt2[SRt2$log2FoldChange<0,]


SCt1<-SCt1$X
SCt2<-SCt2$X
SRt1<-SRt1$X
SRt2<-SRt2$X
par(mar=c(0,0,0,0),mfrow=c(1,1),bty="n")
plot.new()
venn.plot<-venn.diagram(list(SCt1 ,SRt2, SCt2,SRt1), height = 100, width = 300,NULL, fill=c("green", "blue", "yellow","red"), alpha=c(0.2, 0.2,0.2,0.2), cex = 2, cat.fontface=4, category.names=c("Spist CBASS T1", "Spist RSS T2","Spist CBASS T2","Spist RSS T1"), main="Gene DEG for 27 vs 29.5°C")
grid.draw(venn.plot)


####make venn plot with list of  GENE deg for 27 vs 32 in Spist across CBASS,t1,T2,RSS, t1,t2
#UP
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/")

SCt1<-read.csv("Spist_CBASS_32T1_vs_27T1_temp_and_gen_model.csv",header=T)
SCt2<-read.csv("Spist_CBASS_32T3_vs_27T3_temp_and_gen_model.csv",header=T)
SRt1<-read.csv("Spist_RSS_32T1_vs_27T1_temp_and_gen_model.csv",header=T)
SRt2<-read.csv("Spist_RSS_32T3_vs_27T3_temp_and_gen_model.csv",header=T)

SCt1<-SCt1[SCt1$log2FoldChange>0,]
SCt2<-SCt2[SCt2$log2FoldChange>0,]
SRt1<-SRt1[SRt1$log2FoldChange>0,]
SRt2<-SRt2[SRt2$log2FoldChange>0,]

SCt1<-SCt1$X
SCt2<-SCt2$X
SRt1<-SRt1$X
SRt2<-SRt2$X
par(mar=c(0,0,0,0),mfrow=c(1,1),bty="n")
plot.new()
venn.plot<-venn.diagram(list(SCt1 ,SRt2, SCt2,SRt1), height = 100, width = 300,NULL, fill=c("green", "blue", "yellow","red"), alpha=c(0.2, 0.2,0.2,0.2), cex = 2, cat.fontface=4, category.names=c("Spist CBASS T1", "Spist RSS T2","Spist CBASS T2","Spist RSS T1"), main="Gene DEG for 27 vs 32°C")
grid.draw(venn.plot)

####make venn plot with list of  GENE deg for 27 vs 32 in Spist across CBASS,t1,T2,RSS, t1,t2
#DOWN
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/")

SCt1<-read.csv("Spist_CBASS_32T1_vs_27T1_temp_and_gen_model.csv",header=T)
SCt2<-read.csv("Spist_CBASS_32T3_vs_27T3_temp_and_gen_model.csv",header=T)
SRt1<-read.csv("Spist_RSS_32T1_vs_27T1_temp_and_gen_model.csv",header=T)
SRt2<-read.csv("Spist_RSS_32T3_vs_27T3_temp_and_gen_model.csv",header=T)

SCt1<-SCt1[SCt1$log2FoldChange<0,]
SCt2<-SCt2[SCt2$log2FoldChange<0,]
SRt1<-SRt1[SRt1$log2FoldChange<0,]
SRt2<-SRt2[SRt2$log2FoldChange<0,]

SCt1<-SCt1$X
SCt2<-SCt2$X
SRt1<-SRt1$X
SRt2<-SRt2$X
par(mar=c(0,0,0,0),mfrow=c(1,1),bty="n")
plot.new()
venn.plot<-venn.diagram(list(SCt1 ,SRt2, SCt2,SRt1), height = 100, width = 300,NULL, fill=c("green", "blue", "yellow","red"), alpha=c(0.2, 0.2,0.2,0.2), cex = 2, cat.fontface=4, category.names=c("Spist CBASS T1", "Spist RSS T2","Spist CBASS T2","Spist RSS T1"), main="Gene DEG for 27 vs 32°C")
grid.draw(venn.plot)


####make venn plot with list of  GENE deg for 27 vs 34.5 in Spist across CBASS,t1,T2,RSS, t1,t2
#UP
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/")

SCt1<-read.csv("Spist_CBASS_34T1_vs_27T1_temp_and_gen_model.csv",header=T)
SCt2<-read.csv("Spist_CBASS_34T3_vs_27T3_temp_and_gen_model.csv",header=T)
SRt1<-read.csv("Spist_RSS_34T1_vs_27T1_temp_and_gen_model.csv",header=T)
SRt2<-read.csv("Spist_RSS_34T3_vs_27T3_temp_and_gen_model.csv",header=T)

SCt1<-SCt1[SCt1$log2FoldChange>0,]
SCt2<-SCt2[SCt2$log2FoldChange>0,]
SRt1<-SRt1[SRt1$log2FoldChange>0,]
SRt2<-SRt2[SRt2$log2FoldChange>0,]

SCt1<-SCt1$X
SCt2<-SCt2$X
SRt1<-SRt1$X
SRt2<-SRt2$X

par(mar=c(0,0,0,0),mfrow=c(1,1),bty="n")
plot.new()
venn.plot<-venn.diagram(list(SCt1 ,SRt2, SCt2,SRt1), height = 100, width = 300,NULL, fill=c("brown1", "brown3", "brown2","brown4"), alpha=c(0.2, 0.2,0.2,0.2), cex = 2, cat.fontface=4, category.names=c("Spist CBASS T1", "Spist RSS T2","Spist CBASS T2","Spist RSS T1"), main="Gene DEG for 27 vs 34.5°C")
grid.draw(venn.plot)

#DOWN
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/")

SCt1<-read.csv("Spist_CBASS_34T1_vs_27T1_temp_and_gen_model.csv",header=T)
SCt2<-read.csv("Spist_CBASS_34T3_vs_27T3_temp_and_gen_model.csv",header=T)
SRt1<-read.csv("Spist_RSS_34T1_vs_27T1_temp_and_gen_model.csv",header=T)
SRt2<-read.csv("Spist_RSS_34T3_vs_27T3_temp_and_gen_model.csv",header=T)

SCt1<-SCt1[SCt1$log2FoldChange<0,]
SCt2<-SCt2[SCt2$log2FoldChange<0,]
SRt1<-SRt1[SRt1$log2FoldChange<0,]
SRt2<-SRt2[SRt2$log2FoldChange<0,]

SCt1<-SCt1$X
SCt2<-SCt2$X
SRt1<-SRt1$X
SRt2<-SRt2$X

par(mar=c(0,0,0,0),mfrow=c(1,1),bty="n")
plot.new()
venn.plot<-venn.diagram(list(SCt1 ,SRt2, SCt2,SRt1), height = 100, width = 300,NULL, fill=c("brown1", "brown3", "brown2","brown4"), alpha=c(0.2, 0.2,0.2,0.2), cex = 2, cat.fontface=4, category.names=c("Spist CBASS T1", "Spist RSS T2","Spist CBASS T2","Spist RSS T1"), main="Gene DEG for 27 vs 34.5°C")
grid.draw(venn.plot)

