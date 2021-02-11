######## Fig 3 all 

library("ggplot2")



#### barplot A

a<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/02_CBASS_vs_RSS_T1-T2/Percentage_gene_change_with_the_experiements.csv",header=T)


g1<-ggplot(data=a, aes(x=temp, y=per, fill=exp))  + theme(plot.title = element_text(size = 8, face = "bold")) +
  xlab("") + xlab("Temperature comparisons (°C)") + ylab("% of DEGs in genomes")   + geom_bar(stat="identity", position=position_dodge()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_manual("Experiments", values = c("green4","green3","green2","green1","brown4","brown3","brown2","brown1")) + ylim(0, 40)

##### 

#####
a<-read.delim("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/02_CBASS_vs_RSS_T1-T2/01_Paper_Versions/04_Version/01_Fig/Plot_resilience/Resilience_percentage_replacement.txt",header=T)

a<-a[1:12,]
g2<-ggplot(data=a, aes(x=temp, y=resilient., group=paste(a$species,a$exp,sep="_"))) + geom_line(aes(color=paste(a$species,a$exp,sep="_"))) + geom_point(aes(color=paste(a$species,a$exp,sep="_"))) + labs(x="Temperature comparisons (°C)", y = "% of DEGs at T1 returned to baseline expression at T2") + geom_text(data=a,aes(label=paste("(",genenb1,")"),hjust=-0.3, vjust=0),size=3) +scale_color_manual(values=c("#008b00","#9aff9a","#8b2323","#ff4040"))

####
a<-read.delim("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/02_CBASS_vs_RSS_T1-T2/01_Paper_Versions/04_Version/01_Fig/Plot_resilience/Common_acute_chronic_percentage_replacement.txt",header=T)

g3<-ggplot(data=a, aes(x=temp, y=Per, group=Exp)) + ylim(0, 17.5) + geom_line(aes(color=Exp)) + geom_point(aes(color=Exp)) + labs(x="Temperature comparisons (°C)", y = "% of shared DEGs between the short- and long-term heat stress") + scale_color_manual(values=c("#008b00","#9aff9a","#8b2323","#ff4040"))

####

###### nb gene only found at T2 Fig. S5

a<-read.delim("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/02_CBASS_vs_RSS_T1-T2/01_Paper_Versions/04_Version/01_Fig/Plot_resilience/Resilience_percentage_replacement.txt",header=T)

a<-a[1:12,]

p<-ggplot(data=a, aes(x=temp, y=recoverygenes, group=paste(a$species,a$exp,sep="_"))) + geom_line(aes(color=paste(a$species,a$exp,sep="_"))) + geom_point(aes(color=paste(a$species,a$exp,sep="_"))) + labs(x="Temperature comparisons (°C)", y = "nb recovery genes (only affected at T2)") + geom_text(data=a,aes(label=paste("(",recovery.,")"),hjust=-0.3, vjust=0),size=3)

p+scale_color_manual(values=c("#008b00","#9aff9a","#8b2323","#ff4040"))





