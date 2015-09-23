setwd("C:/Users/millimanej/Desktop/JHU_DaSH/Alzheimers_Dash/")

library(slam)
library(gplots)
library(RColorBrewer)
load("biogrid2.Rdata")

A<-biogrid2[,c(1,2)] 
colnames(A)<-c("Gene1","Gene2")  
T<-table(A) 
Td<-dist(as.matrix(T))

setwd("C:/Users/millimanej/Desktop/Alzheimers_Dash/")

nhgri<-read.table("NHGRI.signGene.list",header=T,sep="\t", stringsAsFactors = F, quote = "")
mouse<-read.table("normHMD_HumanPhenotype.txt", header=F, quote="")

sum(colSums(table(nhgri))>4)
(levels(as.factor(nhgri$trait)))

t<-data.frame(ID=letters, mm=0, hs=0)
t[sample(1:26, 5),c(2:3)]<-1
head(mouse)
mouse<-mouse[,-c(1,8)]
mouse$V6<-gsub("MGI:", "", mouse$V6)
sum(!duplicated(mouse$V7))
colSums(table(mouse[,c("V2","V7")]))

cv<-read.table("clinvar_genes_phenos_updated_grouped.txt",header=T, sep="\t", quote="")
colSums(table(cv[,c(2,5)]))

colnames(cv)

mgi<-read.table("../MGI_PhenoGenoMP.rpt", quote="", comment.char = "", sep="\t")
mgi<-mgi[,c("V6","V4")]
colnames(mgi)<-c("MGI","Pheno")
mgi_hu<-merge(mgi, mouse, by.x = "MGI", by.y="V6", all.y=T,all.x=T, suffixes = c("MGI", "EI"))
head(mgi_hu)
head(mgi)
head(mouse)

load("phenotype_by_phenotype_fishers.Rda")
summary(-log10(as.vector(reshape[,-1])[!is.na(as.vector(reshape[,-1]))]))
hist(-log10(as.vector(reshape[,-1])[!is.na(as.vector(reshape[,-1]))]), breaks=1000)
reshape[1:10,1:10]
class(reshape[1,10])
rownames(reshape)<-reshape$Var1
reshape<-reshape[,-1]
reshape[is.na(reshape)]<-1
reshape.transform<--log10(reshape)
reshape.transform[reshape.transform==Inf]<-300

load("phenotype_by_phenotype_fishers_oddsratio.Rda")
hist(log2(as.vector(reshape2[,-1])[!is.na(as.vector(reshape[,-1]))]), breaks=100)
hist(ylim=c(0,500),(as.vector(reshape2[,-1])[!is.na(as.vector(reshape[,-1]))]), breaks=100)
reshape[1:10,1:10]
class(reshape[1,10])
rownames(reshape)<-reshape$Var1
reshape<-reshape[,-1]

png("Fishers_heatmap.png", height = 6000, width = 12000, units="px")
heatmap.2(as.matrix(reshape), trace="none",
          dendrogram ="none", Rowv = F, Colv = F, labRow = "", labCol="")
dev.off()

png("Fishers_heatmap2.png", height = 6000, width = 12000, units="px")
par(cex=2.5)
heatmap.2(as.matrix(reshape.transform), trace="none",
          dendrogram ="none", breaks=seq(0,100,by=25), Rowv = F, Colv = F, labRow = "", labCol="", col=colorRampPalette(c("white","red"))((n=4)))
dev.off()

png("Fishers_test_distribution.png", height = 800, width = 1200, units="px")
par(cex=2)
hist(as.vector(as.matrix(reshape.transform)), breaks=100, ylim=c(0,10000), xlab="-log10(p-value)", main="Mouse phenotypes to human disease test results")
dev.off()

png("Fishers_pvalue_QQ.png", height=800, width = 1200, units="px")
qqnorm(as.vector(as.matrix(reshape)), pch=19)
dev.off()

