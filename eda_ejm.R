library(slam)

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
hist(-log10(as.vector(reshape[,-1])[!is.na(as.vector(reshape[,-1]))]), breaks=1000)
reshape[1:10,1:10]
class(reshape[1,10])
rownames(reshape)<-reshape$Var1
reshape<-reshape[,-1]
sparsity(as.matrix(reshape))
