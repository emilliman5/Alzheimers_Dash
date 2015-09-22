setwd("/Users/doehun/Desktop/Alzheimers_Dash")
install.packages("reshape2")
require(reshape2)

#===============================================

Clinvar<-read.delim("clinvar_genes_phenos_updated_grouped.txt")
Global<-read.delim("mgi_Phenotype.txt" ,header=FALSE )
Disease<-read.delim("NHGRI.signGene.list")

Clinvar<-Clinvar[,c(2,5)]
colnames(Clinvar)<-colnames(Disease)
Disease<-rbind(Disease,Clinvar)

data<-load("test_Contingency.Rda",envir = parent.frame(), verbose = FALSE)
data1<-read.delim("mgi_Phenotype.txt", header=FALSE) 
length(unique(data1))

size<-length(unique(c(Disease$gene, Global$V2)))
x3[,"Neither"] <- list(rep(size, length(x3$Var1)))
x3$Neither<-x3$Neither - x3$Overlap - x3$DiseaseOnly - x3$PhenotypeOnly

x3.1<-x3[!x3$Neither<0,]

pvalue<-apply(as.matrix(x3.1[,3:6]), 1, function(x) 
  fisher.test(matrix(unlist(x), ncol=2))$p.value)
estimate<-apply(as.matrix(x3.1[,3:6]), 1, function(x) 
  fisher.test(matrix(unlist(x), ncol=2))$estimate)

x3.1[,"pvalue"]<-pvalue
x3.1[,"oddsratio"]<-estimate

head(x3.1)

#table<-table(x3.1[,c("Var1", "Var2", "pvalue")])
#table2<-table(x3.1[,c("Var1", "Var2", "oddsratio")])

reshape<-reshape(x3.1[,c("Var1","Var2","pvalue")], idvar = "Var1", timevar = "Var2", direction = "wide")
reshape2<-reshape(x3.1[,c("Var1","Var2","oddsratio")], idvar = "Var1", timevar = "Var2", direction = "wide")


save(reshape,file = "phenotype_by_phenotype_fishers_pvalue.Rda")
save(reshape2,file = "phenotype_by_phenotype_fishers_oddsratio.Rda")

heatmap(reshape)

class(reshape[2,2])
hist(-log10(as.vector(reshape[,-1])[!is.na(as.vector(reshape[,1]))]), breaks=100)
heatmap.2(-log10(as.matrix(reshape2), trace="none"))