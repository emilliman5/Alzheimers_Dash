library(tm)
library(stringi)
library(sqldf)

load("pval_table.Rda")

text<-read.delim("VOC_MammalianPhenotype.rpt",header=FALSE)
text2<-removeWords(removePunctuation(removeNumbers(stri_trans_tolower(text[,2]))),stopwords("SMART"))
text[,3]<-removeWords(removePunctuation(removeNumbers(stri_trans_tolower(text[,3]))),stopwords("SMART"))


t2<-unlist(strsplit(text2," "))
t3<-unlist(strsplit(text3," "))

t2a<-as.data.frame(table(t2))
t2b<-t2a[t2a$Freq>1,] 

t3a<-as.data.frame(table(t3))
t3b<-t3a[t3a$Freq>1,] 

save(t2b,file="WordListMammalianPhenotype.Rdata")
save(t3b,file="WordListObservableMorpPhysBehave.Rdata")

tests<-read.delim("Term.txt",header=FALSE)
test<-tests[1,1]
text$ID<-test
test1<-grepl(test,text[,3])
test2<-text[test1,c(1,4)]
for (i in 2:20){
test<-tests[i,1]
text$ID<-test
test1<-grepl(test,text[,3])
test3<-text[test1,c(1,4)]
test2<-rbind(test2,test3)
}

x4<-x3.1

subX<-sqldf('select b.*, a.ID from test2 a left join x4 b on a.V1=b.Var2')

save(x4,file="x4.Rdata")
save(subX,file="subX.Rdata")


rownames(tests)<-tests$V1

x<-"ALL"
tests2<-c(x,as.character(tests[,1]))

colnames(subX)<-c("Human Disease","M. Phenotype","Overlap","DiseaseOnly","PhenotypeOnly","Neither","pvalue","oddsratio")




