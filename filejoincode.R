library(sqldf)

Global<-read.delim("normHMD_HumanPhenotype.txt",header=FALSE )
Disease<-read.delim("NHGRI.signGene.list")

Global<-Global[,c(2:7)]
colnames(Global)<-c("HGS","HGN","MGN","MGS","M","Phenotype")

Mergetable<-sqldf('select a.*,b.* from Global a join Disease b on a.HGN=b.gene')

D<-unique(Disease[,2])
G<-unique(Global[,6])
x<-expand.grid(D[!is.na(D)],G[!is.na(G)])
x1<-sqldf('select a.Var1, a.Var2, count(b.HGS) as Overlap 
        from x a inner join Mergetable b on a.Var1=b.trait and a.Var2=b.Phenotype
          group by a.Var1, a.Var2')

x2<-sqldf('select a.Var1, a.Var2, Overlap, count(b.HGS) as DiseaseOnly 
        from x1 a inner join Mergetable b on a.Var1=b.trait where a.Var2<>b.Phenotype
          group by a.Var1, a.Var2, a.Overlap')

x3<-sqldf('select a.Var1, a.Var2, Overlap, DiseaseOnly, count(b.HGS) as PhenotypeOnly 
        from x2 a inner join Mergetable b on a.Var2=b.Phenotype where a.Var1<>b.trait
          group by a.Var1, a.Var2, a.Overlap, a.DiseaseOnly')