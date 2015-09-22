library(readr)

# read in the MGI Alzheimer disease allele data
mgi_allele <- read_csv("MGI_alzheimer.csv")

# change column name to replace spaces in name
names(mgi_allele) <- make.names(names(mgi_allele), unique=TRUE)

# get only targeted mutation
mgi_target <- subset(mgi_allele, Allele.Type == "Targeted")
# get allele information
mgi_target[,2]

# get unique gene name from allele name
mgi_gene <- unique(sapply(strsplit(mgi_target[,2], "<"), "[", 1))



library(biomaRt)

# set Biomart human gene set
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# set Biomart mouse gene set
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# map mouse gene symbol to human gene symbol
hs_gene <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mgi_gene ,mart = mouse, attributesL = c("hgnc_symbol","chromosome_name", "start_position"), martL = human, uniqueRows=T)
