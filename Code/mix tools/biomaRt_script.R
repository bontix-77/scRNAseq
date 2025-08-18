
library(biomaRt)
# to see al the datasets availble you can use listEnsembl() listDatasets(), to see the attributes contained in a dataset listAttributes(mart_obj)
#to find an attibute containing a substring inside a mart_obj use: grep(substring,mart_bobj@attributes. It returns the positio indexes.

mart  <- useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")
see <- getBM(attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name"),
filters = "external_gene_name",
values = "NFKB1  ",
mart=mart)

prova <- matrix("",ncol = 2, nrow=5)
prova[,1] <- "NFKB"
prova[,2] <- sapply(prova[,1],function(x){

mart  <- useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")
see <- getBM(attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name"),
filters = "external_gene_name",
values = x,
mart=mart)
return (see[1])


})
