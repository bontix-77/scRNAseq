
library(biomaRt)
# to see al the datasets availble you can use listEnsembl() listDatasets(), to see the attributes contained in a dataset listAttributes(mart_obj)
#to find an attibute containing a substring inside a mart_obj use: grep(substring,mart_bobj@attributes. It returns the positio indexes.

mart  <- useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")
see <- getBM(attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name"),
filters = "external_gene_name",
values = "NFKB1  ",
mart=mart)
#example of a for loop to convert from gene symbols to ensembl ID
prova <- matrix("",ncol = 2, nrow=5)
prova[,1] <- c("NFKB1","ACTB","NFKB2","PPARG","PPARA")
for (i in 1:length(prova[,1])){
     gene <- prova[i,1]
     out <- getBM(attributes = "ensembl_gene_id",
     filters = "external_gene_name",
     values =gene,
     mart=mart)
     prova[i,2] <- out$ensembl_gene_id
}



# to convert from Enseml ID to gene symbols

prova <- matrix("",ncol = 2, nrow=5)
prova[,1] <- c("ENSG00000109320" ,"ENSG00000075624" ,"ENSG00000077150", "ENSG00000132170" ,"ENSG00000186951")
for (i in 1:length(prova[,1])){
     gene <- prova[i,1]
     out <- getBM(attributes = "external_gene_name",
     filters = "ensembl_gene_id",
     values =gene,
     mart=mart)
     prova[i,2] <- out$external_gene_name
}
