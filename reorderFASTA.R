setwd(snakemake@params[["wd"]])

library(Biostrings)

fasta<-readAAStringSet(unlist(snakemake@input))
fasta<-fasta[!duplicated(seq(fasta))]

fasta.names<-strsplit(names(fasta),";")
fasta.names<-unlist(lapply(fasta.names,function(x){x[2]}))
fasta.names<-paste0("pg|",fasta.names,"|")

names(fasta)<-fasta.names

writeXStringSet(fasta,unlist(snakemake@output))