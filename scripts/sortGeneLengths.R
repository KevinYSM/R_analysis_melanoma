gene_lengths<-read.table(("gene_lengths/human_gene_lengths.txt"), sep="", header=TRUE) #reads gene counts files and converts to table

row.names(gene_lengths)<-gene_lengths[[1]]
gene_lengths[[1]]<-NULL

gene_lengths<-as.data.frame.matrix(gene_lengths)
gene_lengths<-gene_lengths[order(row.names(gene_lengths)),]
rm(list=ls())
