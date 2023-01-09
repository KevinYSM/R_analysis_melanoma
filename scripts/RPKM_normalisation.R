
gene_names<-rownames(gene_counts_df)
num_rows<-length(rownames(gene_counts_df))
for (i in (1:num_rows)){
  gene<-gene_names[i]
  mean_length=gene_lengths[gene,][[1]]
  
  gene_counts_df[i,]<-gene_counts_df[i,]*mean_length
 
}



