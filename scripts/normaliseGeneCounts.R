#Normalise gene_counts_df

scaled_gene_counts_df=scale (gene_counts_df,center=FALSE, scale=colSums(gene_counts_df)/1000000)
