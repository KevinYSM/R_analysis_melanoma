get_RPKM_normalised_data<-function(){
  gene_counts_df<-open_gene_counts()
  gene_lengths<-get_gene_lengths()
  gene_names<-rownames(gene_counts_df)
  num_rows<-length(rownames(gene_counts_df))
  num_cols<-length(colnames(gene_counts_df))
  
  normalised_gene_counts_df<-gene_counts_df
  #need to get counts in each row
  
  
  #1.  Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
  perMillionScalingFactor<-colSums(gene_counts_df)/1000000
  
  #2. Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
  
  for (i in (1:num_cols)){
    normalised_gene_counts_df[,i]<-gene_counts_df[,i]/perMillionScalingFactor[i]
  }
  
  
  #3. Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM."
  for (j in (1:num_rows)){
    gene<-gene_names[j]
    longest_isoform_length=gene_lengths[gene,][[3]]
    normalised_gene_counts_df[j,]<-normalised_gene_counts_df[j,]/(longest_isoform_length/1000)
    
  }
  
  #Remove N/A rows
  normalised_gene_counts_df<-na.omit(normalised_gene_counts_df)
  return(normalised_gene_counts_df)
}



