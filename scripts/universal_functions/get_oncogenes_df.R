get_oncogenes_df <- function() {
  oncogene_data<-read.csv("scripts/melanoma_genes_heatmap/Census_allSun Jan 15 04_37_51 2023.csv")
  
  num_cols<-length(colnames(oncogene_data))
  num_rows<-length(rownames(oncogene_data))
  oncogenes<-list()
  for (i in (1:num_rows)){
    tumour_types<-oncogene_data[i,10]
    if (grepl("melanoma",tumour_types,fixed = TRUE)){
      oncogenes<-append(oncogenes,oncogene_data[i,1])
    }
  }
  normalised_gene_counts_df<-get_RPKM_normalised_data()
  
  filtered_gene_counts_df<-subset(normalised_gene_counts_df, rownames(normalised_gene_counts_df) %in% oncogenes)
  return(filtered_gene_counts_df)
}

