get_filtered_genes <- function(min_samples=5, top_genes_expression=15000, top_genes_coefficient=10000) {
  if (!exists("RPKM_data")){
    RPKM_data<-get_RPKM_normalised_data()
  }

  
  
  #FILTER GENES EXPRESSED IN LESS THAN X samples
  filtered_expression_df_1<-RPKM_data[rowSums(RPKM_data > 0)>=min_samples,] 
  
  #KEEP TOP Y genes based on gene expression
  MEAN_df<-data.matrix(apply(filtered_expression_df_1,1,mean))
  MEAN_df<-data.matrix(sort(MEAN_df[,1], decreasing=TRUE),)
  top_N_genes<-data.matrix(head(rownames(MEAN_df),top_genes_expression))
  filtered_expression_df_2<-subset(filtered_expression_df_1, rownames(filtered_expression_df_1) %in% top_N_genes)
  
  #Filter top Z most impactful genes from coefficient of variation DF.
  coefficient_variation_df<-(apply(filtered_expression_df_2,1,sd)/rowMeans(filtered_expression_df_2))
  top_Y_genes<-names(head(sort(abs(coefficient_variation_df),decreasing=TRUE),n=top_genes_coefficient))
  
  return(top_Y_genes)
  
}

get_filtered_rpkm <- function(min_samples=5, top_genes_expression=15000, top_genes_coefficient=10000) {
  if (!exists("RPKM_data")){
    RPKM_data<-get_RPKM_normalised_data()
  }
  
  
  #FILTER GENES EXPRESSED IN LESS THAN X samples
  filtered_expression_df_1<-RPKM_data[rowSums(RPKM_data > 0)>=min_samples,] 
  
  #KEEP TOP Y genes based on gene expression
  MEAN_df<-data.matrix(apply(filtered_expression_df_1,1,mean))
  MEAN_df<-data.matrix(sort(MEAN_df[,1], decreasing=TRUE),)
  top_N_genes<-data.matrix(head(rownames(MEAN_df),top_genes_expression))
  filtered_expression_df_2<-subset(filtered_expression_df_1, rownames(filtered_expression_df_1) %in% top_N_genes)
  
  #Filter top Z most impactful genes from coefficient of variation DF.
  coefficient_variation_df<-(apply(filtered_expression_df_2,1,sd)/rowMeans(filtered_expression_df_2))
  top_Y_genes<-names(head(sort(abs(coefficient_variation_df),decreasing=TRUE),n=top_genes_coefficient))
  
  #Subset RPKM values
  filtered_RPKM_values<-subset(RPKM_data, rownames(RPKM_data) %in% top_Y_genes)
  return(filtered_RPKM_values)
  
}
