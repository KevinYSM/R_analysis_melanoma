#data from https://cancer.sanger.ac.uk/census
#Cancer Gene Census data. Genes that have been causally implicated in cancer
get_oncogenes_count_df <- function() {
  oncogene_data<-read.csv("data/COSMIC_data/cancer_gene_census_data/Census_allSun Jan 15 04_37_51 2023.csv")
  
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

get_census_genes_list <- function(gene_role) {
  #Gene_role = "oncogene" or "TSG"
  #Filters for TSG, Oncogene, and Tier 1 cancer census genes
  oncogene_data<-read.csv("data/COSMIC_data/cancer_gene_census_data/Census_allSun Jan 15 04_37_51 2023.csv")
  num_cols<-length(colnames(oncogene_data))
  num_rows<-length(rownames(oncogene_data))
  oncogenes<-list()
  
  #Filter by:
  # 1. Tier ==1
  # 2. Include if Role.in.cancer includes  "TSG" or "oncogene"
  
  
  #1 
  filter_1<-oncogene_data[oncogene_data$Tier=='1',]
  
  #2
  filter_2<-filter_1[grepl(gene_role,filter_1$Role.in.Cancer),]
 

  #Perhaps as a separate analysis/additional information:
  # 1. Tumor.types.somatic includes "melanoma"
  oncogenes<-filter_2[,1]
  return(oncogenes)
}
filtered_census_genes<-as.list(get_oncogenes_list())

