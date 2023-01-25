library("maftools")
#Open COSMIC cancer mutations file
cancer_mutations_data<-read.table("data/COSMIC_data/cancer_mutation_census_data/cmc_export.tsv", sep="\t", header=TRUE)

#Get maf files paths
maf_files<-list.files("data/vcf_data/", include.dirs = TRUE, full.names=TRUE) #gets names of all files in raw/ directory
maf_files<-maf_files[-1]

laml<-(read.maf(maf_files[2]))
#Filter Cosmic data based on oncogenes
ONCOGENE_list<-get_census_genes_list("oncogene")
cancer_mutations_data_filtered<-subset(cancer_mutations_data, cancer_mutations_data$GENE_NAME %in% ONCOGENE_list)


##Join mutation AA and gene name columns in COSMIC
cancer_mutations_data_filtered$joined_gene_aa<-paste0(cancer_mutations_data_filtered$GENE_NAME,":", cancer_mutations_data_filtered$Mutation.AA)


#Get TSG list
TSG_list<-get_census_genes_list("TSG")

#filter_maf<-function(maf_file){
  laml<-(read.maf(maf_file))
  
  #Filter maf data based on TSGs
  laml_tsg_filtered<-subsetMaf(laml, genes=TSG_list)@data
  
  
  #Filter MAF based on oncogene AA variations
  ##Join mutation AA and gene name columns in MAF
  maf_table<-laml@data
  maf_table$joined_gene_aa<-paste0(maf_table$Hugo_Symbol,":",maf_table$HGVSp_Short)
  
  
  #Filter MAF data based on cosmic mutation data
  laml_oncogenes_filtered<-subset(maf_table, maf_table$joined_gene_aa %in% cancer_mutations_data_filtered$joined_gene_aa)
  
  
#}








