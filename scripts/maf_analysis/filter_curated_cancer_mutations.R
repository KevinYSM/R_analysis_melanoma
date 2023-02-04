library("maftools")
#Open COSMIC cancer mutations file

if (!exists("cancer_mutations_data")){
  cancer_mutations_data<-read.table("data/COSMIC_data/cancer_mutation_census_data/cmc_export.tsv", sep="\t", header=TRUE)
}

#Get maf files paths
maf_files<-list.files("data/vcf_data/", include.dirs = TRUE, full.names=TRUE) #gets names of all files in raw/ directory
maf_files<-maf_files[-1]


#Get TSG, Oncogene list and combined list
TSG_list<-get_census_genes_list("TSG")
ONCOGENE_list<-get_census_genes_list("oncogene")
ALL_GENE_list<-sort(unique(c(TSG_list,ONCOGENE_list)), decreasing=TRUE)

#Filter Cosmic data based on oncogenes
filtered_cosmic_oncogenes<-subset(cancer_mutations_data, cancer_mutations_data$GENE_NAME %in% ONCOGENE_list)


##Join mutation AA and gene name columns in COSMIC
filtered_cosmic_oncogenes$joined_gene_aa<-paste0(filtered_cosmic_oncogenes$GENE_NAME,":", filtered_cosmic_oncogenes$Mutation.AA)






update_gene_sample_counts<-function(maf_list=maf_files){
  #Count mutated genes found across samples
  gene_sample_counts<-data.frame(matrix(nrow=length(ALL_GENE_list),ncol=3))
  colnames(gene_sample_counts)<-c("Hugo_Symbol","Samples Present", "Samples List")
  gene_sample_counts[,1]<-ALL_GENE_list
  gene_sample_counts[,2]<-0
  
  maf_index<-0
  for (maf in maf_list){
    filtered_maf<-filter_maf(maf)
    print(maf)

    for (index in (1:length(filtered_maf))){
      gene<-filtered_maf$Hugo_Symbol[index]
      sample_name<-substr(basename(maf),1,6)
    
      if (gene %in% gene_sample_counts$Hugo_Symbol)
        {
        if (is.na(gene_sample_counts[index,3])){
          gene_sample_counts[which(gene_sample_counts$Hugo_Symbol==gene),3]<-sample_name
          gene_sample_counts[which(gene_sample_counts$Hugo_Symbol==gene),2]<-1
        }
        else{
          gene_sample_counts[which(gene_sample_counts$Hugo_Symbol==gene),3]<-paste(gene_sample_counts[which(gene_sample_counts$Hugo_Symbol==gene),3], sample_name, sep=",")
          gene_sample_counts[which(gene_sample_counts$Hugo_Symbol==gene),2]<-gene_sample_counts[which(gene_sample_counts$Hugo_Symbol==gene),2]+1
        }
      }
    }
  }
  return(gene_sample_counts)
}
gene_sample_counts<-update_gene_sample_counts(maf_files)

#Returns a filtered maf for a single sample
filter_maf<-function(maf_file_path){
  laml<-(read.maf(maf_file_path))
  raw_maf<-laml@data
  raw_maf$joined_gene_aa<-paste0(raw_maf$Hugo_Symbol,":",raw_maf$HGVSp_Short)
  
  #Filter MAF data based on TSGs
  maf_tsg_filtered<-subsetMaf(laml, genes=TSG_list)@data
  maf_tsg_filtered$joined_gene_aa<-paste0(maf_tsg_filtered$Hugo_Symbol,":",maf_tsg_filtered$HGVSp_Short)
 
  ##remove TSG columns with rs in $dpSNP_rs
  maf_tsg_filtered<-maf_tsg_filtered[! grep("rs", maf_tsg_filtered$dbSNP_RS)]
  
  #Filter MAF based on oncogene AA variations
  ##Join mutation AA and gene name columns in MAF
  laml_oncogenes_filtered<-laml@data
  laml_oncogenes_filtered$joined_gene_aa<-paste0(laml_oncogenes_filtered$Hugo_Symbol,":",laml_oncogenes_filtered$HGVSp_Short)
  
  
  #Filter MAF data based on cosmic mutation data
  maf_oncogenes_filtered<-subset(raw_maf, raw_maf$joined_gene_aa %in% filtered_cosmic_oncogenes$joined_gene_aa)

  # Join maf_tsg_filtered with maf_oncogenes_filtered
  joined_filtered_maf<-unique(rbind(maf_oncogenes_filtered, maf_tsg_filtered))
  oncogene_data<-read.csv("data/COSMIC_data/cancer_gene_census_data/Census_allSun Jan 15 04_37_51 2023.csv")
  
  #add column for role in cancer
  final_filtered_maf<-merge(joined_filtered_maf, unique(oncogene_data[,c("Gene.Symbol","Role.in.Cancer")]),by.x="Hugo_Symbol", by.y="Gene.Symbol",all.x=T,all.y=F)
  return(final_filtered_maf)
  }







