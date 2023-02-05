library("maftools")
#Open COSMIC cancer mutations file

if (!exists("cancer_mutations_data")){
  cancer_mutations_data<-read.table("data/COSMIC_data/cancer_mutation_census_data/cmc_export.tsv", sep="\t", header=TRUE)
}

#Get maf files paths
#maf_files<- as.list(Sys.glob("N:\\local\\proj\\bioinformatics_project\\data\\processed\\vcf_processing_workflow\\*.filtered.ann.maf")) #gets names of all files in raw/ directory

  
maf_files<-as.list(Sys.glob("N:\\local\\proj\\bioinformatics_project\\data\\processed\\vcf_processing_workflow\\*mutect2.filtered.ann.maf"))



#Get gene counts
#gene_sample_counts<-update_gene_sample_counts(maf_files)


#Filter all genes present in less than 2 samples
#filtered_gene_list<-subset(gene_sample_counts, gene_sample_counts$`Samples Present`>1)$'Hugo_Symbol'


##Merge and filter all mafs
if (!exists(laml_filtered)){
  laml_filtered<-merge_mafs(lapply(maf_files,function(x) subsetMaf(read.maf(x),query="FILTER == 'PASS'", includeSyn=F)))
}

#Keep FILTER=PASS
df_filtered_maf<-subset(laml_filtered@data,laml_filtered@data$FILTER=="PASS")

#Remove existion_variation contains an RS value
TSG_list<-get_census_genes_list("TSG")
#If contains RS, remove if TSG

df_filtered_maf<-df_filtered_maf[!(grepl("rs", df_filtered_maf$Existing_variation) & !(df_filtered_maf$Hugo_Symbol %in% ONCOGENE_list)),]


#Filter if not an oncogene or TSG
df_filtered_maf<-subset(df_filtered_maf, df_filtered_maf$Hugo_Symbol %in% ALL_GENE_list)

# df_filtered_maf<-subset(df_filtered_maf, 
#                         (!grepl("rs",df_filtered_maf$Existing_variation)
#                         | df_filtered_maf$Hugo_Symbol %in% c("BRAF","KRAS","NRAS","HRAS")) & 
#                           (df_filtered_maf$Variant_Classification != "Splice_Site" | 
#                              df_filtered_maf$Hugo_Symbol %in% TSG_list))


#laml<-merge_mafs(lapply(maf_files, read.maf))
df_filtered_maf$Tumor_Sample_Barcode<-as.character(df_filtered_maf$Tumor_Sample_Barcode)
df_filtered_maf$Tumor_Sample_Barcode<-substr((df_filtered_maf$Tumor_Sample_Barcode),1,nchar(df_filtered_maf$Tumor_Sample_Barcode)-17)

df_filtered_maf$Tumor_Sample_Barcode<-gsub("PCB-","",df_filtered_maf$Tumor_Sample_Barcode)
df_filtered_maf$Tumor_Sample_Barcode<-gsub("PCB","",df_filtered_maf$Tumor_Sample_Barcode)
df_filtered_maf$Tumor_Sample_Barcode<-gsub("PDX","",df_filtered_maf$Tumor_Sample_Barcode)
df_filtered_maf$Tumor_Sample_Barcode<-gsub("_vs_","_",df_filtered_maf$Tumor_Sample_Barcode)
df_filtered_maf$Tumor_Sample_Barcode<-gsub("-","",df_filtered_maf$Tumor_Sample_Barcode)
df_filtered_maf$Tumor_Sample_Barcode<-substr(df_filtered_maf$Tumor_Sample_Barcode,1,2)
oncoplot(maf=read.maf(df_filtered_maf), minMut=3, top=5, removeNonMutated = F,showTumorSampleBarcodes = T)
plotmafSummary(maf = read.maf(df_filtered_maf), rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

update_gene_sample_counts<-function(maf_list=maf_files){
  #Count mutated genes found across samples
  gene_sample_counts<-data.frame(matrix(nrow=length(ALL_GENE_list),ncol=3))
  colnames(gene_sample_counts)<-c("Hugo_Symbol","Samples Present", "Samples List")
  gene_sample_counts[,1]<-ALL_GENE_list
  gene_sample_counts[,2]<-0
  
  maf_index<-0
  for (maf in maf_list){
    filtered_maf<-filter_maf(maf)
    
    for (index in (1:length(filtered_maf))){
      gene<-filtered_maf$Hugo_Symbol[index]
      sample_name<-substr(basename(maf),1,6)
      
      if (gene %in% unique(gene_sample_counts$Hugo_Symbol))
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



#Get gene counts in df_filtered_maf

res <- list() 
for (gene in ALL_GENE_list){
res[[gene]] <- length(unique(df_filtered_maf$Tumor_Sample_Barcode[df_filtered_maf$Hugo_Symbol==gene])) 
}
res<-as.data.frame(res)
res<-names(res)[res>3]
oncogene_data<-read.csv("data/COSMIC_data/cancer_gene_census_data/Census_allSun Jan 15 04_37_51 2023.csv")
present_oncogenes_exome<-subset(oncogene_data, oncogene_data$Gene.Symbol %in% res)
present_oncogenes_exome<-present_oncogenes_exome[grep(pattern = "melanoma", present_oncogenes_exome$Tumour.Types.Somatic., invert=T),]
present_oncogenes_exome<-present_oncogenes_exome[grep(pattern = "melanoma", present_oncogenes_exome$Tumour.Types.Germline., invert=T),]
present_oncogenes_exome<-present_oncogenes_exome[grep(pattern = "skin", present_oncogenes_exome$Tumour.Types.Germline., invert=T),]
present_oncogenes_exome<-present_oncogenes_exome[grep(pattern = "multiple other", present_oncogenes_exome$Tumour.Types.Germline., invert=T),]

subset(cosmic, Hugo_symbol==rownames(res))

