rm(list=ls())

#Open Genes Count Data, convert to dataframe
counts_files<-list.files("RNA_gene_counts_data/raw/", include.dirs=TRUE) #gets names of all files in raw/ directory
gene_counts_df <-data.frame(matrix(nrow=33126))
cols<-list()


for (file in counts_files){
  cols<-append(cols,substr(file,1,6))
  counts_column<-read.table(paste("RNA_gene_counts_data/raw/",file, sep=""), sep="\t")[2] #reads gene counts files and converts to table
  
  
  gene_counts_df[substr(file,1,6)]<-counts_column
}
rows=read.table(paste("RNA_gene_counts_data/raw/",file, sep=""), sep="\t")[1]

rownames(gene_counts_df)<-rows[[1]]

gene_counts_df[1]<-NULL
#Remove first 5 rows
gene_counts_df<-head(gene_counts_df,-5)



#Open Gene Lengths data, convert to data frame and sort
gene_lengths<-read.table(("gene_lengths/human_gene_lengths.txt"), sep="", header=TRUE) #reads gene counts files and converts to table

row.names(gene_lengths)<-gene_lengths[[1]]
gene_lengths[[1]]<-NULL

gene_lengths<-as.data.frame.matrix(gene_lengths)
gene_lengths<-gene_lengths[order(row.names(gene_lengths)),]

#Filter gene lengths data based on known cancer causing genes
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


#Perform RPKM normalisation
gene_names<-rownames(gene_counts_df)
num_rows<-length(rownames(gene_counts_df))
num_cols<-length(colnames(gene_counts_df))

normalised_gene_counts_df<-gene_counts_df
##need to get counts in each row


##1.  Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
perMillionScalingFactor<-colSums(gene_counts_df)/1000000

##2. Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)

for (i in (1:num_cols)){
  normalised_gene_counts_df[,i]<-gene_counts_df[,i]/perMillionScalingFactor[i]
}

oncogenes_normalised_gene_counts_df<-normalised_gene_counts_df[(row.names(normalised_gene_counts_df) %in% oncogenes),]


num_rows<-length(rownames(oncogenes_normalised_gene_counts_df))
num_cols<-length(colnames(oncogenes_normalised_gene_counts_df))
##3. Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM."
for (j in (1:num_rows)){
  gene<-gene_names[j]
  
  longest_isoform_length=gene_lengths[gene,][[3]]
  oncogenes_normalised_gene_counts_df[j,]<-oncogenes_normalised_gene_counts_df[j,]/(longest_isoform_length/1000)
 
  
  
}

##Remove N/A rows
oncogenes_normalised_gene_counts_df<-na.omit(oncogenes_normalised_gene_counts_df)


#Generate Heatmap 
library("pheatmap")
pheatmap(cor(oncogenes_normalised_gene_counts_df))
heatmap(cor(oncogenes_normalised_gene_counts_df))


