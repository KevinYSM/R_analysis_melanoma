library(ggfortify)
library(cluster)
library("pheatmap")
library("readxl")
#Get RPKM Values
if (!exists("RPKM_data")){
  RPKM_data<-get_RPKM_normalised_data()
}
#Filter RPKM Values
if (!exists("filtered_genes")){
  filtered_genes<-get_filtered_genes()
}

if (!exists("filtered_RPKM")){
  filtered_RPKM<-get_filtered_rpkm()
}


#Get log2RPKM
if (!exists("log2_RPKM")){
  log2_filtered_RPKM<-log2(filtered_RPKM+1)
}

if (!exists("immune_related_genes")){
  immune_related_genes<-as.data.frame(read_excel("data/immune_related_genes.xlsx",trim_ws=TRUE))
}


if (!exists("immune_gene_filtered_RPKM")){
  immune_gene_filtered_RPKM<-as.data.frame(log2_filtered_RPKM[rownames(log2_filtered_RPKM) %in% immune_related_genes$immune_related_genes,])
}

if (!exists("chemokine_gene_filtered_RPKM")){
  chemokine_gene_filtered_RPKM<-as.data.frame(log2_filtered_RPKM[rownames(log2_filtered_RPKM) %in% immune_related_genes$chemokines,])
}

if (!exists("interleukin_gene_filtered_RPKM")){
  interleukin_gene_filtered_RPKM<-log2_filtered_RPKM[rownames(log2_filtered_RPKM) %in% immune_related_genes$interleukins,]
}


#Generate metadata
heatmap_metadata<-as.data.frame(read_excel("data/sample_info_v3.xlsx", trim_ws=TRUE))

#Remove Leading 0s
heatmap_metadata$`PCB code`<- sub("^0+", "", heatmap_metadata$`PCB code`)
#Add PCB-
heatmap_metadata$`PCB code`<-paste0("PCB-", heatmap_metadata$`PCB code`)

heatmap_metadata$`PCB code`<-gsub(pattern="PCB-9", replacement="PCB-09",heatmap_metadata$`PCB code`)

heatmap_metadata$`PCB code`<-substr(heatmap_metadata$`PCB code`,1,nchar(heatmap_metadata$`PCB code`)-1)

rownames(heatmap_metadata)<-heatmap_metadata$`PCB code`

heatmap_metadata<-heatmap_metadata[,c("Gender", "Tumor site","Treatment pre or post (none have previous other treatment)")]
colnames(heatmap_metadata)<-c("Gender","Tumor Site","Treatment Status")
pheatmap(cor(immune_gene_filtered_RPKM), border_color = NA, annotation_col = heatmap_metadata)
pheatmap(cor(chemokine_gene_filtered_RPKM), border_color = NA, annotation_col = heatmap_metadata)
pheatmap(cor(interleukin_gene_filtered_RPKM), border_color = NA, annotation_col = heatmap_metadata)

pheatmap(immune_gene_filtered_RPKM, border_color = NA, annotation_col = heatmap_metadata)
pheatmap(chemokine_gene_filtered_RPKM, border_color = NA, annotation_col = heatmap_metadata)
pheatmap(interleukin_gene_filtered_RPKM, border_color = NA, annotation_col = heatmap_metadata)
pheatmap(unique(rbind(rbind(interleukin_gene_filtered_RPKM, immune_gene_filtered_RPKM),chemokine_gene_filtered_RPKM),border_color = NA, annotation_col = heatmap_metadata))

cor_16_36_genes=c("B2M", "HLA-DQA1","HLA-DQB1","CD274","CXCL10","CXCL11","CXCL9","IL32")
pheatmap(unique(rbind(rbind(interleukin_gene_filtered_RPKM, immune_gene_filtered_RPKM),chemokine_gene_filtered_RPKM)[cor_16_36_genes,],border_color = NA, annotation_col = heatmap_metadata))

