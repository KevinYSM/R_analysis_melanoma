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


#Generate Heatmap 
heatmap_metadata<-as.data.frame(read_excel("data/sample_info_v3.xlsx", trim_ws=TRUE))
#Remove Leading 0s
heatmap_metadata$`PCB code`<- sub("^0+", "", heatmap_metadata$`PCB code`)
#Add PCB-
heatmap_metadata$`PCB code`<-paste0("PCB-", heatmap_metadata$`PCB code`)

heatmap_metadata$`PCB code`<-gsub(pattern="PCB-9", replacement="PCB-09",heatmap_metadata$`PCB code`)

heatmap_metadata$`PCB code`<-substr(heatmap_metadata$`PCB code`,1,nchar(heatmap_metadata$`PCB code`)-1)

rownames(heatmap_metadata)<-heatmap_metadata$`PCB code`

heatmap_metadata<-heatmap_metadata[,c("Gender", "Tumor site","Treatment pre or post (none have previous other treatment)", "TIDE Values")]
colnames(heatmap_metadata)<-c("Gender","Tumor Site","Treatment Status", "TIDE Values")
pheatmap(cor(log2_filtered_RPKM), border_color = NA, annotation_col = heatmap_metadata$`Tumor Site`)








