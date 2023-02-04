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

#get sample metadata
heatmap_metadata<-as.data.frame(read_excel("data/sample_info_v3.xlsx", trim_ws=TRUE))

#Remove Leading 0s
heatmap_metadata$`PCB code`<- sub("^0+", "", heatmap_metadata$`PCB code`)
#Add PCB-
heatmap_metadata$`PCB code`<-paste0("PCB-", heatmap_metadata$`PCB code`)

heatmap_metadata$`PCB code`<-gsub(pattern="PCB-9", replacement="PCB-09",heatmap_metadata$`PCB code`)

heatmap_metadata$`PCB code`<-substr(heatmap_metadata$`PCB code`,1,nchar(heatmap_metadata$`PCB code`)-1)

colnames(heatmap_metadata)<-make.names(colnames(heatmap_metadata))
rownames(heatmap_metadata)<-heatmap_metadata$`PCB.code`

#Perform PCA analysis 
immune_log2_filtered_rpkm<-t(log2_filtered_RPKM[rownames(log2_filtered_RPKM) %in% immune_related_genes$immune_related_genes,])
chemokine_log2_filtered_rpkm<-t(log2_filtered_RPKM[rownames(log2_filtered_RPKM) %in% immune_related_genes$chemokines,])
interleukin_log2_filtered_rpkm<-t(log2_filtered_RPKM[rownames(log2_filtered_RPKM) %in% immune_related_genes$interleukins,])
all_immune_genes_log2_filtered_rpkm<-cbind(cbind(immune_log2_filtered_rpkm,chemokine_log2_filtered_rpkm),interleukin_log2_filtered_rpkm)

pca_res<-prcomp((immune_log2_filtered_rpkm), scale. = TRUE)
a<-merge(immune_log2_filtered_rpkm, heatmap_metadata, by="row.names")
autoplot(pca_res, data=a, label=TRUE, colour="Gender")
autoplot(pca_res, data=a, label=TRUE, colour="TIDE.Values")
autoplot(pca_res, data=a, label=TRUE, colour="TIDE.Prediction")
autoplot(pca_res, data=a, label=TRUE, colour="Treatment.before.surgery")
autoplot(pca_res, data=a, label=TRUE, colour="Cell.line.successfully.generated.from.PDX")
autoplot(pca_res, data=a, label=TRUE, colour="Treatment.pre.or.post..none.have.previous.other.treatment.")

colnames(heatmap_metadata)

pca_res<-prcomp((chemokine_log2_filtered_rpkm), scale. = TRUE)
a<-merge(chemokine_log2_filtered_rpkm, heatmap_metadata, by="row.names")
autoplot(pca_res, data=a, label=TRUE, colour="Gender")
autoplot(pca_res, data=a, label=TRUE, colour="TIDE.Values")
autoplot(pca_res, data=a, label=TRUE, colour="TIDE.Prediction")
autoplot(pca_res, data=a, label=TRUE, colour="Treatment.before.surgery")
autoplot(pca_res, data=a, label=TRUE, colour="Cell.line.successfully.generated.from.PDX")
autoplot(pca_res, data=a, label=TRUE, colour="Treatment.pre.or.post..none.have.previous.other.treatment.")

pca_res<-prcomp((interleukin_log2_filtered_rpkm), scale. = TRUE)
a<-merge(interleukin_log2_filtered_rpkm, heatmap_metadata, by="row.names")
autoplot(pca_res, data=a, label=TRUE, colour="Gender")
autoplot(pca_res, data=a, label=TRUE, colour="TIDE.Values")
autoplot(pca_res, data=a, label=TRUE, colour="TIDE.Prediction")
autoplot(pca_res, data=a, label=TRUE, colour="Treatment.before.surgery")
autoplot(pca_res, data=a, label=TRUE, colour="Cell.line.successfully.generated.from.PDX")
autoplot(pca_res, data=a, label=TRUE, colour="Treatment.pre.or.post..none.have.previous.other.treatment.")

pca_res<-prcomp((all_immune_genes_log2_filtered_rpkm), scale. = TRUE)
a<-merge(all_immune_genes_log2_filtered_rpkm, heatmap_metadata, by="row.names")
autoplot(pca_res, data=a, label=TRUE, colour="Gender")
autoplot(pca_res, data=a, label=TRUE, colour="TIDE.Values")
autoplot(pca_res, data=a, label=TRUE, colour="TIDE.Prediction")
autoplot(pca_res, data=a, label=TRUE, colour="Treatment.before.surgery")
autoplot(pca_res, data=a, label=TRUE, colour="Cell.line.successfully.generated.from.PDX")
autoplot(pca_res, data=a, label=TRUE, colour="Treatment.pre.or.post..none.have.previous.other.treatment.")


g<-autoplot(pca_res, data=t(log2_filtered_RPKM), label=TRUE)



library("readxl")
sample_data<-as.data.frame(read_excel("data/sample_info_v3.xlsx", trim_ws=TRUE))
sample_data$`PCB code`<-substr(sample_data$`PCB code`,1,nchar(sample_data$`PCB code`)-1)
#Remove Leading 0s
sample_data$`PCB code`<- sub("^0+", "", sample_data$`PCB code`)
#Add PCB-
sample_data$`PCB code`<-paste0("PCB-", sample_data$`PCB code`)

sample_data$`PCB code`<-gsub(pattern="PCB-9", replacement="PCB-09",sample_data$`PCB code`)



rownames(sample_data)<-sample_data$`PCB code`
sample_data<-subset(sample_data, sample_data$RNAseq=="Y")
merged_pca_sample_data<-merge(sample_pca, sample_data, by.x="sample_id", by.y="PCB code")


#Clean up plots
g + theme_minimal()


sample_pca<-as.data.frame(pca_res$x)
sample_pca$sample_id<-rownames(sample_pca)
ggplot(merged_pca_sample_data,aes(x=PC1,y=PC2, color=`Gender`))+geom_text(aes(label=sample_id), size=1.5) + theme_minimal()


#Looking for enriched pathways that are driving variation along PC1. 
write.table(names(sort(pca_res$rotation[,1], decreasing=TRUE)), quote=F, row.names=F, col.names=F)

#MSigDB -> Investigate 
#Perform Clustering

#Generate Heatmap 

