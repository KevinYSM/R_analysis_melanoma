library(ggfortify)
library(cluster)
library("pheatmap")
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

#Perform PCA analysis 
pca_res<-prcomp(t(log2_filtered_RPKM), scale. = TRUE)
autoplot(pca_res, data=t(log2_filtered_RPKM), label=TRUE)


#Looking for enriched pathways that are driving variation along PC1. 
write.table(names(sort(pca_res$rotation[,1], decreasing=TRUE)), quote=F, row.names=F, col.names=F)

#MSigDB -> Investigate 
#Perform Clustering

#Generate Heatmap 

