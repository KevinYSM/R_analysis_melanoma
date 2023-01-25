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


#Generate Heatmap 

pheatmap(cor(log2_filtered_RPKM))
heatmap(cor(log2_filtered_RPKM))

#Clustering
plot(pheatmap(cor(log2_filtered_RPKM))$tree_row)

