library(factoextra)
library(cluster)
oncogenes_normalised_df<-b
scaled_normalised_gene_counts_df<-scale(oncogenes_normalised_df)

#define linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
ac <- function(x) {
  agnes(scaled_normalised_gene_counts_df, method = x)$ac
}

#calculate agglomerative coefficient for each clustering linkage method
sapply(m, ac)


#ward produced the highest agglomerative coefficient

#perform hierarchical clustering using Ward's minimum variance
clust <- agnes(scaled_normalised_gene_counts_df, method = "ward")

#produce dendrogram
pltree(clust, cex = 0.6, hang = -1, main = "Dendrogram") 



#calculate gap statistic for each number of clusters (up to 10 clusters)
gap_stat <- clusGap(scaled_normalised_gene_counts_df, FUN = hcut, nstart = 25, K.max = 60, B = 50)

#produce plot of clusters vs. gap statistic
fviz_gap_stat(gap_stat)
