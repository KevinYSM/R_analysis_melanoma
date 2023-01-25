RPKM_gene_counts_df<-get_RPKM_normalised_data()


#FILTER GENES EXPRESSED IN LESS THAN X samples
X=5
filtered_expression_df_1<-RPKM_gene_counts_df[rowSums(RPKM_gene_counts_df > 0)>=X,] 

#KEEP TOP N genes based on gene expression
N=15000

MEAN_df<-data.matrix(apply(filtered_expression_df_1,1,mean))
MEAN_df<-data.matrix(sort(MEAN_df[,1], decreasing=TRUE),)
top_N_genes<-data.matrix(head(rownames(MEAN_df),N))


filtered_expression_df_2<-subset(filtered_expression_df_1, rownames(filtered_expression_df_1) %in% top_N_genes)




#Filter top Y most impactful genes from coefficient of variation DF.
Y=10000
coefficient_variation_df<-(apply(filtered_expression_df_2,1,sd)/rowMeans(filtered_expression_df_2))


top_Y_genes<-names(head(sort(abs(coefficient_variation_df),decreasing=TRUE),n=Y))

#Subset RPKM values
filtered_RPKM_values<-subset(RPKM_gene_counts_df, rownames(RPKM_gene_counts_df) %in% top_Y_genes)



hist(coefficient_variation_df, breaks=30)
head(sort(coefficient_variation_df, decreasing=TRUE), 50)


#

#log2RPKM+1 data

log2_RPKM_data<-log2(filtered_expression_df_2+0.00001)

#Coefficient of Variation
a<-apply(m,1,sd)
b<-apply(m,1,mean)
c<-a/b #returns single column matrix


#TO get RPKM data (gene_counts_df["CDH10",])/((3436/1000)*(colSums(gene_counts_df)/1000000))


#DETERMINE OPTIMUM FILTRATION VALUES
library("plot3D")
N_values<-seq(15000,20000,1000)
X_values<-seq(3,10)
filtration_results<-data.frame(matrix(ncol=3, nrow=length(N_values)*length(X_values)))
colnames(filtration_results)<-c("Top X Genes","Less than Y samples","Gene Count")

MEAN_df<-data.matrix(apply(RPKM_gene_counts_df,1,mean))
MEAN_df<-data.matrix(sort(MEAN_df[,1], decreasing=TRUE),)
index=1
for (N in N_values){
  top_N_genes<-data.matrix(head(rownames(MEAN_df),N))
  filtered_expression_df_1<-subset(RPKM_gene_counts_df, rownames(RPKM_gene_counts_df) %in% top_N_genes)
  for (X in X_values){
    filtered_expression_df_2<-filtered_expression_df_1[rowSums(filtered_expression_df_1 > 0)>=X,] 
    filtration_results[index,]<-c(N,X,length(rownames(filtered_expression_df_2)))
    index<-index+1
  }
}

lines3D(filtration_results$`Top X Genes`,filtration_results$`Less than Y samples`,filtration_results$`Gene Count`, main="Gene Filtration Results", bty="b2", theta=10, phi=0)

