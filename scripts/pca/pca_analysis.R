#https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html

<<<<<<< HEAD
library("ggfortify")
gene_counts_df<-open_gene_counts()
oncogenes_df<-get_oncogenes_df()
pca_res<-prcomp(oncogenes_df, scale. = TRUE)
autoplot(pca_res)
autoplot(pca_res, data=oncogenes_df,  label=TRUE,label.size = 3, loadings=TRUE)


#with eigen vectors
autoplot(pca_res, data=oncogenes_df,  label=TRUE,label.size = 3, loadings=TRUE)


#plotting K-means
set.seed(1)
autoplot(kmeans(oncogenes_df, 5), data = oncogenes_df)



rm(list=ls())

=======
library(ggfortify)
# GET FILTERED GENE VALUES
filtered_RPKM_values

# LOG TRANSFORM GENE VALUES
log_transformed_sample_counts<-t(log2(filtered_RPKM_values+1))


# PERFORM PCA

pca_res<-prcomp(log_transformed_sample_counts, scale. = TRUE)
autoplot(pca_res, data=log_transformed_sample_counts)
>>>>>>> bc169e2257b3591c736c2aaef48f98e3ed221a44
