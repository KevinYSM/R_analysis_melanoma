#https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html

library("ggfortify")
gene_counts_df<-get_RPKM_normalised_data()
pca_res<-prcomp(gene_counts_df, scale. = TRUE)
autoplot(pca_res, data=gene_counts_df, colour='col')
