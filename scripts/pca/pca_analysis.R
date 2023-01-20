#https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html

library(ggfortify)
# GET FILTERED GENE VALUES
filtered_RPKM_values

# LOG TRANSFORM GENE VALUES
log_transformed_sample_counts<-t(log2(filtered_RPKM_values+1))


# PERFORM PCA

pca_res<-prcomp(log_transformed_sample_counts, scale. = TRUE)
autoplot(pca_res, data=log_transformed_sample_counts)
