set.seed(123)  # Set seed for reproducibility

library(pheatmap)
data <- matrix(rnorm(100, 0, 10), nrow = 10, ncol = 10)           # Create example data
colnames(data) <- paste0("col", 1:10)                             # Column names
rownames(data) <- paste0("row", 1:10)                             # Row names
cor(data)
p<-pheatmap( cor(data))
p

#Principle component analysis
prcomp(data)
princomp(data)
