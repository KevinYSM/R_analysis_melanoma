library(rnaSeqPanCanClassifier)
classify_output<-classify("data/RNA_gene_counts_data_new", genome="hg38", outdir="classify_output_new")

run_tsne(inhouse=classify_output$inhouse, refset=classify_output$refset, outdir="tsne_output")
