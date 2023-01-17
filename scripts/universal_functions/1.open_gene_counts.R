open_gene_counts <-function(){
  counts_files<-list.files("RNA_gene_counts_data/raw/", include.dirs=TRUE) #gets names of all files in raw/ directory
  
  
  gene_counts_df <-data.frame(matrix(nrow=33126))
  cols<-list()
  for (file in counts_files){
    cols<-append(cols,substr(file,1,6))
    counts_column<-read.table(paste("RNA_gene_counts_data/raw/",file, sep=""), sep="\t")[2] #reads gene counts files and converts to table
    gene_counts_df[substr(file,1,6)]<-counts_column
  }
  rows=read.table(paste("RNA_gene_counts_data/raw/",file, sep=""), sep="\t")[1]
  
  rownames(gene_counts_df)<-rows[[1]]
  
  gene_counts_df[1]<-NULL
  #Remove first 5 rows
  gene_counts_df<-head(gene_counts_df,-5)
  #colnames(gene_counts_df)<-cols
  
  #setnames(gene_counts_df,NULL,"gene")
  #colnames(gene_counts_df)[0]<-"genes"
  #colnames(gene_counts_df)[1]<-"gene_id"
  #test=read.table(paste("raw/",counts_files[9],sep=""),sep="\t")[2]
  #gene_counts_df<-data.matrix(gene_counts_df)
  
  #return(gene_counts_df)
  
  
  
  
}

