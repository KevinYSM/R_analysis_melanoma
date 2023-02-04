library(ggplot2)
mouse_data<-maf_files<-list.files("data/disambiguate_txts/", include.dirs = TRUE, full.names=TRUE) #gets names of all files in raw/ directory

counts_df<-data.frame(matrix(ncol=4))
colnames(counts_df)<-c("sample","human_reads","mouse_reads","ambiguous_reads")
for (counts in mouse_data){
  
  mouse_read<-data.frame(read.table(counts, sep="\t"))
  colnames(mouse_read)<-c("sample","human_reads","mouse_reads","ambiguous_reads")

  mouse_read[1,]<-NA
  mouse_read<-na.omit(mouse_read)
  counts_df<-rbind(counts_df,mouse_read)
}
counts_df<-na.omit(counts_df)
counts_df$sample<-substr(counts_df$sample,1,6)



bar_plot_df<-data.frame(matrix(ncol=3))
colnames(bar_plot_df)<-c("sample","read_type","read_count")
for (counts in mouse_data){
  sample_id<-substr(basename(counts),1,6)
  bar_plot_df<-rbind(bar_plot_df,c(sample_id,"human",as.numeric(counts_df$human_reads[which(counts_df$sample==sample_id)])))
  bar_plot_df<-rbind(bar_plot_df,c(sample_id,"mouse",as.numeric(counts_df$mouse_reads[which(counts_df$sample==sample_id)])))
  bar_plot_df<-rbind(bar_plot_df,c(sample_id,"ambiguous",as.numeric(counts_df$ambiguous_reads[which(counts_df$sample==sample_id)])))

}
bar_plot_df<-na.omit(bar_plot_df)
bar_plot_df$read_count<-as.numeric(bar_plot_df$read_count)
bar_plot_df$sample<-substr(bar_plot_df$sample,5,6) 
pdf(file= "test", width=3, height=5)
bar_plot_df$read_type<-factor(bar_plot_df$read_type, levels=c("ambiguous","mouse","human"))
ggplot(bar_plot_df, aes(fill=read_type, y=read_count, x=sample)) + 
  geom_bar(position="fill", stat="identity", inherit.aes = T)
dev.off()
#position: fill, dodge, stack
ggplot(counts_df, aes(x=sample, y=human_reads))+geom_bar(stat="identity")
