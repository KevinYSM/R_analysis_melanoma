
if (!exists("RPKM_data")){
  RPKM_data<-get_RPKM_normalised_data()
}

TIDE_input<-RPKM_data-rowMeans(RPKM_data)
previous_immunotherapy_list<-c("PCB-30","PCB-31","PCB-54")

write.table(TIDE_input[,previous_immunotherapy_list],file="data/tide_data/tide_input_treated", row.names=T, col.names=T, sep="\t")
write.table(TIDE_input[,!(colnames(TIDE_input) %in% previous_immunotherapy_list)], file="data/tide_data/tide_input_untreated", row.names=T, col.names=T, sep="\t")
