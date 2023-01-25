library("maftools")
maf_files<-list.files("data/vcf_data/", include.dirs = TRUE, full.names=TRUE) #gets names of all files in raw/ directory
###???laml.maf=system.file('extdata', "data/vcf_data/PCB-09-PDX_S1_L004.disambiguatedSpeciesA.sorted.MarkDuplicates.split.recal.pass2.bam.haplotypecaller.filtered.ann.maf*", package='maftools')
maf_files<-maf_files[-1]

#Merge MAF files
laml<-merge_mafs(head(maf_files,1))
single_maf<-read.maf(maf_files[1])


#Basic Summary of file
laml

#Show sample summary
getSampleSummary(laml)

#Show gene summary
getGeneSummary(laml)

#Show clinical data associated with samples
getClinicalData(laml)

#show all fields in MAF
getFields(laml)

#Writes maf summary to an output file with basename laml
write.mafSummary(maf=laml, basename='laml')


#Plot MAF summary. Displays number of variants in each sample as a stacked bar plot 
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)


#Better representation of maf file
oncoplot(maf=laml, top=10)


#Transition and Transversions. titv function classifies SNPs into 
#Transitions and Transversions and returns a list of summarized tables 
#in various ways. Summarized data can also be visualized as a boxplot 
#showing overall distribution of six different conversions and as a 
#stacked barplot showing fraction of conversions in each sample.
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)

#plot titv summary
plotTiTv(res = laml.titv)

#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
lollipopPlot(
  maf = laml,
  gene = 'DNMT3A',
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
  labelPos = 882
)

#General protein domains can be drawn with the function plotProtein
plotProtein(gene = "TP53", refSeqID = "NM_000546")

#Rainfall Plots
brca <- system.file("extdata", "brca.maf.gz", package = "maftools")
brca = read.maf(maf = brca, verbose = FALSE)
rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.4)
