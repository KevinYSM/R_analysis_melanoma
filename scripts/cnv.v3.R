library(GenomicRanges)
library(tidyverse)
library(ggrepel)
library(WriteXLS)
library(biomaRt) #
library(GenomicFeatures) #
library(Rsubread) #
library(readxl)
library(circlize)
library(ChIPpeakAnno) #

# Functions for reading and formatting CNV data --------------------------------
assemble_annotation <- function(gtf_path){
  ncbi_gtf <- flattenGTF(gtf_path)

  colnames(ncbi_gtf) <- c("gene_id","chromosome","start","end","strand")

  ncbi <- ncbi_gtf

  ncbi <- ncbi[!grepl("chrUn_",ncbi$chromosome),]
  ncbi <- ncbi[!grepl("_alt",ncbi$chromosome),]
  ncbi <- ncbi[!grepl("_random",ncbi$chromosome),]
 
  ncbi$strand <- NULL

  map2 <- ncbi
  colnames(map2)[colnames(map2)=="gene_id"] <- "hgnc_symbol"
  
  return(map2)
}

process_cnvs <- function(fname, cgc_path, map){

  sname <- gsub(basename(fname),pattern=".p.value.txt",replacement = "",fixed = T)
  sname <- str_split_fixed(sname,"\\.",2)[,1]

  cnv <- read.table(fname,sep="\t",header=T)

  gr <- GRanges(seqnames = map$chromosome,ranges = IRanges(start=map$start,end=map$end))
  values(gr) <- DataFrame(hgnc_symbol=map$hgnc_symbol)
  
  seg <- data.frame(ID=sname,
                    chrom=paste0("chr",cnv$chr),
                    loc.start=cnv$start,
                    loc.end=cnv$end,
                    num.mark=cnv$end-cnv$start,
                    seg.mean=cnv$copy.number)

  cnv.gr <- cnv[,1:3]
  cnv.gr$chr <- paste0("chr",cnv.gr$chr)
  cnv.gr <- GRanges(cnv.gr)
  values(cnv.gr) <- cnv[,4:ncol(cnv)]

  cgc <- read.table(cgc_path, sep = ",", header = T)

  values(cnv.gr)$hgnc_symbol <- ""
  glist <- list()
  for (i in 1:length(cnv.gr)){
    gene_list <- subsetByOverlaps(gr, cnv.gr[i])
    cnv.gr.df <- as.data.frame(gene_list)

    values(cnv.gr)[i,"hgnc_symbol"] <- gsub(paste0(
      unique(cnv.gr.df$hgnc_symbol),collapse=", "),pattern="^, ",replacement = "")
    glist[[i]] <- cnv.gr.df$hgnc_symbol
  }

  glist.cgc <- list()
  for (i in 1:length(glist)){
    tmp <- data.frame(Gene.Symbol=glist[[i]])
    glist.cgc[[i]] <- merge(tmp,cgc,all.x=F,all.y=F,by="Gene.Symbol")
  }

  glist.cgc.annot <- list()
  for (i in 1:length(glist.cgc)){
    if (nrow(glist.cgc[[i]])>0){
      glist.cgc.annot[[i]] <- cbind(glist.cgc[[i]],as.data.frame(cnv.gr[i]))
      glist.cgc.annot[[i]]$n_genes_in_region <- length(unlist(strsplit(
        unique(glist.cgc.annot[[i]]$hgnc_symbol),", ")))
      glist.cgc.annot[[i]]$hgnc_symbol <- NULL
    }
  }

  if(length(glist.cgc.annot)>0){
    glist.cgc.annot <- glist.cgc.annot[which(unlist(lapply(glist.cgc.annot,function(x){!is.null(x)})))]
    glist.cgc.annot <- unique(do.call("rbind",glist.cgc.annot))
    

    glist.cgc.annot.tsg <- glist.cgc.annot[grepl(glist.cgc.annot[["Role.in.Cancer"]],pattern="TSG") &
                                             (glist.cgc.annot[["copy.number"]] < 2),][,
                                                                                      c("Gene.Symbol","Tier","copy.number","n_genes_in_region",
                                                                                        "Role.in.Cancer","Mutation.Types","Tumour.Types.Somatic.")]
    
    glist.cgc.annot.tsg <- glist.cgc.annot.tsg[order(glist.cgc.annot.tsg$n_genes,decreasing = F),]

    glist.cgc.annot.oncogene <- glist.cgc.annot[grepl(glist.cgc.annot[["Role.in.Cancer"]],pattern="oncogene") &
                                                  glist.cgc.annot[["copy.number"]] > 2,][,
                                                                                         c("Gene.Symbol","Tier","copy.number","n_genes_in_region",
                                                                                           "Role.in.Cancer","Mutation.Types","Tumour.Types.Somatic.")]#,
    
    glist.cgc.annot.oncogene <- glist.cgc.annot.oncogene[order(glist.cgc.annot.oncogene$n_genes,decreasing = F),]

    cnv.cgc.excel <- list(Oncogenes=glist.cgc.annot.oncogene,
                          Tumor_suppressors=glist.cgc.annot.tsg,
                          All_CGC_in_regions=glist.cgc.annot[,c("Gene.Symbol","Tier","copy.number","n_genes_in_region",
                                                                "Role.in.Cancer","Mutation.Types","Tumour.Types.Somatic.")])#,
    
  } else {
    cnv.cgc.excel <- list()
    glist.cgc.annot <- data.frame()
    glist.cgc.annot.tsg <- data.frame()
    glist.cgc.annot.oncogene <- data.frame()
  }

  return(list(cnv.cgc.excel=cnv.cgc.excel,
              glist.cgc.annot=glist.cgc.annot,
              glist.cgc.annot.tsg=glist.cgc.annot.tsg,
              glist.cgc.annot.oncogene=glist.cgc.annot.oncogene,
              cnv=cnv,
              seg=seg))
}

get_cnvs <- function(fnames, cgc_path, map){
  out <- lapply(fnames,function(x){process_cnvs(x,cgc_path = cgc_path, map = map)})
  names(out) <- gsub(basename(fnames), pattern=".p.value.txt", replacement = "", 
                     fixed = T)
  names(out) <- str_split_fixed(names(out), "\\.", 2)[,1]
  
  subcategs.tab <- list()
  categs.tab <- list()
  categs <- names(out[[1]])
  for (categ in categs){
    for (sname in names(out)){
      subcategs <- names(out[[sname]][[categ]])
      if (length(subcategs)>0){
        subcategs.tab[[sname]] <- list()
        if (class(out[[sname]][[categ]])=="list"){
          for (subcateg in subcategs){
            if (nrow(out[[sname]][[categ]][[subcateg]])>0){
              subcategs.tab[[sname]][[subcateg]] <- as.data.frame(
                out[[sname]][[categ]][[subcateg]],stringsAsFactors=F)
              subcategs.tab[[sname]][[subcateg]]$sname <- sname
            }
          }
        } else if (class(out[[sname]][[categ]])=="data.frame"){
          if (nrow(out[[sname]][[categ]])>0){
            subcategs.tab[[sname]] <- as.data.frame(out[[sname]][[categ]],
                                                    stringsAsFactors=F)
            subcategs.tab[[sname]]$sname <- sname
          }
        }
      }
    }
    if (class(out[[1]][[categ]])=="list"){
      subcategs <- names(out[[1]][[categ]])
      subcategs.tab.merged <- list()
      for (subcateg in subcategs){
        subcategs.tab.merged[[subcateg]] <- do.call("rbind",lapply(
          subcategs.tab,function(x){x[[subcateg]]}))
        rownames(subcategs.tab.merged[[subcateg]]) <- NULL
      }
    } else if (class(out[[1]][[categ]])=="data.frame"){
      subcategs.tab.merged <- do.call("rbind",subcategs.tab)
      rownames(subcategs.tab.merged) <- NULL
    }
    categs.tab[[categ]] <- subcategs.tab.merged
  }
  
  snames_vs <- str_split_fixed(gsub(
    basename(fnames),pattern=".p.value.txt",replacement = "",fixed = T),"\\.",2)[,1][
      grepl(basename(dirname(fnames)),pattern="_vs_")]
  for (categ in names(categs.tab[["cnv.cgc.excel"]])){
    categs.tab[["cnv.cgc.excel"]][[categ]]$vs_normal <- F
    categs.tab[["cnv.cgc.excel"]][[categ]]$vs_normal[
      categs.tab[["cnv.cgc.excel"]][[categ]]$sname %in% snames_vs] <- T
  }
  
  return(categs.tab)
}

read_rpkm <- function(fnames){
  library(DESeq2)
  
  snames <- str_split_fixed(basename(fnames),"_",2)[,1]
  print(snames)
  
  sampleTable <- data.frame(snames=snames,
                            fnames=fnames,
                            condition="same")
  
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "/",
                                    design = ~ 1)
  
  cnts <- as.data.frame(counts(dds))
  counts.norm <- cgtools::normalize_counts.rpkm(cnts,assembly="hg38")
  
  return(counts.norm)
}

merge_cnv_rna <- function(counts.norm, categs.tab){
  counts.norm.long <- reshape2::melt(counts.norm)
  colnames(counts.norm.long) <- c("Gene.Symbol","sname","RPKM")
  
  cnv <- categs.tab[["cnv.cgc.excel"]]$All_CGC_in_regions
  cnv$sname <- str_split_fixed(cnv$sname,"_vs",2)[,1]
  cnv$sname <- gsub(pattern="PCB-",replacement = "",cnv$sname)
  cnv$sname <- gsub(pattern="PCB",replacement = "",cnv$sname)
  cnv$sname <- str_split_fixed(cnv$sname,"-",2)[,1]
  cnv$sname <- paste0("PCB-",cnv$sname)
  
  cnv <- cnv[cnv$Gene.Symbol %in% rownames(counts.norm),
             c("Gene.Symbol","copy.number","n_genes_in_region",
               "Role.in.Cancer","sname")]
  
  cnv.rna.merged <- merge(cnv, counts.norm.long, by=c("Gene.Symbol","sname"))
  
  return(cnv.rna.merged)
}

get_genome_cytoband <- function(path_genome_fa_sizes){
  lengths <- read.table(path_genome_fa_sizes, header=F)
  lengths <- setNames(lengths$V2,lengths$V1)
  genome_cytoband <- as.data.frame(lengths)
  genome_cytoband$V1 <- rownames(genome_cytoband)
  rownames(genome_cytoband) <- NULL
  genome_cytoband$V2 <- 0
  genome_cytoband$V4 <- ""
  genome_cytoband$V5 <- "gneg"
  genome_cytoband <- genome_cytoband[,c("V1","V2","lengths","V4","V5")]
  colnames(genome_cytoband) <- paste0("V",1:5)
  genome_cytoband <- genome_cytoband[genome_cytoband$V1!="chrM",]
  genome_cytoband <- genome_cytoband[!grepl("chrUn_",genome_cytoband$V1),]
  genome_cytoband <- genome_cytoband[!grepl("_alt",genome_cytoband$V1),]
  genome_cytoband <- genome_cytoband[!grepl("_random",genome_cytoband$V1),]
  genome_cytoband <- genome_cytoband[!grepl("HLA",genome_cytoband$V1),]
  genome_cytoband <- genome_cytoband[!grepl("EBV",genome_cytoband$V1),]
  
  return(genome_cytoband)
}

# Functions for plotting -------------------------------------------------------
format_seg_for_plot <- function(seg){
  seg_mean_met_df <- as.data.frame(seg,stringsAsFactors=F)
  seg_mean_met_df$ID <- NULL
  seg_mean_met_df$strand <- NULL
  colnames(seg_mean_met_df)[colnames(seg_mean_met_df)=="seg.mean"] <- "mean"
  seg_mean_met_df$num.mark <- NULL
  seg_mean_met_df$mean <- sapply(seg_mean_met_df$mean,function(x){
    ifelse(x==0,-4,log2(x/2))})
  return(seg_mean_met_df)
}

match_regions <- function(seg_mean_df,seg_mean_met_df){
  s <- seg_mean_df[seg_mean_df$seqnames!="chrX",c("seqnames","start","end")]
  s$mean <- NA
  seg_mean_met_df_gr <- GRanges(seg_mean_met_df)
  for (i in 1:nrow(s)){
    gr <- GRanges(s[i,])
    olaps <- findOverlaps(gr,seg_mean_met_df_gr)
    gr_idx <- queryHits(olaps)
    seg_mean_met_df_gr_idx <- subjectHits(olaps)
    s[i,]$mean <- mean(seg_mean_met_df_gr[seg_mean_met_df_gr_idx]$mean)
  }
  seg_mean_met_df <- s
  return(seg_mean_met_df)
}

get_bed_all <- function(map){
  gr <- GRanges(seqnames = map$chromosome, ranges = IRanges(start=map$start, 
                                                            end=map$end))
  values(gr) <- DataFrame(hgnc_symbol = map$hgnc_symbol)
  bed_all <- as.data.frame(gr)
  
  return(bed_all)
}

get_paired_bed <- function(seg, id_pr, id_pdx, bed_all, genome_cytoband, 
                           cnv_genes, cnv.rna.merged){
  
  seg.pdx <- format_seg_for_plot(seg[seg$ID==id_pdx,])
  seg.pr <- format_seg_for_plot(seg[seg$ID==id_pr,])
  
  cnames <- colnames(seg.pr)
  
  colnames(seg.pdx) <- c("seqnames","start","end","mean")
  colnames(seg.pr) <- c("seqnames","start","end","mean")
  
  uni <- rbind(seg.pdx, seg.pr)
  seg.pdx.matched <- match_regions(uni, seg.pdx)
  seg.pr.matched <- match_regions(uni, seg.pr)
  
  colnames(seg.pr.matched) <- cnames
  colnames(seg.pdx.matched) <- cnames
  
  seg.pr <- seg.pr.matched
  seg.pdx <- seg.pdx.matched
  seg.pr$mean[is.na(seg.pr$mean)] <- 0
  seg.pdx$mean[is.na(seg.pdx$mean)] <- 0
  
  cnv_genes <- cnv_genes[cnv_genes$sname %in% c(id_pdx,id_pr),]
  cnv_genes.symbols <- intersect(cnv_genes$Gene.Symbol[cnv_genes$sname==id_pdx],
                                 cnv_genes$Gene.Symbol[cnv_genes$sname==id_pr])
  cnv_genes <- cnv_genes[cnv_genes$Gene.Symbol %in% cnv_genes.symbols,]
  cnv_genes$color <- "black"
    cnv_genes$color[grepl("oncogene",cnv_genes$Role.in.Cancer)] <- "red"
      cnv_genes$color[grepl("TSG",cnv_genes$Role.in.Cancer)] <- "blue"
        cnv_genes$color[cnv_genes$Role.in.Cancer=="fusion"] <- "green"
          
        bed_all <- unique(bed_all[bed_all$hgnc_symbol %in% cnv_genes.symbols,
                                  c("seqnames","start","end","hgnc_symbol")])
        
        tmp <- list()
        for (gene in unique(bed_all$hgnc_symbol)){
          idx <- bed_all$hgnc_symbol==gene
          tmp[[gene]] <- data.frame(seqnames=unique(bed_all$seqnames[idx]),
                                    start=min(bed_all$start[idx]),
                                    end=max(bed_all$end[idx]),
                                    hgnc_symbol=gene)
        }
        bed_all <- do.call("rbind",tmp)
        bed_all$seqnames <- paste0("chr",bed_all$seqnames)
        bed_all <- bed_all[order(bed_all$seqnames),]
        
        bed_all <- merge(bed_all,unique(cnv_genes[,c("Gene.Symbol","color")]),
                         by.x="hgnc_symbol",by.y="Gene.Symbol")
        bed_all <- bed_all[,c("seqnames","start","end","hgnc_symbol","color")]
        bed_all <- bed_all[bed_all$color!="green",]
        
        bed_all <- merge(bed_all,unique(cnv.rna.merged[
          cnv.rna.merged$sname %in% c(id_pdx, id_pr), 
          c("Gene.Symbol","RPKM","sname")]), by.x="hgnc_symbol",by.y="Gene.Symbol")
        
        bed_all <- bed_all[,c("seqnames","start","end","hgnc_symbol","color","RPKM","sname")]
        bed_all$color <- "black"
          
        return(list(bed_all = bed_all,
                    seg.pr = seg.pr,
                    seg.pdx = seg.pdx))
}

get_bed <- function(seg, id_pr, bed_all, genome_cytoband, cnv_genes, cnv.rna.merged){
  seg.pr <- format_seg_for_plot(seg[seg$ID==id_pr,])
  
  cnames <- colnames(seg.pr)
  
  colnames(seg.pr) <- c("seqnames","start","end","mean")
  
  seg.pr$mean[is.na(seg.pr$mean)] <- 0
  
  cnv_genes <- cnv_genes[cnv_genes$sname %in% c(id_pr),]
  cnv_genes.symbols <- cnv_genes$Gene.Symbol[cnv_genes$sname==id_pr]
  cnv_genes <- cnv_genes[cnv_genes$Gene.Symbol %in% cnv_genes.symbols,]
  cnv_genes$color <- "black"
    cnv_genes$color[grepl("oncogene",cnv_genes$Role.in.Cancer)] <- "red"
      cnv_genes$color[grepl("TSG",cnv_genes$Role.in.Cancer)] <- "blue"
        cnv_genes$color[cnv_genes$Role.in.Cancer=="fusion"] <- "green"
          
        bed_all <- unique(bed_all[bed_all$hgnc_symbol %in% cnv_genes.symbols,
                                  c("seqnames","start","end","hgnc_symbol")])
        tmp <- list()
        for (gene in unique(bed_all$hgnc_symbol)){
          idx <- bed_all$hgnc_symbol==gene
          tmp[[gene]] <- data.frame(seqnames=unique(bed_all$seqnames[idx]),
                                    start=min(bed_all$start[idx]),
                                    end=max(bed_all$end[idx]),
                                    hgnc_symbol=gene)
        }
        bed_all <- do.call("rbind",tmp)
        bed_all <- bed_all[order(bed_all$seqnames),]
        
        bed_all <- merge(bed_all,unique(cnv_genes[,c("Gene.Symbol","color")]),
                         by.x="hgnc_symbol",by.y="Gene.Symbol")
        bed_all <- bed_all[,c("seqnames","start","end","hgnc_symbol","color")]
        bed_all <- bed_all[bed_all$color!="green",]
        
        bed_all <- merge(bed_all,unique(cnv.rna.merged[
          cnv.rna.merged$sname %in% c(id_pr), 
          c("Gene.Symbol","RPKM","sname")]), by.x="hgnc_symbol",by.y="Gene.Symbol")
        
        bed_all <- bed_all[,c("seqnames","start","end","hgnc_symbol","color","RPKM","sname")]
        bed_all$color <- "black"
          
        return(list(bed_all = bed_all,
                    seg.pr = seg.pr))
}

plot_cnv_circular <- function(bed_all, seg.pr, seg.pdx=NULL, genome_cytoband){
  circos.clear()
  circos.par(cell.padding = c(0, 0, 0, 0), gap.after = 0, start.degree = 0,
             track.margin = c(0,0))
  circos.initializeWithIdeogram(plotType = NULL, cytoband = genome_cytoband[genome_cytoband$V1!="chrX",])
  circos.genomicLabels(unique(bed_all[,! colnames(bed_all) %in% c("sname","RPKM")]), labels.column = 4, side = "outside",
                       cex = 0.5,
                       padding=0,
                       line_lwd=0.25,
                       niceFacing = T,
                       col = unique(bed_all[,! colnames(bed_all) %in% c("sname","RPKM")])$color,
                       line_col = unique(bed_all[,! colnames(bed_all) %in% c("sname","RPKM")])$color)
  circos.track(ylim = c(0, 10), panel.fun = function(x, y) {
    chr = gsub(CELL_META$sector.index,pattern="chr",replacement = "")
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.lines(c(xlim[1]+10^7/2, xlim[2]-10^7/2),c(0,0),col = "black",lwd = 0.25)
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.5, col = "black",
                facing = "inside", niceFacing = TRUE)
  }, track.height = 0.05, bg.border = NA)
  
  if (!is.null(seg.pdx)){
    circos.genomicTrack(seg.pdx,
                        panel.fun = function(region, value, ...) {
                          circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                                             col = ifelse(value[[1]] > 0, "red", "blue"), ...,border = NA,track.margin=NA,
                                             cell.padding=NA)
                          circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = NA,border = NA)
                        },
                        ylim=c(min(seg.pdx$mean),max(seg.pdx$mean)),bg.border=NA,track.height=0.1,
                        track.margin=c(0,0))
  }
  
  
  circos.genomicTrack(seg.pr,
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                                           col = ifelse(value[[1]] > 0, "red", "blue"), ...,border = NA,track.margin=NA,
                                           cell.padding=NA)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = NA,border = NA)
                      },
                      ylim=c(min(seg.pr$mean),max(seg.pr$mean)),bg.border=NA,track.height=0.1,track.margin=c(0,0))
  
}

plot_cnv_circular_scatter <- function(bed_all, id_pr, id_pdx = NULL, 
                                      genome_cytoband, cnv.rna.merged){
  dfgg <- bed_all
  dfgg$seqnames <- as.character(dfgg$seqnames)
  dfgg <- dfgg[!(grepl("X",dfgg$seqnames) | grepl("Y",dfgg$seqnames)),]
  dfgg$chr_nr <- as.numeric(gsub("chr","",dfgg$seqnames))
  if (is.null(id_pdx)){
    dfgg <- merge(dfgg[dfgg$sname==id_pr,],
                  cnv.rna.merged[cnv.rna.merged$sname==id_pr,
                                 c("Gene.Symbol","copy.number",
                                   "Role.in.Cancer","n_genes_in_region")],
                  by.x="hgnc_symbol",by.y="Gene.Symbol")
  } else {
    dfgg <- rbind(merge(dfgg[dfgg$sname==id_pr,],
                        cnv.rna.merged[cnv.rna.merged$sname==id_pr,
                                       c("Gene.Symbol","copy.number",
                                         "Role.in.Cancer","n_genes_in_region")],
                        by.x="hgnc_symbol",by.y="Gene.Symbol"),
                  merge(dfgg[dfgg$sname==id_pdx,],
                        cnv.rna.merged[cnv.rna.merged$sname==id_pdx,
                                       c("Gene.Symbol","copy.number",
                                         "Role.in.Cancer","n_genes_in_region")],
                        by.x="hgnc_symbol",by.y="Gene.Symbol"))
  }
  
  len <- genome_cytoband[genome_cytoband$V1!="chrX" & genome_cytoband$V1!="chrY",]
  len <- len[order(as.numeric(gsub("chr","",len$V1))),]
  len$chr_nr <- as.numeric(gsub("chr","",len$V1))
  len$cumulative <- cumsum(as.numeric(len$V3))
  
  dfgg$linear_coord <- NA
  for (i in 1:nrow(dfgg)){
    dfgg$linear_coord[i] <- (dfgg[i,]$start/10^6 + sum(len[len$chr_nr<dfgg[i,]$chr_nr,]$V3)/10^6)*10^6
  }
  
  g <- ggplot(dfgg, aes(x=linear_coord, y=log2(RPKM+1), fill=Role.in.Cancer, 
                        label=hgnc_symbol, color=Role.in.Cancer,
                        size=copy.number)) +
    geom_point() + guides(fill="none") +
    coord_polar() +
    scale_x_continuous(limits = c(0, max(len$cumulative)), breaks = c(0,len$cumulative),
                       labels = c(gsub("chr","",len$V1),"")) +
    geom_text_repel() + facet_wrap(~ sname) + theme_minimal() +
    theme(plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()) +
    xlab(NULL) + ylab("log2(RPKM + 1)")
  
  return(g)
}

# Find copy number changes matching known cancer genes and write to file -------
gtf_path <- "N:\\local\\reference\\igenomes\\Homo_sapiens\\NCBI\\GRCh38\\Annotation\\Genes\\genes.gtf"
cgc_path <- "data/COSMIC_data/cancer_gene_census_data/cancer_gene_census.csv"
fai_path <- "N:\\local\\reference\\igenomes\\Homo_sapiens\\GATK\\GRCh38\\Sequence\\WholeGenomeFasta\\Homo_sapiens_assembly38.fasta.fai" 
outdir <- "data/cnv_output/"
dir.create(outdir,showWarnings = F,recursive = T)
fnames.controlfreec <- Sys.glob("N:\\local\\proj\\bioinformatics_project\\scripts\\nextflow\\sarek_workflow\\sarek_final\\results\\variant_calling\\controlfreec\\*\\*.p.value.txt")

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
map <- assemble_annotation(gtf_path = gtf_path)

categs.tab <- get_cnvs(fnames.controlfreec, cgc_path = cgc_path, map = map)

WriteXLS(categs.tab[["cnv.cgc.excel"]],ExcelFileName = paste0(outdir,"cnv.cancer_genes_in_regions.xlsx"),
         col.names = T,row.names = F,SheetNames = names(categs.tab[["cnv.cgc.excel"]]),
         AutoFilter = T,AdjWidth = T)

for (sname in unique(categs.tab[["glist.cgc.annot"]]$sname)){
  df <- categs.tab[["glist.cgc.annot"]][
    grepl(categs.tab[["glist.cgc.annot"]][["Role.in.Cancer"]],pattern="oncogene") & categs.tab[["glist.cgc.annot"]]$sname == sname,]
  genes_highlight <- head(df[with(df, order(-copy.number, n_genes_in_region)), ]$Gene.Symbol,n=5)

  g <- ggplot(df,aes(x=n_genes_in_region,y=copy.number,
                     label=ifelse(Gene.Symbol %in% genes_highlight & copy.number > 2 & n_genes_in_region < 100,Gene.Symbol,""))) +
    geom_point() + geom_text_repel(max.overlaps = 100) + theme_classic() + ylab("Absolute copy number") + xlab("Number of genes in region")
  ggsave(g,filename = paste0(outdir,sname,".cnv.cancer_genes_in_regions.oncogenes.pdf"),width = 4,height = 4)


  df <- categs.tab[["glist.cgc.annot"]][
    grepl(categs.tab[["glist.cgc.annot"]][["Role.in.Cancer"]],pattern="TSG") & categs.tab[["glist.cgc.annot"]]$sname == sname,]
  genes_highlight <- head(df[with(df, order(copy.number, n_genes_in_region)), ]$Gene.Symbol,n=5)

  g <- ggplot(df,aes(x=n_genes_in_region,y=copy.number,
                     label=ifelse(Gene.Symbol %in% genes_highlight & copy.number < 2 & n_genes_in_region < 100,Gene.Symbol,""))) +
    geom_point() + geom_text_repel(max.overlaps = 100) + theme_classic() + ylab("Absolute copy number") + xlab("Number of genes in region")
  ggsave(g,filename = paste0(outdir,sname,".cnv.cancer_genes_in_regions.tsg.pdf"),width = 4,height = 4)

}

write.table(categs.tab$seg[,colnames(categs.tab$seg)!="sname"],
            file=paste0(outdir,"cnv.merged.seg"), sep="\t", col.names = T,
            row.names = F,quote = F)

# Plot copy number vs RNA ------------------------------------------------------

# Change these three lines to read in your own rpkm data instead

if (!exists("RPKM_data")){
  RPKM_data<-get_RPKM_normalised_data()
}
counts.norm <- RPKM_data

cnv.rna.merged <- merge_cnv_rna(counts.norm, categs.tab)

# Copy number vs RNA scatter
pdf(file=paste0(outdir,"cnv_vs_rna_cgc.pdf"),width = 17,height = 12)
ggplot(cnv.rna.merged, aes(x=copy.number,y=log2(RPKM+1),
                          fill=Gene.Symbol,color=Role.in.Cancer,
                          alpha=-log10(n_genes_in_region),
                          label=Gene.Symbol)) +
  geom_point() + guides(fill="none") + theme(legend.position="none") +
  geom_text_repel() + facet_wrap(~ sname) + theme_minimal() + xlab("Copy number")
dev.off()

# Copy number vs RNA circular scatter
genome_cytoband <- get_genome_cytoband(fai_path)

cnv_genes <- categs.tab$cnv.cgc.excel$All_CGC_in_regions
cnv_genes$sname <- str_split_fixed(cnv_genes$sname, "_vs", 2)[,1]
cnv_genes$sname <- gsub(pattern="PCB-", replacement = "", cnv_genes$sname)
cnv_genes$sname <- gsub(pattern="PCB", replacement = "", cnv_genes$sname)
cnv_genes$sname <- str_split_fixed(cnv_genes$sname,"-",2)[,1]
cnv_genes$sname <- paste0("PCB-",cnv_genes$sname)

cnv_genes.matched_rna <- cnv_genes[cnv_genes$sname %in% intersect(
  cnv_genes$sname,cnv.rna.merged$sname),]

seg <- read.table(file = paste0(outdir,"cnv.merged.seg"), header=T)
seg$ID <- str_split_fixed(seg$ID,"_vs",2)[,1]
seg$ID <- gsub(pattern="-PDX", replacement = "", seg$ID)
seg$ID <- gsub(pattern="PCB-", replacement = "", seg$ID)
seg$ID <- gsub(pattern="PCB", replacement = "", seg$ID)
seg$ID<-paste0("PCB-",seg$ID)
seg.matched_rna <- seg[seg$ID %in% intersect(seg$ID, cnv.rna.merged$sname),]

bed_all <- get_bed_all(map)

bed <- list()
for (sname in unique(seg.matched_rna$ID)){
  bed[[sname]] <- get_bed(seg = seg.matched_rna, 
                          id_pr = sname, 
                          bed_all = bed_all, 
                          genome_cytoband = genome_cytoband, 
                          cnv_genes = cnv_genes.matched_rna, 
                          cnv.rna.merged = cnv.rna.merged)
}

for (indv in names(bed)){
  g <- plot_cnv_circular_scatter(bed_all = bed[[indv]]$bed_all, 
                            id_pr = indv,
                            id_pdx = NULL, 
                            genome_cytoband = genome_cytoband, 
                            cnv.rna.merged = cnv.rna.merged)
  ggsave(plot = g, 
         filename = paste0(outdir, "/", indv, ".cnv.circular_scatter.pdf"),
         width = 10, height = 6)
}

# Write session info ------------------------------------------------------------
s <- sessionInfo()
saveRDS(s,paste0(outdir,"/sessionInfo.rda"))

