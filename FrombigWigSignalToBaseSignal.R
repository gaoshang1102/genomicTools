library(rtracklayer)
library(ChIPpeakAnno)
generate_signal_basewise <- function(data_path, genome_region_bedFile_path, output_path){
 genome_region_bedFile <- read.table(genome_region_bedFile_path, sep = "\t", stringsAsFactors = F)
 chr_name <- paste0("chr", c(1:22, "X", "Y"))
 bigWig_data <- import(data_path, format = "BigWig")
 bigWig_data <- data.frame(chr = as.character(seqnames(bigWig_data)),
                      start = start(bigWig_data),
                      end = end(bigWig_data),
                      score = score(bigWig_data),
                      stringsAsFactors=F)
 gene_bigWig_data <- list()
 for(i in chr_name){
  gene_bigWig_data_chr <- list()
  genome_region_bedFile_chr <- genome_region_bedFile[genome_region_bedFile$chr == i, ]
  bigWig_data_chr <- bigWig_data[bigWig_data$chr == i, ]
  bigWig_data_score <- rep(bigWig_data_chr$score, bigWig_data_chr$end - bigWig_data_chr$start + 1)
  for(j in 1:nrow(genome_region_bedFile_chr)){
    if(genome_region_bedFile_chr$strand[j] == "+"){
      gene_bigWig_data_chr[[genome_region_bedFile_chr$ensembl_ID[j]]] <- bigWig_data_score[genome_region_bedFile_chr$start[j]:genome_region_bedFile_chr$end[j]]
    }else{
      gene_bigWig_data_chr[[genome_region_bedFile_chr$ensembl_ID[j]]] <- rev(bigWig_data_score[genome_region_bedFile_chr$start[j]:genome_region_bedFile_chr$end[j]])
    }
  }
  gene_bigWig_data[[i]] <- as.data.frame(do.call(rbind, gene_bigWig_data_chr))
 }
  gene_bigWig_data <- do.call(rbind, gene_bigWig_data)
  gene_bigWig_data <- gene_bigWig_data[genome_region_bedFile$rowname, ]
  print(sum(rownames(gene_bigWig_data) == genome_region_bedFile$rowname))
  save(gene_bigWig_data, file = output_path)
  return(gene_bigWig_data)
}
