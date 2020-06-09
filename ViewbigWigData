generate_signal_basewise <- function(data_path, genome_region_bedFile){
 <- import(data_path, format = "BigWig")
H3K4me3 <- data.frame(chr = as.character(seqnames(H3K4me3)),
                      start = start(H3K4me3),
                      end = end(H3K4me3),
                      score = score(H3K4me3),
                      stringsAsFactors=F)

gene_H3K4me3 <- list()
for(i in chr_name){
  gene_H3K4me3_chr <- list()
  hg38_promoter_chr <- hg38_promoter[hg38_promoter$chr == i, ]
  H3K4me3_chr <- H3K4me3[H3K4me3$chr == i, ]
  H3K4me3_score <- rep(H3K4me3_chr$score, H3K4me3_chr$end - H3K4me3_chr$start + 1)
  for(j in 1:nrow(hg38_promoter_chr)){
    if(hg38_promoter_chr$strand[j] == "+"){
      gene_H3K4me3_chr[[hg38_promoter_chr$ensembl_ID[j]]] <- H3K4me3_score[hg38_promoter_chr$start[j]:hg38_promoter_chr$end[j]]
    }else{
      gene_H3K4me3_chr[[hg38_promoter_chr$ensembl_ID[j]]] <- rev(H3K4me3_score[hg38_promoter_chr$start[j]:hg38_promoter_chr$end[j]])
    }
  }
  gene_H3K4me3[[i]] <- as.data.frame(do.call(rbind, gene_H3K4me3_chr))
}
gene_H3K4me3 <- do.call(rbind, gene_H3K4me3)
gene_H3K4me3 <- gene_H3K4me3[hg38_promoter$rowname, ]
sum(rownames(gene_H3K4me3) == hg38_promoter$rowname)
save(gene_H3K4me3, file = "/SingleCell/Shang/Data/HepG2/H3K4me3_2500bp.RData")
}
