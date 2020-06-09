view_individual_gene_importance <- function(importance_matrix, geneName, epi_type = "ALL"){
  tmp <- matrix(unname(unlist(importance_matrix[geneName, ])), nrow = 2500, byrow = F)
  df <- data.frame(site =seq(1, 2500))
  df$DNase <- tmp[, 1]
  df$H3K4me1 <- tmp[, 2]
  df$H3K4me3 <- tmp[, 3]
  df$H3K9me3 <- tmp[, 4]
  df$H3K27me3 <- tmp[, 5]
  df$H3K36me3 <- tmp[, 6]
  df$WGBS <- tmp[, 7]
  df$DNase <- loess(df$DNase ~ df$site, span=0.1)$fitted
  df$H3K4me1 <-  loess(df$H3K4me1 ~ df$site, span=0.1)$fitted
  df$H3K4me3 <-  loess(df$H3K4me3 ~ df$site, span=0.1)$fitted
  df$H3K9me3 <-  loess(df$H3K9me3 ~ df$site, span=0.1)$fitted
  df$H3K27me3 <-  loess(df$H3K27me3 ~ df$site, span=0.1)$fitted
  df$H3K36me3 <-  loess(df$H3K36me3 ~ df$site, span=0.1)$fitted
  df$WGBS <-  loess(df$WGBS ~ df$site, span=0.1)$fitted
  df <- gather(df, "epi", "importance", 2:8)
  df$epi <- factor(df$epi, levels = c("DNase", "H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3", "WGBS"))
  if(epi_type != "ALL"){
    colors <- gg_color_hue(7)
    epis <- c("DNase", "H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3", "WGBS")
    color_ <- vector()
    for(epi in epi_type){
      color_ <- c(color_, which(epis == epi))
    }
    color_ <- colors[sort(color_)]
    df_sub <- df[df$epi %in% epi_type, ]
    plots <- ggplot(df_sub, aes(site, importance, color = epi)) +
      geom_point(color = color_) + theme_bw()
    return(plots)
  }else{
    plots <- ggplot(df, aes(site, importance, color = epi)) +
      geom_point() + theme_bw()
    return(plots)
  }
}
