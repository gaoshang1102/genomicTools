View_signal <- function(data, geneName){
  df <- data.frame(sites = seq(1, 2500))
  df$signal <- unlist(data[geneName, ])
  ggplot(df, aes(sites, signal)) + geom_point() + theme_bw()
}
