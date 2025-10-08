library(ggplot2)
library(ggrepel)
#' @param feature2show this should be the name of the column representing the features ID in your dataset
#
doVolcano = function(feature2show="Gene", resi, NAME, Pcut=0.05, logFCcut=1, Annotate_features=TRUE) {
  selectedCols <- c(feature2show, "logFC", "P.Value", "adj.P.Val")
  toExport <- resi[, selectedCols]
  colnames(toExport)[1] = "Feature"
  toExport <- toExport[order(toExport$P.Value),]

  tmp <- as.data.frame(toExport) %>%
    filter(!is.na(adj.P.Val)) %>%
    mutate(sig = ifelse(adj.P.Val<Pcut&abs(logFC)>logFCcut,
                        ifelse(logFC > 0, "up", "down"), "nonSig"))
  tmpSignif <- tmp %>%
    filter(P.Value < quantile(tmp$P.Value, 0.001))
  
  p = ggplot(data = tmp, aes(logFC, -log10(P.Value), col = sig)) + 
    geom_point() + 
    xlab("Log2 Fold Change") +
    geom_hline(yintercept = -log10(Pcut), lty = 2) +
    geom_vline(xintercept = c(-logFCcut, logFCcut), lty = 2) + 
    scale_colour_manual(values = c("down"="blue", "nonSig"="black", "up"="red")) + 
    theme(legend.position="none") +
    ggtitle(NAME)

  if (Annotate_features) {
    p = p + geom_text_repel(data = tmpSignif, aes(logFC, 
                              -log(P.Value), label = Feature))
  }
  return(p)
}

