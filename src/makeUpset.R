makeUpset = function(rslt.list, log2FoldChange_thresh, padj_thresh, OUTDIR, h=3.5, w=6) {
  feature.list = lapply(names(rslt.list), function(i) {
    up = as.data.frame(rslt.list[[i]]) %>% 
      filter(logFC>log2FoldChange_thresh&adj.P.Val<padj_thresh)
    down = as.data.frame(rslt.list[[i]]) %>% 
      filter(logFC<log2FoldChange_thresh&adj.P.Val<padj_thresh)
    list(up$Protein_IDs, down$Protein_IDs)
  })  %>% setNames(names(rslt.list))
  
  feature.list_UP = c(feature.list[[1]][1], feature.list[[2]][1], feature.list[[3]][1]) %>% 
    setNames(c(paste0(names(feature.list)[1], "_UP"), 
               paste0(names(feature.list)[2], "_UP"), 
               paste0(names(feature.list)[3], "_UP")))
  feature.list_DOWN = c(feature.list[[1]][2], feature.list[[2]][2], feature.list[[3]][2]) %>% 
    setNames(c(paste0(names(feature.list)[1], "_DN"), 
               paste0(names(feature.list)[2], "_DN"), 
               paste0(names(feature.list)[3], "_DN")))
  
  library(ComplexHeatmap)
  pdf(paste0(OUTDIR, "/Sig.Diff.features_upset.pdf"), h=h, w=w)
  m = make_comb_mat(c(feature.list_UP, feature.list_DOWN))
  par(mar = c(4, 8, 2, 2)) # par(mar = c(bottom, left, top, right))
  print(UpSet(m,
        top_annotation = upset_top_annotation(m, add_numbers = TRUE),
        right_annotation = upset_right_annotation(m, add_numbers = TRUE),
        comb_order = rev(order(comb_size(m))),
        set_order = c(paste0(names(feature.list)[1], "_UP"), 
                      paste0(names(feature.list)[2], "_UP"), 
                      paste0(names(feature.list)[3], "_UP"),
                      paste0(names(feature.list)[1], "_DN"), 
                      paste0(names(feature.list)[2], "_DN"), 
                      paste0(names(feature.list)[3], "_DN")),
        column_title = paste0("Significant proteins\n(adj.P.Val<", padj_thresh, ")"),
        ))

  invisible(dev.off()) 
  return(list(feature.list_UP=feature.list_UP, feature.list_DOWN=feature.list_DOWN))
}

# log2FoldChange_thresh = 0
# padj_thresh = 0.2
# protein.list = makeUpset(log2FoldChange_thresh, padj_thresh)
