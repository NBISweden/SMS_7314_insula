#' @importFrom ComplexHeatmap make_comb_mat UpSet
#' @param SampleNames this should be the same names in the same order as in pathway_rslt.list.
#' @param qval.thresh qval is less conservative than adjusted p value in fdr.
#
pathwayUpset = function(pathway_rslt.list, SampleNames, qval.thresh=0.05){
  print("Preparation...")
  pathway.list.up = lapply(pathway_rslt.list, function(r) {
    r %>% as.data.frame() %>% filter(qvalue<qval.thresh) %>%
      filter(NES>0) %>% .[, "ID", drop=T] 
  }) %>%
    setNames(paste0(SampleNames, "_UP"))
  pathway.list.down = lapply(pathway_rslt.list, function(r) {
    r %>% as.data.frame() %>% filter(qvalue<qval.thresh) %>%
      filter(NES<0) %>% .[, "ID", drop=T] 
  }) %>%
    setNames(paste0(SampleNames, "_DOWN"))
    print("Draw upset plot...")   
    #library(ComplexHeatmap)
    m = make_comb_mat(c(pathway.list.up, pathway.list.down))
    p=UpSet(m,
            top_annotation = upset_top_annotation(m, add_numbers = TRUE),
            right_annotation = upset_right_annotation(m, add_numbers = TRUE),
            comb_order = rev(order(comb_size(m))),
            set_order = c(paste0(SampleNames, "_UP"),
                          paste0(SampleNames, "_DOWN")),
            column_title = paste0("Significant Pathways\n(qvalue<",
                                  qval.thresh, ")"))
    
    return(list(plot=p, pathway.list.up = pathway.list.up, 
                pathway.list.down=pathway.list.down))  
}