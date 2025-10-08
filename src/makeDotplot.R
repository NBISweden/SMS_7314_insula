#' @param rslt pathway analysis result list
#
makeDotplot_ORA_KEGG = function(rslt, index.list=1:2, 
                       p.adjust.thresh=0.1) {
  
  print("Combine the ORA result to a long table...")
  if (is_empty(index.list)) {
    rslt1 = rslt
  } else {
    rslt1 = rslt[index.list]
  }
  table.long = lapply(seq_along(rslt1), function(i) {
    up = rslt1[[i]]$up@result %>% mutate(Regulation = "UP",
                                    Count_reg=Count,
                                    Condition=names(rslt1)[i]) 
    down = rslt1[[i]]$down@result %>% mutate(Regulation = "DOWN",
                                    Count_reg=-Count,
                                    Condition=names(rslt1)[i])
    r = rbind(up, down) %>% 
      dplyr::select(Condition, Regulation, Count_reg, p.adjust, everything()) %>%
      filter(p.adjust<p.adjust.thresh)
    r
  }) %>% Reduce(full_join, .)
  
  print("Sort the pathways to shared and unique...")
  pathway.list_up = lapply(rslt1, function(r) {
    r$up@result %>% filter(p.adjust<p.adjust.thresh) %>% 
      .[, "Description"]
  })
  pathway.list_down = lapply(rslt1, function(r) {
    r$down@result %>% filter(p.adjust<p.adjust.thresh) %>% 
      .[, "Description"]
  })
  pathway.list = mapply(function(x,y) c(x,y),
                        pathway.list_up, pathway.list_down) %>%
    setNames(., names(pathway.list_up))
  pathways.sorted = c(Reduce(intersect, pathway.list_up), 
                      Reduce(intersect, pathway.list_down),
                      Reduce(intersect, pathway.list),
                      setdiff(pathway.list[[1]], 
                              pathway.list[[2]]),
                      table.long$Description) %>%
    unique()
  
  table.long$Description = factor(table.long$Description,
                              levels = rev(pathways.sorted))
  
  p = ggplot(data=table.long, mapping=aes_string(x="Count_reg", y="Description")) +
    geom_point(shape=21, aes(size=p.adjust, fill=Regulation)) +
    facet_wrap(~Condition, ncol=3) + 
    # ggtitle("ORA_KEGG") +
    xlab("Count") +
    ylab("") + labs(size="p.adjust") +
    scale_fill_manual(values=c("DOWN"="#6D6D6D", "UP"="#B32B16")) +
    theme_bw(base_size = 15) + 
    scale_size(trans="reverse") +
    guides(alpha="none", size=guide_legend(override.aes=list(shape=21))) +
    theme(strip.text.x = element_text(size = 14)) +
    labs(caption=paste0("p.adjust.cutoff = ", p.adjust.thresh, ", pAdjustMethod = fdr"))
  return(p)
}