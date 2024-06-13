suppressPackageStartupMessages({
  library(cowplot)
})
checksignificance = function(x) {
  x$sig = as.factor(x$adj.P.Val<0.2 & abs(x$logFC)>0)
  x$reg = factor("ns", levels = c("up", "down", "ns"))
  x$reg[x$adj.P.Val<0.2 & x$logFC>0] = "up"
  x$reg[x$adj.P.Val<0.2 & x$logFC<0] = "down"
  x
}

addLM = function(pl) {
  vars = gsub("~", "", as.character(pl$mapping))
  names(vars) = names(pl$mapping)
  pllm = lm(formula = pl$data[, vars["x"]]~pl$data[, vars["y"]])
  pl + geom_smooth(method="lm", formula="y~x") +
    annotate("text", x=Inf, y=-Inf, hjust=1, vjust=0,
             label=paste("r2=", signif(summary(pllm)$r.squared, 2)))
}

deplots = function(list1, list2, name1, name2, Chisq_test=T) {
  pcatheme = theme(legend.key.size = unit(x = 6, "pt"),
                   axis.title = element_text(size = 8),
                   axis.text = element_text(size = 8),
                   legend.title = element_text(size = 8),
                   legend.text = element_text(size = 8),
                   plot.title = element_text(size = 8))
  list1 = checksignificance(list1)
  list2 = checksignificance(list2)
  both = merge(x=list1, y=list2, by="Protein_IDs",
               suffixes = c("_1", "_2"))
  depoints = addLM(ggplot(both, aes(x=logFC_1, y=logFC_2))+
                     geom_point(aes(shape=sig_1, color=sig_2))) +
    labs(x=paste("log2FC_", name1), y=paste("log2FC_",name2),
         shape=paste("Sig_", name1), color=paste("Sig_", name2)) +
    pcatheme
  
  debars = ggplot(both[both$reg_2!="ns"|both$reg_1!="ns", ],
                  aes(x=reg_2, fill=reg_1)) + geom_bar()+ pcatheme +
    labs(x=paste("Reg_", name2), fill=paste("Reg_", name1))
  
  if (Chisq_test) {
    DE_chi = chisq.test(table(both[,paste0("reg_1")], both[, paste0("reg_2")]))
    
    if (DE_chi$p.value==0) {
      DE_chi$p.value = ">2.2e-16"
    } else {
      DE_chi$p.value = paste0("<", signif(DE_chi$p.value,2))
    }
    
    comparison = merge(as.data.frame(DE_chi$observed),
                       as.data.frame(as.table(DE_chi$expected)),
                       suffixes=c("Obs", "Exp"),
                       by=c("Var1", "Var2"))
    comparison$Var2 = factor(comparison$Var2, levels=c("ns", "down", "up"))
    colnames(comparison) = c(name1, name2, colnames(comparison)[3:ncol(comparison)])
    comparison$RoundFreqExp = round(comparison$FreqExp)
    comparison$Label = paste0(comparison$FreqObs,
                              "\n(",comparison$RoundFreqExp, ")")
    comparison$Difference = comparison$FreqObs-comparison$RoundFreqExp
    
    dechi = ggplot(comparison, aes_string(x=name1, y=name2)) +
      geom_tile(aes(fill=Difference)) +
      geom_text(aes(label=Label), size=3) +
      scale_fill_gradient2(high="#c87b7b", low="#7b81c8", midpoint=0) +
      labs(title=paste("Chi square test \n p", DE_chi$p.value),
           fill="Difference \n obs-exp") +
      pcatheme + 
      theme(legend.position = "bottom",
            legend.margin=margin(t=0, b=0, l=0, r=0, unit="cm"),
            legend.key.width=unit(15, "pt"))
    
    return(plot_grid(
      textGrob(paste0(name1, " and ", name2)),
      plot_grid(depoints,
                plot_grid(ncol=1, rel_heights=c(3,1),
                          debars+theme(legend.position="none"),
                          get_legend(debars)),
                dechi, ncol=3, rel_widths=c(2,0.7,1.2)), ncol=1,
      rel_heights = c(1,5)))     
  } else {
    return(plot_grid(
      textGrob(paste0(name1, " and ", name2)),
      plot_grid(depoints,
                plot_grid(ncol=1, rel_heights=c(3,1),
                          debars+theme(legend.position="none"),
                          get_legend(debars)),
                NULL, ncol=3, rel_widths=c(2,0.45, 1)), ncol=1,
      rel_heights = c(1,5))) 
  }
}