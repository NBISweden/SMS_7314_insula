convertCoreEnrichment2Symbol = function(GSEA_rslt) {
  core_enrichment_geneSymbol = sapply(GSEA_rslt$core_enrichment, function(i) {
    symbol = strsplit(i, "/") %>% as.data.frame() %>% 
      setNames(., c("CoreEnrichment")) %>%
      mutate(symbol = mapIds(org.Hs.eg.db, keys = CoreEnrichment, 
                    column = "SYMBOL", keytype = "ENTREZID", 
                    multiVals = "first"))
    symbol$symbol %>% as.character() %>%
      paste(., collapse = "/")
  })
  mutate(GSEA_rslt, core_enrichment_geneSymbol = core_enrichment_geneSymbol)
}