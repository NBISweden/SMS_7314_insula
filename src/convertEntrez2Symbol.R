#' @param Entrez e.g. "5581/989/375"
#' @param mapping_list e.g.:       
#'  6843        2162       93487 
#' "VAMP1"     "F13A1" "MAPK1IP1L" 
#
convertEntrez2Symbol = function(Entrez, mapping_list=NULL) {
  suppressMessages(library(AnnotationDbi))
  suppressMessages(library(org.Hs.eg.db))
  geneSymbol = sapply(Entrez, function(k) {
    symbol = strsplit(k, "/") %>% as.data.frame() %>% 
      setNames(., c("Entrez0"))
    #symbol
    # Entrez0
    # 1    7416
    # 2    1327
    if (length(mapping_list)==0) {
       symbol = symbol %>%
        mutate(symbol = AnnotationDbi::mapIds(org.Hs.eg.db, keys = Entrez0, 
                               column = "SYMBOL", keytype = "ENTREZID", 
                               multiVals = "first"))
      symbol$symbol %>% as.character() %>%
        paste(., collapse = "/")     
    } else {
      new_mapping_list = setNames(names(mapping_list), 
                                  as.character(mapping_list))
      new_mapping_list[symbol$Entrez] %>% as.character() %>%
        paste(., collapse = "/") 
    }
  })
  if (length(geneSymbol) != length(Entrez)) {
    print("Be careful some Entrez does not have gene gene symbol...")
  }
  return(geneSymbol)
}