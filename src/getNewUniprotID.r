protein_list="A0A024R6I7"
getNewUniprotID = function(protein_list) {
  library(UniProt.ws)
  uniprotID = sapply(protein_list, function(protein) {
    #print(protein)
    a = suppressMessages(rba_uniprot_uniparc_search(accession = protein))
    j=NULL
    for (uniparc in unique(names(a))) {
      for (i in 1:length(a[[uniparc]]$dbReference)) {
        if (length(j)==0) {
          is_Swiss_Prot_present = function(a) {
            tryCatch(a[[uniparc]]$dbReference[[i]]$type=="UniProtKB/Swiss-Prot",
                     error = function(e) {
                       FALSE
                     })
          }
          is_Swiss_Prot_presentr = is_Swiss_Prot_present(a)
          if (length(is_Swiss_Prot_presentr)==0) {is_Swiss_Prot_presentr=FALSE}
          
          
          is_Swiss_Prot_iso_present = function(a) {
            tryCatch(a[[uniparc]]$dbReference[[i]]$type=="UniProtKB/Swiss-Prot protein isoforms",
                     error = function(e) {
                       FALSE
                     })
          }
          is_Swiss_Prot_iso_presentr = is_Swiss_Prot_iso_present(a)
          if (length(is_Swiss_Prot_iso_presentr)==0) {is_Swiss_Prot_iso_presentr=FALSE}
          
          
          is_human1 = function(a){
            tryCatch(a[[uniparc]]$dbReference[[i]]$property[[1]]$value=="9606",
                     error = function(e) {
                       FALSE
                     })      
          }
          is_human1r = is_human1(a)
          if (length(is_human1r)==0) {is_human1r=FALSE}
          
          is_human2 = function(a){
            tryCatch(a[[uniparc]]$dbReference[[i]]$property[[2]]$value=="9606",
                     error = function(e) {
                       FALSE
                     })      
          }
          is_human2r = is_human2(a)
          if (length(is_human2r)==0) {is_human2r=F}
          
          is_TrEMBL_present = function(a) {
            tryCatch(a[[uniparc]]$dbReference[[i]]$type=="UniProtKB/TrEMBL",
                     error = function(e) {
                       FALSE
                     })
          }
          is_TrEMBL_presentr = is_TrEMBL_present(a)
          if (length(is_TrEMBL_presentr)==0) {is_TrEMBL_presentr=F}
          
          
          if (is_Swiss_Prot_presentr) {
            #print(i)
            if (is_human1r) {
              j = a[[uniparc]]$dbReference[[i]]$id
              isValid = queryUniProt(
                query = c(paste0("accession:", j), "organism_id:9606"),
                fields = c("accession", "id", "reviewed"),
                collapse = " AND "
              )
              if (nrow(isValid)==0) {j = NULL}
            } else if (is_human2r){
              j = a[[uniparc]]$dbReference[[i]]$id
              isValid = queryUniProt(
                query = c(paste0("accession:", j), "organism_id:9606"),
                fields = c("accession", "id", "reviewed"),
                collapse = " AND "
              )
              if (nrow(isValid)==0) {j = NULL}
            } #else {j = a[[names(a)[1]]]$dbReference[[i]]$id}
            if (length(j)>0) {break}
          } else if (is_Swiss_Prot_iso_presentr) {
            if (is_human1r) {
              j = a[[uniparc]]$dbReference[[i]]$id
              isValid = queryUniProt(
                query = c(paste0("accession:", j), "organism_id:9606"),
                fields = c("accession", "id", "reviewed"),
                collapse = " AND "
              )
              if (nrow(isValid)==0) {j = NULL}
            } else if (is_human2r){
              j = a[[uniparc]]$dbReference[[i]]$id
              isValid = queryUniProt(
                query = c(paste0("accession:", j), "organism_id:9606"),
                fields = c("accession", "id", "reviewed"),
                collapse = " AND "
              )
              if (nrow(isValid)==0) {j = NULL}
            } 
            if (length(j)>0) {break}
          } else if (is_TrEMBL_presentr) {
            if (is_human1r) {
              j = a[[uniparc]]$dbReference[[i]]$id
              isValid = queryUniProt(
                query = c(paste0("accession:", j), "organism_id:9606"),
                fields = c("accession", "id", "reviewed"),
                collapse = " AND "
              )
              if (nrow(isValid)==0) {j = NULL}
              
            } else if (is_human2r){
              j = a[[uniparc]]$dbReference[[i]]$id
              isValid = queryUniProt(
                query = c(paste0("accession:", j), "organism_id:9606"),
                fields = c("accession", "id", "reviewed"),
                collapse = " AND "
              )
              if (nrow(isValid)==0) {j = NULL}
            } 
          } 
        }
      }
    }
    if (length(j)==0) {j=protein}
    j
  })
  return(uniprotID)
}
