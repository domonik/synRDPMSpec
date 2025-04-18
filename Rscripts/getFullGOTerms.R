library("GO.db")
library(tidyr)
library(parallel)

symbolTable <- snakemake@input[["go_table"]]
goTerms <- read.table(symbolTable, sep="\t", header=TRUE)
numCores <- snakemake@threads

# Serialize keys into lists for parallel use
GOCCANCESTOR_list <- as.list(GOCCANCESTOR)
GOMFANCESTOR_list <- as.list(GOMFANCESTOR)
GOBPANCESTOR_list <- as.list(GOBPANCESTOR)

checkfct <- function(term, GOCC, GOMF, GOBP) {
  if (term %in% names(GOCC)) {
    return(GOCC)
  } else if (term %in% names(GOMF)) {
    return(GOMF)
  } else if (term %in% names(GOBP)) {
    return(GOBP)
  } else {
    return(1)
  }
}

parentfct2 <- function(term, GOCC, GOMF, GOBP) {
  anc <- checkfct(term, GOCC, GOMF, GOBP)
  if (is.numeric(anc)) {
    return(term)
  }
  terms <- unlist(anc[[term]], use.names = FALSE)
  terms <- terms[-length(terms)]
  terms <- append(terms, term)
  terms <- paste(terms, collapse = ",")
  return(terms)
}


cl <- makeCluster(numCores)
clusterExport(cl, varlist = c("checkfct", "parentfct2", "GOCCANCESTOR_list", "GOMFANCESTOR_list", "GOBPANCESTOR_list"))
clusterEvalQ(cl, library(parallel))  # Add any required libraries here
clusterEvalQ(cl, library("GO.db"))  # Add any required libraries here

# Parallelized computation
goTerms$allTerms <- parSapply(cl, goTerms$GOTerm, function(term) {
  parentfct2(term, GOCCANCESTOR_list, GOMFANCESTOR_list, GOBPANCESTOR_list)
})
print(goTerms$allTerms[1])
goTerms <- separate_rows(goTerms, allTerms, sep=",")
goTerms <- goTerms[c("old_locus_tag", "GOTerm")]
write.table(goTerms , file = snakemake@output[["all_terms"]], sep="\t")

print(goTerms)
