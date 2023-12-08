library("GO.db")
library(tidyr)

symbolTable <- snakemake@input[["go_table"]]
goTerms <- read.table(symbolTable, sep="\t", header=TRUE)


checkfct <- function(term) {
  if (term %in% keys(GOCCANCESTOR)) {
    return(GOCCANCESTOR)
  } else if (term %in% keys(GOMFANCESTOR)) {
    return(GOMFANCESTOR)
  } else if(term %in% keys(GOBPANCESTOR)) {
      return(GOBPANCESTOR)
  } else { return(1)}
}

parentfct2 <- function(term) {
  anc <- checkfct(term)
  if (is.numeric(anc)) {
       return(term)

  }
  terms <- unlist(as.list(anc[term]), use.names = FALSE)
  terms <- terms[-length(terms)]
  terms <- append(terms, term)
  terms <- paste(terms, collapse=",")
  return(terms)
}
print(head(goTerms))
print(dim(goTerms))
goTerms$allTerms <- sapply(goTerms$GOTerm, parentfct2)
print(goTerms$allTerms[1])
goTerms <- separate_rows(goTerms, allTerms, sep=",")
goTerms <- goTerms[c("old_locus_tag", "GOTerm")]
write.table(goTerms , file = snakemake@output[["all_terms"]], sep="\t")

print(goTerms)
