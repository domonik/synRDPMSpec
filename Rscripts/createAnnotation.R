
library(AnnotationForge)
library(tidyr)

createAnnotation <- function(goTable, symbolTable, outDir, genus="Bulbasaurus", species="phylloxyron"){
    print("FOO")
    mapping_data <- read.table(goTable, sep="\t", header=TRUE)

    if(!("EVIDENCE" %in% colnames(mapping_data))){
        mapping_data$EVIDENCE <- "experimental"
    }
    colnames(mapping_data) <- c("GID", "GO", "EVIDENCE")
    symbol_data <- read.table(symbolTable, sep="\t", header=TRUE)
    colnames(symbol_data) <- c("GID", "SYMBOL")
    dir.create(outDir, recursive=TRUE)
    rv <- makeOrgPackage(
        gene_info=symbol_data,
        go=mapping_data,
        version="0.1",
        maintainer="DR <schdruzzi@gmail.com>",
        author="DR <schdruzzi@gmail.com>",
        outputDir = snakemake@output[["annotation_db"]],
        tax_id = "42",
        genus = genus,
        species=species,
        goTable = "go"
    )
    return(rv)
}

print(snakemake@input[["go_terms"]])
print(snakemake@input[["symbols"]])
print(snakemake@output[["annotation_db"]])
print(snakemake@params[["species"]])
print(snakemake@params[["genus"]])
package <- createAnnotation(
    goTable=snakemake@input[["go_terms"]],
    symbolTable=snakemake@input[["symbols"]],
    outDir=snakemake@output[["annotation_db"]],
    genus=snakemake@params[["genus"]],
    species=snakemake@params[["species"]]
)
install.packages(package, repos=NULL)
library(basename(package), character.only = TRUE)
cat(NULL, file=snakemake@output[["finished_file"]])