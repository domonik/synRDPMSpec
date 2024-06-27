suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
})


ranked_file <- snakemake@input[["file"]]
ranked_table <- read.table(ranked_file, sep="\t", header=TRUE)
ranked_table <- na.omit(ranked_table)
print(head(ranked_table))
my_list <- list()
experiment_str <- snakemake@wildcards[["experiment"]]
if (snakemake@wildcards[["experiment"]] %in% c("egf_8min", "egf_2min")) {
  organism <- "hsa"
} else {
  organism <- "mmu"
}

for (nr in 1:6) {
  up <- ranked_table$"ANOSIM.R" >= 0.5 & grepl(paste0("FR", nr), ranked_table$"position.strongest.shift")
  up <- ranked_table[up,]
  entrezid_entries <- head(select(org.Mm.eg.db, keys = keys(org.Mm.eg.db), columns = "UNIPROT"))
  overlap_entrezid <- any(as.character(ranked_table$"PG.ProteinGroups") %in% keys(org.Mm.eg.db, keytype = "UNIPROT"))

  ekeggup <- enrichKEGG(gene = as.character(up$"PG.ProteinGroups"),
                        universe = as.character(ranked_table$"PG.ProteinGroups"),
                        organism = organism,
                        keyType="uniprot",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
  )
  summary <- data.frame(ekeggup)
  if (dim(summary)[1] == 0) {
    df <- data.frame(matrix(ncol = 9, nrow = 0))
    x <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
    colnames(df) <- x
    summary <- df
  }
  my_list[[nr]] <- summary
}

up <- ranked_table$"ANOSIM.R" >= 0.5
up <- ranked_table[up,]

ekeggup <- enrichKEGG(gene = as.character(up$"PG.ProteinGroups"),
                      universe = as.character(ranked_table$"PG.ProteinGroups"),
                      organism = organism,
                      keyType = "uniprot",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
)
summary <- data.frame(ekeggup)
if (dim(summary)[1] == 0) {
  df <- data.frame(matrix(ncol = 9, nrow = 0))
  x <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  colnames(df) <- x
  summary <- df
}



dir.create(snakemake@output[["enriched"]])
file_path <- file.path(snakemake@output[["enriched"]], "AllFractions.tsv")
write.table(summary, file = file_path, row.names = FALSE, sep = "\t") # write for all fractions
for (nr in 1:6) {
  summary <- my_list[[nr]]
  file_path <- file.path(snakemake@output[["enriched"]], paste0("Fraction", nr, ".tsv"))
  write.table(summary, file = file_path, row.names = FALSE, sep = "\t")


}

