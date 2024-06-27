suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
})


ranked_file <- snakemake@input[["file"]]
ranked_table <- read.table(ranked_file, sep="\t", header=TRUE)
ranked_table <- na.omit(ranked_table)
print(head(ranked_table))
my_list <- list()
experiment_str <- snakemake@wildcards[["experiment"]]

if (grepl("egf", experiment_str) & grepl("min", experiment_str)) {
  organism <- org.Hs.eg.db
} else {
  organism <- org.Mm.eg.db
}
for (nr in 1:6) {
  up <- ranked_table$"ANOSIM.R" >= 0.5 & grepl(paste0("FR", nr), ranked_table$"position.strongest.shift")
  up <- ranked_table[up,]
  entrezid_entries <- head(select(org.Mm.eg.db, keys = keys(org.Mm.eg.db), columns = "UNIPROT"))
  overlap_entrezid <- any(as.character(ranked_table$"PG.ProteinGroups") %in% keys(org.Mm.eg.db, keytype = "UNIPROT"))

  egoBPup <- enrichGO(gene = as.character(up$"PG.ProteinGroups"),
                      universe = as.character(ranked_table$"PG.ProteinGroups"),
                      OrgDb = organism,
                      keyType = "UNIPROT",
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      minGSSize = 3,
                      maxGSSize = 5000,
                      readable = FALSE)
  summary <- data.frame(egoBPup)
  if (dim(summary)[1] == 0) {
    df <- data.frame(matrix(ncol = 10, nrow = 0))
    x <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
    colnames(df) <- x
    summary <- df
  }
  my_list[[nr]] <- summary
}
dir.create(snakemake@output[["enriched"]])
for (nr in 1:6) {
  summary <- my_list[[nr]]
  file_path <- file.path(snakemake@output[["enriched"]], paste0("Fraction", nr, ".tsv"))
  write.table(summary, file = file_path, row.names = FALSE, sep = "\t")


}
