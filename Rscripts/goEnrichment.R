suppressPackageStartupMessages({
  library(clusterProfiler)
})
package <- list.files(snakemake@input[["annotation_db"]])[1]
library(basename(package), character.only = TRUE)

ranked_file <- snakemake@input[["ranked_file"]]
ranked_table <- read.table(ranked_file, sep="\t", header=TRUE)
rownames(ranked_table) <- ranked_table$old_locus_tag
ranked_table <- na.omit(ranked_table)

up <- ranked_table$Rank <= 200
up <- ranked_table[up, ]



minGSSize <- snakemake@config[["minGSSize"]]
maxGSSize <- snakemake@config[["maxGSSize"]]

egoBPup <- enrichGO(gene = as.character(rownames(up)),
                    universe = as.character(rownames(ranked_table)),
                    OrgDb = basename(package),
                    keyType = "GID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    readable = FALSE)
summary <- data.frame(egoBPup )
if (dim(summary)[1] == 0){
  df <- data.frame(matrix(ncol = 10, nrow = 0))
  x <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  colnames(df) <- x
  summary <- df
}
write.table(summary , file = snakemake@output[["enriched"]], row.names=FALSE, sep="\t")


