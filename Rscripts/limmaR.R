library(limma)


ranked_file <- snakemake@input[["tsv"]]

data <- read.table(ranked_file, sep="\t", header=TRUE)
print(head(data))

row.names(data) <- data$RAPDORid
data <- data[ , -1]
Conditions <- factor(c(rep("CTRL", 4), rep("EGF", 4)))


design <- model.matrix(~ 0 + Conditions)
print(head(data))
colnames(design) <- levels(Conditions)

# Contrast matrix
contrast.matrix <- makeContrasts(EGF - CTRL, levels=design)
fit <- lmFit(data, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

print("_______")
results <- topTable(fit2, adjust="BH", number=nrow(data))
write.table(results , file = snakemake@output[["tsv"]], row.names=TRUE, sep="\t")
