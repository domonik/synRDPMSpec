library(BiocManager)
library(curl)

BiocManager::install(c("mzR", "SummarizedExperiment", "GenomeInfoDb", "GenomicRanges", "MSnbase", "QFeatures", "DAPAR", "PSMatch"))

library(DAPAR)
library(MSnbase)
library(Biobase)

# Load the file
file <- snakemake@input[["tsv"]]

# Read the data
df <- read.csv(file, sep = "\t", row.names = 1)
colnames(df) <- gsub("\\.", "_", colnames(df))
assay_data <- as.matrix(df)
feature_data <- data.frame(featureNames = rownames(df))
rownames(feature_data) <- rownames(df)


# Convert data to matrix

# Create phenotypic data with a 'Condition' column
conditions <- colnames(df)
split_parts <- strsplit(conditions, "\\_")
last_parts <- sapply(split_parts, function(x) tail(x, 1))# Use column names as conditions
modified_conditions <- sapply(split_parts, function(x) {
  # Extract head (all parts except the last one)
  head_part <- paste(head(x, -1), collapse = "_")
  # Remove the last character from the head
  if (nchar(head_part) > 0) {
    head_part <- substr(head_part, 1, nchar(head_part) - 1)
  }
  # Extract the last part
  tail_part <- tail(x, 1)
  # Combine modified head and tail
  paste(head_part, tail_part, sep = "_")
})

pheno_data <- data.frame(Condition = last_parts, Bio.rep = seq(1, length(last_parts)), Sample.name = colnames(df) )
rownames(pheno_data) <- colnames(df)
new_pheno_data <- data.frame(Condition = modified_conditions, Bio.rep = seq(1, length(last_parts)), Sample.name = colnames(df))
rownames(new_pheno_data) <- colnames(df)# Row names must match column names of `df`
df["Protein"] <- rownames(df)

# Create an ExpressionSet
msnset <- createMSnset(df, pheno_data, indExpData = seq(1, 120), colnameForID = "Protein",logData=F, pep_prot_data = "protein", software="maxquant")
level <- 'peptide'
pattern <- "Quantified"
type <- "AtLeastOneCond"
percent <- FALSE
op <- ">="
th <- 2
indices <- GetIndices_MetacellFiltering(msnset, level, pattern, type, percent, op, th)
bla <- MetaCellFiltering(msnset, indices, "keep")$new# Perform normalization using DAPAR
norm_data <- wrapper.normalizeD(obj = msnset, method = "vsn", conds = pheno_data$Condition, type="within conditions")# View normalized counts
norm_data_df <- as.data.frame(exprs(norm_data))
write.table(norm_data_df , file = snakemake@output[["normalized"]], sep="\t")


for (condition in unique(new_pheno_data$Condition)){
  pdata <- new_pheno_data[new_pheno_data$Condition == condition,]
  samples <- pdata$Sample.name
  sub_df <- norm_data_df[samples]
  sub_df["Protein"] <- rownames(sub_df)
  msnset_impute <- createMSnset(sub_df, pdata, indExpData = seq(1, 3), colnameForID = "Protein",logData=F, pep_prot_data = "protein", software="maxquant")
  level <- 'peptide'
  pattern <- "Quantified"
  type <- "AtLeastOneCond"
  percent <- FALSE
  op <- ">="
  th <- 2
  indices <- GetIndices_MetacellFiltering(msnset_impute, level, pattern, type, percent, op, th)
  msnset_impute_filtered <- MetaCellFiltering(msnset_impute, indices, "keep")$new
  imputed_data <- wrapper.impute.KNN(msnset_impute_filtered, 15)
  imputed_counts <- exprs(imputed_data)
  imputed_sub_df <- as.data.frame(imputed_counts)
  norm_data_df[rownames(imputed_sub_df), colnames(imputed_sub_df)] <- imputed_sub_df


}

write.table(norm_data_df , file = snakemake@output[["imputed"]], sep="\t")

