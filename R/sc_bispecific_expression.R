source("scbsabs_coexpress.R")
source("scbsabs_plot.R")

library(ggplot2)
library(dplyr)
library(Seurat)

sc_bispecific_expression <- function(so, ident, geneA, geneB) {
  
# Step 1: Check gene presence
geneA_present <- geneA %in% rownames(so@assays$RNA$counts)
geneB_present <- geneB %in% rownames(so@assays$RNA$counts)

# If both genes are absent, stop and display message
if (!geneA_present && !geneB_present) {
  stop("Both genes are absent in the Seurat object.")
}

# Step 2: Check identity column
if (!ident %in% colnames(so@meta.data)) {
  stop("Identity column '", ident, "' not found in given surat object.")
}

if (!geneA_present || !geneB_present) {
  missing_genes <- c()
  if (!geneA_present) missing_genes <- c(missing_genes, geneA)
  if (!geneB_present) missing_genes <- c(missing_genes, geneB)
  
  warning("The given gene not found in the seurat dataset: ", paste(missing_genes, collapse = ", "))
}

#so_name <- deparse(substitute(so))
so_name=rlang::expr_text(substitute(so))
variable_name <- paste0(so_name, "_boolean")

#assign(variable_name, make_boolean_matrix(so, ident), envir = .GlobalEnv)

if (!exists(variable_name, envir = .GlobalEnv)){
  print("generating boolean matrix")
  boolean_matrix <- make_boolean_matrix(so, ident)
  assign(variable_name, boolean_matrix, envir = .GlobalEnv)
  #assign(variable_name, make_boolean_matrix(so, ident), envir = .GlobalEnv)
  print("Boolean matrix generated.....")
}


message("Calculating percentage of cells expressing given gene(s).....")

result_df <- coexpression_calculation(get(variable_name), geneA, geneB)
output_name <- paste0("Coexpress_", so_name, "_", geneA, "_", geneB)
#assign(output_name, result_df, envir = .GlobalEnv)
write.csv(result_df,paste0(output_name,".csv"),row.names = F,quote=F)

message("Coexpression matrix saved to ", output_name, ".csv.")

message("Generating plots to visualize expression pattern of gene(s).....")

plot_seurat_data(so, ident, geneA, geneB,geneA_present,geneB_present)

}










