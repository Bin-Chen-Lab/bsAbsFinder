#' @import ggplot2
#' @import dplyr
#' @import Seurat
#' @import rlang

plot_seurat_data <- function(so, ident, geneA, geneB,geneA_present,geneB_present) {
  # Subset Seurat object to remove cells with NA in the identity column
  if (ident %in% colnames(so@meta.data)) {
    so[[ident]][is.na(so[[ident]])] <- "Unknown"
    Idents(so) <- so@meta.data[,ident]  # Set the identity column
  } else {
    stop("Identity column '", ident, "' not found in metadata.")
  }
  # Step 3: Plotting
  #file=rlang::expr_text(substitute(so))
  #pdf_file <- paste0("Coexpress_",file,"_",geneA, "_", geneB, ".pdf")
  pdf_file <- paste0("Coexpress_",geneA, "_", geneB, ".pdf")
  #wd= getwd()
  pdf(pdf_file, width = 6, height = 6)
  P0 <- DimPlot(so, label = TRUE)
  print(P0)
  if (geneA_present) {
    P1 <- FeaturePlot(so, features = geneA)
    print(P1)
    D1 <- DotPlot(so, features = geneA) + RotatedAxis() + coord_flip() +
      theme(axis.text.x = element_text(size = 5, angle = 45),
            axis.text.y = element_text(size = 5),
            legend.title = element_text(size = 5),
            legend.text = element_text(size = 5))
    print(D1)
    #return(c(D1, P1))
  } else {
    P1 <- NULL  # Set to NULL if gene A is not present
    D1 <- NULL
  }

  if (geneB_present) {
    P2 <- FeaturePlot(so, features = geneB)
    print(P2)
    D2 <- DotPlot(so, features = geneB) + RotatedAxis() + coord_flip() +
      theme(axis.text.x = element_text(size = 5, angle = 45),
            axis.text.y = element_text(size = 5),
            legend.title = element_text(size = 5),
            legend.text = element_text(size = 5))
    print(D2)
  } else {
    P2 <- NULL  # Set to NULL if gene B is not present
    D2 <- NULL
  }

  if (geneA_present & geneB_present) {
    D3 <- DotPlot(so, features = c(geneA, geneB)) + RotatedAxis() + coord_flip() +
      theme(axis.text.x = element_text(size = 5, angle = 45),
            axis.text.y = element_text(size = 5),
            legend.title = element_text(size = 5),
            legend.text = element_text(size = 5))
    print(D3)
  } else {
    D3 <- NULL
  }

  dev.off()
  message("Plots saved to ", pdf_file, ".")
}

make_boolean_matrix <- function(so) {
  counts=so@assays$RNA$counts
  rownames(counts)=rownames(so@assays$RNA)
  colnames(counts)=colnames(so@assays$RNA)
  subset_counts <- counts[1:nrow(counts), 1:ncol(counts)]
  binary_df <- as.data.frame(subset_counts > 0)
  binary_df <- as.data.frame(t(binary_df))
  binary_df <- binary_df[order(rownames(binary_df)), ]
  return(binary_df)
}



calculate_percentages2 <- function(df, geneA, geneB)
{
  if (geneA %in% colnames(df)) {
    num_A <- sum(df[[geneA]] == TRUE)
    A_names <- if (num_A > 0) {
      rownames(df[df[[geneA]], , drop = FALSE])
    } else {
      NULL
    }

  }else
  {
    num_A=0
    num_AB=0
  }

  if (geneB %in% colnames(df)) {
    num_B <- sum(df[[geneB]] == TRUE)
    B_names <- if (num_B > 0) {
      rownames(df[df[[geneB]], , drop = FALSE])
    } else {
      NULL
    }

  }else
  {
    num_B=0
    num_AB=0
  }

  num_AB <- ifelse(num_A > 0 && num_B > 0, length(intersect(A_names, B_names)), 0)

  total <- nrow(df)

  percent_A <- round((num_A / total) * 100, 2)
  percent_B <- round((num_B / total) * 100, 2)
  percent_AB <- round((num_AB / total) * 100, 2)

  return(c(percent_A = percent_A, percent_B = percent_B, percent_AB = percent_AB))
}

coexpression_calculation<- function(binary_df, geneA, geneB)
{
  result_df <- data.frame(
    Cell_Type = character(),
    Percent_Expressing_A = numeric(),
    Percent_Expressing_B = numeric(),
    Percent_Coexpressing_A_B = numeric(),
    stringsAsFactors = FALSE
  )

  for (ctype in unique(binary_df)) {

    subset=binary_df[binary_df$Ident==ctype,]
    res=calculate_percentages2(subset, geneA, geneB)
    new_row <- data.frame(
      Cell_Type = ctype,
      Percent_Expressing_A = res[["percent_A"]],
      Percent_Expressing_B = res[["percent_B"]],
      Percent_Coexpressing_A_B = res[["percent_AB"]],
      stringsAsFactors = FALSE
    )

    # Bind the new row to result_df
    result_df <- rbind(result_df, new_row)

  }

  return(result_df)
}
#' @title sc_bispecific_expression
#'
#' @description
#' This function calculates the % of cells expressing given markers in seurat obect.
#'
#'
#' @param so It should be loaded SuratObject in enviornment
#' @param ident Identity column of the seurat obect
#' @param geneA character string specifying first marker of a pair
#' @param geneB character string specifying second marker of a pair
#'
#' @return dataframe showing the percentage of cell expressing marker A , marker B and both marker A and marker B
#' @return PDF containing diplot of seuratobject, featureplot showing marker A and marker B and dotplot showing expression of marker A and marker B
#'
#' @usage
#'
#'  sc_bispecific_expression(so = HCC_GSE151530 ,ident = "Celltype_MajorLineage",geneA = "GPC3",geneB="MUC13")
#'
#'
#' @export
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

if (!exists(variable_name, envir = .GlobalEnv)){
  print("generating boolean matrix")
  boolean_matrix <- make_boolean_matrix(so)
  assign(variable_name, boolean_matrix, envir = .GlobalEnv)
  print("Boolean matrix generated.....")
}

bool_df=get(variable_name)
Ident_col <- as.data.frame(so@meta.data[, ident, drop = FALSE])
Ident_col <- Ident_col[order(rownames(Ident_col)), , drop = FALSE]
Ident_col[,1]=factor(Ident_col[,1])
bool_df <- cbind(bool_df, Ident_col)
colnames(bool_df)[which(colnames(bool_df)==ident)] <- "Ident"
bool_df <- bool_df[!is.na(bool_df$Ident), ]
#return(bool_df)
message("Calculating percentage of cells expressing given gene(s).....")

result_df <- coexpression_calculation(bool_df, geneA, geneB)
colnames(result_df)=c(ident,paste0("percent", "_", geneA), paste0("percent", "_", geneB), paste0("percent", "_", geneA, "_", geneB))
output_name <- paste0("Coexpress_", so_name, "_", geneA, "_", geneB)
write.csv(result_df,paste0(output_name,".csv"),row.names = F,quote=F)

message("Coexpression matrix saved to ", output_name, ".csv.")

message("Generating plots to visualize expression pattern of gene(s).....")

plot_seurat_data(so, ident, geneA, geneB,geneA_present,geneB_present)

}










