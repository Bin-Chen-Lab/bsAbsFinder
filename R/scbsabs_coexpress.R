make_boolean_matrix <- function(so, ident) {
  counts=so@assays$RNA$counts
  rownames(counts)=rownames(so@assays$RNA)
  colnames(counts)=colnames(so@assays$RNA)
  subset_counts <- counts[1:nrow(counts), 1:ncol(counts)]
  binary_df <- as.data.frame(subset_counts > 0)
  binary_df <- as.data.frame(t(binary_df))
  binary_df <- binary_df[order(rownames(binary_df)), ]
  Ident_col <- as.data.frame(so@meta.data[, ident, drop = FALSE])
  Ident_col <- Ident_col[order(rownames(Ident_col)), , drop = FALSE]
  binary_df <- cbind(binary_df, Ident_col)
  colnames(binary_df)[ncol(binary_df)] <- "Ident"
  binary_df <- binary_df[!is.na(binary_df$Ident), ]
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
  # geneA_present <- geneA %in% colnames(binary_df)
  # geneB_present <- geneB %in% colnames(binary_df)
  # 
  # if (!geneA_present && !geneB_present) {
  #   #stop("Both genes are absent in the dataset.")
  #   message("Both genes ", geneA, " and ", geneB, " are absent in the dataset.")
  #   return(NULL)
  #   
  # } else {
  #   
  #   if (!geneA_present || !geneB_present) {
  #     missing_genes <- c()
  #     if (!geneA_present) missing_genes <- c(missing_genes, geneA)
  #     if (!geneB_present) missing_genes <- c(missing_genes, geneB)
  #     
  #     warning("The following gene(s) are absent in the dataset: ", paste(missing_genes, collapse = ", "))
  #   }
  #   
    result_df <- data.frame(
      Cell_Type = character(),
      Percent_Expressing_A = numeric(),
      Percent_Expressing_B = numeric(),
      Percent_Coexpressing_A_B = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (ctype in unique(binary_df$Ident)) {
      
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
  






