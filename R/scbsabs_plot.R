library(ggplot2)
library(dplyr)
library(Seurat)
#library(ggthemes)
library(rlang)


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




