#' @export
show_singlecell_data <- function() {
  #csv_path <- system.file("data", "BSABS_SINGLE_CELL_DATA.csv", package = "scBSABS")
  #csv_path <- "data/BSABS_SINGLE_CELL_DATA.csv"
  #seurat_data <- read.csv(csv_path)
  #seurat_data = BSABS_SINGLE_CELL_DATA
  cat("Available Single cell Seurat Data from TISCH:\n")
  cat("----------------------------------------------------\n")

  BSABS_SINGLE_CELL_DATA[,c(1:2)]

}
