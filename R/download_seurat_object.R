#' Download Seurat object based on object name
#'
#' @description
#' Before running this function,please run.\cr
#' SeuratData=bsabsfinder::show_singlecell_data()\cr\cr
#' This function downloads the Seurat object corresponding to the given SEURAT_OBJECT from show_singlecell_data().
#'
#' @param object_name A character string specifying the Seurat object name. Case-sensitive. Just use as given in 'SEURAT_OBJECT' of show_singlecell_data().
#' @param save_dir By default it will be saved in working directory.
#
#' @return Seurat object will be downloaded in working directory and also loaded in environment.
#' @usage
#'
#' download_seurat_object(object_name = "HCC_GSE151530",save_dir = getwd())
#'
#'
#' @export
download_seurat_object <- function(object_name, save_dir = getwd()) {

  # Check if the object already exists in the environment
  if (exists(object_name, inherits = FALSE) && class(get(object_name)) == "Seurat") {
    cat("Seurat object", object_name, "is already loaded in the environment. No need to download.\n")
    #return(get(object_name))
  }

  # Construct the file path for the local RDS file
  local_file <- file.path(save_dir, paste0(object_name, ".RDS"))

  #cat("Local file path:", local_file, "\n")

  # Check if the RDS file exists in the given directory
  if (file.exists(local_file)) {
    cat("Local RDS file found. Attempting to read Seurat object from the local file.\n")
    so <- readRDS(local_file)
    assign(object_name, so, envir = .GlobalEnv)
    return(so)
  }

  # Read the CSV file containing Seurat object information
  seurat_data <- BSABS_SINGLE_CELL_DATA

  # Check if the specified object name exists in the data
  if (object_name %in% seurat_data$SEURAT_OBJECT) {
    # Get the corresponding URL for the specified object name
    url <- seurat_data$DOWNLOAD_LINK[seurat_data$SEURAT_OBJECT == object_name]

    # Print the URL for debugging
    cat("Attempting to download Seurat object", object_name, "from URL:", url, "\n")

    # Download the Seurat object from the URL
    tryCatch({
      # Download the file from URL
      download.file(url, local_file, method = "auto", quiet = TRUE)

      # Read the downloaded RDS file
      so <- readRDS(local_file)

      # Load the Seurat object into the current environment with the specified name
      assign(object_name, so, envir = .GlobalEnv)

      cat("Seurat object", object_name, "downloaded successfully.\n")
      return(so)
    }, error = function(e) {
      cat("Error downloading Seurat object", object_name, ":", conditionMessage(e), "\n")
    })
  } else {
    cat("Error: Seurat object", object_name, "not found in the list.\n")
    return(NULL)
  }
}


