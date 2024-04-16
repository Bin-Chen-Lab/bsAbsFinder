#' Download Seurat object based on object name
#'
#' This function downloads the Seurat object corresponding to the given object name.
#'
#' @param object_name A character string specifying the Seurat object name.
#' @param csv_path A character string specifying the path to the CSV file containing Seurat object information.
#' @param save_dir A character string specifying the directory to save the downloaded Seurat object.
#' @importFrom utils tryCatch
#' @export

download_seurat_object <- function(object_name, csv_path = "data/BSABS_SINGLE_CELL_DATA.csv", save_dir = getwd()) {
  #cat("Function started with object name:", object_name, "\n")
  
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
  seurat_data <- read.csv(csv_path)
  
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

#download_seurat_object("CESC_GSE168652", save_dir = "/mnt/ufs18/rs-013/chenlab/Bispecific/")
