source("show_singlecell_data.R")
source("download_seurat_object.R")
source("sc_bispecific_expression.R")


show_singlecell_data() # This function will show avialble single cell datasets 

create_directory <- function(main_dir_path, subdir_name) {
  subdir_path <- file.path(main_dir_path, subdir_name)
  
  # Check if directory exists, and create it if it doesn't
  if (!dir.exists(subdir_path)) {
    dir.create(subdir_path, recursive = TRUE)
  }
}

create_directory(getwd(),"PRAD_GSE141445")
create_directory(getwd(),"VITAL_ORGANS")

setwd("PRAD_GSE141445")

download_seurat_object(object_name = "PRAD_GSE141445",save_dir = getwd()) 

Idents(PRAD_GSE141445)=PRAD_GSE141445$Celltype_Malignancy
DimPlot(PRAD_GSE141445,label=T)
FeaturePlot(PRAD_GSE141445,features=c('OR51E2','ILDR1','GJB1','SLC27A2','HPN','MS4A8','GPR160','SLC9A2'))
DotPlot(PRAD_GSE141445,features=c('OR51E2','ILDR1','GJB1','SLC27A2','HPN','MS4A8','GPR160','SLC9A2'))+
  theme(axis.text.x=element_text(size=10, angle = 45), axis.text.y=element_text(size=10))


#sc_bispecific_expression()  will calculate the percentage of cells expressing gene A , gene B and both
#It also plots the featureplot and dotplot of given pair 

sc_bispecific_expression(so = PRAD_GSE141445 ,ident = "Celltype_Malignancy",geneA = "OR51E2",geneB="ILDR1")
sc_bispecific_expression(so = PRAD_GSE141445 ,ident = "Celltype_Malignancy",geneA = "HPN",geneB="GJB1")

pair1_prad=read.csv("Coexpress_PRAD_GSE141445_HPN_GJB1.csv")
pair2_prad=read.csv("Coexpress_PRAD_GSE141445_OR51E2_ILDR1.csv")


setwd("../VITAL_ORGANS")

download_seurat_object(object_name = "VITAL_ORGANS",save_dir = getwd())

Idents(VITAL_ORGANS)=VITAL_ORGANS$Organ
DimPlot(VITAL_ORGANS,label=T)
FeaturePlot(VITAL_ORGANS,c('OR51E2','ILDR1','GJB1','SLC27A2','HPN','MS4A8','GPR160','SLC9A2'))
sc_bispecific_expression(so = VITAL_ORGANS ,ident = "Organ",geneA = "HPN",geneB="GJB1")
sc_bispecific_expression(so = VITAL_ORGANS ,ident = "Organ",geneA = "OR51E2",geneB="ILDR1")

pair1_vo=read.csv("Coexpress_VITAL_ORGANS_HPN_GJB1.csv")
pair2_vo=read.csv("Coexpress_VITAL_ORGANS_OR51E2_ILDR1.csv")





