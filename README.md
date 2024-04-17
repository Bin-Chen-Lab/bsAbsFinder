# **Bispecific Antibody Targets Identification using Bulk RNA Sequencing Data**

The R package `bsAbsFinder` enlists the bispecific antibody marker pairs based on OCTAD bulk RNA sequencing data. 

Please refer Sample Code Folder for the code.

## **Requirements**
You need to download following from data folder:
- [basabsfinder_0.0.0.9001.tar.gz](https://chenlab-data-public.s3.amazonaws.com/BISPECIFIC_ANTIBODY/bsAbsFinder_installation/basabsfinder_0.0.0.9001.tar.gz)
- [octad.counts.and.tpm.h5](https://chenlab-data-public.s3.amazonaws.com/octad/octad.counts.and.tpm.h5)

## **Installation**

```r
install.packages('BiocManager')
BiocManager::install(c("DESeq2","edgeR","EDASeq","RUVSeq","EnsDb.Hsapiens.v86"))

library(devtools)
install_github("Bin-Chen-Lab/octad.db")
install_github("Bin-Chen-Lab/octad")

install.packages("~/Downloads/basabsfinder_0.0.0.9001.tar.gz", repos=NULL, type='source')
devtools::install_github('Lionir/bsAbsFinder')

install.packages(c('cluster','dplyr','ggplot2','ggpubr'))
```
## **Input**

To use the package, provide the name of a particular cancer (case) and its normal tissue (control).

```r
case=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&sample.type == 'primary') 
case_id=case$sample.id #select cases
control=subset(phenoDF,sample.type=='normal'&biopsy.site=='LIVER')
control_id=control$sample.id
```
Just replace the 'liver hepatocellular carcinoma' and 'LIVER' in above code with 
other cancer and its corresponding healthy tissue. 
Following code snippet displays the number of samples for various cancers and normal tissues. Select from that list.

```r
total_cancer_count=as.data.frame(table(phenoDF$cancer))
total_cancer_count=total_cancer_count[total_cancer_count$Var1!="normal",]
total_cancer_count=total_cancer_count[order(total_cancer_count$Freq,decreasing = T),]

total_normal=phenoDF[phenoDF$sample.type=="normal",c(1:3)]
total_normal_count=as.data.frame(table(total_normal$biopsy.site))
total_normal_count=total_normal_count[order(total_normal_count$Var1),]
```
The code then extracts the expression matrix for those samples, performs DE analysis, and identifies marker pairs.

## **Output**

The generated output consists of a candidate marker pair table, featuring pair scores and p-values. The pair score is derived from the marker genes' or antigens' ability to distinguish between case and control samples. Marker pairs include instances where:
- Both markers show high expression in case samples
- Both markers show high expression in control samples
- Either marker shows high expression in either case or control sample

Given the focus on identifying candidates for bispecific antibodies, interest lies specifically in pairs where both markers are highly expressed in case samples. Therefore, the output table includes a 'case_greater' column employing Boolean logic.

For further refinement, the generated table should be filtered based on desired cutoffs for pair score, p-value, and specifically selecting those pairs where the 'case_greater' column reads 'TRUE_TRUE,' indicating higher expression of both markers in case samples than control samples. Also, both markers should have zero to low expression in healthy vital organs.

## **Visualization**

Total 4 plots are generated (Manuscript Fig 2):
- Volcano plot for all marker pairs based on pair score
- Marker frequency among top pairs
- Expression pattern of top markers in case, control, and healthy vital organs
- Expression pattern of top pair in case, control, and healthy vital organs

## **scBSABS** : Single Cell Data for Validation of identified marker pair genes
We have provided some single cell datasets to visulaize the expression patterns of identified marker genes. 

```
# To enlist the available single cell dataset

show_singlecell_data()


#  To download the single cell data. Object_name can be selected from available datasets

download_seurat_object(object_name = "HCC_GSE151530",save_dir = getwd())


# To visualize the expression pattern and calculate the percentage of cells expressing given pair 
  
  sc_bispecific_expression(so = HCC_GSE151530 ,ident = "Celltype_Malignancy",geneA = "GPC3",geneB="MUC13")
  
  # For a given pair this function calculates the percetage of cells expressing gene A , gene B or both.
  # It also plots the Featureplot and Dotplot for given pair.
  # Percentage table for above call will be stored as "Coexpress_HCC_GSE151530_GPC3_MUC13.csv" and 
  # plots will be stored "Coexpress_GPC3_MUC13.pdf"

```
  

  
  











