# **Bispecific Antibody Targets Identification using Bulk RNA Sequencing Data**

The R package `bsAbsFinder` enlists the bispecific antibody marker pairs based on OCTAD bulk RNA sequencing data. 

## **Requirements**
You need to download following from data folder:

- [octad.counts.and.tpm.h5](https://chenlab-data-public.s3.amazonaws.com/octad/octad.counts.and.tpm.h5)

## **Installation**

```r

library(devtools)
install_github("Bin-Chen-Lab/octad.db")
install_github("Bin-Chen-Lab/octad")

devtools::install_github('shreya1704/bsAbsFinder')
library(bsabsfinder)

```
## **Input**

To use the package, provide the name of a particular cancer (case) and its normal tissue (control).

Following code snippet displays the number of samples for various cancers and normal tissues. Select from that list.
```r
total_cancer_count=bsabsfinder::total_cancer_count
total_normal_count=bsabsfinder::total_normal_count
```
Just replace the 'liver hepatocellular carcinoma' and 'LIVER' in below code with other cancer and its corresponding healthy tissue.<br><br>
The input is case sensitive so please add it accordingly:<br>
cancer.type as shown in total_cancer_count and 
normal.tissue as shown in total_normal_count.<br><br>
Also, specify the location of octad.counts.and.tpm.h5 and output file name.

```r
bulk_DE_surface_antigen(
    cancer.type = 'liver hepatocellular carcinoma',
    normal.tissue = 'LIVER',
    octad_counts_data_path = "~/Downloads/octad.counts.and.tpm.h5",
    output_file_name = "liver_output"
)
```

The function then extracts the expression matrix for those samples, performs DE analysis, and identifies marker pairs.

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

Please refer Example_Code Folder for all figures.

## **scBSABS** : Single Cell Data for Validation of identified marker pair genes
We have provided some single cell datasets to visulaize the expression patterns of identified marker genes. 

```
# To enlist the available single cell dataset

show_singlecell_data()


#  To download the single cell data. Object_name can be selected from available datasets

download_seurat_object(object_name = "HCC_GSE151530",save_dir = getwd())


# To visualize the expression pattern and calculate the percentage of cells expressing given pair 
  
  sc_bispecific_expression(so = HCC_GSE151530 ,ident = "Celltype_MajorLineage",geneA = "GPC3",geneB="MUC13")
  
  # For a given pair this function calculates the percentage of cells expressing gene A , gene B or both.
  # It also plots the Featureplot and Dotplot for given pair.
  # Percentage table for above call will be stored as "Coexpress_HCC_GSE151530_GPC3_MUC13.csv" and 
  # plots will be stored "Coexpress_GPC3_MUC13.pdf"

```

## **Manuscript Single Cell Figures : Data and Code**
- [FIG4](https://chenlab-data-public.s3.us-west-2.amazonaws.com/BISPECIFIC_ANTIBODY/Single_cell_data/FIG4_VITAL_ORGANS.zip)
- [FIG5](https://chenlab-data-public.s3.us-west-2.amazonaws.com/BISPECIFIC_ANTIBODY/Single_cell_data/FIG5_HCC.zip)
- [FIG6](https://chenlab-data-public.s3.us-west-2.amazonaws.com/BISPECIFIC_ANTIBODY/Single_cell_data/FIG6AB_COEXPRESSION_PLOT.zip)

## **Please cite**
Evgenii Chekalin, Shreya Paithankar, Rama Shankar, Jing Xing, Wenfeng Xu, Bin Chen, Computational discovery of co-expressed antigens as dual targeting candidates for cancer therapy through bulk, single-cell, and spatial transcriptomics, Bioinformatics Advances, Volume 4, Issue 1, 2024, vbae096, https://doi.org/10.1093/bioadv/vbae096

### For more transcriptomics-based tools, check out http://apps.octad.org/   
  











