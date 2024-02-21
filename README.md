# **Bispecific Antibody Markers Identification using Bulk RNA Sequencing Data**

The R package `bsAbsFinder` enlists the bispecific antibody marker pairs based on OCTAD bulk RNA sequencing data. 
Kindly refer provided Rmarkdown for code. 

## **Requirements**
You need to download following from data folder:
- basabsfinder_0.0.0.9001.tar.gz 
- octad.counts.and.tpm.h5

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

To use the package, provide the name of a particular cancer (case) and its normal tissue (control). The code then extracts the expression matrix for those samples, performs DE analysis, and identifies marker pairs.

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

## **Example : Bispecific Antibody Identification for HCC**

```
HCC_primary=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&sample.type == 'primary') #select HCC data
case_id=HCC_primary$sample.id #select cases
Healthy=subset(phenoDF,sample.type=='normal'&biopsy.site=='LIVER')#select Normal liver data
control_id=Healthy$sample.id

cases=loadOctadCounts(case_id,type='tpm',file='~/Downloads/octad.counts.and.tpm.h5')
cases=as.data.frame(cases)
controls=loadOctadCounts(control_id ,type='tpm',file='~/Downloads/octad.counts.and.tpm.h5') #Windows
controls=as.data.frame(controls)

#final data
hcc_with_liver=cbind(cases,controls)

#convert ensg to hgnc and select surface-expressed genes according to  compartments.jensenlab.org
hcc_with_liver=ensg_to_hgnc(hcc_with_liver,select_surface=TRUE)
hcc_with_liver_2=ensg_to_hgnc(hcc_with_liver,select_surface=FALSE)

#create phenotype vector
phenotype_vector=as.factor(c(rep('case',ncol(cases)),rep('control',ncol(controls))))

################################
#perform DE to filter out non-significant genes to speed up the computation
################################
annotation=data.frame(sample=c(colnames(cases),colnames(controls)),phenotype=c(rep('cancer',length(colnames(cases))),rep('control',length(colnames(controls)))))
annotation$phenotype=as.factor(annotation$phenotype)
expression=DGEList(counts=as.matrix(hcc_with_liver),group=annotation$phenotype)
dim(expression)
keep <- rowSums(cpm(expression)>100) >= 2
expression <- expression[keep,]
dim(expression)
expression$samples$lib.size <- colSums(expression$counts)
#expression$samples
expression<- calcNormFactors(expression)
expression_disp <- estimateCommonDisp(expression, verbose=T)
expression_disp <- estimateTagwiseDisp(expression_disp)
DE <- exactTest(expression_disp, pair=c(1,2)) # compare groups 1 and 2
DE=DE$table
DE$padj=p.adjust(DE$PValue,method='BH')
DE=subset(DE,padj<0.01&abs(logFC)>1.3)#Eugene's original cutoff
DE=subset(DE,padj<0.01&abs(logFC)>1)#Shreya's cut off to include CD3 as per reviewer

head(DE) #list of DEs
#filter out only surface-expressed DE genes. Just to speed up. 
hcc_with_liver=hcc_with_liver[row.names(hcc_with_liver)%in%row.names(DE),]

dataframe_for_computation=as.data.frame(t(hcc_with_liver)) #plug, will fix asap
small_res=compute_bsabs(antigene_1=colnames(dataframe_for_computation),data_input=dataframe_for_computation,pheno_input=phenotype_vector)
head(small_res)

p=plot_bsabs(small_res,label='case',pval_cut_off=0.01,pair_score_cut_off=quantile(small_res$pair_score,.99))

backup_res=small_res

#subset result table to filter out only top pairs:
small_res=subset(small_res,pair_score>quantile(small_res$pair_score,.99)&case_greater=='TRUE_TRUE'&p.adj<0.01)

#order and filter top-20
small_res=small_res[order(small_res$pair_score,decreasing = T),][1:20,]
marker_list=unique(c(small_res$antigen_1,small_res$antigen_2))

healthy_tissues=subset(phenoDF,sample.type=='normal')  
healthy_tissues=subset(healthy_tissues,grepl('BRAIN',biopsy.site)|biopsy.site=='LIVER'|biopsy.site=='LUNG'|grepl('HEART',biopsy.site)| grepl('KIDNEY',biopsy.site))

healthy_tissues <- healthy_tissues %>%mutate(biopsy.site = ifelse(grepl("BRAIN", biopsy.site), "BRAIN", biopsy.site))
healthy_tissues <- healthy_tissues %>%mutate(biopsy.site = ifelse(grepl("HEART", biopsy.site), "HEART", biopsy.site))
healthy_tissues <- healthy_tissues %>%mutate(biopsy.site = ifelse(grepl("KIDNEY", biopsy.site), "KIDNEY", biopsy.site))

healthy_tissues_expr=loadOctadCounts(healthy_tissues$sample.id ,type='tpm',file='~/Downloads/octad.counts.and.tpm.h5') #Windows
#healthy_tissues=loadOctadCounts(healthy_tissues$sample.id ,type='tpm',file='~/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')#Unix
healthy_tissues_expr=as.data.frame(healthy_tissues_expr)

#annotate & select only genes from the result table
healthy_tissues_expr=ensg_to_hgnc(healthy_tissues_expr,select_surface=FALSE)
healthy_tissues_expr=healthy_tissues_expr[row.names(healthy_tissues_expr)%in%marker_list,]

healthy_tissues_expr=healthy_tissues_expr[order(rownames(healthy_tissues_expr)),]

healthy_tissues_expr=as.data.frame(t(healthy_tissues_expr))

healthy_tissues_expr$Sample=healthy_tissues$biopsy.site[match(rownames(healthy_tissues_expr),healthy_tissues$sample.id)]

hcc_with_liver=hcc_with_liver[row.names(hcc_with_liver)%in%marker_list,]

hcc=hcc_with_liver[,colnames(hcc_with_liver)%in%case_id]

hcc=hcc[order(rownames(hcc)),]

hcc=as.data.frame(t(hcc))

hcc$Sample="HCC"

colnames(hcc)==colnames(healthy_tissues_expr)

marker_expr=rbind(hcc,healthy_tissues_expr)

table(marker_expr$Sample)

sample_order <- c("HCC", "LIVER", "BRAIN", "HEART", "LUNG","KIDNEY")
marker_expr$Sample <- factor(marker_expr$Sample, levels = sample_order)

pdf("~/Downloads/BSAB_marker.pdf")
{
  for (i in 1:18) {
    gg <- ggplot(marker_expr, aes(x = Sample, y = marker_expr[[i]], fill = Sample)) +
      geom_boxplot() +
      scale_fill_manual(values = c("red", "blue",rep("gray", length(sample_order) - 2))) +
      labs(x = "Biopsy Site", y = names(marker_expr)[i], title = paste("Boxplot of", names(marker_expr)[i], "by Biopsy Site")) +
      theme_minimal() +
      guides(fill = FALSE)
    
    print(gg)
      
  }
}

dev.off()

###############################################################

#Validation purpose

hcc_with_liver_2=ensg_to_hgnc(hcc_with_liver,select_surface=FALSE)

exp=hcc_with_liver_2[rownames(hcc_with_liver_2)%in%c('GPC3','MUC13','CD3D','CD3E','CD3G','MSLN'),]

exp=exp[order(rownames(exp)),]

exp=as.data.frame(t(exp))

exp$Sample=ifelse(rownames(exp)%in%Healthy$sample.id,"Normal_Liver","HCC")

marker_list=colnames(exp)[1:6]

pdf("~/Downloads/BSAB_HCC_validation.pdf")
{
  for (i in 1:length(marker_list)) {
    gg <- ggplot(exp, aes(x = Sample, y = exp[[i]], fill = Sample)) +
      geom_boxplot() +
      scale_fill_manual(values = c("red", "blue",guide = "none")) +
      labs(x = "Biopsy Site", y = names(expr)[i], title = paste("Marker expression by biopsy site: ", names(exp)[i])) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      #guides(fill = none)+
      stat_compare_means(method = "t.test",label="p.signif") #ref.group = "HCC"
     
    print(gg)
    
  }
}

dev.off()
```
