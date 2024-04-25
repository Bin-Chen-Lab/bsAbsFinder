#' @import edgeR
#' @import bsabsfinder
#' @import cluster
#' @import dplyr
#' @import ggplot2
#' @import ggpubr
#' @export



bulk_DE_surface_antigen=function(cancer.type, normal.tissue, octad_counts_data_path){

  case=subset(phenoDF,cancer==cancer.type & sample.type == 'primary')
  case_id=case$sample.id #select cases
  control=subset(phenoDF, biopsy.site== normal.tissue & sample.type=='normal')
  control_id=control$sample.id
case_expr=loadOctadCounts(case_id,type='tpm',file=octad_counts_data_path)
case_expr=as.data.frame(case_expr)
control_expr=loadOctadCounts(control_id ,type='tpm',file=octad_counts_data_path)
control_expr=as.data.frame(control_expr)

#final data
case_with_control_expr=cbind(case_expr,control_expr)

#convert ensg to hgnc and select surface-expressed genes according to  compartments.jensenlab.org
case_with_control_expr=ensg_to_hgnc(case_with_control_expr,select_surface=TRUE)


phenotype_vector= as.factor(c(rep('case',ncol(case_expr)),rep('control',ncol(control_expr))))

################################
#perform DE to filter out non-significant genes to speed up the computation
################################
annotation=data.frame(sample=c(colnames(case_expr),colnames(control_expr)),phenotype=c(rep('cancer',length(colnames(case_expr))),rep('control',length(colnames(control_expr)))))
annotation$phenotype=as.factor(annotation$phenotype)
expression=DGEList(counts=as.matrix(case_with_control_expr),group=annotation$phenotype)
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

DE1=subset(DE,padj<0.01&abs(logFC)>1) # The cutoff criteria can be changed

#filter out only surface-expressed DE genes. Just to speed up.
case_with_control_expr=case_with_control_expr[row.names(case_with_control_expr)%in%row.names(DE1),]
dataframe_for_computation=as.data.frame(t(case_with_control_expr))

#assign(dataframe_for_computation, as.data.frame(t(case_with_control_expr)), envir = .GlobalEnv)
#this takes for a while
small_res=compute_bsabs(antigene_1=colnames(dataframe_for_computation),data_input=dataframe_for_computation,pheno_input=phenotype_vector)
assign("Antigen_pair_table", small_res,  envir = .GlobalEnv)
write.csv(small_res, "Antigen_pairs_calculation.csv")


result_table=small_res
result_table=subset(result_table,pair_score>quantile(result_table$pair_score,.90))
result_table=subset(result_table,pair_score>2)
result_table=subset(result_table,case_greater=='TRUE_TRUE')
result_table=result_table[order(result_table$pair_score,decreasing = T),]
result_table$antigen_1=gsub('-','_',result_table$antigen_1)
result_table$antigen_2=gsub('-','_',result_table$antigen_2)
boxplot=table(c(result_table$antigen_1,result_table$antigen_2))/nrow(result_table)
boxplot=boxplot[order(boxplot,decreasing = T)][1:30]


pdf("BSAB_antigen_output.pdf", width = 10, height = 4)
{
  plot_bsabs(small_res,label='case',pval_cut_off=0.01,pair_score_cut_off=quantile(small_res$pair_score,.99))
  print(barplot(boxplot,las=2, main='Density of genes across bispecific pairs',cex.names=0.9,col=c('black','white'),ylim=c(0,0.5)))

}

dev.off()
}
