#' @title bulk_DE_surface_antigen
#'
#' @description
#' This function enlists the bispecific antibody target pairs for given cancer.
#' The Pairtable.csv and related plots will be saved as pdf in working directory.
#'
#' @param cancer.type A character string specifying the cancer name.It should be selected from total_cancer_count. Case-sensitive.
#' @param normal.tissue A character string specifying the corresponding normal tissue for the cancer.It should be selected from normal.tissue.Case-sensitive.
#' @param octad_counts_data_path path for octad.counts.and.tpm.h5
#' @param output_file_name A character string specifying the desired outfile file name
#'
#' @return dataframe with maximum 30 bispecific antibody target pairs
#' @return Barplot pdf showing frequency of top markers
#' @return Boxplot pdf showing expression pattern of top markers across cancer, normal tissue and normal vital organs
#' @return Volcanoplot pdf for each pair
#'
#' @usage
#' bulk_DE_surface_antigen(
#' cancer.type = 'liver hepatocellular carcinoma',
#' normal.tissue = 'LIVER',
#' octad_counts_data_path = "~/Downloads/octad.counts.and.tpm.h5",
#' output_file_name = "liver_output"
#' )
#' @export
#' @import edgeR
#' @import bsabsfinder
#' @import cluster
#' @import dplyr
#' @import ggplot2
#' @import ggpubr

bulk_DE_surface_antigen=function(cancer.type, normal.tissue, octad_counts_data_path,output_file_name){

case=subset(phenoDF,cancer==cancer.type & sample.type == 'primary')
case_id=case$sample.id #select cases

control <- data.frame()

# Loop through each tissue type in normal.tissue and subset phenoDF
if (length(normal.tissue) == 1) {
  control = subset(phenoDF, biopsy.site == normal.tissue & sample.type == 'normal')
} else if (length(normal.tissue) > 1) {
  for (tissue in normal.tissue) {
  subset_df <- subset(phenoDF, biopsy.site == tissue & sample.type == 'normal')
  control <- rbind(control, subset_df) # Combine subsets
  }
  # control_list <- lapply(normal.tissue, function(tissue) {
  #   subset(phenoDF, biopsy.site == tissue & sample.type == 'normal')
  # })
  # # Combine subsets into a single data frame
  # control <- do.call(rbind, control_list)
  #control = subset(phenoDF, biopsy.site %in% normal.tissue & sample.type == 'normal')
}else {
  stop("normal.tissue must be either a single tissue type or a multiple of tissue.")
}

control_id = control$sample.id
#control=subset(phenoDF, biopsy.site== normal.tissue & sample.type=='normal')
#control_id=control$sample.id
#

case_expr=loadOctadCounts(case_id,type='tpm',file=octad_counts_data_path)
case_expr=as.data.frame(case_expr)
control_expr=loadOctadCounts(control_id ,type='tpm',file=octad_counts_data_path)
control_expr=as.data.frame(control_expr)

#final data
case_with_control_expr=cbind(case_expr,control_expr)

#convert ensg to hgnc and select surface-expressed genes according to  compartments.jensenlab.org
case_with_control_expr=ensg_to_hgnc(case_with_control_expr,select_surface=TRUE)

rownames(case_with_control_expr) <- gsub("HLA-", "HLA_", rownames(case_with_control_expr))
colnames(healthy_tissues_expr)<- gsub("HLA-", "HLA_", colnames(healthy_tissues_expr))

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
#assign("Antigen_pair_table", small_res,  envir = .GlobalEnv)
#write.csv(small_res, "Antigen_pairs_calculation.csv")


result_table=small_res
result_table=subset(result_table,pair_score>quantile(result_table$pair_score,.90))
result_table=subset(result_table,pair_score>2)
result_table=subset(result_table,case_greater=='TRUE_TRUE')
result_table=result_table[order(result_table$pair_score,decreasing = T),]
result_table$antigen_1=gsub('-','_',result_table$antigen_1)
result_table$antigen_2=gsub('-','_',result_table$antigen_2)
boxplot=table(c(result_table$antigen_1,result_table$antigen_2))/nrow(result_table)
boxplot=boxplot[order(boxplot,decreasing = T)][1:30]


#output_file_name
plot_name1=paste0(output_file_name, "_BSAB_allpairs_overview.pdf")
plot_name2=paste0(output_file_name, "_BSAB_marker_density.pdf")
plot_name3=paste0(output_file_name, "_BSAB_marker_boxplot.pdf")
plot_name4=paste0(output_file_name, "_BSAB_markerpairs.pdf")
table_name <- paste0(output_file_name, "_BSABS_top30_pairtable.csv")


print("Generating All Marker Pair Overview Plot")

pdf(plot_name1, width = 10, height = 4)
{
  plot_bsabs(small_res,label='case',pval_cut_off=0.01,pair_score_cut_off=quantile(small_res$pair_score,.99))

}
dev.off()


print("Generating Marker Pair Table")

assign("Marker_pair_table_all", result_table,  envir = .GlobalEnv)

Marker_pair_table_top30=data.frame()

if(nrow(result_table)>30)

{
  Marker_pair_table_top30=result_table[c(1:30),]

} else

{
  Marker_pair_table_top30=result_table

}


write.csv(Marker_pair_table_top30, table_name ,row.names = FALSE)

assign("Marker_pair_table_top30", Marker_pair_table_top30,  envir = .GlobalEnv)


marker_list1=c(Marker_pair_table_top30$antigen_1,Marker_pair_table_top30$antigen_2)

marker_counts <- table(marker_list1)
marker_counts_df <- as.data.frame(marker_counts)
names(marker_counts_df) <- c("Marker", "Frequency")
marker_counts_df <- marker_counts_df[order(-marker_counts_df$Frequency), ]
marker_counts_df$Marker <- factor(marker_counts_df$Marker, levels = marker_counts_df$Marker)
proxy_var <- as.numeric(marker_counts_df$Frequency)

gg_freq<-ggplot(marker_counts_df, aes(x= Marker, y=Frequency, fill=proxy_var)) +
  geom_bar(stat="identity")+
  scale_y_continuous(breaks = seq(0, max(marker_counts_df$Frequency)+2, by = 2), limits = c(0, max(marker_counts_df$Frequency)))+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  labs(x = "Marker", y = "Frequency", title = "Marker Frequency among top 30 BSAB pairs")+
  scale_fill_gradient(low = "black",high = "black", guide = "none")


pdf(plot_name2, width = 10, height = 5)
{
  #plot_bsabs(small_res,label='case',pval_cut_off=0.01,pair_score_cut_off=quantile(small_res$pair_score,.99))
  print(barplot(boxplot,las=2, main='Density of markers across all bispecific pairs',cex.names=0.9,col=c('black','white'),ylim=c(0,0.5)))
  #print(gg_freq)
}

dev.off()

marker_list=unique(marker_list1)

print("Top Markers are:")
print(marker_list)

print("Generating Marker boxplots")
data_a=case_with_control_expr[rownames(case_with_control_expr)%in%marker_list,]
data_a=as.data.frame(t(data_a))
data_a=data_a[,order(colnames(data_a))]
data_a$Sample=ifelse(rownames(data_a)%in%case_id,"CASE","CONTROL")

data_b=healthy_tissues_expr[,colnames(healthy_tissues_expr)%in%marker_list]
data_b$Sample=healthy_tissues_expr$Sample

#colnames(data_a)==colnames(data_b)

common=length(intersect(rownames(data_a),rownames(data_b)))

if (common>0)
{
  com_samp=intersect(rownames(data_a),rownames(data_b))
  data_b=data_b[!rownames(data_b)%in%com_samp,]

}

marker_expr=rbind(data_a,data_b)

assign("Marker_Expression", marker_expr,  envir = .GlobalEnv)

plot_marker_boxplot(marker_expr,marker_list,plot_name3)

print("Generating Markerpair volacano plots")

plot_markerpair_volcano(marker_expr,Marker_pair_table_top30,plot_name4)

}



