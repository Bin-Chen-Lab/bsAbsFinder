library(edgeR)
library(bsabsfinder)
library(cluster)
library(dplyr)
library(ggplot2)
library(ggpubr)

case=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&sample.type == 'primary') 
case_id=case$sample.id #select cases
control=subset(phenoDF,sample.type=='normal'&biopsy.site=='LIVER')
control_id=control$sample.id

case_expr=loadOctadCounts(case_id,type='tpm',file='~/Downloads/octad.counts.and.tpm.h5')
case_expr=as.data.frame(case_expr)
control_expr=loadOctadCounts(control_id ,type='tpm',file='~/Downloads/octad.counts.and.tpm.h5') 
control_expr=as.data.frame(control_expr)

#final data
case_with_control_expr=cbind(case_expr,control_expr)

#convert ensg to hgnc and select surface-expressed genes according to  compartments.jensenlab.org
case_with_control_expr=ensg_to_hgnc(case_with_control_expr,select_surface=TRUE)


phenotype_vector=as.factor(c(rep('case',ncol(case_expr)),rep('control',ncol(control_expr))))

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

DE1=subset(DE,padj<0.01&abs(logFC)>1.3) # The cutoff criteria can be changed 

#filter out only surface-expressed DE genes. Just to speed up. 
case_with_control_expr=case_with_control_expr[row.names(case_with_control_expr)%in%row.names(DE1),]
dataframe_for_computation=as.data.frame(t(case_with_control_expr)) 

#this takes for a while
small_res=compute_bsabs(antigene_1=colnames(dataframe_for_computation),data_input=dataframe_for_computation,pheno_input=phenotype_vector)
head(small_res)

result_table=small_res
result_table=subset(result_table,pair_score>quantile(result_table$pair_score,.90))
result_table=subset(result_table,pair_score>2)
result_table=subset(result_table,case_greater=='TRUE_TRUE')
result_table=result_table[order(result_table$pair_score,decreasing = T),]
result_table$antigen_1=gsub('-','_',result_table$antigen_1)
result_table$antigen_2=gsub('-','_',result_table$antigen_2)
boxplot=table(c(result_table$antigen_1,result_table$antigen_2))/nrow(result_table)
boxplot=boxplot[order(boxplot,decreasing = T)][1:30]


pdf("BSAB_FIG2C.pdf", width = 10, height = 4)
{
  print(barplot(boxplot,las=2, main='Density of genes across bispecific pairs',cex.names=0.9,col=c('black','white'),ylim=c(0,0.5)))
  
}

dev.off()


# FIG 2D
plot_bsabs(small_res,label='case',pval_cut_off=0.01,pair_score_cut_off=quantile(small_res$pair_score,.99)) 

#subset result table to filter out only top pairs:
small_res=small_res[small_res$case_greater=="TRUE_TRUE",]

#ordering as per pair score , highest score should be at top
small_res=small_res[order(small_res$pair_score,decreasing = T),] 

#Subsetting top 20 pairs
small_res_selective=small_res[c(1:20),]

# unique marker genes in top 20 pairs
marker_list=unique(c(small_res_selective$antigen_1,small_res_selective$antigen_2))
marker_list


# PLOT FIG 2E with STATS

healthy_tissues=subset(phenoDF,sample.type=='normal')  
healthy_tissues=subset(healthy_tissues,grepl('BRAIN',biopsy.site)|biopsy.site=='LUNG'|grepl('HEART',biopsy.site)| grepl('KIDNEY',biopsy.site))

healthy_tissues <- healthy_tissues %>%mutate(biopsy.site = ifelse(grepl("BRAIN", biopsy.site), "BRAIN", biopsy.site))
healthy_tissues <- healthy_tissues %>%mutate(biopsy.site = ifelse(grepl("HEART", biopsy.site), "HEART", biopsy.site))
healthy_tissues <- healthy_tissues %>%mutate(biopsy.site = ifelse(grepl("KIDNEY", biopsy.site), "KIDNEY", biopsy.site))


healthy_tissues_expr=loadOctadCounts(healthy_tissues$sample.id ,type='tpm',file='~/Downloads/octad.counts.and.tpm.h5') #Windows
#healthy_tissues=loadOctadCounts(healthy_tissues$sample.id ,type='tpm',file='~/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')#Unix
healthy_tissues_expr=as.data.frame(healthy_tissues_expr)

healthy_tissues_expr=ensg_to_hgnc(healthy_tissues_expr,select_surface=FALSE)
healthy_tissues_expr=healthy_tissues_expr[row.names(healthy_tissues_expr)%in%marker_list,]

healthy_tissues_expr=healthy_tissues_expr[order(rownames(healthy_tissues_expr)),]

healthy_tissues_expr=as.data.frame(t(healthy_tissues_expr))

healthy_tissues_expr$Sample=healthy_tissues$biopsy.site[match(rownames(healthy_tissues_expr),healthy_tissues$sample.id)]

case_with_control_expr2=case_with_control_expr[row.names(case_with_control_expr)%in%marker_list,]

case_with_control_expr2=case_with_control_expr2[order(rownames(case_with_control_expr2)),]

case_with_control_expr2=as.data.frame(t(case_with_control_expr2))

case_with_control_expr2$Sample=ifelse(rownames(case_with_control_expr2)%in%case_id,"HCC","LIVER")

#table(case_with_control_expr2$Sample)

colnames(case_with_control_expr2)==colnames(healthy_tissues_expr)

marker_expr=rbind(case_with_control_expr2,healthy_tissues_expr)

table(marker_expr$Sample)

sample_order <- c("HCC","LIVER", "LUNG","KIDNEY", "BRAIN", "HEART")

marker_expr$Sample <- factor(marker_expr$Sample, levels = sample_order)

pdf("~/Desktop/BSAB_LIVER_marker.pdf", width = 10, height = 4)
{
  for (i in 1:length(marker_list)) {
    gg <- ggplot(marker_expr, aes(x = Sample, y = marker_expr[[i]], fill = Sample)) +
      geom_boxplot() +
      scale_fill_manual(values = c("red", "blue",rep("gray", length(sample_order) - 2)),guide = "none") +
      labs(x = "Biopsy Site", y = names(marker_expr)[i], title = paste("Marker expression by biopsy site: ", names(marker_expr)[i])) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      #guides(fill = none)+
      stat_compare_means(method = "t.test", ref.group = "HCC", label="p.signif")
    
    
    print(gg)
    
  }
}

dev.off()

# Marker Frequency plot 

markers <- c(small_res_selective$antigen_1,small_res_selective$antigen_2)

# Count the occurrences of each antigen in the combined vector
marker_counts <- table(markers)

# Convert the result to a data frame
marker_counts_df <- as.data.frame(marker_counts)
names(marker_counts_df) <- c("Marker", "Frequency")

marker_counts_df <- marker_counts_df[order(-marker_counts_df$Frequency), ]

marker_counts_df$Marker <- factor(marker_counts_df$Marker, levels = marker_counts_df$Marker)

proxy_var <- as.numeric(marker_counts_df$Frequency)

pdf("~/Desktop/BSAB_HCC_marker_frequency.pdf",width = 10, height = 4)

{
  gg_freq<-ggplot(marker_counts_df, aes(x= Marker, y=Frequency, fill=proxy_var)) +
    geom_bar(stat="identity")+
    scale_y_continuous(breaks = seq(0, max(marker_counts_df$Frequency)+2, by = 2), limits = c(0, max(marker_counts_df$Frequency)))+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    labs(x = "Marker", y = "Frequency", title = "Marker Frequency among top 20 BSAB pairs")+
    scale_fill_gradient(low = "black",high = "black", guide = "none") 
  
  print(gg_freq)
}

dev.off()

##############################################################################################

marker_expr$Category=NA

unique(marker_expr$Sample)

marker_expr$Category[marker_expr$Sample=="HCC"]="Case_HCC"
marker_expr$Category[marker_expr$Sample=="LIVER"]="Control_Liver_Normal"
marker_expr$Category[marker_expr$Sample%in%c("HEART","BRAIN", "LUNG","KIDNEY")]="Healthy_Normal"


pdf("~/Desktop/BSAB_HCC_marker_pair.pdf") # PLOT FIG 2A and 2B

{
  
  markerpair1=ggplot(marker_expr, aes(x = PLVAP, y = GPC3 , color = Category)) +
    geom_point() +
    labs(x = "PLVAP", y = "GPC3", title = "Scatter Plot with Categorical Color") +
    scale_color_manual(values = c("Case_HCC" = "red", "Control_Liver_Normal" = "blue", "Healthy_Normal" = "gray")) +
    theme_minimal()
  
  markerpair2=ggplot(marker_expr, aes(x = GPC3, y = MUC13 , color = Category)) +
    geom_point() +
    labs(x = "GPC3", y = "MUC13", title = "Scatter Plot with Categorical Color") +
    scale_color_manual(values = c("Case_HCC" = "red", "Control_Liver_Normal" = "blue", "Healthy_Normal" = "gray")) +
    #scale_shape_manual(values = c(1, 2,4), guide = "none") + 
    theme_minimal() 
  
  
  print(markerpair1)
  print(markerpair2)
  
}

dev.off()

save(list=c('case_expr','control_expr','case_with_control_expr','case_with_control_expr2','DE','DE1','healthy_tissues_expr','small_res','small_res_selective','marker_expr','healthy_tissues','marker_counts'),file="BSAB_HCC.RData")




