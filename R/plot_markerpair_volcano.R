#' @import ggplot2
#' @import ggrepel

plot_markerpair_volcano=function(marker_expr,Marker_pair_table_top30,plot_name)
{

marker_expr$Category[marker_expr$Sample=="CASE"]="Case"
marker_expr$Category[marker_expr$Sample=="CONTROL"]="Control"
marker_expr$Category[marker_expr$Sample%in%c("HEART","BRAIN", "LUNG","KIDNEY","LIVER")]="Healthy_Normal"


pdf(plot_name) # PLOT FIG 2A and 2B
{
for (i in 1:30)

{
  subset=Marker_pair_table_top30[i,]
  g1=subset$antigen_1
  g2=subset$antigen_2
  gg_pair=ggplot(marker_expr, aes_string(x = g1, y = g2 , color = "Category")) +
    geom_point() +
    labs(x = g1, y = g2, title = paste0("Pair_",g1,"_",g2)) +
    scale_color_manual(values = c("Case" = "red", "Control" = "blue", "Healthy_Normal" = "gray")) +
    theme_minimal()
    print(gg_pair)

}

}

dev.off()
}


