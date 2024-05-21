#' @import ggplot2
#' @import ggrepel

plot_marker_boxplot=function(marker_expr,marker_list,plot_name)
{
unique_levels <- unique(marker_expr$Sample)
other_levels <- unique_levels[!unique_levels %in% c("CASE", "CONTROL")]
new_levels <- c("CASE", "CONTROL", sort(other_levels))

marker_expr$Sample <- factor(marker_expr$Sample, levels = new_levels)

pdf(plot_name, height = 5,width=8)

{

  for (marker in marker_list) {
    gg <- ggplot(marker_expr, aes(x = Sample, y = !!sym(marker), fill = Sample)) +
      geom_boxplot() +
      scale_fill_manual(values = c("red", "blue", rep("gray", length(unique(marker_expr$Sample)) - 2)), guide = "none") +
      labs(x = "Samples", y = marker, title = paste("Marker expression across samples: ", marker)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      stat_compare_means(method = "t.test", ref.group = "CASE", label = "p.signif")

    # Print the current plot
    print(gg)
  }
}
dev.off()
}





