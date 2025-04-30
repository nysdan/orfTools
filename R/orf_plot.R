#' Plot genemap
#' @import ggplot2
#' @import dplyr
#' @import gggenes
#' @param orf_table orf_table
#' @return ggplot2 object
#' @export
#'
plot_orfs <- function(orf_table) {
  orf_plot_data <- orf_table %>%
    mutate(y = case_when(
      strand == "+" ~ as.numeric(frame),
      strand == "-" ~ -as.numeric(frame)
    ))

  p <- ggplot(orf_plot_data, aes(
    xmin = start, xmax = end,
    y = y, forward = (strand == "+"), fill = strand, label = ORF
  )) +
    geom_gene_arrow(arrowhead_height = unit(2, "mm"), arrowhead_width = unit(1, "mm")) +
    geom_gene_label(align = "left", rejoin = TRUE, size = 2.5) +
    scale_fill_manual(values = c("+" = "skyblue", "-" = "tomato")) +
    scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
    labs(title = "ORFs and Overlaps Map", x = "Position (bp)", y = "Frame/Strand") +
    theme_genes() +
    theme(text = element_text(size=12))

  return(p)
}
