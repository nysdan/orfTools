#' Plot a gene map of ORFs
#'
#' Dibuja un mapa genico de los ORFs de un genoma con `gggenes`, separando
#' visualmente por hebra (arriba/abajo del eje) y por marco de lectura
#' (distancia al eje).
#'
#' @import ggplot2
#' @import dplyr
#' @import gggenes
#' @param orf_table A data frame con las columnas `start`, `end`, `strand`,
#'   `frame`, `ORF`, como el devuelto por [find_all_orfs()].
#' @return Un objeto ggplot2.
#' @export
plot_orfs <- function(orf_table) {
  if (nrow(orf_table) == 0) {
    stop("orf_table no tiene filas para graficar.", call. = FALSE)
  }

  orf_plot_data <- mutate(orf_table, y = case_when(
    strand == "+" & frame == 0 ~ 0.1,
    strand == "-" & frame == 0 ~ -0.1,
    strand == "+" & frame != 0 ~ as.numeric(frame),
    strand == "-" & frame != 0 ~ -as.numeric(frame))
  )

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
