#' Detect and classify overlaps between ORFs
#'
#' Compara todos los pares de ORFs de una tabla y reporta los que se solapan
#' en coordenadas genomicas, indicando si comparten hebra/frame y si el
#' solapamiento es total (uno anidado dentro del otro) o parcial (terminal).
#'
#' @importFrom utils combn
#' @import dplyr
#' @import tibble
#' @param orf_table A data frame con (al menos) las columnas `ORF`, `start`,
#'   `end`, `strand`, `frame`, como el devuelto por [find_all_orfs()].
#' @return Un tibble con una fila por par de ORFs solapados (`ORF1`, `ORF2`,
#'   sus coordenadas, `Overlap_bp`, `Same_Strand`, `Same_Frame`,
#'   `Overlap_Type`). Vacio (0 filas) si hay menos de 2 ORFs o ninguno se
#'   solapa.
#' @export
detect_and_classify_overlaps <- function(orf_table) {
  empty_result <- function() {
    tibble(
      ORF1 = character(), ORF2 = character(),
      End1 = numeric(), Start2 = numeric(), End2 = numeric(), Start1 = numeric(),
      Overlap_bp = numeric(), Same_Strand = logical(), Same_Frame = logical(),
      Overlap_Type = character()
    )
  }

  n <- nrow(orf_table)
  if (n < 2) return(empty_result())

  pairs <- combn(n, 2)
  orf1 <- orf_table[pairs[1, ], ]
  orf2 <- orf_table[pairs[2, ], ]

  overlap_bp <- pmin(orf1$end, orf2$end) - pmax(orf1$start, orf2$start) + 1
  keep <- overlap_bp > 0
  if (!any(keep)) return(empty_result())

  orf1 <- orf1[keep, ]
  orf2 <- orf2[keep, ]
  overlap_bp <- overlap_bp[keep]

  overlap_type <- case_when(
    orf1$start <= orf2$start & orf1$end >= orf2$end ~ "Embedded overlap",
    orf2$start <= orf1$start & orf2$end >= orf1$end ~ "Embedded overlap",
    TRUE ~ "Terminal overlap"
  )

  tibble(
    ORF1 = orf1$ORF,
    ORF2 = orf2$ORF,
    End1 = orf1$end,
    Start2 = orf2$start,
    End2 = orf2$end,
    Start1 = orf1$start,
    Overlap_bp = overlap_bp,
    Same_Strand = orf1$strand == orf2$strand,
    Same_Frame = orf1$frame == orf2$frame,
    Overlap_Type = overlap_type
  )
}
