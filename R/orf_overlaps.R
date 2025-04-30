#' Detect overlaps
#' @importFrom dplyr case_when bind_rows
#' @import tibble
#' @param orf_table orf_table
#' @return A data frame with predicted overlapping ORFs
#' @export
detect_and_classify_overlaps <- function(orf_table) {
  overlaps <- list()

  for (i in 1:(nrow(orf_table)-1)) {
    for (j in (i+1):nrow(orf_table)) {
      orf1 <- orf_table[i, ]
      orf2 <- orf_table[j, ]

      # Verificar si se superponen
      overlap_bp <- min(orf1$end, orf2$end) - max(orf1$start, orf2$start) + 1

      if (overlap_bp > 0) {  # Si hay solapamiento
        same_strand <- (orf1$strand == orf2$strand)
        same_frame <- (orf1$frame == orf2$frame)

        overlap_type <- case_when(
          orf1$start <= orf2$start & orf1$end >= orf2$end ~ "Embedded overlap",
          orf2$start <= orf1$start & orf2$end >= orf1$end ~ "Embedded overlap",
          TRUE ~ "Terminal overlap"
        )

        overlaps[[length(overlaps) + 1]] <- tibble(
          ORF1 = orf1$ORF,
          ORF2 = orf2$ORF,
          End1 = orf1$end,
          Start2 = orf2$start,
          End2 = orf2$end,
          Start1 = orf1$start,
          Overlap_bp = overlap_bp,
          Same_Strand = same_strand,
          Same_Frame = same_frame,
          Overlap_Type = overlap_type
        )
      }
    }
  }

  overlaps_df <- bind_rows(overlaps)
  return(overlaps_df)
}
