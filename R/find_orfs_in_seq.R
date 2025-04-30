#' Find seqs of  ORFs in a genome
#' @param seq A character object
#' @param strand strand
#' @param frame_offset frame offset
#' @param min_length Minimal ORF length
#' @return A data frame with predicted ORFs
#' @export
find_orfs_in_seq <- function(seq, strand, frame_offset, min_length = 30) {
  start_codon <- "ATG"
  stop_codons <- c("TAA", "TAG", "TGA")
  orfs <- list()

  seq_len <- nchar(seq)
  i <- frame_offset + 1

  while (i <= (seq_len - 2)) {
    codon <- substr(seq, i, i+2)
    if (codon == start_codon) {
      for (j in seq(i + 3, seq_len-2, by = 3)) {
        stop_codon <- substr(seq, j, j+2)
        if (stop_codon %in% stop_codons) {
          orf_length <- j + 3 - i
          if (orf_length >= min_length) {
            frame_value <- ifelse(strand == "+", frame_offset, -frame_offset)
            orfs[[length(orfs) + 1]] <- list(
              start = i, end = j+2, strand = strand, frame = frame_value, length = orf_length
            )
          }
          break
        }
      }
      i <- j + 3
    } else {
      i <- i + 3
    }
  }

  return(orfs)
}
