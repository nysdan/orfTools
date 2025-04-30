#' Find all ORFs in a genome
#' @importFrom stringr str_count
#' @importFrom Biostrings reverseComplement DNAString readDNAStringSet
#' @import tibble
#' @param fasta_file Path to a FASTA file
#' @return A data frame with predicted ORFs
#' @export
#'


find_all_orfs <- function(fasta_file) {
  seqs <- readDNAStringSet(fasta_file)
  results <- list()

  for (i in seq_along(seqs)) {
    genome_id <- names(seqs)[i]
    seq <- toupper(as.character(seqs[[i]]))
    rev_seq <- as.character(reverseComplement(seqs[[i]]))
    seq_len <- nchar(seq)

    for (frame in 0:2) {
      orfs <- find_orfs_in_seq(seq, "+", frame)
      for (orf in orfs) {
        orf_id <- paste0(genome_id, "_ORF_+", orf$start, "_", orf$end)
        results[[length(results) + 1]] <- tibble(
          genome = genome_id,
          ORF = orf_id,
          start = orf$start,
          end = orf$end,
          strand = "+",
          frame = frame,
          length = orf$length
        )
      }
    }

    for (frame in 0:2) {
      orfs <- find_orfs_in_seq(rev_seq, "-", frame)
      for (orf in orfs) {
        real_start <- seq_len - orf$end + 1
        real_end <- seq_len - orf$start + 1
        orf_id <- paste0(genome_id, "_ORF_-", real_start, "_", real_end)
        results[[length(results) + 1]] <- tibble(
          genome = genome_id,
          ORF = orf_id,
          start = real_start,
          end = real_end,
          strand = "-",
          frame = frame,
          length = orf$length
        )
      }
    }
  }

  final_df <- bind_rows(results)
  return(final_df)
}
