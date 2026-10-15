#' Find, translate and tabulate all ORFs in a genome
#'
#' Escanea uno o mas registros de un archivo FASTA en los 6 marcos de lectura
#' (3 en la hebra +, 3 en la reversa complementaria) usando
#' [find_orfs_in_seq()], y traduce cada ORF encontrado a proteina.
#'
#' @importFrom Biostrings readDNAStringSet reverseComplement subseq translate
#' @import tibble
#' @import dplyr
#' @param fasta_file Path a un archivo FASTA (una o mas secuencias).
#' @param min_length Longitud minima del ORF en nucleotidos, se pasa
#'   directamente a [find_orfs_in_seq()]. Por defecto 30 (10 codones
#'   incluyendo el stop).
#' @param keep_stop_codon Si es `FALSE` (default), elimina el `*` final que
#'   deja [Biostrings::translate()] para el codon de stop -- conveniente para
#'   herramientas rio abajo (p. ej. DeepTMHMM) que esperan una secuencia
#'   proteica sin caracteres de stop.
#' @details Codones con bases ambiguas IUPAC (N, R, Y, etc.) se resuelven via
#'   `translate(..., if.fuzzy.codon = "solve")`: si todas las combinaciones
#'   posibles del codon ambiguo traducen al mismo aminoacido, se usa ese
#'   aminoacido; si no, se reporta `X`. Antes de este fix, un codon ambiguo
#'   hacia fallar translate() y abortaba el procesamiento del genoma completo.
#' @return Un tibble con una fila por ORF: `genome`, `ORF` (id que codifica
#'   genoma, hebra y coordenadas genomicas de inicio/fin), `start`, `end`
#'   (1-based; para la hebra `-` ya expresadas en coordenadas de la hebra +),
#'   `strand`, `frame` (0, 1 o 2), `length` (en codones, stop incluido) y
#'   `aa_sequence`.
#' @export
find_all_orfs <- function(fasta_file, min_length = 30, keep_stop_codon = FALSE) {
  seqs <- readDNAStringSet(fasta_file)

  build_row <- function(orf, genome_id, start, end, strand, frame, dna_seq) {
    orf_seq <- subseq(dna_seq, start = start, end = end)
    if (strand == "-") orf_seq <- reverseComplement(orf_seq)

    aa_seq <- as.character(translate(orf_seq, if.fuzzy.codon = "solve"))
    if (!keep_stop_codon) aa_seq <- sub("\\*$", "", aa_seq)

    tibble(
      genome = genome_id,
      ORF = paste0(genome_id, "_ORF_", strand, start, "_", end),
      start = start,
      end = end,
      strand = strand,
      frame = frame,
      length = orf$length,
      aa_sequence = aa_seq
    )
  }

  genome_results <- vector("list", length(seqs))

  for (i in seq_along(seqs)) {
    genome_id <- names(seqs)[i]
    dna_seq   <- seqs[[i]]
    seq_len   <- nchar(as.character(dna_seq))
    fwd_seq   <- as.character(dna_seq)
    rev_seq   <- as.character(reverseComplement(dna_seq))

    rows <- list()

    for (frame in 0:2) {
      for (orf in find_orfs_in_seq(fwd_seq, "+", frame, min_length)) {
        rows[[length(rows) + 1L]] <- build_row(
          orf, genome_id, orf$start, orf$end, "+", frame, dna_seq
        )
      }
    }

    for (frame in 0:2) {
      for (orf in find_orfs_in_seq(rev_seq, "-", frame, min_length)) {
        real_start <- seq_len - orf$end + 1L
        real_end   <- seq_len - orf$start + 1L
        rows[[length(rows) + 1L]] <- build_row(
          orf, genome_id, real_start, real_end, "-", frame, dna_seq
        )
      }
    }

    genome_results[[i]] <- if (length(rows) > 0) bind_rows(rows) else NULL
  }

  bind_rows(genome_results)
}
