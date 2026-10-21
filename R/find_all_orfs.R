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
#'   `strand`, `frame` (0, 1 o 2; en la hebra `+` es `(start - 1) %% 3`, en la
#'   `-` se mide desde el extremo 3' del genoma, `(longitud - end) %% 3`, y no
#'   es comparable con el de la `+` ni entre genomas de distinta longitud; para
#'   marcos relativos usar `Relative_Frame` de
#'   [detect_and_classify_overlaps()]), `length` (en codones, stop incluido) y
#'   `aa_sequence`.
#' @export
find_all_orfs <- function(fasta_file, min_length = 30, keep_stop_codon = FALSE) {
  seqs <- readDNAStringSet(fasta_file)

  start_codons <- c("ATG", "GTG", "TTG")
  stop_codons  <- c("TAA", "TAG", "TGA")

  # Coordenadas de los ORFs de una hebra/frame, vectorizado. Misma logica que
  # find_orfs_in_seq(): todos los starts en fase, cada uno con el primer stop
  # en fase que le sigue. Devuelve start/end sobre `seq_chr` (hebra ya orientada).
  orf_coords <- function(seq_chr, frame) {
    n <- nchar(seq_chr)
    if (n - 2L < frame + 1L) return(NULL)
    codon_starts <- seq.int(frame + 1L, n - 2L, by = 3L)
    codons <- substring(seq_chr, codon_starts, codon_starts + 2L)
    k <- which(codons %in% start_codons)
    stop_pos <- which(codons %in% stop_codons)
    if (length(k) == 0L || length(stop_pos) == 0L) return(NULL)
    idx <- findInterval(k, stop_pos) + 1L
    ok <- idx <= length(stop_pos)
    k <- k[ok]; idx <- idx[ok]
    i <- codon_starts[k]
    j <- codon_starts[stop_pos[idx]]
    keep <- (j + 3L - i) >= min_length
    if (!any(keep)) return(NULL)
    list(start = i[keep], end = j[keep] + 2L)
  }

  genome_results <- vector("list", length(seqs))
  orf_dna <- vector("list", length(seqs))

  for (g in seq_along(seqs)) {
    genome_id <- names(seqs)[g]
    dna_seq   <- seqs[[g]]
    contig_len <- length(dna_seq)
    strand_seq <- list(
      "+" = as.character(dna_seq),
      "-" = as.character(reverseComplement(dna_seq))
    )

    parts <- list()
    dna   <- character()
    for (strand in c("+", "-")) {
      for (frame in 0:2) {
        o <- orf_coords(strand_seq[[strand]], frame)
        if (is.null(o)) next
        # ORF ya en orientacion 5'->3' de su hebra (para "-" es el rev. complementario)
        dna <- c(dna, substring(strand_seq[[strand]], o$start, o$end))
        if (strand == "-") {
          real_start <- contig_len - o$end + 1L
          real_end   <- contig_len - o$start + 1L
        } else {
          real_start <- o$start
          real_end   <- o$end
        }
        parts[[length(parts) + 1L]] <- tibble(
          genome = genome_id,
          ORF = paste0(genome_id, "_ORF_", strand, real_start, "_", real_end),
          start = real_start,
          end = real_end,
          strand = strand,
          frame = frame,
          length = (o$end - o$start + 1L) / 3
        )
      }
    }
    genome_results[[g]] <- if (length(parts) > 0) bind_rows(parts) else NULL
    orf_dna[[g]] <- dna
  }

  res <- bind_rows(genome_results)
  if (nrow(res) == 0L) return(res)

  # Traduccion en un solo bloque: translate(if.fuzzy.codon = "solve") reconstruye
  # el codigo genetico en cada llamada (~0.1 s), asi que llamarlo por ORF era el
  # cuello de botella.
  aa <- as.character(translate(Biostrings::DNAStringSet(unlist(orf_dna)), if.fuzzy.codon = "solve"))
  if (!keep_stop_codon) aa <- sub("\\*$", "", aa)
  res$aa_sequence <- unname(aa)
  res
}
