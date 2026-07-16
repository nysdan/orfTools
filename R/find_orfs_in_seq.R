#' Find ORFs within a single reading frame
#'
#' Escanea una secuencia de nucleotidos en un frame de lectura fijo (0, 1 o 2)
#' y devuelve, para cada codon de inicio (ATG/GTG/TTG), el primer codon de
#' stop en fase que le sigue. Una vez reportado un ORF, la busqueda continua
#' despues de su stop: ORFs anidados dentro de uno ya reportado (mismo marco)
#' no se listan por separado, ya que se asume traduccion desde el primer start
#' disponible.
#'
#' @param seq A character string con la secuencia de nucleotidos (una hebra,
#'   5' -> 3'). Para la hebra reversa, pasar el reverso complementario ya
#'   calculado (ver [find_all_orfs()]).
#' @param strand `"+"` o `"-"`. Se usa unicamente para etiquetar los ORFs
#'   devueltos; la funcion no calcula el reverso complementario.
#' @param frame_offset Offset del marco de lectura (0, 1 o 2) relativo al
#'   inicio de `seq`.
#' @param min_length Longitud minima del ORF en nucleotidos, incluyendo los
#'   codones de inicio y de stop.
#' @return Una lista de ORFs; cada uno es una lista con `start`, `end`
#'   (coordenadas 1-based, inclusive, dentro de `seq`), `strand`, `frame`,
#'   `length` (en codones, incluyendo el stop) y `sequence` (secuencia
#'   nucleotidica del ORF).
#' @export
find_orfs_in_seq <- function(seq, strand, frame_offset, min_length = 30) {
  start_codons <- c("ATG", "GTG", "TTG")
  stop_codons  <- c("TAA", "TAG", "TGA")

  seq_len <- nchar(seq)
  codon_starts <- seq.int(frame_offset + 1L, seq_len - 2L, by = 3L)
  if (length(codon_starts) == 0) return(list())

  codons <- substring(seq, codon_starts, codon_starts + 2L)
  is_start <- codons %in% start_codons
  stop_positions <- which(codons %in% stop_codons)

  orfs <- list()
  k <- 1L
  n_codons <- length(codons)

  while (k <= n_codons) {
    if (!is_start[k]) {
      k <- k + 1L
      next
    }

    # siguiente codon de stop en fase, si existe (busqueda binaria sobre
    # las posiciones de stop, ya ordenadas)
    idx <- findInterval(k, stop_positions) + 1L
    if (idx > length(stop_positions)) break  # no hay mas stops en este marco

    stop_k <- stop_positions[idx]
    i <- codon_starts[k]
    j <- codon_starts[stop_k]
    orf_length <- j + 3L - i

    if (orf_length >= min_length) {
      orfs[[length(orfs) + 1L]] <- list(
        start = i,
        end = j + 2L,
        strand = strand,
        frame = paste0(strand, frame_offset),
        length = orf_length / 3,
        sequence = substr(seq, i, j + 2L)
      )
    }

    k <- stop_k + 1L
  }

  orfs
}
