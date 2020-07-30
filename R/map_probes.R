
#' Map 450k/EPIC probes to the reference genome
#'
#' @param probes Probe names (450k or EPIC).
#' @param path directory were bisulfite-converted genomes are stored (see details).
#' @param array Array used (450k or EPIC)
#' @param chromosomes Which chromosomes to include (by default: chr1:22, X, Y and M).
#' @param min_width Minimum probe length to map (starting from the 3'-end of the probe).
#' @param max_width Maximum probe length to map (starting from the 3'-end of the probe).
#' @param step_size Map probe lengths from min_width to max_width in these steps.
#' @param allow_mismatch Allow a mismatch in matching? (TRUE/FALSE)
#' @param allow_INDEL Allow an INDEL in matching? (TRUE/FALSE)
#' @param cores Number of cores to use (default = 1).
#' @param verbose Verbose (TRUE/FALSE).
#' @return A data frame with one row for each match and .. columns 
#'     \item{Probe}{Probe ID}
#'     \item{chr, start, end, strand}{chromosome positions}
#'     \item{next_base}{base preceding the 'C' base for type I probes}
#'     \item{mismatch_pos}{position of mismatch (bp from 3'end of probe), NA if exact match}
#'     \item{Type2}{Type: II, I_Methylated or I_Unmethylated}
#'     \item{width}{width (in basepairs) of the match}
#'     \item{channel}{predicted color channel for type I probes}
#' @export
#' @examples
#' # Map some probes to all chromosomes, allowing no mismatch/INDELs
#' probes <- c("cg19883843", "cg24299847", "cg05542153")
#' \dontrun{
#' matches <- DNAmCrosshyb::map_probes(probes,
#'            path = "~/Desktop/genome_bs/hg19",
#'            chromosomes = "all",
#'            min_width = 30,
#'            max_width = 50,
#'            step_size = 5,
#'            allow_mismatch = FALSE,
#'            allow_INDEL = FALSE,
#'            cores = 1
#'             )
#' }
map_probes <- function(probes, path,
                       array = "450k",
                       chromosomes = "all",
                       min_width = 15,
                       max_width = 50,
                       step_size = 5,
                       allow_mismatch = FALSE,
                       allow_INDEL = FALSE,
                       cores = 1,
                       verbose = TRUE) {
  
  # Some checks
  if(min_width < 10 || max_width > 50) stop("widths should be >=10 and <=50.")
  if(cores < 1 || !is.numeric(cores)) cores <- 1
  if(step_size < 1 || !is.numeric(step_size)) stop("step_size should be a positive integer.")
  #if(!genome_build %in% c("hg19", "hg38")) stop("genome_build should be either hg19 or hg38.")
  if(allow_INDEL) stop("INDELs are not incorporated yet..")
  
  # Load appropriate probe sequences -> UPDATE LATER!
  if(array == "450k") {
    probe_sequences <- DNAmCrosshyb:::probe_sequences_450k
    anno <- DNAmCrosshyb:::anno450k
  } else if(array == "EPIC") {
    probe_sequences <- DNAmCrosshyb:::probe_sequences_EPIC
    anno <- DNAmCrosshyb:::annoEPIC
  }
  
  ## Check if probes are in probe sequences
  check <- mean(probes %in% probe_sequences$Probe)
  if(check < 1) {
    message("Not all probes found, are you sure these probes are present on the specified array?")
  }
  
  widths <- seq(from = min_width, to = max_width, by = step_size)
  if(widths[length(widths)] != max_width) widths <- c(widths, max_width)
  
  probe_sequences <- probe_sequences %>%
    dplyr::filter(Probe %in% probes) %>%
    dplyr::mutate(j = 1:nrow(.))
  
  ## Create stringset
  stringset <- Biostrings::DNAStringSet(probe_sequences$Probe_sequence)
  #stringset <- Biostrings::reverseComplement(stringset)
  ids <- probe_sequences$j
  
  ## Run for all chromosomes
  if(length(chromosomes) == 1 && chromosomes == "all") {
    chromosomes <- c(1:22, "X", "Y", "M")
  } else {
    chromosomes <- as.character(chromosomes)
  }
  
  if(cores > 1) {
    matches <- BiocParallel::bplapply(widths,
                             .f = map_per_width_dict,
                             ids = ids,
                             stringset = stringset,
                             anno = anno,
                             path = path,
                             allow_mismatch = allow_mismatch,
                             allow_INDEL = allow_INDEL,
                             chromosomes = chromosomes,
                             probe_sequences = probe_sequences,
                             BPPARAM = BiocParallel::MulticoreParam(nworkers = nr_cores)
    )
    matches <- do.call(rbind, matches)
  } else {
    matches <- purrr::map_df(widths,
                             .f = map_per_width_dict,
                             ids = ids,
                             stringset = stringset,
                             anno = anno,
                             path = path,
                             allow_mismatch = allow_mismatch,
                             allow_INDEL = allow_INDEL,
                             chromosomes = chromosomes,
                             probe_sequences = probe_sequences,
                             verbose = verbose
    )
  }
  if(length(widths) > 1) {
    # Select max match length at each position
    matches <- select_max_width(matches)
  }
  matches <- annotate_color_channel(matches, anno = anno)
  
  # matches <- matches %>%
  #   dplyr::left_join(probe_sequences[,c("Probe", "i", "j")], by = "j") %>%
  #   dplyr::select(-c("j", "k")) %>%
  #   dplyr::select(Probe, i, dplyr::everything())
  matches
  
}

map_per_width_dict <- function(width, ids, stringset, anno, path, allow_mismatch, allow_INDEL, chromosomes,
                               probe_sequences, verbose) {
  if(verbose) message("Width: ", width)
  stringset <- Biostrings::subseq(stringset, start = 50 - width + 1, end = 50)
  stringset <- Biostrings::reverseComplement(stringset)
  if(!allow_mismatch && !allow_INDEL) {
    stringset_dict <- Biostrings::PDict(stringset)
  } else if(allow_mismatch && !allow_INDEL) {
    stringset_dict <- Biostrings::PDict(stringset, max.mismatch = 1)
  } else {
    stringset_dict <- stringset
  }
  
  if(verbose) message("Mapping probes to chromosome ..")
  matches_all <- purrr::map_df(chromosomes,
                               .f = map_per_chr_dict,
                               ids = ids,
                               stringset_dict = stringset_dict,
                               stringset = stringset,
                               anno = anno,
                               path = path,
                               allow_mismatch = allow_mismatch,
                               allow_INDEL = allow_INDEL,
                               probe_sequences = probe_sequences,
                               width = width,
                               verbose = verbose
                          
  )
  if(nrow(matches_all) > 0) {
    matches_all$width <- width
  }
  matches_all
}

map_per_chr_dict <- function(chr, ids, stringset, stringset_dict, anno, path, allow_mismatch, allow_INDEL, probe_sequences, width, verbose) {
  if(verbose) message(chr)
  forward <- readRDS(sprintf("%s/forward/chr%s.rds", path, chr ))
  matches_forward <- get_matches_strand("forward", chr = chr, genome_bs = forward, stringset = stringset,
                                        stringset_dict = stringset_dict, allow_mismatch = allow_mismatch,
                                        allow_INDEL = allow_INDEL, width = width)
  
  reverse <- readRDS(sprintf("%s/reverse/chr%s.rds", path, chr ))
  matches_reverse <- get_matches_strand("reverse", chr = chr, genome_bs = reverse, stringset = stringset,
                                        stringset_dict = stringset_dict, allow_mismatch = allow_mismatch,
                                        allow_INDEL = allow_INDEL, width = width)
  
  forward_complement <- Biostrings::reverseComplement(forward)
  matches_forward_complement <- get_matches_strand("forward_complement", chr = chr, genome_bs = forward_complement, stringset = stringset,
                                                   stringset_dict = stringset_dict, allow_mismatch = allow_mismatch,
                                                   allow_INDEL = allow_INDEL, width = width)
  
  reverse_complement <- Biostrings::reverseComplement(reverse)
  matches_reverse_complement <- get_matches_strand("reverse_complement", chr = chr, genome_bs = reverse_complement, stringset = stringset,
                                                   stringset_dict = stringset_dict, allow_mismatch = allow_mismatch,
                                                   allow_INDEL = allow_INDEL, width = width)
  
  matches <- dplyr::bind_rows(matches_forward, matches_reverse, matches_forward_complement, matches_reverse_complement)
  if(nrow(matches) > 0) {
    matches <- clean_matches(matches, anno = anno, probe_sequences = probe_sequences)
  } else {
    matches <- tibble::tibble(Probe = character(), chr = character(), start = double(), end = double(),
                              next_base = character(), mismatch_pos = double(), Type2 = character())
  }
  
  matches
}

get_matches_strand <- function(strand, chr, genome_bs, stringset, stringset_dict,
                               allow_mismatch, allow_INDEL, width) {
  # Run
  matches <- match_pattern_dict(
    stringset = stringset_dict,
    subject = genome_bs, chr = chr,
    allow_INDEL = allow_INDEL, allow_mismatch = allow_mismatch)
  
  if(nrow(matches) > 0) {
    matches$strand <- strand
    # Add next_base
    matches <- add_next_base(matches, genome_bs = genome_bs)
    matches$mismatch_pos <- NA
  }
  if(nrow(matches) > 0 && allow_mismatch) {
    matches$mismatch_pos <- unlist(purrr::map(unique(matches$j), get_mismatch_position,
                                              matches = matches, genome_bs = genome_bs, probes = stringset, width = width))
  }
  if(nrow(matches) > 0 && strand %in% c("reverse", "forward_complement")) {
    start <- matches$start
    end <- matches$end
    matches$start <- length(genome_bs) - end + 1
    matches$end <- length(genome_bs) - start + 1
  }
  matches
}

match_pattern_dict <- function(subject, stringset,
                               strand,
                               chr, allow_mismatch, allow_INDEL) {
  if(allow_INDEL) {
    ## TODO
    ##  mtch <- Biostrings::matchPDict(stringset, subject,
    #   with.indels = TRUE, max.mismatch = 1, fixed = "pattern")
  } else if(allow_mismatch) {
    mtch <- Biostrings::matchPDict(stringset, subject,
                                   with.indels = FALSE, max.mismatch = 1, fixed = "pattern")
  } else {
    mtch <- Biostrings::matchPDict(stringset, subject,
                                   with.indels = FALSE, fixed = "pattern")
  }
  if(sum(IRanges::elementNROWS(mtch) > 0)) {
    positions <- create_indices(mtch, chr = chr)
  } else {
    positions <- tibble::tibble(j = numeric(), k = numeric(), chr = character(), start = double(), end = double())
  }
  positions
}

create_indices <- function(match, chr, k = TRUE) {
  index_matches <- IRanges::elementNROWS(match)
  index_matches <- tibble::tibble(j = 1:length(index_matches), nr_matches = index_matches)
  index_matches <- index_matches %>% dplyr::filter(nr_matches > 0)
  new_index <- rep(index_matches$j, times = index_matches$nr_matches)
  starts <- unlist(Biostrings::startIndex(match))
  ends <- unlist(Biostrings::endIndex(match))
  if(k) {
    index_matches2 <- tibble::tibble(j = new_index, k = 1:length(new_index), chr = chr, start = starts, end = ends)
  } else {
    index_matches2 <- tibble::tibble(j = new_index, chr = chr, start = starts, end = ends)
  }
  index_matches2
}


get_mismatch_position <- function(index, matches, genome_bs, probes, width) {
  positions.tmp <- matches %>% dplyr::filter(j == index)
  
  # Get sequences for matches
  seqs <- get_seqs(genome_bs, positions.tmp$start, positions.tmp$end)
  check <- matrix(!Biostrings::hasLetterAt(seqs,
                                           probes[[index]],
                                           at = 1:width, fixed = FALSE),
                  nrow = length(seqs),
                  ncol = stringr::str_length(probes[[index]]))
  mismatch_pos <- apply(check, 1, FUN = which)
  no_mismatch <- !apply(check, 1, FUN = any)
  mismatch_pos[no_mismatch] <- NA
  mismatch_pos <- unlist(mismatch_pos)
  mismatch_pos
}

get_seqs <- function(dnastring, starts, ends) {
  seqs <- purrr::map2(.x = starts, .y = ends, .f = get_seq, dnastring = dnastring)
  seqs <- Biostrings::DNAStringSet(seqs)
  seqs
}

get_seq <- function(start, end, dnastring) {
  seq <- dnastring[start:end]
  seq
}

add_next_base <- function(matches, genome_bs) {
  seq <- as.character(genome_bs[matches$start-1])
  matches$next_base <- stringr::str_split(seq, pattern = "")[[1]]
  matches
}


## Get 'cleaned' matches (i.e. for type II probes multiple sequences exist (number of R bases)^2)), select one
clean_matches <- function(matches, anno, probe_sequences) {
  
  matches <- matches %>%
    dplyr::left_join(probe_sequences[,c("j", "i", "ProbeSeq", "Probe")], by = "j") %>%
    dplyr::left_join(anno[,c("Name", "Type")], by = c("Probe" = "Name"))
  
  
  matches <- matches %>% dplyr::mutate(
    Type2 = dplyr::case_when(
      Type == "I" & ProbeSeq == "ProbeSeqA" ~ "I_Unmethylated",
      Type == "I" & ProbeSeq == "ProbeSeqB" ~ "I_Methylated",
      TRUE ~ "II"
    )
  )
  # Add ID
  matches <- matches %>%
    dplyr::mutate(probe_start_end_strand = paste(Probe, Type2, start, end, strand, sep = "_"))
  
  # Check
  matches <- matches %>%
    dplyr::distinct(probe_start_end_strand, .keep_all = TRUE) %>%
    dplyr::select(Probe, chr, start, end, strand, next_base, mismatch_pos, Type2)
  
  matches
}

select_max_width <- function(matches) {
  # Add ID
  matches <- matches %>%
    dplyr::mutate(id = dplyr::case_when(
      strand %in% c("forward", "reverse_complement") ~ paste(Probe, Type2, strand, chr, start, sep = "_"),
      strand %in% c("reverse", "forward_complement") ~ paste(Probe, Type2, strand, chr, end, sep = "_")
    )) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(width == max(width)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-id)

  matches
}

annotate_color_channel <- function(matches, anno) {
  matches <- matches %>%
    dplyr::left_join(anno[,c("Name", "Color")],
                     by = c("Probe" = "Name")) %>%
    dplyr::mutate(predicted_color = ifelse(next_base %in% c("A", "T"), "Red",
                                           ifelse(next_base %in% c("C", "G"), "Grn", NA)),
                  channel = ifelse(Color == predicted_color & Type2 %in% c("I_Methylated", "I_Unmethylated"), "IB", 
                                   ifelse(Color != predicted_color & Type2 %in% c("I_Methylated", "I_Unmethylated"), "OOB", NA))
                  ) %>%
    dplyr::select(-c("Color", "predicted_color"))
}

