#' Calculate 3'-subsequence overlap between a set of probes
#'
#' @param probes Probes to test (450k/EPIC)
#' @return A data frame with one row for each match and .. columns 
#'   \item{Probe_Index}{Probe ID}
#'   \item{Bead_Index}{Bead ID: Probe ID with '_methylated' or '_unmethylated' appended for type I probes}
#'   \item{Probe_Target}{Probe ID}
#'   \item{Bead_Target}{Bead ID: Probe ID with '_methylated' or '_unmethylated' appended for type I probes}
#'   \item{overlap}{overlap in bp}
#'   \item{mismatch_pos}{position of mismatch (bp from 3'end of probe), NA if exact match}
#' @export
#' @example 
#' # Get overlaps for set of probes 
#' library(DNAmCrosshyb)
#' library(dplyr)
#' # Run 
#' overlaps <- get_probe_overlaps(c("cg09994391", "cg19403339", "cg01370437", "cg18002896", "cg11613875", "cg21249376", 
#'                                "cg23074747", "cg15793563", "cg20307896", "cg00521048", "cg00801568", "cg16517021", "cg12078510",
#'                                 "cg22383472", "cg09022230", "cg05990720", "cg14026202", "cg14363787"))
#' # Show number of >=10bp overlaps per probe 
#' overlaps %>% dplyr::group_by(Bead_Index) %>% summarize(n = sum(overlap >= 10))

get_probe_overlaps <- function(probes) {
  
  # Check if probes are 450k or EPIC 
  if(all(probes %in% DNAmCrosshyb:::probes_450k)) {
    array <- "450k"
  } else if(all(probes %in% DNAmCrosshyb:::probes_epic)) {
    array <- "EPIC"
  } else if(all(probes %in% c(DNAmCrosshyb:::probes_epic, DNAmCrosshyb:::probes_450k))) {
    array <- "450k/EPIC"
  } else {
    stop("Not all probes are recognized, are they not from the 450k or EPIC array?")
  }
  
  if(array == "450k") {
    anno <- DNAmCrosshyb:::anno450k
  } else if(array == "EPIC") {
    anno <- DNAmCrosshyb:::annoEPIC
  } else if(array == "450k/EPIC") {
    anno_450k <- DNAmCrosshyb:::anno450k
    anno_EPIC <- DNAmCrosshyb:::annoEPIC
    anno_EPIC <- anno_EPIC %>% dplyr::filter(!Name %in% anno_450k$Name)
    anno <- dplyr::bind_rows(anno_450k, anno_EPIC)
    rm(anno_450k, anno_EPIC)
  }
  
  typeI <- anno %>% dplyr::filter(Name %in% probes, Type == "I")
  typeII <- anno %>% dplyr::filter(Name %in% probes, Type == "II")
  seqs <- tibble::tibble(Probe = c(typeI$Name, typeI$Name, typeII$Name), Bead = c(paste0(typeI$Name, "_unmethylated"), paste0(typeI$Name, "_methylated"), typeII$Name ),
                 seq = c(typeI$ProbeSeqA, typeI$ProbeSeqB, typeII$ProbeSeqA))
  seqs$i <- 1:nrow(seqs)
  rm(typeI, typeII)
  
  message("Calculating overlaps..")
  overlaps <- purrr::map_df(seqs$i, .f = get_probe_overlap, seqs = seqs)
  message("Done!")
  overlaps
}

get_probe_overlap <- function(i_, seqs) {
  seq <- IRanges::reverse(Biostrings::DNAString(seqs %>% dplyr::filter(i == i_) %$% seq))
  seqs_ <- IRanges::reverse(Biostrings::DNAStringSet(seqs %>% dplyr::filter(i != i_) %$% seq))
  check <- matrix(Biostrings::hasLetterAt(seqs_,
                                          seq,
                                          at = 1:50, fixed = FALSE),
                  nrow = length(seqs_),
                  ncol = stringr::str_length(seq))
  sums <- rowSums(check)
  mismatch_pos <- t(unlist(apply(!check, 1, FUN = function(x) which(x)[1:2])))
  mismatch_pos <- tibble::tibble(Probe_Index = seqs %>% dplyr::filter(i == i_) %$% Probe,
                         Bead_Index = seqs %>% dplyr::filter(i == i_) %$% Bead,
                         Probe_Target = seqs %>% dplyr::filter(i != i_) %$% Probe,
                         Bead_Target = seqs %>% dplyr::filter(i != i_) %$% Bead,
                         overlap = mismatch_pos[,2],
                         mismatch_pos = mismatch_pos[,1]
                         )
  mismatch_pos  <- mismatch_pos %>% dplyr::mutate(overlap = ifelse(is.na(overlap), 50, overlap))
  mismatch_pos
}
