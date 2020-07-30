
#' Check if matches overlap with repeat sequences
#'
#' @param matches Matches (generated using the 'map_probes' function)
#' @param genome_build Genome build (hg19 or hg38)
#' @param min_overlap Minimum overlap with repeat. Defaults to "any".
#' @export 
find_repeat_overlaps <- function(matches, genome_build, min_overlap = "any") {
  stopifnot(genome_build %in% c("hg19", "hg38"))
  if(genome_build == "hg19") {
    genome <- BSgenome.Hsapiens.UCSC.hg19.masked::BSgenome.Hsapiens.UCSC.hg19.masked
  } else if(genome_build == "hg38") {
    genome <- BSgenome.Hsapiens.UCSC.hg38.masked::BSgenome.Hsapiens.UCSC.hg38.masked
  } 
  chromosomes <- unique(matches$chr)
  matches <- purrr::map_df(chromosomes, .f = find_repeat_overlaps_perchr, matches = matches, 
                           genome = genome, min_overlap)
  matches
}

find_repeat_overlaps_perchr <- function(chromosome, matches, genome, min_overlap) {
  matches <- matches %>% dplyr::filter(chr == chromosome)
  repeat_locs <- genome[[paste0("chr",chromosome)]]@masks$RM
  repeat_locs <- GenomicRanges::GRanges(seqnames = paste0("chr",chromosome),
                                        ranges = IRanges::IRanges(start = BiocGenerics::start(repeat_locs), 
                                                                  end = BiocGenerics::end(repeat_locs)))
  match_locs <- GenomicRanges::GRanges(seqnames = paste0("chr",chromosome),
                                       ranges = IRanges::IRanges(start = matches$start, end = matches$end),
                                       strand = "*")
  if(min_overlap == "any") {
    overlaps <- GenomicRanges::findOverlaps(match_locs, repeat_locs)
  } else {
    overlaps <- GenomicRanges::findOverlaps(match_locs, repeat_locs, minoverlap = min_overlap)
  }
  
  matches <- matches %>% 
    dplyr::mutate(i = 1:nrow(.)) %>%
    dplyr::mutate(repeat_overlap = ifelse(i %in% S4Vectors::queryHits(overlaps), TRUE, FALSE)) %>%
    dplyr::select(-i)
  matches
}