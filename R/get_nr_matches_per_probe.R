#' Get number of matches per probe
#'
#' @param matches Matches (generated using the `map_probes` function).
#' @export 
get_nr_matches_per_probe <- function(matches) {
  
  matches <- matches %>%
    dplyr::mutate(id = dplyr::case_when(
      strand %in% c("forward", "reverse_complement") ~ paste(Probe, strand, chr, start, sep = "_"),
      strand %in% c("reverse", "forward_complement") ~ paste(Probe, strand, chr, end, sep = "_")
    )) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(width == max(width)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!duplicated(id)) %>%
    dplyr::select(Probe, width) %>%
    dplyr::group_by(Probe, width) %>%
    dplyr::summarize(n = n()) %>%
    tidyr::spread(key = width, value = n) 
  
  matches[is.na(matches)] <- 0
  for(i in ncol(matches):3) {
    matches[[i-1]] <- matches[[i]] + matches[[i-1]]
  }
  colnames(matches)[2:ncol(matches)] <- paste0("bp", colnames(matches)[2:ncol(matches)])
  matches
}

  
