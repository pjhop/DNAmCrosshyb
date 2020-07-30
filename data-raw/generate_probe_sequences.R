# Create a table containing all possible probe sequences for the 450k/EPIC arrays.
# For Type I probes there are two sequences per site (methylated/unmethylated),
# For Type II probes there are up to 8 (2^3) possible sequences per site,
# since the probe can contain up to 3 'R' bases (which hybridize to either a
# 'A' or a 'C' base.


# Libraries ---------------------------------------------------------------------
library(stringr)
library(tidyverse)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Functions ---------------------------------------------------------------------

# Get all the possible probe sequences for a type II probe (2^n)
permute <- function(string) {
  locs <- stringr::str_locate_all(string, "R")[[1]][,'start']
  if(length(locs) == 0) {
    return(string)
  } else {
    string <- stringr::str_split(string, pattern = "")[[1]]
    combs <- t(expand.grid(rep(list(c('A', 'G')), length(locs))))
    strings <- purrr::map_chr(1:ncol(combs), .f = replace, string = string, locations = locs, combinations = combs)
    return(strings)
  }
}

replace <- function(i, string, locations, combinations) {
  string[locations] <- combinations[,i]
  stringr::str_c(string, collapse = "")
}

# Create tibble with all probe sequences ---------------------------------------------------------------------

## Annotations
# EPIC
data(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annoEPIC <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))

# 450k
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_450k <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

## EPIC
# Create one dataframe with probes + information
annoEPIC_typeII <- annoEPIC %>% dplyr::filter(Type == "II")
typeII_probes <- purrr::map(annoEPIC_typeII$ProbeSeqA, permute) # Takes a minute or two
names(typeII_probes) <- annoEPIC_typeII$Name
lengths <- purrr::map_dbl(typeII_probes, length)
cpgs <- rep(names(typeII_probes), lengths)

probe_sequences_EPIC <- tibble::tibble(Probe = cpgs, Probe_sequence = unlist(typeII_probes), Type = "II", ProbeSeq = "ProbeSeqA")

# Add ProbeSeqA and ProbeSeqB for type I probes
annoEPIC_typeI <- annoEPIC %>% dplyr::filter(Type == "I")
probe_sequences_EPIC <- dplyr::bind_rows(probe_sequences_EPIC, tibble::tibble(Probe = annoEPIC_typeI$Name,
                                                                              Probe_sequence = annoEPIC_typeI$ProbeSeqA, Type = "I", ProbeSeq = "ProbeSeqA"))
probe_sequences_EPIC <- dplyr::bind_rows(probe_sequences_EPIC, tibble::tibble(Probe = annoEPIC_typeI$Name,
                                                                              Probe_sequence = annoEPIC_typeI$ProbeSeqB, Type = "I", ProbeSeq = "ProbeSeqB"))
probe_sequences_EPIC$i <- 1:nrow(probe_sequences_EPIC)

## 450k
# Create one dataframe with probes + information
anno_450k_typeII <- anno_450k %>% dplyr::filter(Type == "II")
typeII_probes <- purrr::map(anno_450k_typeII$ProbeSeqA, permute) # Takes a minute or two
names(typeII_probes) <- anno_450k_typeII$Name
lengths <- purrr::map_dbl(typeII_probes, length)
cpgs <- rep(names(typeII_probes), lengths)

probe_sequences_450k <- tibble::tibble(Probe = cpgs, Probe_sequence = unlist(typeII_probes), Type = "II", ProbeSeq = "ProbeSeqA")

# Add ProbeSeqA and ProbeSeqB for type I probes
anno_450k_typeI <- anno_450k %>% dplyr::filter(Type == "I")
probe_sequences_450k <- dplyr::bind_rows(probe_sequences_450k, tibble::tibble(Probe = anno_450k_typeI$Name,
                                                                              Probe_sequence = anno_450k_typeI$ProbeSeqA, Type = "I", ProbeSeq = "ProbeSeqA"))
probe_sequences_450k <- dplyr::bind_rows(probe_sequences_450k, tibble::tibble(Probe = anno_450k_typeI$Name,
                                                                              Probe_sequence = anno_450k_typeI$ProbeSeqB, Type = "I", ProbeSeq = "ProbeSeqB"))
probe_sequences_450k$i <- 1:nrow(probe_sequences_450k)

# Save files ---------------------------------------------------------------------
readr::write_tsv(probe_sequences_EPIC, "data-raw/probe_sequences_EPIC.txt")
readr::write_tsv(probe_sequences_450k, "data-raw/probe_sequences_450k.txt")