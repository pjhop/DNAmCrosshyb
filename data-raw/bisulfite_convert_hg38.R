# Libraries ---------------------------------------------------------------------
library(tidyverse)
library(stringr)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)

# Bisulfite conversion function ---------------------------------------------------------------------
bisulfite_convert <- function(sequence, nextbase, methylation_status = c("methylated", "unmethylated"),
                              use_Y = FALSE) {
  stopifnot(is.character(sequence))
  methylation_status <- match.arg(methylation_status)

  convert_vec <- c("C" = "T", "G" = "G", "A" = "A", "T" = "T", "N" = ".")
  sequence <- str_split(sequence, pattern = "")[[1]]

  if(use_Y) {
    sequence_bs <- vector("character", length = length(sequence))
    for(i in seq_along(sequence)) {
      if(i == length(sequence) && sequence[i] == 'C') {
        if(nextbase == "G") {
          sequence_bs[i] <- "Y"
        } else {
          sequence_bs[i] <- "T"
        }
      } else if(sequence[[i]] %in% c("T", "A", "G", "N")) {
        sequence_bs[[i]] <- unname(convert_vec[sequence[[i]]])
      } else if(sequence[[i]] == "C" && sequence[[i + 1]] != "G") {
        sequence_bs[[i]] <- "T"
      } else {
        sequence_bs[[i]] <- "Y"
      }
    }
    sequence_bs <- str_c(sequence_bs, collapse = "")
    sequence_bs
  } else if(methylation_status == "methylated") {
    sequence_bs <- vector("character", length = length(sequence))
    for(i in seq_along(sequence)) {
      if(i == length(sequence) && sequence[i] == 'C') {
        if(nextbase == "G") {
          sequence_bs[i] <- "C"
        } else {
          sequence_bs[i] <- "T"
        }
      } else if(sequence[[i]] %in% c("T", "A", "G")) {
        sequence_bs[[i]] <- unname(convert_vec[sequence[[i]]])
      } else if(sequence[[i]] == "C" && sequence[[i+1]] != "G") {
        sequence_bs[[i]] <- "T"
      } else {
        sequence_bs[[i]] <- "C"
      }
    }
    sequence_bs <- str_c(sequence_bs, collapse = "")
    sequence_bs
  } else if(methylation_status == 'unmethylated') {
    sequence_bs <- convert_vec[sequence]
    sequence_bs <- str_c(sequence_bs, collapse = "")
    sequence_bs
  }
}

# Run conversion for all chromosomes ---------------------------------------------------------------------

run <- function(chr, outdir, direction = "forward") {

  print(sprintf("Bisulfite converting %s in the %s direction", chr, direction))
  # Get current chromosome
  genome <- Hsapiens[[chr]]
  if(direction == "forward") {
    genome_bs <- bisulfite_convert(as.character(genome), nextbase="A", use_Y = TRUE)
  } else if(direction == "reverse") {
    genome_bs <- bisulfite_convert(as.character(reverseComplement(genome)), nextbase="A", use_Y = TRUE)
  }
  genome_bs <- DNAString(genome_bs)

  # Save
  saveRDS(genome_bs, file = sprintf("%s%s.rds", outdir, chr))
  return(NULL)
}

## Run in parallel
chromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))

# forward
null <- bplapply(chromosomes, FUN = run,
                 outdir = "../data/genome_bs/hg38/forward/",
                 direction = "forward", BPPARAM = MulticoreParam(workers = 4))
# reverse
null <- bplapply(chromosomes, FUN = run,
                outdir = "../data/genome_bs/hg38/reverse/",
                direction = "reverse", BPPARAM = MulticoreParam(workers = 4))
