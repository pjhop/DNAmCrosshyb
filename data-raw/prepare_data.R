library(tidyverse)

# Probe Sequences
probe_sequences_EPIC <- read_tsv("data-raw/probe_sequences_EPIC.txt")
probe_sequences_450k <- read_tsv("data-raw/probe_sequences_450k.txt")

# Annotation files 
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
anno450k <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
annoEPIC <- data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
anno450k <- anno450k %>%
  dplyr::select(chr, pos, strand, Name, ProbeSeqA, ProbeSeqB, AddressA, AddressB, Type, NextBase, Color)
annoEPIC <- annoEPIC %>%
  dplyr::select(chr, pos, strand, Name, ProbeSeqA, ProbeSeqB, AddressA, AddressB,Type, NextBase, Color)
probes_epic <- annoEPIC$Name
probes_450k <- anno450k$Name

usethis::use_data(probe_sequences_EPIC,
                  probe_sequences_450k,
                  anno450k,
                  annoEPIC,
                  probes_epic,
                  probes_450k,
                  internal = TRUE,
                  overwrite = TRUE)



