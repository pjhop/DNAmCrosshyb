% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_repeat_overlaps.R
\name{find_repeat_overlaps}
\alias{find_repeat_overlaps}
\title{Check if matches overlap with repeat sequences}
\usage{
find_repeat_overlaps(matches, genome_build, min_overlap = "any")
}
\arguments{
\item{matches}{Matches (generated using the 'map_probes' function)}

\item{genome_build}{Genome build (hg19 or hg38)}

\item{min_overlap}{Minimum overlap with repeat. Defaults to "any".}
}
\description{
Check if matches overlap with repeat sequences
}
