
`DNAmCrosshyb` : Collection of functions useful in detecting
cross-reactive probes on Illumina 450k/EPIC DNA methylation arrays.

The following functions are implemented:

  - `map_probes` : Map 450k/EPIC probes to the reference genome,
    allowing for mismatches/INDELs.
  - `get_nr_matches_per_probe` : Get number of matches for per match
    length (on output of `map_probes` function).
  - `find_repeat_overlaps` : Check if matches overlap with repeats (UCSC
    masks) (on output of `map_probes` function).
  - `map_probes_sequence` : Map 450k/EPIC probes to a user-specified DNA
    sequence.  
  - `get_probe_overlaps` : Check if there are overlapping
    3’-subsequences in a set of probes.  
  - `get_OOB` : Get out-of-band normalized beta-values and/or
    intensities.  
  - `locusplot`: Plot a region of p-values.

### Installation

``` r
## Install from repo 
## Note: installation may take a little while, since the package includes
## probe sequences and annotation as internal data.
devtools::install_github("pjhop/DNAmCrosshyb")

## Install from source
install.packages("DNAmCrosshyb.tar.gz", repos=NULL, type="source")

# current version
packageVersion("DNAmCrosshyb")
```

### Mapping probes to the reference genome

#### Map Probes

Map some probes to reference genome (hg19), for each width from 30bp to
50bp (in steps of 5). Bisulfite-converted reference genomes can be
generated using the following scripts:
<https://github.com/pjhop/DNAmCrosshyb/blob/master/data-raw/bisulfite_convert_hg19.R>
and
<https://github.com/pjhop/DNAmCrosshyb/blob/master/data-raw/bisulfite_convert_hg38.R>

Bisulfite-converted genomes in the R .rds file format are available at:
<https://zenodo.org/record/4088020>

``` r
library(DNAmCrosshyb)
probes <- c("cg00005164", "cg12344104", "cg25521682", "cg27660038", "cg20554142", "cg03890998", "cg00947801")
matches <- map_probes(probes,
                      path = "genome_bs/hg19",
                      chromosomes = "all",
                      min_width = 30,
                      max_width = 50,
                      step_size = 5,
                      allow_mismatch = FALSE,
                      allow_INDEL = FALSE,
                      cores = 1
                     )
head(matches %>% data.frame())
```

    ##        Probe chr     start       end  strand next_base mismatch_pos
    ## 1 cg20554142   2 238821095 238821124 forward         A           NA
    ## 2 cg12344104   6  29976418  29976447 reverse         T           NA
    ## 3 cg12344104   6  29912333  29912362 reverse         Y           NA
    ## 4 cg03890998   9  15817371  15817400 reverse         A           NA
    ## 5 cg00005164  10  69575196  69575225 reverse         T           NA
    ## 6 cg25521682  15  23678247  23678276 forward         Y           NA
    ##            Type2 width channel
    ## 1 I_Unmethylated    30      IB
    ## 2             II    30    <NA>
    ## 3             II    30    <NA>
    ## 4 I_Unmethylated    30     OOB
    ## 5             II    30    <NA>
    ## 6             II    30    <NA>

#### Number of matches per probe

``` r
nr_matches <- get_nr_matches_per_probe(matches)
nr_matches
```

    ## # A tibble: 7 x 5
    ## # Groups:   Probe [7]
    ##   Probe       bp30  bp35  bp45  bp50
    ##   <chr>      <int> <int> <int> <int>
    ## 1 cg00005164     4     1     1     1
    ## 2 cg00947801     4     4     2     1
    ## 3 cg03890998     4     2     1     1
    ## 4 cg12344104     3     1     1     1
    ## 5 cg20554142     2     1     1     1
    ## 6 cg25521682     2     1     1     1
    ## 7 cg27660038     2     1     1     1

#### Overlap with repeats

``` r
matches <- find_repeat_overlaps(matches, genome_build = "hg19", min_overlap = "any")
head(matches %>% data.frame())
```

    ##        Probe chr     start       end  strand next_base mismatch_pos
    ## 1 cg20554142   2 238821095 238821124 forward         A           NA
    ## 2 cg12344104   6  29976418  29976447 reverse         T           NA
    ## 3 cg12344104   6  29912333  29912362 reverse         Y           NA
    ## 4 cg00947801   6  28925610  28925644 reverse         T           NA
    ## 5 cg00947801   6  96117896  96117940 reverse         T           NA
    ## 6 cg12344104   6  29693260  29693309 reverse         Y           NA
    ##            Type2 width channel repeat_overlap
    ## 1 I_Unmethylated    30      IB          FALSE
    ## 2             II    30    <NA>          FALSE
    ## 3             II    30    <NA>          FALSE
    ## 4 I_Unmethylated    35      IB           TRUE
    ## 5 I_Unmethylated    45      IB           TRUE
    ## 6             II    50    <NA>          FALSE

### Map probes to a user-specified DNA sequence

#### Assume the repeat is methylated

Map probes to the *C9orf72* hexanucleotide repeat:

``` r
# Map probes to the C9orf72 hexanucleotide repeat
repeat_sequence <- paste(rep("GGCCCC", 10), collapse="")

matches_c9 <- map_probes_sequence(sequence = repeat_sequence, next_base = "G", prev_base = "C",
                            array = "450k", min_width = 10, max_width = 25, allow_indel = FALSE, 
                            allow_mismatch = FALSE, step_size = 1, use_Y = FALSE, methylation_status = "methylated")
head(matches_c9 %>% data.frame())
```

    ##        Probe start end  strand width sbe_site mismatch_pos max_mismatch_pos
    ## 1 cg05921207     1  10 forward    10        C           NA               NA
    ## 2 cg05921207     7  16 forward    10        C           NA               NA
    ## 3 cg05921207    13  22 forward    10        C           NA               NA
    ## 4 cg05921207    19  28 forward    10        C           NA               NA
    ## 5 cg05921207    25  34 forward    10        C           NA               NA
    ## 6 cg05921207    31  40 forward    10        C           NA               NA
    ##   n_mismatch indel_pos width_incl_indel sequence_bs Type2 channel
    ## 1          0        NA               NA  GGTTTCGGTT    II     OOB
    ## 2          0        NA               NA  GGTTTCGGTT    II     OOB
    ## 3          0        NA               NA  GGTTTCGGTT    II     OOB
    ## 4          0        NA               NA  GGTTTCGGTT    II     OOB
    ## 5          0        NA               NA  GGTTTCGGTT    II     OOB
    ## 6          0        NA               NA  GGTTTCGGTT    II     OOB

#### Allow mismatch (\>=6 bp from 3’-end of probe)

``` r
# Map probes to the C9orf72 hexanucleotide repeat
repeat_sequence <- paste(rep("GGCCCC", 10), collapse="")

matches_c9 <- map_probes_sequence(sequence = repeat_sequence, next_base = "G", prev_base = "C",
                            array = "450k", min_width = 10, max_width = 25, allow_indel = FALSE, 
                            allow_mismatch = TRUE, min_distance = 6,
                            step_size = 1, use_Y = FALSE, methylation_status = "methylated")
head(matches_c9 %>% dplyr::filter(n_mismatch > 0) %>% data.frame())
```

    ##        Probe start end  strand width sbe_site mismatch_pos max_mismatch_pos
    ## 1 cg01927625     1  10 forward    10        C            6                6
    ## 2 cg05921207     1  10 forward    10        C            6                6
    ## 3 cg13788061     1  10 forward    10        C            6                6
    ## 4 cg24751886     1  10 forward    10        C            6                6
    ## 5 cg01877196     1  10 forward    10        C            6                6
    ## 6 cg03887614     1  10 forward    10        C            6                6
    ##   n_mismatch indel_pos width_incl_indel sequence_bs Type2 channel
    ## 1          1        NA               NA  GGTTTCGGTT    II     OOB
    ## 2          1        NA               NA  GGTTTCGGTT    II     OOB
    ## 3          1        NA               NA  GGTTTCGGTT    II     OOB
    ## 4          1        NA               NA  GGTTTCGGTT    II     OOB
    ## 5          1        NA               NA  GGTTTCGGTT    II     OOB
    ## 6          1        NA               NA  GGTTTCGGTT    II     OOB

#### Run in parallel

``` r
# Map probes to the C9orf72 hexanucleotide repeat
repeat_sequence <- paste(rep("GGCCCC", 10), collapse="")

matches_c9_2 <- map_probes_sequence(sequence = repeat_sequence, next_base = "G", prev_base = "C",
                            array = "450k", min_width = 10, max_width = 25, allow_indel = FALSE, 
                            allow_mismatch = TRUE, min_distance = 6,
                            step_size = 1, use_Y = FALSE, methylation_status = "methylated",
                            cores = 4)
head(matches_c9_2 %>% data.frame())
```

    ##        Probe start end  strand width sbe_site mismatch_pos max_mismatch_pos
    ## 1 cg05921207     1  10 forward    10        C           NA               NA
    ## 2 cg05921207     7  16 forward    10        C           NA               NA
    ## 3 cg05921207    13  22 forward    10        C           NA               NA
    ## 4 cg05921207    19  28 forward    10        C           NA               NA
    ## 5 cg05921207    25  34 forward    10        C           NA               NA
    ## 6 cg05921207    31  40 forward    10        C           NA               NA
    ##   n_mismatch indel_pos width_incl_indel sequence_bs Type2 channel
    ## 1          0        NA               NA  GGTTTCGGTT    II     OOB
    ## 2          0        NA               NA  GGTTTCGGTT    II     OOB
    ## 3          0        NA               NA  GGTTTCGGTT    II     OOB
    ## 4          0        NA               NA  GGTTTCGGTT    II     OOB
    ## 5          0        NA               NA  GGTTTCGGTT    II     OOB
    ## 6          0        NA               NA  GGTTTCGGTT    II     OOB

### Sequence overlap between a set of probes

The `get_probe_overlaps` can be used to check for overlapping
3’-subsequences in a set of probes. Here we apply this function to
probes that map to the C9 repeat (identified above) and some random
probes

``` r
matches_c9_20bp <- matches_c9 %>% dplyr::filter(width >= 20)
set.seed(10)
random_probes <- sample(DNAmCrosshyb:::anno450k$Name, size = 10)

# Calculate overlaps
overlaps <- get_probe_overlaps(c(unique(matches_c9_20bp$Probe), random_probes))

## Note this also includes overlap between _unmethylated and _methylated bead types of type I probes
overlaps %>% dplyr::arrange(desc(overlaps$overlap)) %>% head(10)
```

    ## # A tibble: 10 x 6
    ##    Probe_Index Bead_Index      Probe_Target Bead_Target     overlap mismatch_pos
    ##    <chr>       <chr>           <chr>        <chr>             <int>        <int>
    ##  1 cg00754896  cg00754896_unm… cg00754896   cg00754896_met…      26            1
    ##  2 cg00754896  cg00754896_met… cg00754896   cg00754896_unm…      26            1
    ##  3 cg09994391  cg09994391_unm… cg01587390   cg01587390_unm…      21            9
    ##  4 cg12419491  cg12419491_unm… cg01587390   cg01587390_unm…      21            7
    ##  5 cg01587390  cg01587390_unm… cg09994391   cg09994391_unm…      21            9
    ##  6 cg01587390  cg01587390_unm… cg12419491   cg12419491_unm…      21            7
    ##  7 cg09994391  cg09994391_met… cg01587390   cg01587390_met…      21            9
    ##  8 cg12419491  cg12419491_met… cg01587390   cg01587390_met…      21            7
    ##  9 cg01587390  cg01587390_met… cg09994391   cg09994391_met…      21            9
    ## 10 cg01587390  cg01587390_met… cg12419491   cg12419491_met…      21            7

Plot:

``` r
ggplot(overlaps %>% filter(Probe_Index != Probe_Target), 
       aes(x = overlap)) +
  geom_histogram() +
  theme_classic() + 
  xlab("Pair-wise probe overlap (bp)") +
  theme(text = element_text(size=13))
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Extracting and normalizing out-of-band (OOB) intensities and betas

Extracting OOB betas from an RGset (*minfi*).

``` r
# For this example, the minfiData package is used.
library(minfiData)

# Load data
baseDir <- system.file("extdata", package="minfiData")
samplesheet <- read.metharray.sheet(baseDir)
```

    ## [1] "/Users/phop2/Library/R/3.5/library/minfiData/extdata/SampleSheet.csv"

``` r
rgset <- read.metharray(samplesheet$Basename, extended=TRUE)

# Get normalized OOB betas (only keep beta-values, other options are 'intensity' and 'both')
oob_betas_normalized <- get_OOB(rgset, normalized = TRUE, keep = "beta")

# Show subset
oob_betas_normalized$beta[1:5,1:5]
```

    ##            5723646052_R02C02 5723646052_R04C01 5723646052_R05C02
    ## cg02004872         0.5519658         0.6333704         0.6377151
    ## cg02050847         0.7833947         0.4767726         0.5807178
    ## cg02233190         0.4110656         0.3269978         0.3708729
    ## cg02494853         0.3098046         0.3583913         0.5464111
    ## cg02842889         0.4197092         0.4690575         0.4869732
    ##            5723646053_R04C02 5723646053_R05C02
    ## cg02004872         0.6122822         0.6137673
    ## cg02050847         0.5693312         0.3288250
    ## cg02233190         0.4085624         0.5167526
    ## cg02494853         0.6951439         0.3906633
    ## cg02842889         0.5948174         0.4203944
