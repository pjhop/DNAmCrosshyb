#' Locusplot
#'
#' @param cpg CpG-site to plot
#' @param stats data.frame containing association p-values. Should contain a 'Probe' and a 'p' column.
#' @param windowsize Size (bp) of the window surrounding the CpG-site that should be plotted.
#' @param threshold Significance threshold. By default will use bonferonni threshold, using the number of probes in the `stats` file.
#' @param significance_line Plot a significane line(TRUE/FALSE)?
#' @param highlight Highlight significant sites?
#' @param text_size Text size.
#' @param point_size Point size.
#' @import ggplot2
#' @export

locusplot <- function(cpg, stats, windowsize, threshold = NULL,
                      significance_line = FALSE, highlight = TRUE,
                      text_size = 13, point_size = 2.5) {
  
  if(all(stats$Probe %in% DNAmCrosshyb:::probes_450k)) {
    array <- "450k"
  } else if(all(stats$Probe %in% DNAmCrosshyb:::probes_epic)) {
    array <- "EPIC"
  } else {
    stop("Not all probes are recognized, are they not from the 450k or EPIC array, or maybe a combination of the two?")
  }
  
  # Annotation 
  if(array == "450k") {
    anno <- DNAmCrosshyb:::anno450k
  } else if(array == "EPIC") {
    anno <- DNAmCrosshyb:::annoEPIC
  }
  
  # If threshold = NULL, use bonferroni
  if(is.null(threshold)) {
    threshold <- 0.05 / nrow(stats)
  }
  
  # Select necessary info from stats + add annotation 
  stats <- stats %>% 
    dplyr::select(Probe, p) %>%
    dplyr::mutate(logp = -log10(p)) %>%
    dplyr::left_join(anno[,c("Name", "chr", "pos")], by = c("Probe" = "Name"))
  
  # Select CpG 
  dmp <- stats %>% 
    dplyr::filter(Probe == cpg)
  
  # Selection region basd on windowsize
  region <- stats %>% 
    dplyr::filter(chr == dmp$chr,
                  pos <= (dmp$pos + round(windowsize/2)),
                  pos >= (dmp$pos - round(windowsize/2))
    ) %>%
    dplyr::arrange(pos) %>%
    dplyr::mutate(Sig = ifelse(p < threshold, TRUE, FALSE))
  
  max <- max(region$logp)
  
  plot <- ggplot(region, aes(x = pos, y = logp, color = if(highlight) Sig else 1)) +
    geom_point(data = region %>% dplyr::filter(Probe != cpg), size = point_size) +
    geom_point(data = region %>% dplyr::filter(Probe == cpg), shape = 9, color = "red", size = point_size) + 
    xlab(region$chr[1]) +
    ylab(expression(-log[10](P))) +
    ggtitle(cpg) +
    scale_x_continuous(breaks = c(dmp$pos - round(windowsize/2) + 0.025 * windowsize, dmp$pos, dmp$pos + round(windowsize/2) - 0.025 * windowsize),
                       labels = c(sprintf("- %s bp", round(windowsize/2)), as.character(dmp$pos), sprintf("+ %s bp", round(windowsize/2))),
                       limits = c(dmp$pos - round(windowsize/2), dmp$pos + round(windowsize/2) )) +
    ylim(0, max) +
    theme_classic()  +
    theme(text = element_text(size = text_size),
          axis.text = element_text(size = text_size),
          legend.position = "none")
  if(highlight) {
    plot <- plot + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))
  }
  if(significance_line) {
    plot <- plot + geom_hline(yintercept = -log10(threshold), linetype = "dashed")
  }
  
  plot 
}

