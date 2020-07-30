#' Get out-of-band betas/intensities.
#'
#' @param rgset RGset (minfi).
#' @param normalized Normalize intensity-values (using the wateRmelon 'nasen()' function)?
#' @param keep What to keep? ('both', 'intensity', 'beta')
#' @export

get_OOB <- function(rgset, 
                    normalized = TRUE,
                    keep = c("both", "beta", "intensity")
                    ) {
  # Checks
  minfi:::.isRGOrStop(rgset)
  keep <- match.arg(keep)
  stopifnot(is.logical(normalized))
  
  # Guess array type (function adapted from minfi)
  array <- guessArrayTypes(nrow(rgset))
  
  if(array == "450k") {
    anno <- DNAmCrosshyb:::anno450k
  } else if(array == "EPIC") {
    anno <- DNAmCrosshyb:::annoEPIC
  } else {
    stop("Unknown array!")
  }
  
  # Annotation
  anno <- anno %>% dplyr::select(Name, Type, AddressA, AddressB, Color)
  anno_typeII <- anno %>% dplyr::filter(Type == "II")
  anno_grn <- anno %>% dplyr::filter(Type == "I", Color == "Grn")
  anno_red <- anno %>% dplyr::filter(Type == "I", Color == "Red")
  
 
  
  # Red
  red <- minfi::getRed(rgset)
  # Green
  green <- minfi::getGreen(rgset)
  rm(rgset); gc()
  
  # Probes on different versions of EPIC array differ
  if(array == "EPIC") {
    intersec_A <- intersect(rownames(green), c(anno_red$AddressA, anno_grn$AddressA))
    intersec_B <- intersect(rownames(green), c(anno_red$AddressB, anno_grn$AddressB))
    intersec_II <- intersect(rownames(green), anno_typeII$AddressA)
    
    anno_grn <- anno_grn %>% dplyr::filter(AddressA %in% intersec_A, AddressB %in% intersec_B)
    anno_red <- anno_red %>% dplyr::filter(AddressA %in% intersec_A, AddressB %in% intersec_B)
    anno_typeII <- anno_typeII %>% dplyr::filter(AddressA %in% intersec_II)
  }
  
  # Select OOB values in red channel
  red_M <- red[anno_grn$AddressB,, drop = FALSE]
  red_U <- red[anno_grn$AddressA,, drop = FALSE]
  rownames(red_M) <- anno_grn$Name
  rownames(red_U) <- anno_grn$Name
  
  # Select OOB values in green channel
  green_M <- green[anno_red$AddressB,, drop = FALSE]
  green_U <- green[anno_red$AddressA,, drop = FALSE]
  rownames(green_M) <- anno_red$Name
  rownames(green_U) <- anno_red$Name
  
  # Workaround: function does not work if no type II probes are included
  # Therefore, some type II probes are added in, which are removed afterwards.
  anno_typeII <- anno_typeII[1:10,]
  M_II <- green[anno_typeII$AddressA,, drop = FALSE]
  U_II <- red[anno_typeII$AddressA,, drop = FALSE]
  rownames(M_II) <- anno_typeII$Name
  rownames(U_II) <- anno_typeII$Name
  rm(red, green)
  
  # Combine M and U 
  M <- rbind(red_M, green_M, M_II)
  U <- rbind(red_U, green_U, U_II)
  onetwo <- c(anno_grn$Type, anno_red$Type, anno_typeII$Type)
  rm(red_M, green_M, M_II, red_U, green_U, U_II)
  
  if(normalized) {
    values <- wateRmelon::nasen(M, U, 
                      onetwo = onetwo, ret2 = if(keep %in% c("both", "intensity")) TRUE else FALSE)
    if(keep == "beta") {
      values <- list(beta = values, total = NULL)
    }
    if(keep == "intensity" || keep == "both") {
      values$total <- values$methylated + values$unmethylated
    } 
    values$methylated <- values$unmethylated <- NULL
    if(keep == "intensity") values$beta <- NULL
  } else {
    values <- list()
    if(keep == "beta" || keep == "both") {
      values$beta <- M / (M + U + 100)
    } 
    if(keep == "intensity" || keep == "both") {
      values$total <- M + U
    } 
  }
  # Remove type II 
  if(keep == "beta" || keep == "both") {
    values$beta <- values$beta[!rownames(values$beta) %in% anno_typeII$Name,, drop = FALSE]
  }
  if(keep == "intensity" || keep == "both") {
    values$total <- values$total[!rownames(values$total) %in% anno_typeII$Name,, drop = FALSE]
  }
  values
}


## Function taken from minfi github (https://github.com/hansenlab/minfi/blob/master/R/read.meth.R)
## To infer whether data is EPIC/450k/27k
guessArrayTypes <- function(nProbes) {
  if (nProbes >= 622000 && nProbes <= 623000) {
    array = "450k"
  } else if (nProbes >= 1050000 && nProbes <= 1053000) {
    # NOTE: "Current EPIC scan type"
    array = "EPIC"
  } else if (nProbes >= 1032000 && nProbes <= 1033000) {
    # NOTE: "Old EPIC scan type"
    array = "EPIC"
  } else if (nProbes >= 54000 && nProbes <= 56000) {
    array = "27k"
  } else {
    array = "Unknown"
    array
  }
}
  