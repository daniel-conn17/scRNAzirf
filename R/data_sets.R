#' Names of transcription factors for ZIRFs vignette
#'
#' This vector contains the vector of transcription factor
#' names used in the ZIRFs vignette.
#' @format A character vector with transcription factor names.
"TF_names"

#' Un-normalized count data for ZIRFs vignette
#'
#' This matrix contains the filtered un-normalized count data
#' used in ZIRFs vignette.
#' @format A matrix with genes along the rows and cells along the
#' columns. As this is a toy data set, there are only 70 genes
#' and 865 cells.
"countDat"

#' Normalization size factors
#'
#' This data.frame contains the normalization size factors for the
#' single cell RNA-seq data in the ZIRFs vignette. The size factors are
#' used in the zero-inflation component of ZIRFs. The size factors
#' were derived using scran.
#' @format A data.frame with one column containing the size factors.
"sizeFactors"

#' Normalized count data for ZIRFs vignette
#'
#' This matrix contains the normalized count data used in ZIRFs vignette.
#' This matrix is arrived at by applying the following transformation
#' to the un-normalized count data: $log(count_{ij}/size_factor_{ij} +)$.
#' @format A matrix with genes along the rows and cells along the
#' columns. As this is a toy data set, there are only 70 genes
#' and 865 cells.
"normDat"

#' An additional set of target genes used in the ZIRFs vignette.
#'
#' This character vector contains additional genes used as targets in the
#' ZIRFs vignette. The vignette provides a toy example of a ZIRFs run.
#' The target genes were selected because in a full run of SCENIC they
#' were found to be active in multiple regulons.
#' @format A character vector of gene names.
"target_names"
