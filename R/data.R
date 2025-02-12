
#' Phase data
#'
#' Example phased SNP data, in binary form, with individuals in columns
#' and markers in rows. The dataset is a subset of a genomic SNP-phased data
#' from the Center for Genomic Resources Netherlands (CGN).
#'
#' @format ## `exPhase`
#' A data frame with 1229 rows and 2458 columns. Each column is one chromosome
#' of a single individual, so there are ploidy*individuals columns.
#'
"exPhase"

#' Individual labels
#'
#' Example individual labels. Anonimized labels related to the `exPhase` dataset.
#' The order of the labels relates to the order of the individuals in the `exPhase` data.
#'
#' @format ## `exInds`
#' A character vector with individual names.
#'
"exInds"

#' Position of markers
#'
#' Example position information related to the `exPhase` dataset. The data.frame contains
#' four columns, "CHROM" (chromosome name), "POS" (genomic position, in bp), "REF" (reference
#' allele) and "ALT" (alternative allele)
#'
#' @format ## `exPos`
#' A data.frame with 4 columns and 1229 rows, one for each marker in the `exPhase` dataset.
#'
"exPos"
