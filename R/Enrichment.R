#' Enrichment-class
#'
#' An S4 class to store enrichment analysis results, including
#' the enrichment results of input and background protein lists.
#'
#' @slot input An object containing enrichment results for the input set
#' @slot background An object containing enrichment results for the background set
#'
#' @export
setClass(
  "Enrichment",
  slots = c(
    input      = "ANY",
    background = "ANY"
  )
)