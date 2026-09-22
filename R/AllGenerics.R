#' @include AllClasses.R
NULL

#' @title Get input
#' @description Accessor for the input slot of an Enrichment object
#' @param object An Enrichment object
#' @return The input data frame
#' @examples
#' obj <- methods::new("Enrichment",
#'   input      = data.frame(Gene = c("TP53", "BRCA1")),
#'   background = data.frame(Gene = c("TP53", "BRCA1", "EGFR"))
#' )
#' getInput(obj)
#' @exportMethod getInput
setGeneric("getInput", function(object) standardGeneric("getInput"))

#' @title Get background
#' @description Accessor for the background slot of an Enrichment object
#' @param object An Enrichment object
#' @return The background data frame
#' @examples
#' obj <- methods::new("Enrichment",
#'   input      = data.frame(Gene = c("TP53", "BRCA1")),
#'   background = data.frame(Gene = c("TP53", "BRCA1", "EGFR"))
#' )
#' getBackground(obj)
#' @exportMethod getBackground
setGeneric("getBackground", function(object) standardGeneric("getBackground"))
