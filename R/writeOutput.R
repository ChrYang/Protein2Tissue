#' Write Enrichment Analysis Output to File
#'
#' @param out An \code{Enrichment} S4 object containing the results of
#'   tissue enrichment analysis.
#' @param type A character string specifying the output file format.
#'   Must be either \code{"csv"} or \code{"tsv"}. Default is \code{"csv"}.
#' @param filename A character string specifying the output file path.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing a file to disk.
#'#' @examples
#' \dontrun{
#' writeOutput(result,type = "csv","~/Enrichment_Results.csv")
#' }
#'
#' @keywords write
#' @noRd
writeOutput <- function(out, type = "csv", filename) {
  
  # --- input validation ---
  stopifnot(
    "'out' must be an Enrichment object" = is(out, "Enrichment"),
    "'type' must be a character string"  = is.character(type),
    "'filename' must be a character string" = is.character(filename)
  )
  
  type <- match.arg(type, choices = c("csv", "tsv"))
  
  # --- write output ---
  if (type == "csv") {
    utils::write.csv(out@input,
                     file      = filename,
                     row.names = FALSE)
    
  } else if (type == "tsv") {
    utils::write.table(out@input,
                       file      = filename,
                       row.names = FALSE,
                       sep       = "\t")
  }
  
  message("Output written to: ", filename)
  invisible(NULL)
}
