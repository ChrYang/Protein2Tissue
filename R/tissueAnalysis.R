#' Tissue Enrichment Analysis
#'
#' Performs tissue enrichment analysis by mapping input proteins and background
#' proteins against a reference expression database. Proteins not found in the
#' database are removed with a warning, and the matched results are returned
#' as an \code{Enrichment} object for downstream analysis.
#'
#' @param input A character vector of protein identifiers to be analysed
#'   (e.g. gene symbols or Ensembl IDs).
#' @param background A character vector of background protein identifiers
#'   representing the universe of proteins. Defaults to
#'   \code{"Default"}.
#' @param typeKey A character string specifying the identifier type for the
#'   input proteins. Must match a column name in the reference database. 
#'   Default is \code{"Gene"}.
#' @param typeKeyBg A character string indicating the identifier type used for
#'   the background proteins. This must be "Gene", "Ensembl", or "EntrezID".
#'   The default is \code{"Ensembl"}.
#' @param database Defaults to \code{ee_tb} (Human Protein Atlas data bundled 
#' with this package). Can be replaced with a custom data frame provided it 
#' contains columns matching \code{typeKey} and \code{typeKeyBg}.  
#'
#'
#' @return An \code{Enrichment} S4 object containing:
#'   \itemize{
#'     \item \code{input} — a data frame of matched input proteins with
#'       reference database annotations
#'     \item \code{background} — a data frame of matched background proteins
#'       with reference database annotations
#'   }
#'
#' @details
#' The function maps input and background protein lists against the internal
#' reference expression database \code{ee_tb}. Proteins not found in the
#' database are removed before analysis. If no proteins are found in the
#' database, the function stops with an informative error message.
#'
#' Note that \code{typeKey} and \code{typeKeyBg} may differ because input
#' proteins and background proteins can use different identifier systems
#' (e.g. gene symbols vs Ensembl IDs). Both are mapped to the same reference
#' database \code{ee_tb} which contains multiple identifier columns.
#'
#' @examples
#' \dontrun{
#' # using gene symbols for input, Ensembl IDs for background
#' result <- tissueAnalysis(
#'     input      = c("TP53", "BRCA1", "EGFR"),
#'     background = c("347", "948", "1050", "2167", "2819", "3625"),
#'     typeKey    = "Gene",
#'     typeKeyBg  = "EntrezID",
#'     database   = NULL
#' )
#'
#' # access results
#' result@input
#' result@background
#' }
#'
#' @seealso \code{\link{enrichmentAnalysis}}, \code{\link{Enrichment-class}}
#'
#' @importFrom dplyr left_join
#' @export
tissueAnalysis <- function(input,
                           background = "Default",
                           typeKey    = "Gene",
                           typeKeyBg  = "Ensembl",
                           database   = NULL) {
    
    if (is.null(database)) {
        database <- get(utils::data("ee_tb",
                                    package  = "Protein2Tissue",
                                    envir    = environment()))
    }
  
  
    if (identical(background, "Default")){
        background <- database[,typeKeyBg]
    }
  
  
    # --- input validation ---
    stopifnot(
        "'input' must be a non-empty vector"                              = length(input) > 0,
        "'background' must be a non-empty vector or left as default"      = length(input) > 0,
        "'typeKey' must be a character string"                            = is.character(typeKey),
        "'typeKeyBg' must be a character string"                          = is.character(typeKeyBg)
    )
    
    stopifnot(
        "'typeKey' not found in database"   = typeKey   %in% colnames(database),
        "'typeKeyBg' not found in database" = typeKeyBg %in% colnames(database)
    )
    
    
    # --- process background ---
    backgroundDf <- .processProteinList(
        proteinList = background,
        typeKey     = typeKeyBg,
        listLabel   = "background",
        database    = database 
    )
    
    resultBg <- dplyr::left_join(backgroundDf, database, by = typeKeyBg)
    
    # --- process input ---
    inputDf <- .processProteinList(
        proteinList = input,
        typeKey     = typeKey,
        listLabel   = "input",
        database    = database 
    )
    
    resultInput <- dplyr::left_join(inputDf, database, by = typeKey)
    
    # --- return Enrichment S4 object ---
    methods::new("Enrichment",
                 input      = resultInput,
                 background = resultBg)
}


