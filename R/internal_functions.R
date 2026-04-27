## internal functions

#' Process and Validate a Protein List Against the Reference Database
#'
#' Internal helper function that validates a protein list against the reference
#' expression database like the default one \code{ee_tb}, removes unmatched proteins, 
#' and returns a cleaned data frame.
#'
#' @param proteinList A character vector of protein identifiers.
#' @param typeKey A character string specifying the identifier column in
#'  the database like the default one \code{ee_tb}.
#' @param listLabel A character string used in messages to identify whether
#'   this is the \code{"input"} or \code{"background"} list.
#'
#' @return A single-column data frame of matched protein identifiers with
#'   column name set to \code{typeKey}.
#'
#' @keywords internal
#' @noRd
.processProteinList <- function(proteinList, typeKey, listLabel, database) {
    
    # convert to clean data frame
    proteinDf           <- data.frame(id = as.character(proteinList),
                                      stringsAsFactors = FALSE)
    colnames(proteinDf) <- typeKey
    
    # remove NAs and duplicates
    proteinDf <- proteinDf[!is.na(proteinDf[, typeKey]), , drop = FALSE]
    proteinDf <- proteinDf[!duplicated(proteinDf[, typeKey]), , drop = FALSE]
    
    # match against reference database
    matched <- proteinDf[, typeKey] %in% database[, typeKey]
    nUnmatched <- sum(!matched)
    
    if (all(matched)) {
        message("All ", listLabel, " proteins found in the database.")
        
    } else if (!any(matched)) {
        stop(
            "No ", listLabel, " proteins found in the database. ",
            "Please check the ", listLabel, " list and the assigned key type '",
            typeKey, "'.",
            call. = FALSE
        )
        
    } else {
        warning(
            nUnmatched, " protein(s) in the ", listLabel,
            " list not found in the database and have been removed.",
            call. = FALSE
        )
        proteinDf <- proteinDf[matched, , drop = FALSE]
    }
    
    proteinDf
}




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
#'
#' @keywords internal
#' @noRd
.writeOutput <- function(out, type = "csv", filename) {
    
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



#' Convert Enrichment Object from Wide to Long Format
#'
#' @param out An \code{Enrichment} S4 object containing wide-format input and
#'   background data frames.
#' @param removeNA A logical value indicating whether to remove \code{NA}
#'   signal values. If \code{FALSE}, rows with \code{NA} signals are retained.
#'   Default is \code{FALSE}.
#' @param inclusion A character vector specifying which signal categories to
#'   include. Valid options are \code{"Tissue enhanced"},
#'   \code{"Group enriched"}, and \code{"Tissue enriched"}. Use \code{"All"}
#'   to include all three categories. Default is \code{"All"}.
#'
#' @return An \code{Enrichment} S4 object with input and background data frames
#'   in long format, containing columns for tissue, signal, and present.
#'
#' @keywords internal
#' @noRd
.wideToLong <- function(out, removeNA = FALSE, inclusion = "All") {
    
    # --- input validation ---
    stopifnot(
        "'out' must be an Enrichment object" = is(out, "Enrichment"),
        "'removeNA' must be logical"         = is.logical(removeNA),
        "'inclusion' must be a character"    = is.character(inclusion)
    )
    
    validSignals <- c("Tissue enhanced", "Group enriched", "Tissue enriched")
    
    # expand "All" to full set
    if (length(inclusion) == 1 && inclusion == "All") {
        inclusion <- validSignals
    }
    
    # validate inclusion values
    invalidSignals <- inclusion[!inclusion %in% validSignals]
    if (length(invalidSignals) > 0) {
        stop(
            "Invalid inclusion value(s): ",
            paste(invalidSignals, collapse = ", "),
            ". Must be one of: ",
            paste(validSignals, collapse = ", "),
            call. = FALSE
        )
    }
    
    # --- extract data ---
    result   <- out@input
    resultBg <- out@background
    
    tissueCols <- names(out@input)[-c(1:6)]
    
    # --- pivot input to long format ---
    dfLong <- result |>
        tidyr::pivot_longer(
            cols      = tidyr::all_of(tissueCols),
            names_to  = "tissue",
            values_to = "signal"
        ) |>
        dplyr::select(-Tissue.specificity) |>
        dplyr::mutate(present = !is.na(signal))
    
    # --- pivot background to long format ---
    dfLongBg <- resultBg |>
        tidyr::pivot_longer(
            cols      = tidyr::all_of(tissueCols),
            names_to  = "tissue",
            values_to = "signal"
        ) |>
        dplyr::select(-Tissue.specificity) |>
        dplyr::mutate(present = !is.na(signal))
    
    # --- filter based on removeNA ---
    if (!removeNA) {
        
        dfLong   <- dfLong   |> dplyr::filter(signal %in% inclusion | is.na(signal))
        dfLongBg <- dfLongBg |> dplyr::filter(signal %in% inclusion | is.na(signal))
        
    } else {
        
        dfLong   <- dfLong   |> dplyr::filter(signal %in% inclusion & present)
        dfLongBg <- dfLongBg |> dplyr::filter(signal %in% inclusion & present)
        
    }
    
    # --- return Enrichment S4 object ---
    methods::new("Enrichment",
                 input      = data.frame(dfLong),
                 background = data.frame(dfLongBg))
}

