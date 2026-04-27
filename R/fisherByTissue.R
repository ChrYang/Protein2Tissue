#' Statistical Test of significance for Tissue Enrichment
#'
#' Performs Fisher's exact test for each tissue to identify tissues that are
#' significantly enriched in the input protein list compared to the background.
#'
#' @param out An \code{Enrichment} S4 object containing the results of
#'   tissue enrichment analysis.
#' @param inclusion A character vector specifying which signal categories to
#'   include in the analysis. Valid options are \code{"Tissue enhanced"},
#'   \code{"Group enriched"}, and \code{"Tissue enriched"}. Use \code{"All"}
#'   to include all three categories. Default is \code{"All"}.
#' @param padj A character string specifying the method for multiple testing
#'   correction. Passed to \code{\link[stats]{p.adjust}}. Common options
#'   include \code{"BH"}, \code{"bonferroni"}, and \code{"none"}.
#'   Default is \code{"none"}.
#' @param secretoryOnly A logical value indicating whether to restrict the
#'   analysis to secretory proteins only (i.e. rows where \code{Is_secreted}
#'   is \code{TRUE}). Default is \code{FALSE}.
#'
#' @return A data frame with one row per tissue containing the following
#'   columns:
#'   \itemize{
#'     \item \code{tissue} — tissue name
#'     \item \code{input_in} — number of input proteins present in the tissue
#'     \item \code{input_total} — total number of input proteins
#'     \item \code{bg_in} — number of background proteins present in the tissue
#'     \item \code{bg_total} — total number of background proteins
#'     \item \code{input_prop} — proportion of input proteins in the tissue
#'     \item \code{bg_prop} — proportion of background proteins in the tissue
#'     \item \code{OR} — odds ratio from Fisher's exact test
#'     \item \code{p_value} — p-value from Fisher's exact test
#'     \item \code{p_adj} — adjusted p-value
#'   }
#'
#' @details
#' For each tissue, a 2x2 contingency table is constructed and Fisher's exact
#' test is applied to test whether the proportion of input proteins expressed
#' in that tissue differs significantly from the background proportion.
#'
#' The contingency table is structured as follows:
#' \preformatted{
#'                  Expressed    Not expressed
#'   Input              a               b
#'   Background         c               d
#' }
#'
#' Multiple testing correction is applied across all tissues using the method
#' specified in \code{padj}.
#'
#' @examples
#' \dontrun{
#' result <- tissueAnalysis(
#'     input      = c("TP53", "BRCA1", "EGFR"),
#'     background = c("347", "948", "1050", "2167", "2819", "3625")
#' )
#'
#' # run Fisher's exact test across all tissues
#' fisherResult <- fisherByTissue(result)
#'
#' # with BH correction and secretory proteins only
#' fisherResult <- fisherByTissue(
#'     result,
#'     inclusion     = c("Tissue enriched", "Tissue enhanced"),
#'     padj          = "BH",
#'     secretoryOnly = TRUE
#' )
#' }
#'
#' @seealso \code{\link{tissueAnalysis}}, \code{\link{heatmapEnrich}},
#'   \code{\link[stats]{fisher.test}}, \code{\link[stats]{p.adjust}}
#'
#' @importFrom stats fisher.test p.adjust
#' @importFrom dplyr filter
#' @export
fisherByTissue <- function(out,
                           inclusion     = "All",
                           padj          = "none",
                           secretoryOnly = FALSE) {
    
    # --- input validation ---
    stopifnot(
        "'out' must be an Enrichment object"    = is(out, "Enrichment"),
        "'inclusion' must be a character"       = is.character(inclusion),
        "'padj' must be a character string"     = is.character(padj),
        "'secretoryOnly' must be a boolean"       = is.logical(secretoryOnly)
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
    
    # validate padj method
    validPadj <- c("holm", "hochberg", "hommel", "bonferroni",
                   "BH", "BY", "fdr", "none")
    if (!padj %in% validPadj) {
        stop(
            "'padj' must be one of: ",
            paste(validPadj, collapse = ", "),
            call. = FALSE
        )
    }
    
    # --- filter by secretory if needed ---
    if (!secretoryOnly) {
        fg <- out@input
        bg <- out@background
    } else {
        if (!"Is_secreted" %in% colnames(out@input)) {
            stop(
                "'Is_secreted' column not found in input. ",
                "Cannot filter by secretory proteins.",
                call. = FALSE
            )
        }
        fg <- out@input      |> dplyr::filter(Is_secreted)
        bg <- out@background |> dplyr::filter(Is_secreted)
    }
    
    nFg <- nrow(fg)
    nBg <- nrow(bg)
    
    # --- get tissue list ---
    outLong    <- .wideToLong(out, removeNA = TRUE, inclusion = inclusion)
    tissueList <- sort(unique(outLong@input$tissue))
    
    # --- run Fisher's exact test per tissue ---
    resultList <- lapply(tissueList, function(tissueName) {
        
        a <- sum(fg[, tissueName] %in% inclusion)
        b <- nFg - a
        c <- sum(bg[, tissueName] %in% inclusion)
        d <- nBg - c
        
        mat  <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
        test <- stats::fisher.test(mat)
        
        data.frame(
            tissue      = tissueName,
            input_in    = a,
            input_total = nFg,
            bg_in       = c,
            bg_total    = nBg,
            input_prop  = a / nFg,
            bg_prop     = c / nBg,
            OR          = unname(test$estimate),
            p_value     = test$p.value,
            stringsAsFactors = FALSE
        )
    })
    
    # --- combine results ---
    output         <- do.call(rbind, resultList)
    output$p_adj   <- stats::p.adjust(output$p_value, method = padj)
    rownames(output) <- NULL
    
    output
}