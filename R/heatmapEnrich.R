#' Heatmap of Tissue Expression Fingerprint
#'
#' Generates a heatmap showing tissue expression signals per gene from an
#' \code{Enrichment} object, with optional filtering by signal category.
#'
#' @param out An \code{Enrichment} S4 object containing the results of
#'   tissue enrichment analysis.
#' @param inclusion A character vector specifying which signal categories to
#'   include. Valid options are \code{"Tissue enhanced"},
#'   \code{"Group enriched"}, and \code{"Tissue enriched"}. Use \code{"All"}
#'   to include all three categories. Default is \code{"All"}.
#'
#' @return A \code{ggplot2} object showing a heatmap of tissue expression
#'   signals per gene, with tissues on the x-axis and genes on the y-axis.
#'
#' @details
#' The function first converts the \code{Enrichment} object to long format
#' using the internal \code{.wideToLong} function, retaining \code{NA} values
#' to show absent signals in the heatmap. Signal categories are colour-coded
#' for easy interpretation.
#'
#' @examples
#' \dontrun{
#' result <- tissueAnalysis(
#'     input      = c("TP53", "BRCA1", "EGFR"),
#'     background = c("347", "948", "1050", "2167", "2819", "3625"),
#'     typeKey    = "Gene",
#'     typeKeyBg  = "Ensembl"
#' )
#'
#' # plot all signal categories
#' heatmapEnrich(result)
#'
#' # plot specific categories only
#' heatmapEnrich(result, inclusion = c("Tissue enriched", "Tissue enhanced"))
#' }
#'
#' @seealso \code{\link{plotProteinType}}, \code{\link{tissueAnalysis}}
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual labs theme
#'   element_text element_blank
#' @export
heatmapEnrich <- function(out, inclusion = "All") {
    
    # --- input validation ---
    stopifnot(
        "'out' must be an Enrichment object"  = is(out, "Enrichment"),
        "'inclusion' must be a character"     = is.character(inclusion)
    )
    
    # --- convert to long format ---
    outLong <- .wideToLong(out, removeNA = FALSE, inclusion = inclusion)
    
    # --- plot ---
    p <- ggplot2::ggplot(
        outLong@input,
        ggplot2::aes(x = tissue, y = Gene, fill = signal)
    ) +
        ggplot2::geom_tile() +
        ggplot2::theme_minimal() +
        ggplot2::scale_fill_manual(
            values = c(
                "Tissue enhanced" = "pink",
                "Group enriched"  = "orange2",
                "Tissue enriched" = "red2",
                "NA"              = "grey90"
            ),
            na.value = "grey90"
        ) +
        ggplot2::labs(
            title = "Tissue expression fingerprint per gene",
            x     = "Tissue",
            y     = "Feature"
        ) +
        ggplot2::theme(
            axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1),
            panel.grid   = ggplot2::element_blank()
        )
    
    p
}