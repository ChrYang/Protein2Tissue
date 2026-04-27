#' Plot Gene-level RNA Tissue Specificity
#'
#' Generates a bar chart showing the distribution of genes across RNA tissue
#' specificity classes from an \code{Enrichment} object.
#'
#' @param out An \code{Enrichment} S4 object containing the results of
#'   tissue enrichment analysis.
#'
#' @return A \code{ggplot2} object showing the number of genes in each
#'   RNA tissue specificity class.
#'
#' @details
#' The function summarises the \code{Tissue.specificity} column from the
#' input slot of the \code{Enrichment} object and plots the count of genes
#' in each specificity class as a horizontal bar chart.
#'
#' @examples
#' \dontrun{
#' result <- tissueAnalysis(
#'     input      = c("TP53", "BRCA1", "EGFR"),
#'     background = c("347", "948", "1050", "2167", "2819", "3625"),
#'     typeKey    = "Gene",
#'     typeKeyBg  = "Ensembl"
#' )
#' plotProteinType(result)
#' }
#'
#' @seealso \code{\link{heatmapEnrich}}, \code{\link{tissueAnalysis}}
#'
#' @importFrom ggplot2 ggplot aes geom_col labs coord_flip
#' @importFrom dplyr group_by summarise
#' @export
plotProteinType <- function(out) {
    
    # --- input validation ---
    stopifnot(
        "'out' must be an Enrichment object" = is(out, "Enrichment")
    )
    
    # --- plot ---
    p <- out@input |>
        dplyr::group_by(Tissue.specificity) |>
        dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
        ggplot2::ggplot(ggplot2::aes(
            x = Tissue.specificity,
            y = n
        )) +
        ggplot2::geom_col(fill = "#4C72B0") +
        ggpubr::theme_classic2() +
        ggplot2::labs(
            title = "Gene-level RNA tissue specificity",
            x     = "Specificity class",
            y     = "Number of genes"
        ) +
        ggplot2::coord_flip()
    
    p
}