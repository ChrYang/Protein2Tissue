#' Plot Fisher's Exact Test Results by Tissue
#'
#' Generates a horizontal grouped bar chart comparing the proportion of input
#' proteins versus background proteins across tissues, with p-values annotated
#' from Fisher's exact test.
#'
#' @param fisherTable A data frame output from \code{\link{fisherByTissue}},
#'   containing columns \code{tissue}, \code{input_prop}, \code{bg_prop},
#'   and \code{p_value}.
#'
#' @return A \code{ggplot2} object showing a horizontal grouped bar chart with:
#'   \itemize{
#'     \item Blue bars representing the proportion of input proteins per tissue
#'     \item Grey bars representing the proportion of background proteins
#'     \item P-values from Fisher's exact test annotated per tissue
#'   }
#'
#' @details
#' The function reshapes the Fisher's exact test result table from wide to long
#' format and plots input and background proportions side by side per tissue.
#' Tissues are ordered by proportion for easier visual comparison. P-values
#' are rounded to 3 decimal places and displayed to the right of the bars.
#'
#' @examples
#' \dontrun{
#' result <- tissueAnalysis(
#'     input      = c("TP53", "BRCA1", "EGFR"),
#'     background = c("347", "948", "1050", "2167", "2819", "3625")
#' )
#'
#' fisherResult <- fisherByTissue(result, padj = "BH")
#'
#' # plot results
#' fisherPlot(fisherResult)
#' }
#'
#' @seealso \code{\link{fisherByTissue}}, \code{\link{tissueAnalysis}}
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_text position_dodge coord_flip
#'   labs scale_fill_manual theme_classic theme element_text
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr select mutate
#' @export
fisherPlot <- function(fisherTable) {
    
    # --- input validation ---
    stopifnot(
        "'fisherTable' must be a data frame" = is.data.frame(fisherTable)
    )
    
    requiredCols <- c("tissue", "input_prop", "bg_prop", "p_value")
    missingCols  <- requiredCols[!requiredCols %in% colnames(fisherTable)]
    if (length(missingCols) > 0) {
        stop(
            "Missing required columns in 'fisherTable': ",
            paste(missingCols, collapse = ", "),
            call. = FALSE
        )
    }
    
    # --- reshape to long format ---
    tableLong <- fisherTable |>
        dplyr::select(tissue, input_prop, bg_prop) |>
        tidyr::pivot_longer(
            cols      = c(input_prop, bg_prop),
            names_to  = "dataset",
            values_to = "proportion"
        ) |>
        dplyr::mutate(proportion = round(as.numeric(proportion), 3))
    
    # --- p-value annotation data frame ---
    pvalDf <- fisherTable |>
        dplyr::select(tissue, p_value) |>
        dplyr::mutate(p_value = round(as.numeric(p_value), 3))
    
    # --- plot ---
    p <- ggplot2::ggplot(
        tableLong,
        ggplot2::aes(
            x    = reorder(tissue, proportion),
            y    = proportion,
            fill = dataset
        )
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge(width = 0.8),
            width    = 0.7
        ) +
        ggplot2::geom_text(
            ggplot2::aes(label = proportion),
            position = ggplot2::position_dodge(width = 0.8),
            hjust    = -0.05,
            size     = 3
        ) +
        ggplot2::geom_text(
            data         = pvalDf,
            ggplot2::aes(
                x     = tissue,
                y     = max(tableLong$proportion, na.rm = TRUE) * 1.15,
                label = paste0("p=", p_value)
            ),
            inherit.aes = FALSE,
            size        = 3
        ) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x        = "Tissue",
            y        = "Proportion of proteins",
            fill     = "",
            title    = "Protein enrichment by tissue\nInput vs Background",
            subtitle = "p-values from Fisher's exact test"
        ) +
        ggplot2::scale_fill_manual(
            values = c(
                input_prop = "steelblue",
                bg_prop    = "grey70"
            ),
            labels = c(
                input_prop = "Input",
                bg_prop    = "Background"
            )
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            legend.position = "top"
        )
    
    p
}