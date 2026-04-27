#' GO enrichment analysis on tissue-specific genes
#'
#' Performs Gene Ontology (GO) enrichment analysis on genes showing
#' tissue-specific or group-enriched expression patterns, using
#' \code{\link[clusterProfiler]{enrichGO}}.
#'
#' @param out An object compatible with \code{.wideToLong()}, typically
#'   a \code{SummarizedExperiment} or similar Bioconductor object containing
#'   gene expression data across tissues.
#' @param type A character string specifying the GO ontology category.
#'   One of \code{"BP"} (Biological Process), \code{"MF"} (Molecular
#'   Function), or \code{"CC"} (Cellular Component). Default is \code{"BP"}.
#' @param tissueTest A character string specifying the tissue to test.
#'   Must match a value in the \code{tissue} column of the long-format data.
#'   Default is \code{""}.
#' @param inclusion A character vector specifying which signal categories
#'   to include. One of \code{"Tissue enhanced"}, \code{"Group enriched"},
#'   \code{"Tissue enriched"}, or \code{"All"} (default) which includes
#'   all three categories.
#'
#' @return An \code{enrichResult} object from
#'   \code{\link[clusterProfiler]{enrichGO}} containing GO enrichment
#'   results for the specified tissue and signal categories.
#'
#' @examples
#' \dontrun{
#'   data(exampleData)
#'   result <- goEnriched(exampleData, type = "BP", tissueTest = "liver")
#'   print(result)
#' }
#'
#' @importFrom dplyr filter
#' @importFrom clusterProfiler enrichGO
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
#' @export
goEnriched <- function(out,
                       type       = "BP",
                       tissueTest = "",
                       inclusion  = "All") {
    
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop("Package 'org.Hs.eg.db' is required for GO enrichment. ",
             "Install it with: BiocManager::install('org.Hs.eg.db')",
             call. = FALSE)
    }
    
    ## ── Input validation ────────────────────────────────────────────────────
    if (missing(out) || is.null(out)) {
        stop("'out' must be provided and cannot be NULL.")
    }
    
    validTypes <- c("BP", "MF", "CC")
    if (!type %in% validTypes) {
        stop("'type' must be one of: ", paste(validTypes, collapse = ", "), ".")
    }
    
    if (!is.character(tissueTest) || length(tissueTest) != 1) {
        stop("'tissueTest' must be a single character string.")
    }
    
    validSignals <- c("Tissue enhanced", "Group enriched", "Tissue enriched")
    
    if (length(inclusion) == 1 && inclusion == "All") {
        inclusion <- validSignals
    } else if (!all(inclusion %in% validSignals)) {
        stop("'inclusion' must be 'All' or a subset of: ",
             paste(validSignals, collapse = ", "), ".")
    }
    
    ## ── Data preparation ────────────────────────────────────────────────────
    outLong2 <- .wideToLong(out, removeNA = TRUE, inclusion = inclusion)
    
    inputLong <- dplyr::filter(
        outLong2@input,
        signal %in% inclusion,
        tissue == tissueTest
    )
    
    backgroundLong <- dplyr::filter(
        outLong2@background,
        signal %in% inclusion,
        tissue == tissueTest
    )
    
    if (nrow(inputLong) == 0) {
        stop("No genes found for tissue '", tissueTest,
             "' with the specified inclusion categories.")
    }
    
    inputExt     <- as.character(unique(inputLong$EntrezID))
    backgroundExt <- as.character(unique(backgroundLong$EntrezID))
    
    ## ── GO enrichment ───────────────────────────────────────────────────────
    egoBP <- clusterProfiler::enrichGO(
        gene          = inputExt,
        universe      = backgroundExt,
        OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = type,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
    )
    
    return(egoBP)
}