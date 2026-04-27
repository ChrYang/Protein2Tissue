#' Network plot of gene-tissue relationships
#'
#' Creates a network visualisation showing gene-tissue relationships
#' from a wide-format expression object.
#'
#' @param out typically a \code{SummarizedExperiment} or similar Bioconductor 
#' object containing gene expression data across tissues.
#'
#' @return A \code{ggplot} object representing the gene-tissue network.
#'   Tissue nodes are displayed as large labelled circles; gene nodes as
#'   small points. Edges are coloured by signal category.
#'
#' @examples
#' \dontrun{
#'   data(exampleData)
#'   p <- networkPlot(exampleData)
#'   print(p)
#' }
#'
#' @importFrom dplyr select filter
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#'   scale_edge_color_manual
#' @importFrom ggplot2 aes theme_void
#'
#' @export
networkPlot <- function(out) {
    
    ## ── Input validation ────────────────────────────────────────────────────
    if (missing(out) || is.null(out)) {
        stop("'out' must be provided and cannot be NULL.")
    }
    
    ## ── Data preparation ────────────────────────────────────────────────────
    outLong2 <- .wideToLong(out, removeNA = TRUE, inclusion = "All")
    dfLong2  <- outLong2@input
    
    edges <- dplyr::select(
        dfLong2,
        from   = "Gene",
        to     = "tissue",
        signal
    )
    
    edges$to <- gsub("\\.", " ", edges$to)
    edges$to <- gsub("1",   "",  edges$to)
    
    nodes      <- data.frame(name = unique(c(edges$from, edges$to)),
                             stringsAsFactors = FALSE)
    nodes$type <- ifelse(nodes$name %in% dfLong2$Gene, "Gene", "Tissue")
    
    ## ── Graph construction ──────────────────────────────────────────────────
    g <- igraph::graph_from_data_frame(edges,
                                       vertices = nodes,
                                       directed = FALSE)
    
    ## ── Colour palette ──────────────────────────────────────────────────────
    edgeColours <- c(
        "Tissue enhanced" = "#1b9e77",
        "Group enriched"  = "#d95f02",
        "Tissue enriched" = "#7570b3"
    )
    
    ## ── Plot ────────────────────────────────────────────────────────────────
    ggraph::ggraph(g, layout = "fr") +
        
        ggraph::geom_edge_link(
            ggplot2::aes(color = signal),
            linewidth = 1,
            alpha     = 0.8
        ) +
        
        ## Tissue hubs: large bold circles
        ggraph::geom_node_point(
            data   = function(nodes) dplyr::filter(nodes, nodes$type == "Tissue"),
            ggplot2::aes(x = x, y = y),
            shape  = 21,
            size   = 14,
            fill   = "white",
            color  = "black",
            stroke = 2
        ) +
        
        ## Tissue labels inside hubs
        ggraph::geom_node_text(
            data     = function(nodes) dplyr::filter(nodes, nodes$type == "Tissue"),
            ggplot2::aes(label = name),
            size     = 3,
            fontface = "bold"
        ) +
        
        ## Gene nodes: small points
        ggraph::geom_node_point(
            data  = function(nodes) dplyr::filter(nodes, nodes$type == "Gene"),
            ggplot2::aes(x = x, y = y),
            size  = 3,
            color = "black"
        ) +
        
        ## Gene labels with repulsion
        ggraph::geom_node_text(
            data  = function(nodes) dplyr::filter(nodes, nodes$type == "Gene"),
            ggplot2::aes(label = name),
            repel = TRUE,
            size  = 3
        ) +
        
        ggraph::scale_edge_color_manual(values = edgeColours) +
        
        ggplot2::theme_void()
}