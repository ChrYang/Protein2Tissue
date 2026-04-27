#' Classify genes by tissue specificity
#'
#' This function classifies genes into tissue specificity categories
#' ("Tissue enriched", "Group enriched", "Tissue enhanced",
#' "Low tissue specificity", or "Not detected") based on relative
#' expression levels across tissues, following Human Protein Atlas–style
#' rules.
#'
#' Optionally, the function can integrate predicted secretory protein
#' annotations from the \pkg{Protein2Tissue} database or derive gene,
#' Ensembl, Entrez ID, and UniProt annotations using \pkg{org.Hs.eg.db}.
#'
#' @param expression A numeric matrix or data frame with genes as rows and
#' tissues as columns. Row names of it must correspond to gene identifiers 
#' defined by \code{typeKey}.
#'
#' @param tissue_groups Optional data frame with two columns named
#' \code{Tissue} and \code{Group}, defining predefined tissue groupings for
#' group-enriched analysis. If not provided, the tissue groups would be 
#' automatically calculated among the top-expressing tissues.
#'
#' @param min_expression Numeric value specifying the minimum expression
#' level required for a gene to be considered detected. Default is \code{1}.
#'
#' @param secretory_analysis Logical indicating whether to include predicted
#' secretory protein annotations from the human protein atlas.
#' Default is \code{FALSE}.
#'
#' @param typeKey Character string specifying the gene identifier type used
#' in \code{expression} row names. Must be one of \code{"Gene"},
#' \code{"Ensembl"}, or \code{"EntrezID"}.
#'
#'
#' @return data.frame with tissue specificity classification
#' @export
#'
#' @details
#' The classification is performed gene by gene using the following rules:
#'
#' \itemize{
#'   \item \strong{Not detected}: Maximum expression across all tissues is
#'   below \code{min_expression}.
#'
#'   \item \strong{Tissue enriched}: Expression in one tissue is at least
#'   4-fold higher than expression in all other tissues.
#'
#'   \item \strong{Group enriched}: Mean expression of a group of tissues
#'   (either user-defined or automatically detected among the top-expressing
#'   tissues) is at least 4-fold higher than the maximum expression in all
#'   remaining tissues.
#'
#'   \item \strong{Tissue enhanced}: Expression in one tissue is at least
#'   4-fold higher than the mean expression across all tissues.
#'
#'   \item \strong{Low tissue specificity}: Not meeting any of the criteria 
#'   above.
#' }
#'
#' When \code{secretory_analysis = TRUE}, gene annotations and predicted
#' secretory status are retrieved from the internal \code{ee_tb} dataset
#' provided by the \pkg{Protein2Tissue} package. Genes not present in this
#' database are removed.
#'
#'
#' @return
#' A data frame containing:
#' \itemize{
#'   \item Gene annotation columns: \code{Gene}, \code{Ensembl},
#'   \code{EntrezID}, and \code{Uniprot}
#'   \item \code{Is_secreted} (if available)
#'   \item \code{Tissue.specificity} classification
#'   \item One column per tissue indicating classification labels
#' }
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom utils data
#'
#' @examples
#' \dontrun{
#' # Example expression matrix (genes x tissues)
#' # Classify tissue specificity using gene symbols
#' res <- classify_genes(
#'   expression = expr,
#'   typeKey = "Gene",
#'   min_expression = 1,
#'   secretory_analysis = FALSE
#' )
#'
#' head(res)
#'
#' # Example with predefined tissue groups
#' tissue_groups <- data.frame(
#'   Tissue = c("Colon_Transverse", "Colon_Sigmoid", "Liver"),
#'   Group  = c("Colon", "Colon", "Liver")
#' )
#'
#' res_grouped <- classify_genes(
#'   expression = expr,
#'   tissue_groups = tissue_groups,
#'   typeKey = "Gene",
#'   min_expression = 1,
#'   secretory_analysis = FALSE
#' )
#' }
#'
#' @export


geneClassification <- function(expression,
                           tissue_groups = NULL,
                           min_expression = 1,
                           secretory_analysis = FALSE,
                           typeKey) {
  
  # -----------------------------
  # Input checks
  # -----------------------------
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("expression must have rownames (genes) and colnames (tissues)")
  }
  
  if (!is.numeric(expression)) {
    expression <- as.matrix(expression)
    storage.mode(expression) <- "numeric"
  }
  
  stopifnot(
    "'min_expression' must be numeric" = is.numeric(min_expression),
    "'typeKey' must be a character string" = is.character(typeKey),
    "'typeKey' must be one of Gene, Ensembl, EntrezID" =
      typeKey %in% c("Gene", "Ensembl", "EntrezID")
  )
  
  if (!is.null(tissue_groups)) {
    stopifnot(
      identical(colnames(tissue_groups), c("Tissue", "Group"))
    )
    tissue_groups <- split(tissue_groups$Tissue, tissue_groups$Group)
  }
  
  genes   <- rownames(expression)
  tissues <- colnames(expression)
  
  # -----------------------------
  # Initialise output
  # -----------------------------
  out <- data.frame(
    matrix(NA_character_,
           nrow = nrow(expression),
           ncol = ncol(expression),
           dimnames = list(genes, tissues)),
    stringsAsFactors = FALSE
  )
  
  out[, typeKey] <- genes
  outnames <- colnames(out)[colnames(out) %in% tissues]
  
  # -----------------------------
  # Annotation / secretory logic
  # -----------------------------
  if (secretory_analysis) {
    
    data("ee_tb", package = "Protein2Tissue", envir = environment())
    database <- ee_tb[, c("Gene", "Ensembl", "EntrezID",
                          "Uniprot", "Is_secreted",
                          "Tissue.specificity")]
    
    orig_n <- nrow(out)
    out <- merge(database, out, by = typeKey)
    
    message(
      paste0(orig_n - nrow(out),
             " feature(s) removed due to missing database annotation.")
    )
    
    keep <- !is.na(out$Gene) & !is.na(out$Ensembl) & !is.na(out$EntrezID)
    message(paste0(sum(!keep), " feature(s) removed due to invalid IDs."))
    
    out <- out[keep, , drop = FALSE]
    
    
  } else {
    
    # Empty placeholders
    out$Is_secreted <- NA
    out$Tissue.specificity <- NA
    out$Gene <- out$Ensembl <- out$EntrezID <- out$Uniprot <- NA
    
    map_ids <- function(keys, column, keytype) {
      suppressMessages(
        AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = keys,
          column = column,
          keytype = keytype,
          multiVals = "first"
        )
      )
    }
    
    if (typeKey == "Gene") {
      out$Gene <- genes
      out$Ensembl  <- map_ids(out$Gene, "ENSEMBL", "SYMBOL")
      out$EntrezID <- map_ids(out$Gene, "ENTREZID", "SYMBOL")
      out$Uniprot  <- map_ids(out$Gene, "UNIPROT", "SYMBOL")
    }
    
    if (typeKey == "Ensembl") {
      out$Ensembl <- genes
      out$Gene     <- map_ids(out$Ensembl, "SYMBOL", "ENSEMBL")
      out$EntrezID <- map_ids(out$Ensembl, "ENTREZID", "ENSEMBL")
      out$Uniprot  <- map_ids(out$Ensembl, "UNIPROT", "ENSEMBL")
    }
    
    if (typeKey == "EntrezID") {
      out$EntrezID <- genes
      out$Gene     <- map_ids(out$EntrezID, "SYMBOL", "ENTREZID")
      out$Ensembl  <- map_ids(out$EntrezID, "ENSEMBL", "ENTREZID")
      out$Uniprot  <- map_ids(out$EntrezID, "UNIPROT", "ENTREZID")
    }
    
    keep <- !is.na(out$Gene) & !is.na(out$Ensembl) & !is.na(out$EntrezID)
    message(paste0(sum(!keep), " feature(s) removed due to invalid IDs."))
    
    out <- out[keep, , drop = FALSE]
    
    out <- out[,c("Gene","Ensembl","EntrezID","Uniprot",
                  "Is_secreted","Tissue.specificity",
                  outnames)]
  }
  
  # Reorder expression to match output
  expression <- expression[match(out[, typeKey], rownames(expression)), , drop = FALSE]
  
  # -----------------------------
  # Main classification loop
  # -----------------------------
  for (i in seq_len(nrow(expression))) {
    
    x <- unlist(expression[i, ])
    max_val <- max(x, na.rm = TRUE)
    
    if (max_val < min_expression) {
      out$Tissue.specificity[i] <- "Not detected"
      next
    }
    
    max_idx <- which.max(x)
    max_tissue <- tissues[max_idx]
    others <- x[-max_idx]
    mean_all <- mean(x)
    
    # Tissue enriched
    if (all(max_val >= 4 * others)) {
      out[i, c("Tissue.specificity", max_tissue)] <- "Tissue enriched"
      next
    }
    
    is_group_enriched <- FALSE # logic showing whether the gene is group enriched
    
    # Group enriched (with or without predefined groups)
    if (is.null(tissue_groups)) {
      
      xs <- sort(x, decreasing = TRUE)
      names(xs) <- tissues[order(x, decreasing = TRUE)]
      
      for (k in 2:min(7, length(xs) - 1)) {
        if (mean(xs[1:k]) >= 4 * max(xs[(k + 1):length(xs)])) {
          out[i, c("Tissue.specificity", names(xs[1:k]))] <- "Group enriched"
          is_group_enriched <- TRUE
          break
        }
      }
      
    } else {
      
      for (grp in names(tissue_groups)) {
        tg <- intersect(tissue_groups[[grp]], tissues)
        if (length(tg) < 2) next
        
        if (mean(x[tg]) >= 4 * max(x[setdiff(tissues, tg)])) {
          out[i, c("Tissue.specificity", tg)] <- "Group enriched"
          is_group_enriched <- TRUE
          break
        }
      }
    }
    
    # Tissue enhanced
    
    if (!is_group_enriched) {
      
      enhanced_tissues <- tissues[x >= 4 * mean_all]
      
      if (length(enhanced_tissues) > 0) {
        out[i, c("Tissue.specificity", enhanced_tissues)] <-
          "Tissue enhanced"
      }
    }
    

  }
  
  out$Tissue.specificity[is.na(out$Tissue.specificity)] <- "Low tissue specificity"
  
  rownames(out) <- 1:nrow(out)
  
  out
}
