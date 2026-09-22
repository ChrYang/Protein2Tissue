#' @include AllGenerics.R
NULL

#' @rdname getInput
#' @export
setMethod("getInput", "Enrichment", function(object) object@input)

#' @rdname getBackground
#' @export
setMethod("getBackground", "Enrichment", function(object) object@background)

#' @exportMethod show
setMethod("show", "Enrichment", function(object) {
    cat("An object of class 'Enrichment'\n")
    cat("--------------------------------\n")
    cat("Input genes:     ", nrow(object@input),      "\n")
    cat("Background genes:", nrow(object@background), "\n")
    cat("\nInput (first 5 rows):\n")
    print(head(object@input, 5))
})