#' Specific tissue Enrichment Analysis 
#'
#'Performs enrichment analysis in the tissue designated by the 
#'parameter \item{tissueToTest} by mapping input genes or proteins against 
#'a reference expression database. Genes or proteins not found in the database 
#'are removed with a warning, and the matched results are returned as a data 
#'frame showing genes or proteins enriched in the corresponding tissue.
#'
#' @param input A character vector of gene or proteins identifiers to be 
#'   analysed (e.g. gene symbols or Ensembl IDs).
#' @param typeKey A character string specifying the identifier type for the
#' input genes or proteins. Must match a column name in the reference database.
#' Default is \code{"Gene"} (gene symbol).
#' @param database Defaults to \code{ee_tb} (Human Protein Atlas data bundled 
#' with this package). Can be replaced with a custom data frame provided it 
#' contains columns matching \code{typeKey}.
#' @param inclusion A character vector specifying which signal categories to
#' include in the analysis. Valid options are \code{"Tissue enhanced"},
#' \code{"Group enriched"}, and \code{"Tissue enriched"}. Use \code{"All"} 
#' to include all three categories. Default is \code{"All"}.
#' @param tissueToTest The tissue in which the input proteins or genes are 
#' tested for enrichment. Must be one of the tissue listed in the \code{database} 
#' provided. Tissues in the default database are c("Adipose tissue","Breast",
#' "Seminal vesicle","Parathyroid gland","Bone marrow","Lung","Lymphoid tissue",
#' "Liver","Testis","Brain","Esophagus","Gallbladder","Intestine","Pancreas",
#' "Ovary","Kidney","Placenta", "Cervix","Skin","Endometrium","Skeletal muscle",
#' "Retina","Tongue","Heart muscle", "Choroid plexus","Fallopian tube",
#' "Epididymis","Salivary gland","Vagina","Smooth muscle","Adrenal gland",
#' "Blood vessel", "Pituitary gland","Stomach","Urinary bladder",
#' "Thyroid gland","Prostate") 
#' @param secretoryOnly A logical value indicating whether to restrict the
#' analysis to secretory proteins only (i.e. rows where \code{Is_secreted} 
#' is \code{TRUE}). Default is \code{FALSE}.
#'
#' @return A data frame of matched input genes with reference database 
#' annotations for the tissue designated by the \code{tissueToTest} 
#' parameter.
#'
#' @details
#' The function maps input gene or proteins lists against the internal
#' reference expression database like the default\code{ee_tb} in the 
#' tissue designated by the parameter \code{tissueToTest}. Proteins 
#' or genes not found in the database are removed before analysis. 
#' If no proteins are found in the database, the function stops with 
#' an informative error message.
#' 
#'
#' @examples
#' \dontrun{
#' # using gene symbols for input, Ensembl IDs for background
#' result2 <- tissueAnalysis(
#'   input          = c("BRCA1", "EGFR"),
#'   typeKey        = "Gene",
#'   database       = NULL,
#'   inclusion      = "All",
#'   secretoryOnly  = FALSE,
#'   tissueToTest   = "Testis"
#' )
#'
#' # access results
#' result2
#' }
#'
#' @seealso \code{\link{enrichmentAnalysis}}, \code{\link{Enrichment-class}}
#'
#' @importFrom dplyr left_join
#' @export
tissueAnalysis2 <- function(input,
                            typeKey        = "Gene",
                            database       = NULL,
                            inclusion      = "All",
                            secretoryOnly  = FALSE,
                            tissueToTest) {
  
  if (is.null(database)) {
    database <- get(utils::data("ee_tb",
                                package  = "Protein2Tissue",
                                envir    = environment()))
  }
  # --- input validation ---
  stopifnot(
    "'input' must be a non-empty vector"      = length(input) > 0,
    "'typeKey' must be a character string"    = is.character(typeKey),
    "'tissueToTest' must be a character string"  = is.character(tissueToTest),
    "'tissueToTest' must be a boolean"  = is.logical(secretoryOnly)
  )
  
  stopifnot(
    "'typeKey' not found in database"   = typeKey   %in% colnames(database)
  )
  
  database_sel <- database[c("Gene","Ensembl","EntrezID","Uniprot",
                             "Is_secreted","Tissue.specificity",
                             tissueToTest)] 
  
  # --- process input ---
  inputDf <- .processProteinList(
    proteinList = input,
    typeKey     = typeKey,
    listLabel   = "input",
    database = database
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
  
  
  resultInput <- dplyr::left_join(inputDf, database_sel, by = typeKey)
  
  col <- resultInput[[tissueToTest]]
  resultInput <- resultInput[col %in% inclusion, ]
  
  if (secretoryOnly) {
    resultInput <- resultInput[resultInput$Is_secreted,]
  }
  
  # --- return results ---
  return(resultInput[,colnames(resultInput) != "RNA.tissue.specificity"])
}


