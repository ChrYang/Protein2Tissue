#' Human Protein Atlas tissue expression database
#'
#' A data frame containing tissue-specific gene expression data
#' derived from the Human Protein Atlas (HPA).
#'
#' @format A data frame with \code{n} rows and the following columns:
#' \describe{
#'   \item{Gene}{HGNC gene symbol}
#'   \item{Ensembl}{Ensembl gene identifier}
#'   \item{EntrezID}{NCBI Entrez ID}
#'   \item{Uniprot}{UNIPROT ID}
#'   \item{Is_secreted}{Whether the protein is a predicted secretory protein  or not}
#'   \item{Tissue.specificity}{expression category, one of \code{"Tissue enhanced"},
#'     \code{"Group enriched"}, or \code{"Tissue enriched"}}
#'     
#'   \item{tissue}{tissue names}
#' }
#'
#' @source Human Protein Atlas \url{https://www.proteinatlas.org}.
#'   Uhlen M et al. (2017). A pathology atlas of the human cancer
#'   transcriptome. \emph{Science}, 357(6352).
#'   \doi{10.1126/science.aan2507}
#'
"ee_tb"



#' GTEx tissue-average TPM matrix
#'
#' A numeric matrix containing tissue-averaged transcript-per-million (TPM)
#' expression values from the GTEx project. Rows correspond to genes
#' (Ensembl gene IDs) and columns correspond to human tissues.
#'
#' @format
#' A numeric matrix with genes (Ensembl gene IDs) as rows and tissues as columns.
#'
#' @source
#' The Genotype-Tissue Expression (GTEx) Project  
#' \url{https://www.gtexportal.org/}  
#' GTEx Consortium (2013). *The Genotype-Tissue Expression (GTEx) project*.
#' *Nature Genetics*, 45(6). \doi{10.1038/ng.2653}
"tissue_tpm_gtex"


#' An example of the data frame the for the 'tissue_groups' parameter of 
#' \code{geneClassification}. 
#'
#' @format An example of the data frame the for the 'tissue_groups' parameter of 
#' \code{geneClassification}. This data frame has two columns named
#' \code{Tissue} and \code{Group}, defining predefined tissue groupings for
#' group-enriched analysis. 
#' 
"tissue_map"


#' Example plasma proteomics dataset of age-associated protein expression
#'
#' A data frame containing differential expression results from a plasma
#' proteomics study. The dataset captures associations between protein
#' expression levels and age and is intended as example input for
#' \code{tissueAnalysis()}.
#'
#' @format
#' A data frame with \code{n} rows and 4 columns:
#' \describe{
#'   \item{UniqueSymbol}{Character. Protein or gene identifier.}
#'   \item{Coefficient.Age}{Numeric. Estimated effect size for association with age.}
#'   \item{p.Age}{Numeric. Raw p-value for the age association.}
#'   \item{q.Age}{Numeric. Adjusted p-value (q-value) controlling for multiple testing.}
#' }
#'
#' @details
#' This dataset is adapted from Lehallier et al. (2019), which profiled
#' age-associated changes in the human plasma proteome.
#'
#' @source
#' Lehallier et al. (2019). *Undulating changes in human plasma proteome
#' profiles across the lifespan*. Nature Medicine.
#' \doi{10.1038/s41591-019-0673-2}
#'
#' @usage
#' plasma_data
#'
#' @keywords datasets
"plasma_data"
