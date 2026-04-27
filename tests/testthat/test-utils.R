## ── integration: full workflow ────────────────────────────────────────────────

test_that("full workflow runs without error", {
    result <- tissueAnalysis(
        input      = c("TP53", "BRCA1", "EGFR"),
        background = background_default,
        typeKey    = "Gene",
        typeKeyBg  = "Ensembl"
    )
    
    fisherResult <- fisherByTissue(result, padj = "BH")
    
    expect_s4_class(result, "Enrichment")
    expect_true(is.data.frame(fisherResult))
    expect_s3_class(fisherPlot(fisherResult), "ggplot")
    expect_s3_class(heatmapEnrich(result), "ggplot")
    expect_s3_class(plotProteinType(result), "ggplot")
    expect_s3_class(networkPlot(result), "ggplot")
})

## ── edge cases ────────────────────────────────────────────────────────────────

test_that("tissueAnalysis handles unknown genes gracefully", {
    expect_warning(
        tissueAnalysis(
            input      = c("TP53", "FAKEGENE999"),
            background = background_default
        )
    )
})

test_that("fisherByTissue handles single inclusion category", {
    result <- tissueAnalysis(
        input      = c("TP53", "BRCA1", "EGFR"),
        background = background_default
    )
    res <- fisherByTissue(result, inclusion = "Tissue enriched")
    expect_true(is.data.frame(res))
})

test_that("plasma_data genes work as tissueAnalysis input", {
    data("plasma_data", package = "Protein2Tissue")
    
    myGenes <- plasma_data |>
        dplyr::filter(Coefficient.Age > 0 & q.Age < 0.05) |>
        dplyr::pull(X) |>
        toupper()
    
    result <- tissueAnalysis(
        input      = myGenes,
        background = background_default,
        typeKey    = "Gene",
        typeKeyBg  = "Ensembl"
    )
    expect_s4_class(result, "Enrichment")
})