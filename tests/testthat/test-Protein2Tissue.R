## ── Shared test fixtures ──────────────────────────────────────────────────────

## small gene list — all known to exist in ee_tb
testGenes <- c("TP53", "BRCA1", "EGFR", "MYC", "PTEN",
               "ALB",  "APOE",  "FN1",  "TTR", "HP")

## unknown genes — should be silently dropped
unknownGenes <- c("FAKEGENE1", "FAKEGENE2")

## run tissueAnalysis once and reuse across tests
result <- tissueAnalysis(
    input      = testGenes,
    background = background_default,
    typeKey    = "Gene",
    typeKeyBg  = "Ensembl"
)

## ── tissueAnalysis ────────────────────────────────────────────────────────────

test_that("tissueAnalysis returns an Enrichment S4 object", {
    expect_s4_class(result, "Enrichment")
})

test_that("tissueAnalysis output slots are data frames", {
    expect_true(is.data.frame(result@input))
    expect_true(is.data.frame(result@background))
})

test_that("tissueAnalysis input slot contains expected genes", {
    expect_true(all(testGenes %in% result@input$Gene))
})

test_that("tissueAnalysis background slot is non-empty", {
    expect_gt(nrow(result@background), 0)
})

test_that("tissueAnalysis errors on empty input", {
    expect_error(
        tissueAnalysis(input = character(0), background = background_default),
        "'input' must be a non-empty vector"
    )
})

test_that("tissueAnalysis errors on empty background", {
    expect_error(
        tissueAnalysis(input = testGenes, background = character(0)),
        "'background' must be a non-empty vector"
    )
})

test_that("tissueAnalysis errors on invalid typeKey", {
    expect_error(
        tissueAnalysis(input = testGenes, typeKey = "NotAColumn"),
        "'typeKey' not found in database"
    )
})

test_that("tissueAnalysis errors on invalid typeKeyBg", {
    expect_error(
        tissueAnalysis(input = testGenes, typeKeyBg = "NotAColumn"),
        "'typeKeyBg' not found in database"
    )
})

test_that("tissueAnalysis accepts a custom database", {
    customDb <- ee_tb[seq_len(100), ]
    expect_s4_class(
        tissueAnalysis(
            input      = testGenes,
            background = background_default,
            database   = customDb
        ),
        "Enrichment"
    )
})

## ── fisherByTissue ────────────────────────────────────────────────────────────

fisherResult <- fisherByTissue(result, padj = "BH")

test_that("fisherByTissue returns a data frame", {
    expect_true(is.data.frame(fisherResult))
})

test_that("fisherByTissue output has expected columns", {
    expectedCols <- c("tissue", "input_in", "input_total",
                      "bg_in", "bg_total", "input_prop",
                      "bg_prop", "OR", "p_value", "p_adj")
    expect_true(all(expectedCols %in% colnames(fisherResult)))
})

test_that("fisherByTissue p_value is between 0 and 1", {
    expect_true(all(fisherResult$p_value >= 0 & fisherResult$p_value <= 1))
})

test_that("fisherByTissue p_adj is between 0 and 1", {
    expect_true(all(fisherResult$p_adj >= 0 & fisherResult$p_adj <= 1))
})

test_that("fisherByTissue OR is non-negative", {
    expect_true(all(fisherResult$OR >= 0, na.rm = TRUE))
})

test_that("fisherByTissue errors on non-Enrichment input", {
    expect_error(
        fisherByTissue(data.frame(x = 1)),
        "'out' must be an Enrichment object"
    )
})

test_that("fisherByTissue errors on invalid inclusion", {
    expect_error(
        fisherByTissue(result, inclusion = "InvalidCategory"),
        "Invalid inclusion value"
    )
})

test_that("fisherByTissue errors on invalid padj method", {
    expect_error(
        fisherByTissue(result, padj = "notamethod"),
        "'padj' must be one of"
    )
})

test_that("fisherByTissue secretoryOnly works without error", {
    expect_true(is.data.frame(
        fisherByTissue(result, secretoryOnly = TRUE)
    ))
})

test_that("fisherByTissue inclusion subset works", {
    res <- fisherByTissue(result, inclusion = c("Tissue enriched"))
    expect_true(is.data.frame(res))
    expect_gt(nrow(res), 0)
})

## ── fisherPlot ────────────────────────────────────────────────────────────────

test_that("fisherPlot returns a ggplot object", {
    p <- fisherPlot(fisherResult)
    expect_s3_class(p, "ggplot")
})

test_that("fisherPlot errors on non-data-frame input", {
    expect_error(
        fisherPlot("not a data frame"),
        "'fisherTable' must be a data frame"
    )
})

test_that("fisherPlot errors on missing required columns", {
    badDf <- data.frame(tissue = "liver", x = 1)
    expect_error(
        fisherPlot(badDf),
        "Missing required columns"
    )
})

## ── plotProteinType ───────────────────────────────────────────────────────────

test_that("plotProteinType returns a ggplot object", {
    p <- plotProteinType(result)
    expect_s3_class(p, "ggplot")
})

test_that("plotProteinType errors on non-Enrichment input", {
    expect_error(
        plotProteinType(data.frame(x = 1)),
        "'out' must be an Enrichment object"
    )
})

## ── heatmapEnrich ─────────────────────────────────────────────────────────────

test_that("heatmapEnrich returns a ggplot object", {
    p <- heatmapEnrich(result)
    expect_s3_class(p, "ggplot")
})

test_that("heatmapEnrich works with specific inclusion", {
    p <- heatmapEnrich(result, inclusion = c("Tissue enriched", "Tissue enhanced"))
    expect_s3_class(p, "ggplot")
})

test_that("heatmapEnrich errors on non-Enrichment input", {
    expect_error(
        heatmapEnrich(data.frame(x = 1)),
        "'out' must be an Enrichment object"
    )
})

test_that("heatmapEnrich errors on invalid inclusion", {
    expect_error(
        heatmapEnrich(result, inclusion = "InvalidCategory"),
        "inclusion"
    )
})

## ── networkPlot ───────────────────────────────────────────────────────────────

test_that("networkPlot returns a ggplot object", {
    p <- networkPlot(result)
    expect_s3_class(p, "ggplot")
})

test_that("networkPlot errors on NULL input", {
    expect_error(
        networkPlot(NULL),
        "'out' must be provided and cannot be NULL"
    )
})

test_that("networkPlot errors on missing input", {
    expect_error(
        networkPlot(),
        "'out' must be provided and cannot be NULL"
    )
})

## ── goEnriched ────────────────────────────────────────────────────────────────

test_that("goEnriched errors on NULL input", {
    expect_error(
        goEnriched(NULL),
        "'out' must be provided and cannot be NULL"
    )
})

test_that("goEnriched errors on invalid type", {
    expect_error(
        goEnriched(result, type = "XX"),
        "'type' must be one of"
    )
})

test_that("goEnriched errors on invalid tissueTest type", {
    expect_error(
        goEnriched(result, tissueTest = c("liver", "kidney")),
        "'tissueTest' must be a single character string"
    )
})

test_that("goEnriched errors on invalid inclusion", {
    expect_error(
        goEnriched(result, inclusion = "InvalidCategory"),
        "'inclusion' must be 'All' or a subset of"
    )
})

test_that("goEnriched errors when no genes found for tissue", {
    expect_error(
        goEnriched(result, tissueTest = "nonexistenttissue"),
        "No genes found for tissue"
    )
})

test_that("goEnriched returns enrichResult for valid tissue", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("clusterProfiler")
    
    res <- goEnriched(result, type = "BP", tissueTest = "liver")
    expect_s4_class(res, "enrichResult")
})