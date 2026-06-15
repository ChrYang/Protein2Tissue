## ── ee_tb ─────────────────────────────────────────────────────────────────────

test_that("ee_tb loads correctly", {
    data("ee_tb", package = "Protein2Tissue")
    expect_true(exists("ee_tb"))
    expect_true(is.data.frame(ee_tb))
    expect_gt(nrow(ee_tb), 0)
})

test_that("ee_tb has required columns", {
    data("ee_tb", package = "Protein2Tissue")
    requiredCols <- c("Gene", "Ensembl", "EntrezID", "Tissue.specificity")
    expect_true(all(requiredCols %in% colnames(ee_tb)))
})

test_that("ee_tb signal column has valid categories only", {
    data("ee_tb", package = "Protein2Tissue")
    validSignals <- c("Tissue enhanced", "Group enriched",
                      "Tissue enriched", "Low tissue specificity",
                      "Not detected", NA)
    expect_true(all(ee_tb$signal %in% validSignals))
})

## ── plasma_data ───────────────────────────────────────────────────────────────

test_that("plasma_data loads correctly", {
    data("plasma_data", package = "Protein2Tissue")
    expect_true(exists("plasma_data"))
    expect_true(is.data.frame(plasma_data))
    expect_gt(nrow(plasma_data), 0)
})

test_that("plasma_data has required limma columns", {
    data("plasma_data", package = "Protein2Tissue")
    requiredCols <- c("UniqueSymbol", "Coefficient.Age", "p.Age", "q.Age")
    expect_true(all(requiredCols %in% colnames(plasma_data)))
})

test_that("plasma_data q values are all between 0 and 1", {
    data("plasma_data", package = "Protein2Tissue")
    expect_true(all(plasma_data$q.Age >= 0 &
                        plasma_data$q.Age <= 1, na.rm = TRUE))
})