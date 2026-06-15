## ── Enrichment S4 class ───────────────────────────────────────────────────────

test_that("Enrichment class can be instantiated", {
    obj <- methods::new("Enrichment",
                        input      = data.frame(Gene = "TP53"),
                        background = data.frame(Gene = "BRCA1"))
    expect_s4_class(obj, "Enrichment")
})

test_that("Enrichment input slot is a data frame", {
    obj <- methods::new("Enrichment",
                        input      = data.frame(Gene = "TP53"),
                        background = data.frame(Gene = "BRCA1"))
    expect_true(is.data.frame(obj@input))
})

test_that("Enrichment background slot is a data frame", {
    obj <- methods::new("Enrichment",
                        input      = data.frame(Gene = "TP53"),
                        background = data.frame(Gene = "BRCA1"))
    expect_true(is.data.frame(obj@background))
})

test_that("isVirtualClass is FALSE for Enrichment", {
    expect_false(methods::isVirtualClass("Enrichment"))
})