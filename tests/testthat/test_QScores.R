# tests/testthat/test-computeQScoreMods.R
# library(testthat)
# library(SpaceTrooper)
#
# test_that("computeQScoreMods identical to computeQScore on CosMx", {
#     sample_dir <- system.file(file.path("extdata", "CosMx_DBKero_Tiny"),
#                                package="SpaceTrooper")
#     skip_if(sample_dir == "" || !dir.exists(sample_dir),
#             "no CosMx example data found in inst/extdata; salto il test")
#
#     spe <- readCosmxSPE(sample_dir)
#
#     spe <- spatialPerCellQC(spe)
#     library(withr)
#     with_seed(1998,
#         qs_orig <- computeQScore(spe)
#     )
#     with_seed(1998,
#         qs_mod  <- computeQScoreMods(spe)
#     )
#
#     expect_equal(qs_mod$quality_score, qs_orig$quality_score, tolerance = 1e-8)
# })
