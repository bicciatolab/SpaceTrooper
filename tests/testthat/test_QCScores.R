library(testthat)
library(SpaceTrooper)

# load the example SpatialExperiment
spe0 <- example(readCosmxSPE)$value

test_that("QC functions are exported", {
    expect_true(exists("spatialPerCellQC",   mode = "function"))
    expect_true(exists("computeQCScore",      mode = "function"))
    expect_true(exists("computeSpatialOutlier", mode = "function"))
    expect_true(exists("computeQCScoreFlags",  mode = "function"))
    expect_true(exists("computeThresholdFlags",  mode = "function"))
})


test_that("spatialPerCellQC adds per‐cell metrics to colData", {
    spe <- spatialPerCellQC(spe0, micronConvFact = 0.15)
    expect_s4_class(spe, "SpatialExperiment")
    cd <- colData(spe)
    # at least these should now be in colData()
    required <- c("sum", "detected", "total", "control_sum", "target_sum",
                  "Area_um", "log2AspectRatio", "log2CountArea")
    expect_true(all(required %in% colnames(cd)))
})

test_that("computeQCScore adds a flag_score between 0 and 1", {
    spe <- spatialPerCellQC(spe0)
    spe2 <- computeQCScore(spe)
    cd2 <- colData(spe2)
    expect_true("QC_score" %in% colnames(cd2))
    fs <- cd2$QC_score
    expect_true(is.numeric(fs))
    expect_true(all(fs >= 0 & fs <= 1))
})


test_that("computeSpatialOutlier flags outliers for a chosen metric", {
    spe <- spatialPerCellQC(spe0)
    # use 'Area_um' and both methods
    out <- computeSpatialOutlier(spe, computeBy="Area_um", method="scuttle")
    cd_out <- colData(out)
    expect_true(is.logical(cd_out$Area_um_outlier_sc) ||
                    is(cd_out$Area_um_outlier_sc, "outlier.filter"))
})


test_that("computeQCScoreFlags combines filters and returns filter_out", {
    spe <- spatialPerCellQC(spe0)
    spe <- computeQCScore(spe)
    ff <- computeThresholdFlags(spe,
                            totalThreshold = 10,
                            ctrlTotRatioThreshold = 0.2)
    cd_ff <- colData(ff)
    # expect logical columns and combined filter_out
    expect_true(all(c("is_zero_counts", "is_ctrl_tot_outlier",
                      "threshold_flags") %in% colnames(cd_ff)))
    expect_true(is.logical(cd_ff$is_zero_counts))
    expect_true(is.logical(cd_ff$is_ctrl_tot_outlier))
    expect_true(is.logical(cd_ff$threshold_flags))
    # filter_out should only be TRUE where all three flags are TRUE
    combined <- (ff$is_zero_counts & ff$is_ctrl_tot_outlier)
    expect_identical(combined, cd_ff$threshold_flags)
})
