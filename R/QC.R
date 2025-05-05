#' spatialPerCellQC
#' @name spatialPerCellQC
#' @rdname spatialPerCellQC
#' @description
#' Computes quality‐control metrics for each cell and adds them to `colData`.
#'
#' @param spe A `SpatialExperiment` object containing spatial data.
#' @param micronConvFact Numeric factor to convert pixels to microns. Default
#'   `0.12`.
#' @param rmZeros logical for removing zero counts cells (default is TRUE).
#' @param negProbList Character vector of patterns to identify negative probes.
#'   Defaults include:
#'   - Nanostring CosMx: `"NegPrb"`, `"Negative"`, `"SystemControl"`
#'   - Xenium: `"NegControlProbe"`, `"NegControlCodeword"`,
#'     `"UnassignedCodeword"`
#'   - MERFISH: `"Blank"`
#'
#' @return A `SpatialExperiment` object with added QC metrics in `colData`.
#'
#' @details
#' Calculates sums and detected counts for control and target probes,
#' computes ratio and count‐area metrics, converts coords to microns for
#' CosMx, and drops zero‐count cells.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom scater addPerCellQC
#' @importFrom S4Vectors cbind.DataFrame
#' @export
#' @examples
#' example(readCosmxSPE)
#' spe <- spatialPerCellQC(spe)
spatialPerCellQC <- function(spe, micronConvFact=0.12, rmZeros=TRUE,negProbList=
    c("NegPrb", "Negative", "SystemControl", "Ms IgG1", "Rb IgG",
    "NegControlProbe", "NegControlCodeword", "UnassignedCodeword", "Blank")) {
    stopifnot(is(object=spe, "SpatialExperiment"))
    idxlist <- lapply(negProbList, function(ng) {
        grep(paste0("^", ng), rownames(spe))
    })
    names(idxlist) <- negProbList
    idxlist <- idxlist[which(lengths(idxlist)!=0)]
    spe <- addPerCellQC(spe, subsets=idxlist)
    idx <- grep("^subsets_.*_sum$", colnames(colData(spe)))
    npc <- npd <- 0
    if ( length(idx) !=0 ) {
        npc <- rowSums(as.matrix(colData(spe)[ , idx, drop=FALSE]))
        # TODO: not robust at all! the +1 is not a really good choice
        npd <- rowSums(as.matrix(colData(spe)[ , idx+1, drop=FALSE]))
    }
    spe$control_sum <- npc
    spe$control_detected <- npd
    spe$target_sum <- spe$sum - npc
    spe$target_detected <- spe$detected - npd
    if(!all(spatialCoordsNames(spe) %in% names(colData(spe)))) {
        # TODO: CHANGE SPE constructor WITH COORDINATES IN COLDATA
        colData(spe) <- cbind.DataFrame(colData(spe), spatialCoords(spe))
    }

    if(metadata(spe)$technology == "Nanostring_CosMx") {
        spnc <- spatialCoords(spe) * micronConvFact
        colnames(spnc) <- gsub("px", "um", spatialCoordsNames(spe))
        colData(spe) <- cbind.DataFrame(colData(spe), spnc)
        spe$Area_um <- spe$Area * (micronConvFact^2)
        spe <- .computeBorderDistanceCosMx(spe)
    }
    if (metadata(spe)$technology == "10X_Xenium") {
        spe$Area_um <- spe$cell_area # standardized across other techs
    }
    if ("AspectRatio" %in% colnames(colData(spe))) {
        spe$log2AspectRatio <- log2(spe$AspectRatio) # not cosmx
    } else { warning("Missing aspect ratio in colData") }
    spe$ctrl_total_ratio <- spe$control_sum/spe$total
    spe$ctrl_total_ratio[which(is.na(spe$ctrl_total_ratio))] <- 0
    spe$CountArea <- spe$sum/spe$Area_um
    spe$log2CountArea <- log2(spe$CountArea)
    if(sum(spe$sum==0) > 0) { # TODO: add a flag argument?
        message("Removing ", dim(spe[,spe$sum==0])[2], " cells with 0 counts!")
        spe <- spe[,!spe$sum==0]
    }
    return(spe)
}

#' .computeBorderDistanceCosMx
#' @name .computeBorderDistanceCosMx
#' @rdname dot-computeBorderDistanceCosMx
#' @description
#' Calculates the minimum distance of each cell to the field‐of‐view border
#' and adds it to `colData`.
#'
#' @param spe A `SpatialExperiment` object with CosMx data.
#' @param xwindim Width of FOV in x (default from `metadata(spe)$fov_dim`).
#' @param ywindim Height of FOV in y (default from `metadata(spe)$fov_dim`).
#'
#' @return A `SpatialExperiment` object with `dist_border` columns in
#' `colData`.
#'
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @keywords internal
.computeBorderDistanceCosMx <- function(spe,
                                    xwindim=metadata(spe)$fov_dim[["xdim"]],
                                    ywindim=metadata(spe)$fov_dim[["ydim"]])
{
    stopifnot(is(spe, "SpatialExperiment"))

        cd <- colData(spe)
        cdf <- left_join(as.data.frame(cd), metadata(spe)$fov_positions, by="fov")
        spcn <- spatialCoordsNames(spe)
        fovpn <- colnames(metadata(spe)$fov_positions)[colnames(
            metadata(spe)$fov_positions)%in%c("x_global_px", "y_global_px")]

    cd$dist_border_x <- pmin(cdf[,spcn[1]] - cdf[,fovpn[1]],
                            (cdf[,fovpn[1]] + xwindim) - cdf[,spcn[1]])

    cd$dist_border_y <- pmin(cdf[,spcn[2]] - cdf[,fovpn[2]],
                            (cdf[,fovpn[2]] + ywindim) - cdf[,spcn[2]])

    cd$dist_border <- pmin(cd$dist_border_x, cd$dist_border_y)
    colData(spe) <- cd
    return(spe)
}

#' computeSpatialOutlier
#' @name computeSpatialOutlier
#' @rdname computeSpatialOutlier
#' @description
#' Computes outliers based on the Area (in micron) of the experiment.
#' It gives the possibility to choose between the medcouple (mc method argument)
#' and the MADs (scuttle method argument).
#'
#' @details
#' The medcouple method is a measure for the skeweness of univariate
#' distribution as described in Hubert M. et al. (2008).
#' In particular, the computed medcouple value must be in a range between -0.6
#' and 0.6 to computed adjusted boxplots and perform the outlier detection.
#' For median absolute deviations (MADs) method we just wrap the isOutlier
#' function in the scuttle package. Please see McCarthy DJ et al (2017)
#' for further details.
#'
#' @param spe a SpatialExperiment object with target_counts, area in micron
#' and log2 of the aspect ratio in the `colData`.
#' @param compute_by character indicating a `colData` column name on which
#' compute the outlier.
#' @param method one of `mc`, `scuttle`, `both`.
#' Use `mc` for medcouple, `scuttle` for median absolute deviations as computed
#' in `scuttle`, `both` for computing both of them.
#' @param mcDoScale logical indicating if the values to compute the medcouple
#' for the outlier detection should be scaled (default is FALSE, as suggested
#' by the original Medcouple authors.). See \link[robustbase]{mc} for further
#' readings.
#' @param scuttleType One of `"both"`, `"lower"`, `"higher"` for scuttle method.
#'
#' @return a SpatialExperiment object with additional column(s) (named as
#' the column name indicated in `column_by` followed by the outlier_sc/mc
#' nomenclature) with the outlier detection as `outlier.filter` logical class
#' object. This allows to store the thresholds as attributes of the column.
#' use attr(,"thresholds") to retrieve them.
#'
#' @export
#' @importFrom robustbase mc adjbox
#' @importFrom e1071 skewness
#' @importFrom scuttle isOutlier outlier.filter
#'
#' @examples
#' example(spatialPerCellQC)
#' spe <- computeSpatialOutlier(spe, compute_by="log2CountArea", method="both")
#' table(spe$log2CountArea_outlier_mc)
#' table(spe$log2CountArea_outlier_sc)
computeSpatialOutlier <- function(spe, compute_by=NULL,
    method=c("mc", "scuttle", "both"), mcDoScale=FALSE,
    scuttleType=c("both", "lower", "higher")) {
    stopifnot(all(is(spe, "SpatialExperiment"), !is.null(compute_by)))
    stopifnot(compute_by %in% names(colData(spe)))
    options(mc_doScale_quiet=TRUE)
    method <- match.arg(method)
    scuttleType <- match.arg(scuttleType)
    cd <- colData(spe)
    cdcol <- cd[[compute_by]]
    mcfl<-scuttlefl<-FALSE
    switch(method, both={ mcfl<-scuttlefl<-TRUE },
            mc={ mcfl<-TRUE }, scuttle={ scuttlefl<-TRUE },
            {stop("Method is not one of allowed methods")} )
    if (mcfl) {
        skw <- e1071::skewness(cdcol, na.rm = TRUE) # NAs arise problems
        if (skw>-1 & skw<1) warning("Distribution is symmetric: ",
                "mc is for asymmetric distributions. Use scater instead.")
        mcval <- robustbase::mc(cdcol, doScale=mcDoScale, na.rm=TRUE)
        if ( any( (mcval <= -0.6), (mcval >= 0.6) ) )
            stop("mc is: ",round(mcval, digits=4),"outliers reqs not satisfied")
        names(cdcol) <- colnames(spe)
        outl <- robustbase::adjbox(cdcol, plot=FALSE)
        outsmc <- rep("NO", dim(cd)[1])
        outsmc[rownames(cd) %in% names(outl$out)] <-
            ifelse(outl$out <= outl$fence[1], "LOW", "HIGH")
        outlier_mc <- scuttle::outlier.filter(outsmc) #using scuttle class
        thrs <- as.numeric(outl$fence)
        names(thrs) <- c("lower", "higher")
        attr(outlier_mc, "thresholds") <- thrs
        cd$outlier_mc <- outlier_mc
        names(cd)[names(cd)=="outlier_mc"] <- paste0(compute_by, "_outlier_mc")
        # TODO: compute distributions in the adjusted boxplots to store in cd
    }
    if (scuttlefl) {
        outssc <- scuttle::isOutlier(cdcol, type=scuttleType)
        sctri <- rep("NO", dim(cd)[1])
        sctri <- ifelse(outssc==TRUE & cdcol<=attr(outssc, "thresholds")[1],
                        "LOW", sctri)
        outlier_sc <-ifelse(outssc==TRUE & cdcol>=attr(outssc, "thresholds")[2],
                        "HIGH", sctri)
        outlier_sc <- scuttle::outlier.filter(outlier_sc)
        attr(outlier_sc, "thresholds") <- attr(outssc, "thresholds")
        cd$outlier_sc <- outlier_sc
        names(cd)[names(cd)=="outlier_sc"] <- paste0(compute_by, "_outlier_sc")
    }
    colData(spe) <- cd
    return(spe)
}

#' computeFixedFlags
#' @name computeFixedFlags
#' @rdname computeFixedFlags
#' @description
#' Compute Flagged cells using fixed thresholds for SpatialExperiment.
#'
#' This function calculates various flags to identify outliers in a
#' `SpatialExperiment` object based on quality control metrics.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param total_threshold A numeric value for the threshold of total counts to
#' identify cells with zero counts. Default is `0`.
#' @param ctrl_tot_ratio_threshold A numeric value for the threshold of
#' control-to-total ratio to flag outliers. Default is `0.1`.
#'
#' @return The `SpatialExperiment` object with added filter flags in `colData`.
#'
#' @details The function flags cells based on zero counts, control-to-total
#' ratio, and `quality_score` to identify potential outliers. It also combines
#' these flags into a single filter flag.
#'
#' @importFrom SummarizedExperiment colData
#' @export
#' @examples
#' example(readCosmxSPE)
#' spe <- spatialPerCellQC(spe)
#' spe <- computeFixedFlags(spe)
#' table(spe$fixed_filter_out)
computeFixedFlags <- function(spe, total_threshold=0,
                            ctrl_tot_ratio_threshold=0.1)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("total" %in% names(colData(spe)))
    stopifnot("ctrl_total_ratio" %in% names(colData(spe)))

    spe$is_zero_counts <- ifelse(spe$total == total_threshold, TRUE, FALSE)
    #flagging cells with probe counts on total counts ratio > 0.1
    spe$is_ctrl_tot_outlier <- ifelse(spe$ctrl_total_ratio >
                                        ctrl_tot_ratio_threshold, TRUE, FALSE)

    spe$fixed_filter_out <- (spe$is_ctrl_tot_outlier &
                                spe$is_zero_counts)
    return(spe)
}

#' computeQScore
#' @name computeQScore
#' @rdname computeQScore
#' @description
#' Compute quality score and automatically define weights for quality score
#' through glm training. This function computes quality score with a formula
#' that depends on the technology.
#' To automatically define the formula coefficient weights, model training
#' is performed through ridge regression.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param verbose logical for having a verbose output. Default is FALSE.
#'
#' @return The `SpatialExperiment` object with added quality score in `colData`.
#' @export
#' @importFrom dplyr case_when filter mutate distinct
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats as.formula model.matrix predict quantile
#' @examples
#' example(spatialPerCellQC)
#' set.seed(1998)
#' spe <- computeQScore(spe)
#' summary(spe$training_status)
#' summary(spe$quality_score)
computeQScore <- function(spe, verbose=FALSE) {
    spe_temp <- computeSpatialOutlier(spe[,spe$total>0],
                        compute_by="log2CountArea", method="both")
    if(attr(spe_temp$log2CountArea_outlier_mc, "thresholds")[1] <
        min(spe_temp$log2CountArea)){
            low_thr <- quantile(spe$log2CountArea, probs = 0.01)
    } else{
        low_thr <- attr(spe_temp$log2CountArea_outlier_mc, "thresholds")[1]
    }
    high_thr <- attr(spe_temp$log2CountArea_outlier_mc, "thresholds")[2]
    spe$log2CountArea_outlier_train <- case_when(spe$total==0 ~ "NO",
        spe$log2CountArea<low_thr ~ "LOW",spe$log2CountArea>high_thr ~ "HIGH",
        TRUE ~ "NO")
    attr(spe$log2CountArea_outlier_train, "thresholds") <-
        attr(spe_temp$log2CountArea_outlier_mc, "thresholds")
    attr(spe$log2CountArea_outlier_train, "thresholds")[1] <- low_thr
    if(metadata(spe)$technology == "Nanostring_CosMx") {
        spe <- computeSpatialOutlier(spe, compute_by="log2AspectRatio",
                                        method="scuttle")
        train_bad <- data.frame(colData(spe)) |>
            filter((log2AspectRatio_outlier_sc == "HIGH" & dist_border < 50) |
                    (log2AspectRatio_outlier_sc == "LOW" & dist_border < 50) |
                    log2CountArea_outlier_train == "LOW") |>
            mutate(qscore_train = 0)
            train_good <- data.frame(colData(spe)) |>
            filter((log2AspectRatio > quantile(log2AspectRatio, probs = 0.25) &
                    log2AspectRatio < quantile(log2AspectRatio, probs = 0.75) &
                    dist_border > 50) |
                    (log2CountArea > quantile(log2CountArea, probs = 0.90) &
                    log2CountArea < quantile(log2CountArea, probs = 0.99))) |>
            mutate(qscore_train=1, is_a_bad_boy=cell_id%in%train_bad$cell_id)
        model_formula <- paste0("~ log2CountArea + I(abs(log2AspectRatio) ",
            "* as.numeric(dist_border<50)) + ",
            " log2CountArea:I(abs(log2AspectRatio)",
            "* as.numeric(dist_border<50))") #for cosmx
    }
    if(any(metadata(spe)$technology %in% c("10X_Xenium", "Vizgen_MERFISH"))) {
        train_bad <- data.frame(colData(spe)) |>
            filter(log2CountArea_outlier_train=="LOW") |> mutate(qscore_train=0)
        train_good <- data.frame(colData(spe)) |>
            filter((log2CountArea > quantile(log2CountArea, probs = 0.90) &
                    log2CountArea < quantile(log2CountArea, probs = 0.99))) |>
            mutate(qscore_train=1, is_a_bad_boy=cell_id%in%train_bad$cell_id)
        model_formula <- "~log2CountArea" # xen and merf
    }
    train_bad <- train_bad |> distinct(cell_id, .keep_all = TRUE)
    if(verbose) message("Chosen low quality examples: ", dim(train_bad)[1])
    train_good <- train_good |> distinct(cell_id, .keep_all = TRUE)
    train_good <- train_good[!train_good$is_a_bad_boy,]
    train_good <- train_good[sample(rownames(train_good), dim(train_bad)[1],
                                replace=FALSE),]
    if(verbose) message("Chosen good quality examples: ", dim(train_good)[1])
    train_good$is_a_bad_boy <- NULL
    train_df <- rbind(train_bad, train_good)
    train_df <- train_df |> distinct(cell_id, .keep_all = TRUE)
    model_matrix <- model.matrix(as.formula(model_formula), data = train_df)
    # set.seed(1998)
    model <- glmnet(x = model_matrix, y = train_df$qscore_train,
                            family = "binomial", lambda = NULL, alpha=0)
    # set.seed(1998)
    ridge_cv <- cv.glmnet(model_matrix, train_df$qscore_train,
                                family = "binomial", alpha = 0, lambda=NULL)
    best_lambda <- ridge_cv$lambda.min
    train_df$doom <- case_when(train_df$qscore_train==0 ~ "BAD",
                                train_df$qscore_train==1 ~ "GOOD")
    cd <- data.frame(colData(spe))
    full_matrix <- model.matrix(as.formula(model_formula), data = cd)
    cd$quality_score <- as.vector(predict(model, s = best_lambda,
                                        newx = full_matrix,
                                        type = "response"))
    spe$quality_score <- cd$quality_score
    train_identity <- rep("TEST", dim(spe)[2])
    spe$training_status <- dplyr::case_when(
        spe$cell_id%in%train_bad$cell_id ~ "BAD",
        spe$cell_id%in%train_good$cell_id ~ "GOOD",
        TRUE ~ train_identity)
    return(spe)
}

#' computeQScoreFlags
#' @name computeQScoreFlags
#' @rdname computeQScoreFlags
#' @description
#' Compute flagged cells based on a manually chosen threshold on quality score
#'
#' This function Compute flagged cells based on a manually chosen threshold on
#' quality score stored in `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param qs_threshold Numeric threshold or quantile for quality score. Default
#'   `0.5`.
#' @param use_qs_quantiles Logical; if `TRUE`, treat `qs_threshold` as a
#'   percentile.
#'
#' @return The `SpatialExperiment` object with added filter flags in `colData`.
#'
#'
#' @importFrom SummarizedExperiment colData
#' @export
#' @examples
#' example(computeQScore)
#' spe <- computeQScoreFlags(spe)
#' table(spe$is_qscore_outlier)
#' # if fixed filters are defined we have an additional column
#' spe <- computeFixedFlags(spe)
#' spe <- computeQScoreFlags(spe)
#' table(spe$fixed_qscore_out)
computeQScoreFlags <- function(spe, qs_threshold=0.5, use_qs_quantiles=FALSE)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("quality_score" %in% names(colData(spe)))

    if(use_qs_quantiles)
    {
            spe$is_qscore_opt_outlier <- ifelse(spe$quality_score <
                                            quantile(spe$quality_score,
                                                    probs=qs_threshold),
                                            TRUE, FALSE)
    } else {
            spe$is_qscore_outlier <- spe$quality_score < qs_threshold

    }

    if("fixed_filter_out" %in% names(colData(spe)))
    {
        spe$fixed_qscore_out <- (spe$is_qscore_outlier & spe$fixed_filter_out)
    }
    return(spe)
}
