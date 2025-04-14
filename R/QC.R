#' spatialPerCellQC
#' @description
#' Perform Per-Cell Quality Control on a SpatialExperiment Object
#'
#' This function calculates quality control metrics for each cell in a
#' `SpatialExperiment` object and adds them to `colData`.
#'
#' @param spe A `SpatialExperiment` object containing spatial transcriptomics
#' data.
#' @param micronConvFact A numeric value for converting pixel dimensions to
#' microns. Default is `0.12`.
#' @param negProbList A character vector of patterns used to identify negative
#' probes.
#' Default values are:
#' - For Nanostring CosMx: `"NegPrb"`, `"Negative"`, `"SystemControl"`
#' - For Xenium: `"NegControlProbe"`, `"NegControlCodeWord"`,
#' `"UnassignedCodeWord"`
#' - For MERFISH: `"Blank"`
#'
#'
#' @return The `SpatialExperiment` object with added quality control metrics in
#' `colData`.
#'
#' @details The function computes several QC metrics, including control probe
#' sums, target probe sums, and the ratio of control probes to the total.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom scater addPerCellQC
#' @importFrom S4Vectors cbind.DataFrame
#' @export
#' @example
#' # TBD
spatialPerCellQC <- function(spe, micronConvFact=0.12,
    negProbList=c("NegPrb", "Negative", "SystemControl", # CosMx
        "NegControlProbe", "NegControlCodeWord", "UnassignedCodeWord", # Xenium
        "Blank" # MERFISH
    ))
{
    ## CHECK EXISTENCE OF POLYGONS/AREAS ETC -> create function for
    ## metrics creation
    stopifnot(is(object=spe, "SpatialExperiment"))
    idxlist <- lapply(negProbList, function(ng){
        grep(paste0("^", ng), rownames(spe))
    })
    names(idxlist) <- negProbList
    idxlist <- idxlist[which(lengths(idxlist)!=0)]

    spe <- addPerCellQC(spe, subsets=idxlist)
    idx <- grep("^subsets_.*_sum$", colnames(colData(spe)))
    npc = npd = 0
    if ( length(idx) !=0 )
    {
        ## TO TEST -> BENEDETTA
        # meglio dataframe, perché rowsums non funziona -> ?
        npc <- rowSums(as.matrix(colData(spe)[ , idx, drop=FALSE])) #sum
        ## getting detected probes as the column suddenly after the sum column
        ## # not robust at all! the +1 is not a really good choice
        npd <- rowSums(as.matrix(colData(spe)[ , idx+1, drop=FALSE])) #detected
    }

    spe$control_sum <- npc
    spe$control_detected <- npd
    spe$target_sum <- spe$sum - npc
    spe$target_detected <- spe$detected - npd

    if(!all(spatialCoordsNames(spe) %in% names(colData(spe))))
    {
        #### CHANGE SPE constructor WITH COORDINATES IN COLDATA #########
        colData(spe) <- cbind.DataFrame(colData(spe), spatialCoords(spe))
    }

    if(metadata(spe)$technology == "Nanostring_CosMx")
    {
        spnc <- spatialCoords(spe) * micronConvFact
        colnames(spnc) <- gsub("px", "um", spatialCoordsNames(spe))
        colData(spe) <- cbind.DataFrame(colData(spe), spnc)
        spe$Area_um <- spe$Area * (micronConvFact^2)
        spe <- computeBorderDistanceCosMx(spe)
    }

    # adding this line so that Area in Xenium has the same name as the others
    # already in um ## move to a specific standardized function (?)
    if (metadata(spe)$technology == "10X_Xenium"){
        spe$Area_um <- spe$cell_area
    }

    #### compute AspectRatio for other technologies ####
    if ("AspectRatio" %in% colnames(colData(spe)))
    {
        spe$log2AspectRatio <- log2(spe$AspectRatio)
    } else {
        warning(paste0("Missing aspect ratio in colData...\n",
                    "NB: This could lead to additional warnings or errors.\n",
                    "Missing AspectRatio can be computed by loading polygons."))
    }
    spe$ctrl_total_ratio <- spe$control_sum/spe$total
    spe$ctrl_total_ratio[which(is.na(spe$ctrl_total_ratio))] <- 0
    # changed to Area um, now it's the same for every technology
    spe$CountArea <- spe$sum/spe$Area_um
    spe$log2CountArea <- log2(spe$CountArea)
    ## adding a flag argument (?)
    message("Removing ", dim(spe[,spe$sum==0])[2], " cells, they have 0 counts")
    spe <- spe[,!spe$sum==0]

    return(spe)
}

#' computeBorderDistanceCosMx
#' @description
#' Compute Distance to FoV Border in SpatialExperiment for CosMx technology.
#'
#' Calculates the minimum distance of each coordinate in a `SpatialExperiment`
#' object to the nearest border of the field of view (FOV) and adds it to
#' `colData`.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param xwindim Width of the FOV in the x-dimension. Defaults to
#' `metadata(spe)$fov_dim[["xdim"]]`.
#' @param ywindim Height of the FOV in the y-dimension. Defaults to
#' `metadata(spe)$fov_dim[["ydim"]]`.
#'
#' @return The `SpatialExperiment` object with added border distance data in
#' `colData`.
#'
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment colData
#' @export
#' @example
#' # TBD
## move to a specific reading/missing metrics function (?)
computeBorderDistanceCosMx <- function(spe,
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
#'
#' @return a SpatialExperiment object with additional column(s) (named as
#' the column name indicated in `column_by` followed by the outlier_sc/mc
#' nomenclature) with the outlier detection as `outlier.filter` logical class
#' object. This allows to store the thresholds as attributes of the column.
#' use attr(,"thresholds") to retrieve them.
#' @export
#' @importFrom robustbase mc adjbox
#' @importFrom e1071 skewness
#' @importFrom scuttle isOutlier outlier.filter
#'
#' @examples
#' # TBD
computeSpatialOutlier <- function(spe, compute_by=NULL,
                                  method=c("mc", "scuttle", "both"),
                                  mcDoScale=FALSE,
                                  scuttleType=c("both", "lower", "higher"))
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(!is.null(compute_by))
    stopifnot(compute_by %in% names(colData(spe)))

    options(mc_doScale_quiet=TRUE)

    method <- match.arg(method)
    scuttleType <- match.arg(scuttleType)
    cd <- colData(spe)
    cdcol <- cd[[compute_by]]
    mcfl=scuttlefl=FALSE
    switch(method,
           both={ mcfl=scuttlefl=TRUE },
           mc={ mcfl=TRUE },
           scuttle={ scuttlefl=TRUE },
           {stop("Method is not one of allowed methods")}
    )

    if (mcfl)
    {
        # NAs impede the following steps
        skw <- e1071::skewness(cdcol, na.rm = TRUE)

        # this was changed because I realized I didn't know well how to interpret
        # skewness, if skw is >-1 and <1 the distribution is symmetric
        if (skw>-1 & skw<1) warning("The distribution is symmetric. ",
                "The medcouple is not suited for symmetric distributions. ",
                "In your case we suggest to use the scater method instead.")
        # NAs impede the following steps
        mcval <- robustbase::mc(cdcol, doScale=mcDoScale, na.rm=TRUE)
        if ( any( (mcval <= -0.6), (mcval >= 0.6) ) )
            stop("Obtained medcouple value is: ", round(mcval, digits=4),
                 "\nIt doesn't meet the needed requirements for outlier",
                 " identification with this method.")
        names(cdcol) <- colnames(spe)
        outl <- robustbase::adjbox(cdcol, plot=FALSE)
        # now outliers are defined as NO, HIGH or LOW, no more TRUE or FALSE,
        # distinction whether they have values < lower thr. or > than higher thr.
        outsmc <- rep("NO", dim(cd)[1])
        outsmc[rownames(cd) %in% names(outl$out)] <-
            ifelse(outl$out <= outl$fence[1], "LOW", "HIGH")
        ## using scuttle outlier.filter class to store fences for each filtering
        ## in the attributes of the class
        outlier_mc <- scuttle::outlier.filter(outsmc)
        thrs <- as.numeric(outl$fence)
        names(thrs) <- c("lower", "higher")
        attr(outlier_mc, "thresholds") <- thrs
        cd$outlier_mc <- outlier_mc
        names(cd)[names(cd)=="outlier_mc"] <- paste0(compute_by, "_outlier_mc")
        # metadata(spe)$outlier_fences[[compute_by]] <- outl$fence ## DEPRECATED

        ## TODO: compute distributions in the adjusted boxplots to store them
        ## in the coldata
    }

    if (scuttlefl)
    {
        outssc <- scuttle::isOutlier(cdcol, type=scuttleType)
        # outliers are defined as NO, LOW or HIGH even here for scuttle method
        sctri <- rep("NO", dim(cd)[1])
        sctri <- ifelse(outssc==TRUE&cdcol <= attr(outssc, "thresholds")[1],
                        "LOW", sctri)
        outlier_sc <- ifelse(outssc==TRUE&cdcol>=attr(outssc, "thresholds")[2],
                        "HIGH", sctri)
        # this is to keep the attributes even when filtering or they disappear
        outlier_sc <- scuttle::outlier.filter(outlier_sc)
        # adding thresholds also for scuttle
        attr(outlier_sc, "thresholds") <- attr(outssc, "thresholds")
        cd$outlier_sc <- outlier_sc
        names(cd)[names(cd)=="outlier_sc"] <- paste0(compute_by, "_outlier_sc")
    }

    colData(spe) <- cd
    return(spe)
}

#' computeFixedFlags
#' @description
#' Compute Flagged cells using fixed thresholds for SpatialExperiment
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
#' # TBD
computeFixedFlags <- function(spe,
                              total_threshold=0,
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
#' @description
#' Compute quality score and automatically define weights for quality score
#' through glm training
#'
#' This function computes quality score with a formula that depends on the
#' technology.
#' To automatically define the formula coefficient weights, model training
#' is performed
#' through ridge regression.
#'
#' `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#'
#' @return The `SpatialExperiment` object with added quality score in `colData`.
#'
#' @details
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter mutate full_join distinct case_when
#' @importFrom glmnet glmnet cv.glmnet
#' @export
#' @examples
#' # TBD
computeQScore <- function(spe=spe) {

    # this is necessary because I found Xenium and Merfish datasets with 0
    #  counts cells having log2CountArea as -Inf
    # something our functions such as medcouple and skewness don't deal with,
    # thus I exclude cells
    # with 0 counts and don't use them for training

    spe_temp <- computeSpatialOutlier(spe[,spe$total>0],
                        compute_by="log2CountArea", method="both")
    if(attr(spe_temp$log2CountArea_outlier_mc, "thresholds")[1] <
        min(spe_temp$log2CountArea)){
            low_thr <- quantile(spe$log2CountArea, probs = 0.01)
    } else{
        low_thr <- attr(spe_temp$log2CountArea_outlier_mc, "thresholds")[1]
    }
    high_thr <- attr(spe_temp$log2CountArea_outlier_mc, "thresholds")[2]

    spe$log2CountArea_outlier_train <- case_when(spe$total == 0 ~ "NO",
                                         spe$log2CountArea < low_thr ~ "LOW",
                                         spe$log2CountArea > high_thr ~ "HIGH",
                                         TRUE ~ "NO")

    attr(spe$log2CountArea_outlier_train, "thresholds") <-
        attr(spe_temp$log2CountArea_outlier_mc, "thresholds")

    attr(spe$log2CountArea_outlier_train, "thresholds")[1] <- low_thr

    print("How many LOW outliers cells will be used for training?")
    print(table(spe$log2CountArea_outlier_train))

    if(metadata(spe)$technology == "Nanostring_CosMx")
    {
        # I selected scuttle because otherwise it would return a warning, since
        # aspect ratio is always normal in every technology, I would turn to scuttle
        # instead
        spe <- computeSpatialOutlier(spe, compute_by="log2AspectRatio", method="scuttle")
        print("How many outliers were found for log2AspectRatio")
        print(table(spe$log2AspectRatio_outlier_sc))

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

        # model formula definition for glmnet training for CosMx technology

        model_formula <- paste0("~ log2CountArea + I(abs(log2AspectRatio) ",
            "* as.numeric(dist_border<50)) + ",
            " log2CountArea:I(abs(log2AspectRatio)",
            "* as.numeric(dist_border<50))")
    }

    if(metadata(spe)$technology == "10X_Xenium" |
       metadata(spe)$technology =="Vizgen_MERFISH") {

        train_bad <- data.frame(colData(spe)) |>
            filter(log2CountArea_outlier_train == "LOW") |>
            mutate(qscore_train = 0)

        train_good <- data.frame(colData(spe)) |>
            filter((log2CountArea > quantile(log2CountArea, probs = 0.90) &
                    log2CountArea < quantile(log2CountArea, probs = 0.99))) |>
            mutate(qscore_train=1, is_a_bad_boy=cell_id%in%train_bad$cell_id)

        # model formula definition for glmnet training for
        # Xenium and Merfish technology
        model_formula <- "~log2CountArea"
    }

    # check of duplicates

    # bad example duplicates removal without any warning to the user

    train_bad <- train_bad |> distinct(cell_id, .keep_all = TRUE)

    message(paste0("Chosen low quality examples: ", dim(train_bad)[1]))

    # good example duplicates removal without any warning to the user

    train_good <- train_good |> distinct(cell_id, .keep_all = TRUE)

    train_good <- train_good[!train_good$is_a_bad_boy,]
    train_good <- train_good[sample(rownames(train_good), dim(train_bad)[1],
                                replace=FALSE),]

    message(paste0("Chosen good quality examples, (should be the same number
              of bad quality examples): ", dim(train_good)[1]))

    # merge into same training dataset

    # I would prefer this first way to merge the two datasets together without
    # creating any duplicate column, however this didn't work with merfish tech,
    # having cell ids as only numbers, that is
    # probably why this is not well handled inside the following functions,
    # I opted for an rbind instead

    #train_df <- train_bad %>%
    #mutate(rn = data.table::rowid(cell_id)) %>%
    #full_join(train_good %>%
    #mutate(rn = data.table::rowid(cell_id))) %>%
    #select(-rn)

    train_good$is_a_bad_boy <- NULL
    train_df <- rbind(train_bad, train_good)

    train_df <- train_df |> distinct(cell_id, .keep_all = TRUE)

    # glmnet training

    model_matrix <- model.matrix(as.formula(model_formula), data = train_df)

    set.seed(1998)
    model <- glmnet(x = model_matrix, y = train_df$qscore_train,
                            family = "binomial", lambda = NULL, alpha=0)
    ## include training into a .trainModel function
    set.seed(1998)
    #https://stackoverflow.com/questions/34677526/set-seed-with-cv-glmnet-paralleled-gives-different-results-in-r
    # it says that if cv.glmnet is run in different days, results will differ
    # and that you need to define nFolds manually in the function call if you
    # want reproducible results, because they are randomly chosen.
    ridge_cv <- cv.glmnet(model_matrix, train_df$qscore_train,
                                  family = "binomial", alpha = 0, lambda=NULL)

    best_lambda <- ridge_cv$lambda.min

    print("Model coefficients for every term used in the formula:")
    print(round(predict(model, s = best_lambda, type="coefficients"),2))

    train_df$doom <- case_when(train_df$qscore_train==0 ~ "BAD",
                               train_df$qscore_train==1 ~ "GOOD")

    cd <- data.frame(colData(spe))

    full_matrix <- model.matrix(as.formula(model_formula), data = cd)

    cd$quality_score <- as.vector(predict(model, s = best_lambda,
                                          newx = full_matrix,
                                          type = "response"))
    spe$quality_score <- cd$quality_score

    # column with information about which cells were used in the training
    # loaded in spe

    # again, preferred code but Merfish cellids giving problems

    #cd <- left_join(cd, train_df, by = "cell_id")
    #cd$doom <- tidyr::replace_na(cd$doom, "TEST")
    #spe$doom <- cd$doom

    # error output: Error in attributes(lst) <- a :
    #'names' attribute [395215] must be the same length as the vector [5242]

    train_identity <- rep("TEST", dim(spe)[2])
    spe$doom <- dplyr::case_when( ## change into something $training_status !
        spe$cell_id%in%train_bad$cell_id ~ "BAD",
        spe$cell_id%in%train_good$cell_id ~ "GOOD",
        TRUE ~ train_identity)

    return(spe)
}

#' computeQscoreFlags
#' @description
#' Compute flagged cells based on a manually chosen threshold on quality score
#'
#' This function Compute flagged cells based on a manually chosen threshold on
#' quality score stored in `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object with spatial transcriptomics data.
#' @param use_qs_quantiles a boolean. If TRUE uses the value specified in
#' qs_threshold as a percentile
#' @param opt a boolean value to set to TRUE if you want to compute flagged cells
#' for optimized quality score
#'
#' @return The `SpatialExperiment` object with added filter flags in `colData`.
#'
#' @details
#'
#' @importFrom SummarizedExperiment colData
#' @export
#' @examples
#' # TBD
#'
computeQscoreFlags <- function(spe, qs_threshold=0.5,
                               use_qs_quantiles=FALSE)
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
        spe$fixed_qscore_out <- (spe$is_qscore_outlier & spe$fixed_filter_out)

    return(spe)
}
