#' plotCellsFovs
#' @name plotCellsFovs
#' @rdname plotCellsFovs
#'
#' @description
#' Plot cell centroids in FoVs
#' Creates a scatter plot with cell centroids arranged in their FoVs as an
#' overlapping grid.
#'
#' @param spe A `SpatialExperiment` object with `fov` in `colData`.
#' @param sample_id Character string identifying which sample to plot.
#'   Default: `unique(spe$sample_id)`.
#' @param point_col Color for the cell centroids. Default: `"firebrick"`.
#' @param numbers_col Color for the FoV labels. Default: `"black"`.
#' @param alpha_numbers Numeric transparency for FoV labels. Default: `0.8`.
#' @param fov_dim numeric with two named dimensions xdim, ydim. (Default is
#' metadata(spe)$fov_dim)
#'
#' @return A `ggplot` object showing cell centroids and FoV boundaries.
#'
#' @importFrom ggplot2 ggplot aes annotate geom_point geom_text ggtitle
#' coord_fixed
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' example(readCosmxSPE)
#' g <- plotCellsFovs(spe)
#' print(g)
plotCellsFovs <- function(spe, sample_id=unique(spe$sample_id),
                        point_col="firebrick", numbers_col="black",
                        alpha_numbers=0.8, fov_dim=metadata(spe)$fov_dim)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("fov" %in% names(colData(spe)))
    stopifnot( all(names(fov_dim) %in% c("xdim","ydim")) )

    spd <- as.data.frame(spatialCoords(spe))
    x_coord <- spatialCoordsNames(spe)[1]
    y_coord <- spatialCoordsNames(spe)[2]
    ggp <- ggplot() +
        geom_point(data=spd, mapping=aes(x=.data[[x_coord]],
                                        y=.data[[y_coord]]),
                    colour=point_col,
                    fill=point_col,
                    size=0.05, alpha=0.8) +
        annotate("rect",
            xmin=metadata(spe)$fov_positions["x_global_px"][ , , drop=TRUE],
            xmax=metadata(spe)$fov_positions["x_global_px"][ , , drop=TRUE] +
                fov_dim[["xdim"]],
            ymin=metadata(spe)$fov_positions["y_global_px"][ , , drop=TRUE],
            ymax=metadata(spe)$fov_positions["y_global_px"][ , , drop=TRUE] +
                fov_dim[["ydim"]],
            alpha=.2, color="black", linewidth=0.2) +
        geom_text(aes(x=metadata(spe)$fov_positions["x_global_px"][,,drop=TRUE]+
                        fov_dim[["xdim"]]/2,
                    y=metadata(spe)$fov_positions["y_global_px"][,,drop=TRUE]+
                        fov_dim[["ydim"]]/2,
                    label=metadata(spe)$fov_positions["fov"][,,drop=TRUE]),
                    color=numbers_col, fontface="bold", alpha=alpha_numbers) +
        ggtitle(sample_id) +
        .fov_image_theme(back.color="white", back.border="white",
                        title.col="black") + ggplot2::coord_fixed()
    return(ggp)
}

#' plotCentroids
#' @name plotCentroids
#' @rdname plotCentroids
#' @description
#' Plot Spatial Coordinates for a SpatialExperiment Object
#' This function generates a ggplot of spatial coordinates from a
#' `SpatialExperiment` object, optionally coloring the points by a specified
#' column in `colData`.
#'
#' @param spe A `SpatialExperiment` object containing spatial
#' transcriptomics data.
#' @param colour_by An optional character string specifying the column in
#' `colData(spe)` to use for coloring the points. If `NULL`, all points will be
#' colored the same.
#' @param colour_log Logical to log-transform the data to enhance visualization
#' (Default is FALSE).
#' @param sample_id A character string specifying the sample identifier to be
#' used as the plot title. (Default is the unique sample ID from `spe`)
#' @param isNegativeProbe A logical value indicating whether to apply a custom
#' color gradient for negative probe data. (Default is `FALSE`)
#' @param palette A vector of colors to be used as a custom palette. For
#' categorical data, this should be a vector of colors with the same length as
#' the number of levels in `colour_by`. For continuous data, this should be a
#' vector of colors used to create a gradient.
#' @param point_col A character string specifying the color of the points when
#' `colour_by` is `NULL`. (Default is `"darkmagenta"`)
#' @param size A numeric value specifying the size of the points. (Default is
#' `0.05`)
#' @param alpha A numeric value specifying the transparency level of the points.
#' (Default is `0.2`)
#' @param aspect_ratio A numeric value specifying the aspect ratio of the plot.
#' (Default is `1`)
#' @return A `ggplot` object representing the spatial coordinates plot of
#' polygon centroids.
#'
#' @import SpatialExperiment
#' @importFrom ggplot2 geom_point aes_string theme_bw theme ggtitle guides
#' guide_legend scale_color_gradient scale_color_manual scale_color_gradientn
#' @importFrom scater plotColData
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot geom_point aes ggtitle theme_bw coord_fixed
#' @importFrom ggplot2 scale_color_manual scale_color_gradientn
#' @export
#'
#' @examples
#' example(readCosmxSPE)
#' g <- plotCentroids(spe, colour_by="Mean.DAPI")
#' print(g)
plotCentroids <- function(spe, colour_by=NULL, colour_log=FALSE,
                        sample_id=unique(spe$sample_id),
                        isNegativeProbe=FALSE, palette=NULL,
                        point_col="darkmagenta", size=0.05, alpha=0.8,
                        aspect_ratio=1)
{
    stopifnot(is(spe, "SpatialExperiment"))
    if(is.null(colour_by)) {
        ggp <- ggplot(data.frame(spatialCoords(spe)),
                    aes(x=.data[[spatialCoordsNames(spe)[1]]],
                        y=.data[[spatialCoordsNames(spe)[2]]])) +
            geom_point(colour=point_col, fill=point_col,
                        size=size, alpha=alpha) + ggplot2::ggtitle(sample_id) +
            ggplot2::theme_bw() + ggplot2::coord_fixed()
    } else {
        if(colour_log) {
            stopifnot(colour_by %in% names(colData(spe)))
            colour_byo <- colour_by
            colour_by <- paste0("log(", colour_byo, ")")
            colData(spe)[[colour_by]] <- log1p(colData(spe)[[colour_byo]])
        }
        ggp <- ggplot(data.frame(colData(spe), spatialCoords(spe)),
                    aes(x=.data[[spatialCoordsNames(spe)[1]]],
                        y=.data[[spatialCoordsNames(spe)[2]]],
                        colour=.data[[colour_by]],
                        fill=.data[[colour_by]])) +
            geom_point(size=size, alpha=alpha) + coord_fixed()
        if(!is(spe[[colour_by]], "factor")) {
            ggp <- ggp + ggplot2::scale_colour_viridis_c() +
                ggplot2::scale_fill_viridis_c() + coord_fixed()
        }

        if(isNegativeProbe) {
            ggp <- ggp + scale_color_gradient(low="white", high="red",
                name=colour_by) + .negative_image_theme()
        } else if(all(!is.null(palette), (palette %in% names(colData(spe))))) {
            palette <- createPaletteFromColData(spe, palette_names=colour_by,
                                                    palette_colors=palette)
            if(is.factor(colData(spe)[[colour_by]])) {
                ggp <- ggp + scale_color_manual(values=palette)
            } else if(is.numeric(colData(spe)[[colour_by]])) {
                ggp <- ggp + scale_color_gradientn(colors=palette)
            }
        }
    }
    ggp <- ggp + ggtitle(sample_id) +
        theme(aspect.ratio=aspect_ratio, plot.title=element_text(hjust=0.5))
    if(!isNegativeProbe) ggp <- ggp + theme_bw()
    return(ggp)
}


#' plotMetricHist
#' @name plotMetricHist
#' @rdname plotMetricHist
#' @description Plot a Histogram for a Given Metric in a SpatialExperiment
#' Object
#'
#' This function generates a histogram for a specified metric in a
#' `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object containing spatial transcriptomics
#' data.
#' @param metric A character string specifying the name of the metric
#' (column in `colData(spe)`) to plot.
#' @param fill_color A character string specifying the fill color of the
#' histogram bars. (Default is `"#69b3a2"`)
#' @param use_fences A character string specifying the name of the column in
#' `colData(spe)` that contains the fence thresholds (typically from an outlier
#' filter). If `NULL`, no fences will be plotted. (Default is `NULL`)
#' @param fences_colors A named character vector specifying the colors to use
#' for the lower and higher fences. The names should be `"lower"` and `"higher"`
#'. (Default is `c("lower"="purple4", "higher"="tomato")`)
#' @param bins An integer specifying the number of bins to use in the histogram.
#' (Default is `30`)
#' @param bin_width A numeric value specifying the width of the bins. If `NULL`,
#' the bin width will be automatically determined based on the `bins` parameter.
#' (Default is `NULL`)
#'
#' @return A `ggplot` object representing the histogram of the specified metric.
#'
#' @importFrom ggplot2 ggplot geom_histogram aes ggtitle theme_bw geom_vline
#' labs scale_colour_manual
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#' @export
#' @examples
#' example(readCosmxSPE)
#' g <- plotMetricHist(spe, metric="Mean.DAPI")
#' print(g)
plotMetricHist <- function(spe, metric, fill_color="#c0c8cf",
        use_fences=NULL, fences_colors=c("lower"="purple4", "higher"="tomato"),
        bins=30, bin_width=NULL)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(metric %in% names(colData(spe)))

    ggp <- ggplot(data=as.data.frame(colData(spe))) +
            geom_histogram(aes(x=.data[[metric]]), fill=fill_color,
                    bins=bins, binwidth=bin_width)
    if (!is.null(use_fences))
    {
        stopifnot(use_fences %in% names(colData(spe)))
        fences <- getFencesOutlier(spe, use_fences, "both", 2)
        fences_labs <- paste0(names(fences), ": ", fences)
        names(fences_colors) <- fences_labs
        ggp <- ggp +
            geom_vline(aes(xintercept=fences[1], color=fences_labs[1])) +
            geom_vline(aes(xintercept=fences[2], color=fences_labs[2])) +
            labs(color=paste0("Fences ", use_fences)) +
            scale_colour_manual(values=fences_colors)
    }
    ggp <- ggp + ggtitle(metric) + theme_bw()

    return(ggp)
}


#' plotPolygons
#' @name plotPolygons
#' @rdname plotPolygons
#'
#' @description Plot polygons from a `SpatialExperiment` object using ggplot2.
#'
#' @param spe A `SpatialExperiment` object with polygon data as an `sf` object.
#' @param colour_by A column in `colData(spe)` for coloring the polygons or a
#' string color in colors(). (Default is "darkgrey")
#' @param colour_log Logical to log-transform the data to enhance visualization
#' (Default is FALSE).
#' @param sample_id Sample ID for plot title. Default is the unique sample ID.
#' @param fill_alpha Transparency level for polygon fill. Default is `1`.
#' @param palette Colors to use if `colour_by` is a factor. Default is `NULL`.
#' @param border_col Color of polygon borders. Default is `"black"`.
#' @param border_alpha Transparency level for borders. Default is `1`.
#' @param border_line_width Width of polygon borders. Default is `0.1`.
#' @param draw_borders Logical; whether to draw borders. Default is `TRUE`.
#' @param poly_column character for the name of the column where the polygons sf
#' are stored (default is "polygons.global")
#' @param bg_color character indicating color for the background
#' (default is "white")
#'
#' @return A `ggplot` object representing the polygon plot of the spatial data.
#' @export
#'
#' @importFrom ggplot2 ggplot geom_sf aes scale_fill_manual scale_fill_viridis_c
#' scale_fill_identity theme_minimal theme element_text margin labs
#' @importFrom sf st_as_sf
#' @importFrom SummarizedExperiment colData
#' @importFrom grDevices colors
#'
#' @examples
#' example(readAndAddPolygonsToSPE)
#' plotPolygons(spe, colour_by="Mean.DAPI")
plotPolygons <- function(spe, colour_by="darkgrey", colour_log=FALSE,
    poly_column="polygons.global", sample_id=unique(spe$sample_id),
    bg_color="white", fill_alpha=1, palette=NULL, border_col=NA,
    border_alpha=1, border_line_width=0.1, draw_borders=TRUE) {
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("polygons" %in% names(colData(spe)))
    df <- data.frame(colData(spe))
    polflag <- FALSE
    if(!is.null(colour_by)) {
        if(colour_by %in% names(colData(spe))) {
            if(colour_log) {
                colour_byo <- colour_by
                colour_by <- paste0("log(", colour_byo, ")")
                df[[colour_by]] <- log1p(df[[colour_byo]])
            }
            polflag <- TRUE
        } else {
            if(!(colour_by %in% grDevices::colors())) {
                warning("colour_by not in known colors nor in colData")
                colour_by <- "darkgrey"
            }}}
    border_params <- if(draw_borders) {
        list(color=border_col, size=border_line_width)
    } else { list(color=NA, size=0) }
    if(polflag) {
        p <- ggplot(df, aes(geometry=.data[[poly_column]],
            fill=.data[[colour_by]])) + geom_sf(alpha=fill_alpha,
            color=border_params$color, size=border_params$size)
    } else {
        p <- ggplot(df, aes(geometry = .data[[poly_column]])) +
            geom_sf(fill=colour_by, alpha=fill_alpha,
                color=border_params$color, size=border_params$size)
    }
    if(!is.null(colour_by) && (is.factor(df[[colour_by]]) ||
                                is.logical(df[[colour_by]]))) {
        if(!is.null(palette)) { p <- p + scale_fill_manual(values=palette) }
    } else if(!is.null(colour_by) && is.numeric(df[[colour_by]])) {
        p <- p + scale_fill_viridis_c(option="D")
    } else { p <- p + scale_fill_identity() }
    p <- p + theme_minimal() + theme(legend.position="right",
        plot.title.position="plot",plot.title=element_text(face="bold",size=14),
            plot.margin=margin(0, 0, 0, 0),
            panel.background=element_rect(fill=bg_color, color=bg_color,size=1),
            panel.border = element_rect(color="black", fill=NA, linewidth=0.1),
            panel.grid.minor=element_blank()
        ) + labs(title=sample_id, fill=colour_by)
    return(p)
}


#' plotZoomFovsMap
#' @name plotZoomFovsMap
#' @rdname plotZoomFovsMap
#' @description
#'
#' Plot Zoomed-in FOVs with Map and Polygons
#'
#' This function generates a plot that shows a map of all fields of view (FOVs)
#' within a `SpatialExperiment` object, alongside a zoomed-in view of the
#' specified FOVs with an overlay of polygons and optional coloring.
#'
#' @param spe A `SpatialExperiment` object containing spatial transcriptomics
#' data.
#' @param fovs A character vector specifying the FOVs to be zoomed in and
#' plotted. Must match values in the `fov` column of `colData(spe)`.
#' @param map_point_col A character string specifying the color of the points
#' in the map. Default is `"darkmagenta"`.
#' @param map_numbers_col A character string specifying the color of the
#' numbers on the map. Default is `"black"`.
#' @param map_alpha_numbers A numeric value specifying the transparency of the
#' numbers on the map. Default is `0.8`.
#' @param title An optional character string specifying the title of the final
#' plot. If `NULL`, no title is added. Default is `NULL`.
#' @param ... Additional arguments passed to `plotPolygons`.
#'
#' @return A combined plot showing a map of all FOVs with zoomed-in views of
#' the specified FOVs and their associated polygons.
#'
#' @details The function first filters the `SpatialExperiment` object to the
#' specified FOVs, generates a plot of the cells for the entire map, then
#' creates a detailed polygon plot of the selected FOVs, and finally combines
#' these into a single side-by-side visualization. If `title` is not `NULL`, it
#' adds a title to the combined plot.
#'
#' @importFrom ggpubr ggarrange annotate_figure
#' @importFrom tmap tmap_grob
#' @export
#'
#' @examples
#' example(readAndAddPolygonsToSPE)
#' plotZoomFovsMap(spe, fovs=16, title="FOV 16")
plotZoomFovsMap <- function(spe, fovs=NULL,
                            map_point_col="darkmagenta",
                            map_numbers_col="black",
                            map_alpha_numbers=0.8,
                            title=NULL, ...) {
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("fov" %in% names(colData(spe)))
    stopifnot(all(fovs %in% spe$fov))
    spefovs <- spe[, spe$fov %in% fovs]
    map <- plotCellsFovs(spefovs, point_col=map_point_col,
        numbers_col=map_numbers_col, alpha_numbers=map_alpha_numbers,
        sample_id=NULL)
    g2 <- plotPolygons(spefovs, sample_id=NULL, ...)
    final_plot <- ggpubr::ggarrange(map, g2, ncol=2)
    if (!is.null(title)) {
        final_plot <- ggpubr::annotate_figure(final_plot,
            top=ggpubr::text_grob(title, face="bold", size=14))
    }
    return(final_plot)
}

#' plotQScoreTerms
#' @name plotQScoreTerms
#' @rdname plotQScoreTerms
#'
#' @description
#' Plots the individual terms that combine into the quality score formula,
#' allowing assessment of each term’s impact on the final score.
#'
#' @param spe A `SpatialExperiment` object with `quality_score` and term
#'   columns in `colData`.
#' @param sample_id Character string for plot title. Must match values in the
#'   `fov` column of `colData(spe)`. Default: `unique(spe$sample_id)`.
#' @param size Numeric point size for the scatter plots. Default: `0.05`.
#' @param alpha Numeric transparency for the scatter plots. Default: `0.2`.
#' @param aspect_ratio Numeric aspect ratio of the plots. Default: `1`.
#' @param custom Logical; if `TRUE`, use custom polygon‐derived metrics.
#'
#' @return A combined plot (via `cowplot::plot_grid`) showing spatial maps
#'   of each QS term.
#'
#' @importFrom scater plotColData
#' @importFrom ggplot2 ggtitle coord_fixed
#' @importFrom cowplot plot_grid
#' @export
#'
#' @examples
#' example(readAndAddPolygonsToSPE)
#' example(spatialPerCellQC)
#' p <- plotQScoreTerms(spe)
#' print(p)
plotQScoreTerms <- function(spe, sample_id=unique(spe$sample_id), size=0.05,
                            alpha=0.8, aspect_ratio=1, custom = FALSE) {
    if(metadata(spe)$technology=="Nanostring_CosMx") {
        if(custom==TRUE) {
            ggp <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                        y=spatialCoordsNames(spe)[2],
                                        colour_by="cust_log2CountArea",
                                        point_size=size, point_alpha=alpha)+
                ggplot2::ggtitle(sample_id)+ .centroid_image_theme() +
                ggplot2::coord_fixed()
            ggp2 <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                        y=spatialCoordsNames(spe)[2],
                                        colour_by="cust_log2AspectRatio",
                                        point_size=size, point_alpha=alpha)+
                ggplot2::ggtitle(sample_id) + .centroid_image_theme() +
                ggplot2::coord_fixed()
        } else {
            ggp <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                        y=spatialCoordsNames(spe)[2],
                                        colour_by="log2CountArea",
                                        point_size=size, point_alpha=alpha)+
                ggtitle(sample_id)+ .centroid_image_theme() + coord_fixed()
            ggp2 <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                        y=spatialCoordsNames(spe)[2],
                                        colour_by="log2AspectRatio",
                                        point_size=size, point_alpha=alpha)+
                ggplot2::ggtitle(sample_id) + .centroid_image_theme() +
                ggplot2::coord_fixed()
        }
        ggp3 <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                    y=spatialCoordsNames(spe)[2],
                                    colour_by="dist_border",
                                    point_size=size, point_alpha=alpha)+
            ggplot2::ggtitle(sample_id) + .centroid_image_theme() +
            ggplot2::coord_fixed()
        ggp <- cowplot::plot_grid(ggp, ggp2, ggp3, ncol = 2)
    } else {
        ## check if column variable is logical to impose our colors
        ggp <- ggp1 <- scater::plotColData(spe, x=spatialCoordsNames(spe)[1],
                                            y=spatialCoordsNames(spe)[2],
                                            colour_by="log2CountArea",
                                            point_size=size, point_alpha=alpha)+
            ggplot2::ggtitle(sample_id)+ .centroid_image_theme() +
            ggplot2::coord_fixed()
    }
    return(ggp)
}



.assign_outlier_color <- function(vals, outlier_mc, high_label, low_label) {
    thr <- round(attr(outlier_mc, "thresholds"), 2)
    dplyr::case_when(
        vals > thr[2] ~ high_label,
        vals < thr[1] ~ low_label,
        TRUE ~ "unflagged"
    )
}

.assign_collapsed_color <- function(is_ctrl, area_col, dapi_col) {
    color <- rep("unflagged", length(is_ctrl))
    color[is_ctrl] <- "ctrl/total ratio > 0.1"
    color[area_col == "> area um higher thr."] <- "> area um higher thr."
    color[area_col == "< area um lower thr."] <- "< area um lower thr."
    color[dapi_col == "> DAPI higher thr."] <- "> DAPI higher thr."
    color[dapi_col == "< DAPI lower thr."] <- "< DAPI lower thr."
    color
}

.make_outlier_plot <- function(polygons, fov, fillvar, pal, title=NULL,
                            leg=FALSE) {
    p <- ggplot2::ggplot() +
        ggplot2::geom_sf(
            data = polygons[polygons$fov %in% fov, ],
            mapping = ggplot2::aes_string(fill=fillvar, color=fillvar),
            lwd=0, show.legend="polygon"
        ) +
        ggplot2::scale_fill_manual(values=pal) +
        ggplot2::scale_color_manual(values=pal)
    if (!leg) p <- p + ggplot2::theme(legend.position="none")
    if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
    p
}


#' qcFlagPlots
#' @name qcFlagPlots
#' @rdname qcFlagPlots
#' @description
#' Plots the flagged cells identified with first filter, based on control count
#' on total count ratio, area in um and DAPI signal.
#'
#' This function generates a plot that shows selected (FOVs)
#' within a `SpatialExperiment` object, with cells flagged in different colors
#' over a light or dark layout chosen by the user.
#'
#' @param spe A `SpatialExperiment` object containing spatial transcriptomics
#' data.
#' @param fov An integer or numeric vector specifying the FOVs to be plotted
#' Must match values in the `fov` column of `colData(spe)`.
#' @param theme A character string among "light" or "dark".
#' @param custom A boolean value. If TRUE, custom polygons derived metrics will
#' be used.
#'
#' @return A panel with multiple plots showing flagged cells for different
#' variables.
#'
#' @importFrom dplyr case_when
#' @importFrom ggplot2 ggplot geom_sf aes scale_fill_manual scale_color_manual
#' @importFrom ggplot2 ggtitle theme
#' @importFrom cowplot get_legend plot_grid
#' @export
#'
#' @examples
#'
#' example(readAndAddPolygonsToSPE)
#' spe <- spatialPerCellQC(spe)
#' spe <- computeThresholdFlags(spe)
#' p <- qcFlagPlots(spe, fov=16, theme="dark")
#' print(p)
qcFlagPlots <- function(spe, fov=unique(spe$fov),
                            theme=c("light","dark"), custom=FALSE) {
    if (!all(c("is_zero_counts", "is_ctrl_tot_outlier") %in%
            names(colData(spe))))
        stop("Fixed thresholds flag cells not found. Run computeFixedFlags().")
    spe <- computeSpatialOutlier(spe, compute_by="Area_um", method="both")
    spe <- computeSpatialOutlier(spe, compute_by="Mean.DAPI", method="both")
    spe$polygons$fixed_flags_color <- dplyr::case_when(
        spe$is_zero_counts ~ "0 counts",
        spe$is_ctrl_tot_outlier ~ "ctrl/total ratio > 0.1",
        TRUE ~ "unflagged"
    )
    if(all(spe$polygons$fixed_flags_color=="unflagged"))
        warning("No 0 counts or control/total ratio > 0.1 found")
    dapi_col <- .assign_outlier_color(spe$Mean.DAPI, spe$Mean.DAPI_outlier_mc,
                                    "> DAPI higher thr.", "< DAPI lower thr.")
    spe$polygons$dapi_outlier_color <- dapi_col
    area_vals <- if (custom) spe$cust_Area_um else spe$Area_um
    area_mc   <- if (custom) spe$cust_Area_um_outlier_mc else
        spe$Area_um_outlier_mc
    area_col  <- .assign_outlier_color(area_vals, area_mc,
        "> area um higher thr.", "< area um lower thr.")
    spe$polygons$area_outlier_color <- area_col
    spe$polygons$collapsed_color <-
        .assign_collapsed_color(spe$is_ctrl_tot_outlier, area_col, dapi_col)
    outlier_palette <- c(
        "unflagged"="#c0c8cf", "ctrl/total ratio > 0.1"="magenta",
        "< area um lower thr."="darkturquoise", "> area um higher thr."="red",
        "< DAPI lower thr."="purple", "> DAPI higher thr."="greenyellow"
    )
    plot_func <- if(theme[1]=="light") .light_theme else .dark_theme
    ggp1 <- .make_outlier_plot(spe$polygons, fov, "fixed_flags_color",
                            outlier_palette, "Control counts ratio") +
        plot_func()
    ggp2 <- .make_outlier_plot(spe$polygons, fov, "area_outlier_color",
                            outlier_palette, "Area in um") + plot_func()
    ggp3 <- .make_outlier_plot(spe$polygons, fov, "dapi_outlier_color",
                            outlier_palette, "Mean DAPI") + plot_func()
    legp <- .make_outlier_plot(spe$polygons, fov, "collapsed_color",
                            outlier_palette, leg=TRUE) + plot_func() +
        ggplot2::theme(legend.title=element_blank())
    ggp4 <- cowplot::get_legend(legp)
    final <- cowplot::plot_grid(ggp1, ggp2, ggp3, ggp4, ncol=2)
    if(theme[1]=="dark")
        final <- final +
            ggplot2::theme(panel.background=element_rect(fill="black"))
    final
}



