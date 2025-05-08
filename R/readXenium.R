#' readXeniumSPE
#' @name readXeniumSPE
#' @rdname readXeniumSPE
#' @title Load data from a 10x Genomics Xenium experiment
#'
#' @description
#' Creates a [`SpatialExperiment`] from an unzipped Xenium Output Bundle
#' directory containing spatial gene expression data.
#'
#' @param dirname `character(1)`
#'   Path to a Xenium Output Bundle directory.
#' @param sample_name `character(1)`
#'   Sample identifier to assign to `sample_id`. Default: `"sample01"`.
#' @param type `character(1)`
#'   One of `"HDF5"` or `"sparse"`; method to read the feature matrix.
#' @param coord_names `character(2)`
#'   Names of X/Y spatial coordinate columns. Default:
#'   `c("x_centroid", "y_centroid")`.
#' @param boundaries_type `character(1)`
#'   One of `"parquet"` or `"csv"`; format of the polygon file.
#' @param compute_missing_metrics `logical(1)`
#'   If `TRUE`, compute area and aspect‐ratio from boundary polygons.
#' @param keep_polygons `logical(1)`
#'   If `TRUE`, append raw polygon geometries to `colData`.
#' @param countsfilepattern `character(1)`
#'   Pattern to locate the feature matrix file. Default:
#'   `"cell_feature_matrix"`.
#' @param metadatafpattern `character(1)`
#'   Pattern to locate the cell metadata file. Default: `"cells"`.
#' @param polygonsfpattern `character(1)`
#'   Pattern to locate the cell boundaries file. Default:
#'   `"cell_boundaries"`.
#' @param polygonsCol `character(1)`
#'   Name of the polygons column to add to `colData`. Default:
#'   `"polygons"`.
#'
#' @details
#' Expects the unzipped bundle to contain an `outs/` folder with:
#' - `cell_feature_matrix.h5` or `cell_feature_matrix/`
#' - `cells.csv.gz`
#'
#' @return A [`SpatialExperiment`] object with assays, `colData`, spatial
#'   coordinates, and `metadata$polygons` & `metadata$technology`.
#'
#' @author Estella Yixing Dong, Dario Righelli
#'
#' @importFrom DropletUtils read10xCounts
#' @importFrom data.table fread
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SummarizedExperiment assays rowData
#' @export
#' @examples
#' xepath <- system.file(
#'   "extdata", "Xenium_small", package = "SpaceTrooper"
#' )
#' (spe <- readXeniumSPE(
#'   dirname = xepath,
#'   keep_polygons = TRUE
#' ))
readXeniumSPE <- function(dirname,
                          sample_name="sample01",
                          type=c("HDF5", "sparse"),
                          coord_names=c("x_centroid", "y_centroid"),
                          boundaries_type=c("parquet", "csv"),
                          compute_missing_metrics=TRUE, keep_polygons=FALSE,
                          countsfilepattern="cell_feature_matrix",
                          metadatafpattern="cells",
                          polygonsfpattern="cell_boundaries",
                          polygonsCol="polygons",
                          txpattern = "transcripts", add_FOVs = TRUE)
{
    stopifnot(file.exists(dirname))
    type <- match.arg(type)
    boundaries_type <- match.arg(boundaries_type)

    # add "outs/" directory if not already included
    if(basename(dirname) != "outs")
    {
        dirbkup <- dirname
        dirname <- file.path(dirname, "outs")
        if (!file.exists(dirname))
        {
            dirname <- dirbkup
        } else {
            warning("automatically detected/added outs dir in the 10x filepath")
        }
    }

    cfm <- paste0(countsfilepattern, switch(type, HDF5=".h5", ""))
    counts <- file.path(dirname, cfm)

    metadata_file <- file.path(dirname, paste0(metadatafpattern, ".csv.gz"))
    pex <- paste0(polygonsfpattern, switch(boundaries_type,
                                           parquet=".parquet",
                                           csv=".csv.gz"))
    pol_file <- list.files(dirname, pex, full.names=TRUE)
    stopifnot(all(file.exists(c(metadata_file, pol_file))))

    # Count matrix + rowData
    sce <- DropletUtils::read10xCounts(counts, col.names=TRUE)

    # Spatial and colData
    cd <- DataFrame(fread(metadata_file, header=TRUE))
    rownames(cd) <- cd$cell_id

    if ( dim(sce)[2] != dim(cd)[1] )
    {
        sce <- sce[, colnames(sce) %in% rownames(cd)]
        cd <- cd[rownames(cd) %in% colnames(sce), ]
    }
    if (compute_missing_metrics)
    {
        message("Computing missing metrics, this could take some time...")
        cd <- computeMissingMetricsXenium(pol_file, cd, keep_polygons,
                                          polygonsCol)
    }
    if (add_FOVs){
        cd <- addFovFromTx(txpattern, cd)
    }
    # construct 'SpatialExperiment'
    spe <- SpatialExperiment::SpatialExperiment(
        sample_id=sample_name,
        assays = assays(sce),
        rowData = rowData(sce),
        colData = cd,
        spatialCoordsNames = coord_names,
        metadata=list(polygons=pol_file, technology="10X_Xenium")
    )
    return(spe)
}

#' computeMissingMetricsXenium
#' @name computeMissingMetricsXenium
#' @rdname computeMissingMetricsXenium
#'
#' @description
#' Compute Missing Metrics for Xenium Data
#'
#' This function computes missing metrics, such as the aspect ratio, from
#' polygon data in a Xenium dataset and optionally appends the polygon data to
#' the resulting `colData`.
#'
#' @param pol_file A character string specifying the file path to the polygon
#' data.
#' @param coldata A `DataFrame` containing the `colData` for the Xenium dataset.
#' @param keep_polygons A logical value indicating whether to keep the polygon
#' data in the resulting `colData`. Default is `FALSE`.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#'
#' @return A `DataFrame` containing the updated `colData` with computed metrics.
#' If `keep_polygons` is `TRUE`, the polygon data is also included.
#'
#' @details The function reads the polygon data from the specified file,
#' computes the aspect ratio for each polygon, and merges these metrics with
#' the provided `colData`. Optionally, the polygon data can be kept in the
#' returned `colData`.
#'
#' @importFrom S4Vectors cbind.DataFrame
#' @export
#'
#' @examples
#' example(readXeniumSPE)
#' colData(spe) <- computeMissingMetricsXenium(metadata(spe)$polygons,
#'     colData(spe), keep_polygons=TRUE)
computeMissingMetricsXenium <- function(pol_file, coldata, keep_polygons=FALSE,
                                polygonsCol="polygons")
{
    stopifnot(file.exists(pol_file))
    polygons <- readPolygonsXenium(pol_file, keepMultiPol=TRUE)
    cd <- coldata
    cd$AspectRatio <- computeAspectRatioFromPolygons(polygons)
    if(keep_polygons) cd <- .addPolygonsToCD(cd, polygons, polygonsCol)
    return(cd)
}

#' addFovFromTx
#' @name addFovFromTx
#' @rdname addFovFromTx
#'
#' @description
#' Add FOV information from transcript file to cell metadata.
#'
#' This function retrieves FOV information from transcript file and appends
#' the data to the resulting `colData`.
#'
#' @param txpattern A character string specifying the pattern to look for in the
#' directory to find the transcript file, only parquet file is supported.
#' data.
#' @param coldata A `DataFrame` containing the `colData` for the Xenium dataset.
#'
#' @return A `DataFrame` containing the updated `colData` with FOV information.
#'
#' @details The function reads the transcript file then groups it by cell_id
#' and merges the FOV information to the cell metadata in `colData`. Only
#' parquet file is supported for this operation
#'
#' @importFrom S4Vectors cbind.DataFrame
#' @importFrom arrow read_parquet
#' @importFrom dplyr group_by select distinct left_join
#' @export
#'
#' @examples
#' example(readXeniumSPE)
#' colData(spe) <- addFovFromTx("transcripts.parquet", colData(spe))
addFovFromTx <- function(txpattern, coldata){
    tx_file <- file.path(dirname, paste0(txpattern, ".parquet"))
    stopifnot(file.exists(tx_file))
    df <- data.frame(coldata)
    tx <- arrow::read_parquet(tx_file)
    if (!"fov_name" %in% colnames(tx)) {
        stop("No fov_name column was found in tx file. \r\n",
             "Rerun readXeniumSPE without adding FOV information.")}
    g_tx <- group_by(tx, cell_id) |> select(cell_id, fov = fov_name) |>
        distinct(cell_id, .keep_all = TRUE)
    df <- left_join(df, g_tx, by="cell_id")
    cd$fov <- df$fov
    return(cd)
}

