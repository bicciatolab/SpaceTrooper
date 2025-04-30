#' readMerfishSPE
#' @name readMerfishSPE
#' @rdname readMerfishSPE
#' @description
#' `readMerfishSPE()` imports MERFISH outputs (counts, metadata, and optionally
#' cell boundary polygons) from a directory and builds a SpatialExperiment
#' object.
#'
#' @param dirname `character(1)`
#'   Path to a folder containing MERFISH output files.
#' @param sample_name `character(1)`
#'   Identifier to assign to the `sample_id` field in the returned object.
#'   Default: `"sample01"`.
#' @param compute_missing_metrics `logical(1)`
#'   If `TRUE`, compute area and aspect‐ratio metrics from the cell boundary
#'   polygons. Default: `TRUE`.
#' @param keep_polygons `logical(1)`
#'   If `TRUE`, attach raw polygon geometries as extra columns in `colData`.
#'   Default: `FALSE`.
#' @param boundaries_type `character(1)`
#'   One of `"HDF5"` or `"parquet"`. If `"HDF5"`, uses a folder of HDF5
#'   polygon files; if `"parquet"`, reads a single Parquet file of boundaries.
#' @param countmatfpattern `character(1)`
#'   Regex passed to `list.files()` to find the count matrix CSV. Default:
#'   `"cell_by_gene.csv"`.
#' @param metadatafpattern `character(1)`
#'   Pattern to find the cell metadata CSV. Default: `"cell_metadata.csv"`.
#' @param polygonsfpattern `character(1)`
#'   Pattern to find the cell boundaries file. Default:
#'   `"cell_boundaries.parquet"`.
#' @param coord_names `character(2)`
#'   Names of the columns in `colData` that store X/Y spatial coordinates.
#'   Default: `c("center_x", "center_y")`.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#'
#' @return A `SpatialExperiment` object with:
#'   - `assays$counts`: gene × cell count matrix
#'   - `colData`: per‐cell metadata (including computed metrics)
#'   - spatial coordinates named by `coord_names`
#'   - `metadata$polygons`: path to the boundaries file
#'   - `metadata$technology`: `"Vizgen_MERFISH"`.
#' @author Estella Yixing Dong, Dario Righelli
#' @export
#' @importFrom data.table fread
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr left_join
#' @importFrom SpatialExperiment SpatialExperiment
#' @examples
#' path <- system.file("extdata", "Merfish_Tiny",
#'                     package = "SpaceTrooper")
#' spe <- readMerfishSPE(
#'   dirname = path,
#'   sample_name = "Patient2",
#'   keep_polygons = TRUE,
#'   boundaries_type = "parquet"
#' )
#' spe
readMerfishSPE <- function(dirname,
                           sample_name="sample01",
                           compute_missing_metrics=TRUE, keep_polygons=FALSE,
                           boundaries_type=c("HDF5", "parquet"),
                           countmatfpattern = "cell_by_gene.csv",
                           metadatafpattern = "cell_metadata.csv",
                           polygonsfpattern = "cell_boundaries.parquet",
                           coord_names = c("center_x", "center_y"),
                           polygonsCol="polygons")
{
    countmat_file <- list.files(dirname, countmatfpattern, full.names=TRUE)
    metadata_file <- list.files(dirname, metadatafpattern, full.names=TRUE)
    pol_file <- list.files(dirname, polygonsfpattern, full.names=TRUE) #check if parquet

    # stopifnot(all(file.exists(countmat_file), file.exists(metadata_file),
    #               file.exists(fovpos_file), file.exists(pol_file)))

    # Read in
    countmat <- data.table::fread(countmat_file)
    names(countmat)[names(countmat) %in% c("V1", "cell")] <- "cell_id"

    metadata <- data.table::fread(metadata_file) # cell metadata
    names(metadata)[names(metadata) %in% c("EntityID", "V1")] <- "cell_id"

    # Count matrix
    countmat <- left_join(countmat, metadata[, "cell_id"], by = "cell_id")
    cn <- countmat$cell_id
    counts <- subset(countmat, select = -cell_id)

    features <- colnames(counts)
    counts <- t(as.matrix(counts))
    rownames(counts) <- features
    colnames(counts) <- as.character(cn)

    # rowData (does not exist)
    # To be associated to the tx file
    # use readSparseCSV sparseArray from harve pege

    # colData
    colData <- left_join(metadata, countmat[, "cell_id"], by = "cell_id")
    rownames(colData) <- colData$cell_id
    cd <- subset(colData, select = c(2,1,3:dim(colData)[2]))
    if (compute_missing_metrics)
    {
        message("Computing missing metrics, this could take a while...")
        cd <- computeMissingMetricsMerfish(pol_file, colData, boundaries_type,
                                           keep_polygons, polygonsCol)
    }

    spe <- SpatialExperiment::SpatialExperiment(
        sample_id=sample_name,
        assays = list(counts = counts),
        # rowData = rowData,
        colData = cd, # this must be cd, not colData, otherwise it will not store the computed
        # missing metrics
        spatialCoordsNames = coord_names,
        metadata=list(polygons=pol_file, technology="Vizgen_MERFISH")
    )
    colnames(spe) <- spe$cell_id
    return(spe)

}

#' computeMissingMetricsMerfish
#' @name computeMissingMetricsMerfish
#' @rdname computeMissingMetricsMerfish
#' @description
#' `computeMissingMetricsMerfish()` takes cell metadata and boundary
#' polygons, calculates per‐cell area and aspect‐ratio, and optionally
#' appends the raw polygon geometries.
#'
#' @param pol_file `character`
#' Path (or vector of paths) to polygon files (HDF5 or Parquet).
#' @param coldata `DataFrame` or `data.frame`
#' Cell metadata with at least a `cell_id` column.
#' @param boundaries_type `character(1)`
#'   One of `"HDF5"` or `"parquet"`—passed on to
#'   `readPolygonsMerfish()`.
#' @param keep_polygons `logical(1)`
#'   If `TRUE`, cbinds the raw polygon `sf` columns onto `coldata`.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#'
#' @return A `DataFrame` (or `data.frame`) with:
#'   - all columns of `coldata`
#'   - `um_area`: area of each cell’s polygon
#'   - `AspectRatio`: width/height aspect ratio
#'   - (optionally) the polygon geometries
#'
#' @export
#' @importFrom S4Vectors DataFrame
#' @examples
#' example(readMerfishSPE)
#' spe <- computeMissingMetricsMerfish(spe)
computeMissingMetricsMerfish <- function(pol_file, coldata,
                                        boundaries_type=c("HDF5", "parquet"),
                                        keep_polygons=FALSE,
                                        polygonsCol="polygons")
{
    polygons <- readPolygonsMerfish(pol_file, type=boundaries_type)
    cd <- DataFrame(coldata)
    cd$um_area <- computeAreaFromPolygons(polygons)
    cd$AspectRatio <- computeAspectRatioFromPolygons(polygons)

    if (keep_polygons) cd <- .addPolygonsToCD(cd, polygons, polygonsCol)
    return(cd)
}

