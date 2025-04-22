#' readMerfishSPE
#'
#' @description
#'
#'
#' @param dirname
#' @param sample_name
#' @compute_missing_metrics boolean
#' @boundaries_type one of HDF5, parquet.
#' If HDF5 indicate the polygons folder in `polygonspattern` where the HDF5
#' polygons files are stored
#' @param countmatfpattern
#' @param metadatafpattern
#' @param coord_names
#' @param polygonspattern
#' @param fov_dims
#'
#' @return A SpatialExperiment object
#' @export
#'
#' @importFrom data.table fread
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr left_join
#' @importFrom SpatialExperiment SpatialExperiment
#' @examples
#'
#' path <- system.file(
#'   file.path("extdata", "Merfish_Tiny"),
#'   package = "SpaceTrooper")
#'
#' # read the count matrix .h5 file
#' spe <- readMerfishSPE(dirname = path, keep_polygons=TRUE,
#'     boundaries_type="parquet")
#' spe
readMerfishSPE <- function(dirname,
                           sample_name="sample01",
                           compute_missing_metrics=TRUE, keep_polygons=FALSE,
                           boundaries_type=c("HDF5", "parquet"),
                           countmatfpattern = "cell_by_gene.csv",
                           metadatafpattern = "cell_metadata.csv",
                           polygonsfpattern = "cell_boundaries.parquet",
                           coord_names = c("center_x", "center_y"))
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
    colnames(counts) <- cn

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
                                           keep_polygons)
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

computeMissingMetricsMerfish <- function(pol_file, coldata,
                                        boundaries_type=c("HDF5", "parquet"),
                                        keep_polygons=FALSE)
{
    polygons <- readPolygonsMerfish(pol_file, type=boundaries_type)
    cd <- coldata
    cd$um_area <- computeAreaFromPolygons(polygons)
    cd$AspectRatio <- computeAspectRatioFromPolygons(polygons)

    if (keep_polygons) cd <- cbind.DataFrame(cd, polygons)
    return(cd)
}

