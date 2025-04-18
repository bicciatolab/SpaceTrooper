#' readPolygons
#' @description
#'  Read and Validate Polygons from a File
#'
#' This function reads polygon data from a specified file, validates the polygons,
#' and returns them as an `sf` object. It supports multiple file formats and can
#' handle both global and local coordinates.
#'
#' @param polygonsFile A character string specifying the path to the polygon
#' file.
#' @param type A character string specifying the file type. Supported types are
#' `"csv"`, `"parquet"`, and `"h5"`. Default is `"csv"`.
#' @param x A character vector specifying the column names for the x-coordinates
#' in the polygon data. Default is `c("x_global_px", "vertex_x")`.
#' @param y A character vector specifying the column names for the y-coordinates
#' in the polygon data. Default is `c("y_global_px", "vertex_y")`.
#' @param xloc A character string specifying the column name for the local
#' x-coordinates. Default is `"x_local_px"`.
#' @param yloc A character string specifying the column name for the local
#' y-coordinates. Default is `"y_local_px"`.
#' @param keepMultiPol A logical value indicating whether to keep multipolygons
#' during validation. Default is `TRUE`.
#' @param verbose A logical value indicating whether to print additional
#' information during processing. Default is `FALSE`.
#'
#' @return An `sf` object with the loaded and validated polygons.
#'
#' @details The function reads polygon data from the specified file and formats.
#' It validates the polygons and handles both global and local coordinates if
#' provided. If the file type is `"h5"`, the function currently does not handle
#' the data, as this part of the code is not implemented.
#'
#' @importFrom data.table fread
#' @importFrom arrow read_parquet
#' @importFrom sf st_geometry
#' @export
#'
#' @examples
#' # Reading polygon data from a CSV file:
#' # polygons <- readPolygons("~/Downloads/CosMx_data/polygons.csv")
#'
#' # Reading polygon data from a Parquet file with verbose output:
#' # polygons <- readPolygons("~/Downloads/CosMx_data/polygons.parquet",
#' #                         type = "parquet", verbose = TRUE)
readPolygons <- function(polygonsFile, type=c("csv", "parquet", "h5"),
                            x=c("x_global_px", "vertex_x"),
                            y=c("y_global_px", "vertex_y"),
                            xloc="x_local_px", yloc="y_local_px",
                            #micronConvFact=0.12,
                            keepMultiPol=TRUE,
                            verbose=FALSE)
{

    stopifnot(file.exists(polygonsFile))
    type <- match.arg(type)
    # type <- grep("csv", polygonsFile)

    if(type=="h5")
    {
        # polfiles <- list.files(polygonsFolder, pattern=hdf5pattern,
        #                        full.names=TRUE)
        # dfsfl <- lapply(seq_along(polfiles), function(i)
        # {
        #     poll <- readh5polygons(pol_file=polfiles[i])
        #     df <- data.frame(cell_id=paste0("f", i-1, "_c", poll$ids),
        #                      cell_ID=poll$ids,
        #                      fov=i-1, geometry=sf::st_sfc(poll$g))
        #     dfsf <- sf::st_sf(df)
        # })
        # polygons <- do.call(rbind, dfsfl)
    } else {
        spat_obj <- switch(type, csv=fread(polygonsFile),
                                parquet=read_parquet(polygonsFile))
        if (! "cell_id" %in% colnames(spat_obj))
        {
            spat_obj$cell_id <- paste0("f", spat_obj$fov, "_c", spat_obj$cellID)
        }

        spat_obj$cell_id <- as.factor(spat_obj$cell_id)

        polygons <- .createPolygons(spat_obj, x=x, y=y,
                                        polygon_id="cell_id")
        polygons <- .renameGeometry(polygons, "geometry", "global") ## only for cosmx

        if(all(c(xloc, yloc) %in% colnames(spat_obj)))
        {
            idxs <- which(colnames(polygons) %in% c(xloc, yloc))
            polygons <- polygons[,-idxs]
            polygons_loc <- .createPolygons(spat_obj, x=xloc, y=yloc,
                                            polygon_id="cell_id")[,-c(4,5)]
            polygons <- cbind(polygons, polygons_loc$geometry)
            polygons <- .renameGeometry(polygons, "geometry", "local")
        }

        if(verbose) message("Polygons detected: ", dim(polygons)[1])#### otherwise
        #### will print number of columns next to the numer of rows

        polygons <- .checkPolygonsValidity(polygons, keepMultiPol=keepMultiPol,
                                            verbose=verbose)

        rownames(polygons) <- polygons$cell_id #### if not here, then the check
        #### in addPolygonsToSPE cannot be be done

        #### identical cannot work since coordinates are different, but the validity
        #### and the geometries as well should stay the same
        #### needs to be changed or ignored
        # if(!table(st_is_valid(polygons$global))==table(st_is_valid(polygons$global)))
        #     warning("Global and Local geometries are not identical")
        if(verbose) message("Polygons after validity: ", dim(polygons)[1])#### otherwise
        #### will print number of columns next to the numer of rows
        return(polygons)
    }
    ### write polygons as parquet file
}


#' .checkPolygonsValidity
#' @description checks validity on a geometry of `sf` object.
#' It removes multipolygons when `keepMultiPol` is `FALSE`
#'
#' @details
#' In case geometry is NULL validity is checked on the active geometry,
#' otherwise it is checked on the passed geometry without changing the active
#' geometry of the sf object.
#' In case of not valid polygons, these are removed.
#' If keeMultiPol is FALSE, possible detected multipolygons are removed.
#'
#'
#' @param sf An `sf` class object containing the spatial data.
#' @param geometry character for the geometry to check validity, if `NULL`
#' it checks the active geometry (default is `NULL`)
#' @param keepMultiPol logical for keeping/removing moltipolygons, if any
#' (default is `TRUE`, so keeping the multipolygons)
#' @param verbose logical to print verbose output (default is `FALSE`)
#'
#' @return An `sf` object with valid geometries, possibly with multipolygons
#' removed.
#' @keywords internal
#' @importFrom sf st_is_valid st_buffer st_geometry_type
#'
#' @examples
#' # TBD
.checkPolygonsValidity <- function(sf, geometry=NULL, keepMultiPol=TRUE,
                                    verbose=FALSE)
{
    stopifnot(is(sf, "sf"))
    if(is.null(geometry))
    {
        geometry <- .getActiveGeometryName(sf)
    } else {
        act <- .getActiveGeometryName(sf)
        sf <- .setActiveGeometry(sf, geometry)
    }

    sf_tf <- st_is_valid(sf) # to parallelize? how? split sf in multiple sf and parallelize on it?

    if(sum(sf_tf)!=dim(sf)[1]) sf <- st_buffer(sf, dist=0)

    ############ CHECKING MULTIPOLYGONS
    # Subsetting to remove non-polygons
    # cellids <- unlist(apply(sf, 1, function(geom)
    # {

        # just in case of custom in cosmx but present in automatic in xenium
        # if(attr(geom[[geometry]], "class")[2] == "MULTIPOLYGON")
        # {
        #     return(geom$cell_id)
        # }
    # }))
    is_multi <- sf::st_geometry_type(sf$global) == "MULTIPOLYGON"
    ## TRUE is merscope - FALSE is cosmx and xenium
    merscopeFl <- (sum(is_multi) == dim(sf)[1])
    funct <- ifelse(merscopeFl, "lengths", "length")
    ll <- lapply(sf[[geometry]], funct)
    sums <- lapply(ll, sum)
    idx <- which(sums!=1)
    sf$is_multi <- FALSE
    sf$multi_n <- 1
    if(length(idx)!=0)
    {
        sf$is_multi[idx] <- TRUE
        sf$multi_n[idx] <- unlist(sums[idx])

    }
    if( all(merscopeFl, length(idx)!=0) )
    {
        sf$global[-idx] <- st_cast(sf$global[-idx], "POLYGON")
    }

    if(verbose) message("Detected ", sum(sf$is_multi), " multipolygons.")

    if(!keepMultiPol)
    {
        if(verbose) message("Removing ", sum(sf$is_multi), " multipolygons.")
        sf <- sf[!sf$is_multi,]
    }

    if(exists("act")) sf <- .setActiveGeometry(sf, act)

    return(sf)
}
#' readAndAddPolygonsToSPE
#' @description Read and Add Polygons to a SpatialExperiment Object
#'
#' This function reads polygon boundary data based on the technology associated
#' with the provided SpatialExperiment (SPE) object and adds the polygons to the
#' SPE.
#'
#' @param spe A \code{SpatialExperiment} object. The object should contain
#' metadata with the field \code{"technology"}, specifying the technology used
#' (e.g., "Nanostring_CosMx", "Vizgen_MERFISH", or "10X_Xenium").
#' @param keepMultiPol Logical. If \code{TRUE}, multi-polygon features will be
#' kept when reading the boundary data. Defaults to \code{TRUE}.
#' @param boundaries_type Character. Specifies the type of boundary file format
#' to read. Options are \code{"HDF5"} or \code{"parquet"}. Defaults to
#' \code{"HDF5"}.
#'
#' @return A \code{SpatialExperiment} object with the added polygon data.
#'
#' @export
#'
#' @examples
#' # TBD
readAndAddPolygonsToSPE <- function(spe, keepMultiPol=TRUE,
                    boundaries_type=c("HDF5", "parquet"))
{
    boundaries_type<-match.arg(boundaries_type)
    stopifnot("technology" %in% names(metadata(spe)))
    tech <- metadata(spe)$technology
    if(is.null(polygons))
    {
        switch(tech,
               "Nanostring_CosMx"=
                   {
                       polygons <- readPolygonsCosmx(metadata(spe)$polygons)
                   },
               "Vizgen_MERFISH"=
                   {
                       ##### NEED TO HANDLE THE DIFFERENCES BETWEEN HDF5 FILES
                       ##### AND PARQUET, TO PROPAGATE TO READING FUNCTION
                       # ifelse(boundaries_type=="HDF5", merpol=)
                       polygons <- readPolygonsMerfish(polygonsFolder,
                                        keepMultiPol=TRUE, type=boundaries_type)

                   },
               "10X_Xenium"=
                   {
                       polygons <- readPolygonsXenium(pol_file, keepMultiPol=TRUE)
                   },
               stop("Unrecognized technology, please use an SPE from one of ",
                    "Nanostring_CosMx, Vizgen_MERFISH and 10x_Xenium")
        )
    }
    spe <- addPolygonsToSPE(spe, polygons)
    return(spe)
}


#
# addPolygonsToSPE <- function(spe, polygons=NULL)
# {
#     stopifnot(all(is(spe, "SpatialExperiment"), is(polygons, "sf")))
#
#     if (sum(rownames(polygons) == colnames(spe)) != dim(spe)[2])
#     {
#         cd <- data.frame(colData(spe))
#         # polygons <- left_join(polygons, cd[, c("fov", "cellID")],
#         #                       by=c("fov", "cellID"))
#         # cd <- left_join(cd, polygons[ , c("fov", "cellID")],
#         #                 by=c("fov","cellID"))
#         # polygons$cell_id <- paste0("f", polygons$fov, "_c", polygons$cellID)
#         # rownames(cd) <- cd$cell_id
#         # rownames(polygons) <- polygons$cell_id
#         # polygons <- polygons[polygons$cell_id %in% rownames(cd),]
#         cd <- .addPolygonsToCD(cd, polygons)
#         spe <- spe[, spe$cell_id %in% rownames(cd$polygons)]
#     }
#     spe <- spe[, rownames(polygons)]
#     colData(spe)$polygons <- polygons
#     return(spe)
# }

# .addPolygonsToCD <- function(cd, polygons)
# {
#     stopifnot(all(is(cd, "DataFrame"), is(polygons, "sf")))
#     polygons <- left_join(polygons, cd[, c("fov", "cellID")],
#                           by=c("fov", "cellID"))
#     cd <- left_join(cd, polygons[ , c("fov", "cellID")],
#                     by=c("fov","cellID"))
#     polygons$cell_id <- paste0("f", polygons$fov, "_c", polygons$cellID)
#     rownames(cd) <- cd$cell_id
#     rownames(polygons) <- polygons$cell_id
#     cd$polygons <- polygons
#     return(cd)
# }

#' Attach sf polygons to a DataFrame of cell metadata
#'
#' This function enriches a DataFrame (e.g., from colData) with matching polygon geometries.
#'
#' @param cd A DataFrame containing at least `fov` and `cellID` columns.
#' @param polygons An sf object with matching `fov` and `cellID` columns.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#' @return A DataFrame identical to `cd`, but row‑subset to cells present in `polygons` and with a new `polygons` list‑column of sf geometries.
#' @keywords internal
#' @examples
#' # cd2 <- addPolygonsToCD(cd, polygons_sf)
.addPolygonsToCD <- function(cd, polygons, polygonsCol="polygons") {
    stopifnot(inherits(cd, "DataFrame"), inherits(polygons, "sf"))

    if("fov"%in%colnames(cd))
    {
        # Build consistent cell IDs
        cd$cell_id <- paste0("f", cd$fov, "_c", cd$cellID)
        polygons$cell_id <- paste0("f", polygons$fov, "_c", polygons$cellID)
        rownames(cd) <- cd$cell_id
        rownames(polygons) <- polygons$cell_id
    } else {
        rownames(cd) <- cd$cell_id
        rownames(polygons) <- polygons$cell_id
    }

    # Identify matching cells
    common_ids <- intersect(rownames(cd), rownames(polygons))
    if (length(common_ids) == 0) {
        stop("No matching cell IDs between 'cd' and 'polygons'.")
    }
    if (length(common_ids) < nrow(cd)) {
        warning(sprintf(
            "%d/%d cells have polygon data; subsetting to matches.",
            length(common_ids), nrow(cd)
        ))
    }

    # Subset and attach
    cd <- cd[common_ids, , drop = FALSE]
    polygons <- polygons[common_ids, , drop = FALSE]
    cd[[polygonsCol]] <- polygons

    return(cd)
}


#' addPolygonsToSPE
#'
#' @description This function adds polygon data to a `SpatialExperiment` object.
#'
#' @param spe A `SpatialExperiment` object to which polygons will be added.
#' @param polygons An `sf` object containing the polygon data.
#' @param polygonsCol character indicating the name of the polygons column to
#' add into the colData (default is `polygons`).
#' @return The `SpatialExperiment` object with polygons added to the `colData`.
#'
#' @importFrom dplyr left_join
#' @export
#'
#' @examples
#' # spe <- addPolygonsToSPE(spe, polygons)
addPolygonsToSPE <- function(spe, polygons, polygonsCol="polygons") {
    stopifnot(inherits(spe, "SpatialExperiment"), inherits(polygons, "sf"))

    # Extract and enrich colData via DataFrame
    cd <- S4Vectors::DataFrame(colData(spe))
    cd <- .addPolygonsToCD(cd, polygons, polygonsCol)

    # Subset and reorder the SpatialExperiment
    spe <- spe[, cd$cell_id, drop = FALSE]

    # Replace colData with enriched DataFrame
    colData(spe) <- cd

    return(spe)
}


#' .createPolygons
#'
#' @description This internal function creates polygons from a spatial data object.
#'
#' @param spat_obj A data frame or similar object containing spatial data.
#' @param x A character vector specifying the x-coordinates.
#' @param y A character vector specifying the y-coordinates.
#' @param polygon_id A character string specifying the polygon ID.
#'
#' @return An `sf` object containing the created polygons.
#' @keywords internal
#' @importFrom sfheaders sf_polygon
#' @importFrom sf st_as_sf
#'
#' @examples
#' #TBD
.createPolygons <- function(spat_obj, x=NULL, y=NULL, polygon_id=NULL, geometry="Geometry")
{
    if(all(!is.null(x), !is.null(y)))
    {
        polygons <- sfheaders::sf_polygon(obj=spat_obj,
                                      x=x, y=y,
                                      polygon_id=polygon_id, keep=TRUE)
    } else {
        polygons <- sf::st_as_sf(as.data.frame(spat_obj))
        ## check structure of polygons
        # ll <- lapply(polygons$Geometry, lengths)
        # sums <- lapply(ll, sum)
        # idx <- which(sums==1)
        # polygons$is_multi <- FALSE
        # polygons$is_multi[!idx] <- TRUE
        # polygons$multi_n <- 1
        # polygons$multi_n[!idx] <- sums[!idx]
        # polygons[idx, ] <- sf::st_cast(polygons[idx, ], "POLYGON")
    }
    # polygons <- polygons[order(polygons$cell_id),]
    return(polygons)
}
#' readPolygonsCosmx
#'
#' @description This function reads polygon data specific to CosMx technology.
#'
#' @param polygonsFile A character string specifying the file path to the
#' polygon data.
#' @param type A character string specifying the file type ("csv" or "parquet").
#' @param x A character string specifying the x-coordinate column.
#' @param y A character string specifying the y-coordinate column.
#' @param xloc A character string specifying the local x-coordinate column.
#' @param yloc A character string specifying the local y-coordinate column.
#' @param keepMultiPol A logical value indicating whether to keep multipolygons.
#' @param verbose A logical value indicating whether to print additional
#' information.
#'
#' @return An `sf` object containing the CosMx polygon data.
#' @export
#'
#' @examples
#' # Read CosMx polygon data from a CSV file:
#' # polygons <- readPolygonsCosmx("path/to/polygons.csv", type="csv")
readPolygonsCosmx <- function(polygonsFile, type=c("csv", "parquet"),
                              x="x_global_px",
                              y="y_global_px",
                              xloc="x_local_px",
                              yloc="y_local_px",
                              keepMultiPol=TRUE,
                              verbose=FALSE)
{
    type=match.arg(type)
    polygons <- readPolygons(polygonsFile, type=type, x=x, y=y, xloc=xloc,
                    yloc=yloc, keepMultiPol=keepMultiPol, verbose=verbose)
    polygons <- st_cast(polygons, "GEOMETRY")
    mandatory <- c("cell_id", "global", "is_multi", "multi_n")
    cnames <- colnames(polygons)[!colnames(polygons) %in% mandatory]
    polygons <- polygons[,c(mandatory, cnames)]
    polygons <- st_cast(polygons, "GEOMETRY")
    return(polygons)
}

#' readPolygonsXenium
#'
#' @description This function reads polygon data specific to Xenium technology.
#'
#' @param polygonsFile A character string specifying the file path to the
#' polygon data.
#' @param type A character string specifying the file type ("parquet" or "csv").
#' @param x A character string specifying the x-coordinate column.
#' @param y A character string specifying the y-coordinate column.
#' @param keepMultiPol A logical value indicating whether to keep multipolygons.
#' @param verbose A logical value indicating whether to print additional
#' information.
#'
#' @return An `sf` object containing the Xenium polygon data.
#' @export
#'
#' @examples
#' # Read Xenium polygon data from a Parquet file:
#' # polygons <- readPolygonsXenium("path/to/polygons.parquet", type="parquet")
readPolygonsXenium <- function(polygonsFile, type=c("parquet", "csv"),
                   x="vertex_x", y="vertex_y", keepMultiPol=TRUE,
                   verbose=FALSE)
{
    type <- match.arg(type)
    polygons <- readPolygons(polygonsFile=polygonsFile, type=type, x=x, y=y,
        keepMultiPol=keepMultiPol, verbose=verbose)
    mandatory <- c("cell_id", "global", "is_multi", "multi_n")
    cnames <- colnames(polygons)[!colnames(polygons) %in% mandatory]
    polygons <- polygons[,c(mandatory, cnames)]
    return(polygons)
}

#' readPolygonsMerfish
#'
#' @description This function reads polygon data specific to MERFISH technology.
#'
#' @param polygonsFolder A character string specifying the folder containing the
#' polygon data files.
#' @param type A character string specifying the file type ("HDF5" or "parquet").
#' @param hdf5pattern A character string specifying the pattern to match HDF5
#' files.
#' @param keepMultiPol A logical value indicating whether to keep multipolygons.
#' @param z_lev An integer specifying the Z level to filter the data. Default is
#' `3L`.
#' @param zcolumn A character string specifying the column name for the Z index.
#' @param geometry A character string specifying the geometry column name.
#' @param verbose A logical value indicating whether to print additional
#' information.
#'
#' @return An `sf` object containing the MERFISH polygon data.
#' @export
#' @importFrom sf st_sf st_cast st_sfc
#' @importFrom arrow read_parquet
#'
#' @examples
#' # Read MERFISH polygon data from a Parquet file:
#' # polygons <- readPolygonsMerfish("path/to/polygon_folder", type="parquet")
readPolygonsMerfish <- function(polygonsFolder, type=c("HDF5", "parquet"),
                                keepMultiPol=TRUE, hdf5pattern="hdf5",
                                z_lev=3L, zcolumn="ZIndex",
                                geometry="Geometry",
                                verbose=FALSE)
{
    type <- match.arg(type)
    if (type=="HDF5")
    {
        polfiles <- list.files(polygonsFolder, pattern=hdf5pattern,
                                full.names=TRUE)
        dfsfl <- lapply(seq_along(polfiles), function(i)
        {
            poll <- readh5polygons(pol_file=polfiles[i])
            df <- data.frame(cell_id=paste0("f", i-1, "_c", poll$ids),
                             cell_ID=poll$ids,
                             fov=i-1, geometry=sf::st_sfc(poll$g)) ## geometry can be a simple column
            dfsf <- sf::st_sf(df)
        })
        polygons <- do.call(rbind, dfsfl)
        ## check if polygons geometry are polygons
        polygons <- .checkPolygonsValidity(polygons, keepMultiPol=keepMultiPol)
        return(polygons)
    } else { ## case parquet
        polfile <- list.files(polygonsFolder, pattern=type,
                              full.names=TRUE)
        polygons <- arrow::read_parquet(polfile, as_data_frame=TRUE)
        polygons <- polygons[polygons[[zcolumn]]==z_lev,]
        polygons$cell_id <- polygons$EntityID
        polygons <- .createPolygons(polygons)
        polygons <- .renameGeometry(polygons, geometry, "global")
        polygons <- .checkPolygonsValidity(polygons, keepMultiPol=keepMultiPol,
                                        verbose=verbose)

    }
    mandatory <- c("cell_id", "global", "is_multi", "multi_n")
    cnames <- colnames(polygons)[!colnames(polygons) %in% mandatory]
    polygons <- polygons[,c(mandatory, cnames)]
    return(polygons)
}
#' computeCenterFromPolygons
#'
#' @description This function computes the center coordinates on x and y axis
#' from polygon data and adds it to the `colData`. It is necessary only for Merfish.
#'
#' @param polygons An `sf` object containing polygon data.
#' @param coldata A `DataFrame` containing the `colData` to which center coordinates information will be added.
#'
#' @return A `DataFrame` with the added area information.
#' @export
#'
#' @examples
#' # Assuming `polygons` is an sf object and `coldata` is a DataFrame:
#' # coldata <- computeCenterFromPolygons(polygons, coldata)
computeCenterFromPolygons <- function(polygons, coldata)
{
    cd <- coldata
    # cd$Area <- NA
    centroid <- sf::st_centroid(polygons)
    center_x <- lapply(centroid$geometry, function(x) ## get active name of geometry instead of global
    {
        x[1]
    })
    center_y <- lapply(centroid$geometry, function(x) ## get active name of geometry instead of global
    {
        x[2]
    })
    cd$center_x <- unlist(center_x)
    cd$center_y <- unlist(center_y)
    return(cd)
}


#' computeAreaFromPolygons
#'
#' @description This function computes the area from polygon data.
#'
#' @param polygons An `sf` object containing polygon data.
#'
#' @return A `numeric` vector with the area information.
#' @export
#'
#' @examples
#' # Assuming `polygons` is an sf object:
#' # coldata <- computeAreaFromPolygons(polygons, coldata)
computeAreaFromPolygons <- function(polygons)
{
    stopifnot(is(polygons, "sf"))
    area <- sf::st_area(polygons)
    um_area <- unlist(area)
    return(um_area)
}


#' computeAspectRatioFromPolygons
#'
#' @description This function computes the aspect ratio from polygon data.
#'
#' @param polygons An `sf` object containing polygon data.
#'
#' @return A `numeric` vector with the aspect ratio information.
#' @export
#'
#' @examples
#' # Assuming `polygons` is an sf object:
#' # ar <- computeAspectRatioFromPolygons(polygons)
computeAspectRatioFromPolygons <- function(polygons)
{
    stopifnot(all(is(polygons, "sf"), ("is_multi" %in% colnames(polygons))))
    aspRatL <- numeric(nrow(polygons))
    if(any(polygons$is_multi)) {
        aspRatL[which(polygons$is_multi)] <- NA
        warning("Found ", sum(polygons$is_multi), " multi-poligons: returning NA aspect ratio for them.")
    }
    aspRatL[!polygons$is_multi] <- lapply(polygons$global[!polygons$is_multi], function(x) {
        (max(x[[1]][, 2]) - min(x[[1]][, 2]))/(max(x[[1]][, 1]) -
                                                   min(x[[1]][, 1]))
    })
    names(aspRatL) <- polygons$cell_id
    ar <- unlist(aspRatL)
    return(ar)
}

#' readh5polygons
#'
#' @description This function reads polygon data from an HDF5 file.
#'
#' @param pol_file A character string specifying the file path to the HDF5 polygon data.
#'
#' @return A list containing the polygon geometries and their associated cell IDs.
#' @author Lambda Moses
#' @export
#'
#' @examples
#' # Read polygons from an HDF5 file:
#' # polygons <- readh5polygons("path/to/polygons.h5")
readh5polygons <- function(pol_file)
{
    l <- rhdf5::h5dump(pol_file)[[1]]
    cell_ids <- names(l)
    geometries <- lapply(l, function(m) {
        sf::st_polygon(list(t(m[["zIndex_0"]]$p_0$coordinates[,,1])))
    })
    return(list(g=geometries, ids=cell_ids))
}

#' customPolyMetrics
#'
#' @description This function computes centroids, area, area in um, aspect ratio
#' and logged target counts on area ratio for custom polygons. New variables have
#' a prefix cust_ to distinguish them from the previous variables.
#'
#' @param spe A `SpatialExperiment` object with custom polygons stored inside after
#' running addPolygonsToSpe.
#'
#' @return A `SpatialExperiment` with new metrics
#' @author
#' @export
#'
#' @examples
#' #TBD
customPolyMetrics <- function(spe = spe){
    st_geometry(spe$polygons) <- "global"
    centroid_sf <- st_centroid(spe$polygons)
    spe$cust_CenterX_global_px <- unlist(
        lapply(centroid_sf$global, function(x)x[1]))
    spe$cust_CenterY_global_px <- unlist(
        lapply(centroid_sf$global, function(x)x[2]))

    spatialCoordsNames(spe)[1] <- "cust_CenterX_global_px"
    spatialCoordsNames(spe)[2] <- "cust_CenterY_global_px"

    spe$cust_Area <- st_area(spe$polygons)
    spe$cust_Area_um <- st_area(spe$polygons)*(0.12^2)

    custom_Aspect_ratio <- lapply(spe$polygons$global[!spe$polygons$is_multi],
                                  function(x){
        (max(x[[1]][,1]) - min(x[[1]][,1]))/(max(x[[1]][,2]) - min(x[[1]][,2]))
    })

    names(custom_Aspect_ratio) <- spe$polygons$cell_id[!spe$polygons$is_multi]

    spe$cust_AspectRatio <- NA
    posz <- match(names(custom_Aspect_ratio), spe$cell_id)
    spe$cust_AspectRatio[posz] <- unlist(custom_Aspect_ratio)

    spe$cust_log2AspectRatio <- log2(spe$cust_AspectRatio)

    spe$cust_log2CountArea <- log2(spe$target_sum/spe$cust_Area_um)
    return(spe)
}
