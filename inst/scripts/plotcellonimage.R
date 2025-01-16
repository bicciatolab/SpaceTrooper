spe <- readRDS("/Users/inzirio/Downloads/CosMx_data/CosMx_TNBC_146/spe_fov_146_only.rds")

# Load required libraries
library(ggplot2)
library(magick)
library(dplyr)
library(SpatialExperiment)
# Load the image
image_path <- "/Users/inzirio/Downloads/CosMx_data/CosMx_TNBC_146/CellOverlay_F146.jpg"
image <- image_read(image_path)


# Convert the image to a raster for ggplot2
image_raster <- as.raster(image)

dim(image_raster)
# Load your dataset
# Replace 'your_dataset.csv' with the actual dataset file path
dataset <- as.data.frame(colData(spe))

# Prepare the dataset (filter or clean if necessary)
# Make sure the columns 'CenterX_local_px', 'CenterY_local_px', and 'cell_id' exist
dataset <- dataset %>%
    mutate(CenterY_local_px = -CenterY_local_px)  # Flip Y-axis for plotting

# Plot the image and overlay the cell IDs
ggplot() +
    annotation_raster(image_raster, xmin = 0, xmax = 4256, ymin = -4256, ymax = 0) +
    geom_text(data = dataset, aes(x = CenterX_local_px, y = CenterY_local_px, label = cell_id),
              color = "red", size = 2) +
    theme_void() +
    coord_fixed() +
    ggtitle("Cell Overlay with IDs")


plotCellDoubl <- function(spe, cell_id=NULL, sample_id=unique(spe$sample_id), image_id)
{
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot("neighbors" %in% names(colData(spe)))

    message("it could take some time... ")
    image_raster <- imgRaster(spe, sample_id, image_id)
    data <- as.data.frame(colData(spe))

    data <- data %>%
        mutate(CenterY_local_px = -CenterY_local_px)

    if(!is.null(cell_id))
    {
        ## substitute this with recalling the neighbors from the colPair
        # neighb <- colPair(spe)[from(colPair(spe))==as.integer(cell_id),]
        neighb <- as.integer(strsplit(spe$neighbors[spe$cell_id==cell_id], ",")[[1]])
        txt <- geom_text(
            data = data,
            aes(
                x = CenterX_local_px,
                y = CenterY_local_px,
                label = cell_id,
                color = case_when(
                    cell_id == !!cell_id ~ "yellow", # Correct comparison with function argument
                    cell_id %in% neighb ~ "green",  # Cellule in cidx in verde
                    TRUE ~ "red"  # Tutte le altre in rosso
                )
            ),
            size = 2
        )
    } else {
        txt <- geom_text(data = dataset, aes(x = CenterX_local_px, y = CenterY_local_px, label = cell_id),
                         color = "red", size = 2)
    }
    ggplot() +
        annotation_raster(image_raster, xmin = 0, xmax = 4256, ymin = -4256, ymax = 0) +
        txt +
        scale_color_identity() +  # Usa i colori specificati direttamente
        theme_void() +
        coord_fixed() +
        ggtitle("Cell Overlay with IDs")
}



# generalization
library(sf)
pols <- spe$polygons
idxs <- apply(as.data.frame(colData(spe)), 1, function(r){
    which(st_touches(pols, pols[pols$cell_id==r$cell_id,], sparse=FALSE))
})
idxss<-lapply(idxs, function(x){if(length(x)!=0) {return(x)} else {return(NA)}})


spe$neighbors <- sapply(idxss, function(x) if (is.null(x) || all(is.na(x))) NA else paste(x, collapse = ","))

spe <- addImg(spe, "sample01", "fov146", imageSource="/Users/inzirio/Downloads/CosMx_data/CosMx_TNBC_146/CellOverlay_F146.jpg", scaleFactor=1, load=TRUE)

plotCellNeigh(spe, image_id="fov146")

plotCellNeigh(spe, cell_id=1, image_id="fov146")
