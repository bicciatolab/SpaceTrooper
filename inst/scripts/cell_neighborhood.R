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

#
# library(sf)
# pols <- spe$polygons
# pols <- pols[pols$cell_id==227,]
# # cell_227_polygon <- pols[pols$cell_id == 227, ]  # Filtra la cellula 227
#
# ggplot() +
#     annotation_raster(image_raster, xmin = 0, xmax = 4256, ymin = -4256, ymax = 0) +
#     geom_text(
#         data = dataset,
#         aes(
#             x = CenterX_local_px,
#             y = CenterY_local_px,
#             label = cell_id,
#             color = factor(cell_id == 227)  # Colore condizionale
#         ),
#         size = 2
#     ) +
#     scale_color_manual(
#         values = c("FALSE" = "red", "TRUE" = "yellow"),  # Associa colori condizionali
#         guide = "none"  # Nasconde la legenda dei colori
#     ) +
#     theme_void() +
#     coord_fixed() +
#     ggtitle("Cell Overlay with IDs")
#

# library(sf)
# t <- st_touches(spe$polygons, pols, sparse=FALSE)
# cidx <- which(t)

plotCellNeigh <- function(spe, cell_id=NULL, sample_id=unique(spe$sample_id), image_id)
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

################## fin qui tutto ok

# constructing cell graph
library(igraph)

names(idxss) <- lapply(strsplit(names(idxss), "_"), function(x){ x[2]})
idxssna <- lapply(idxss, function(x) {
    if (all(is.na(x))) {
        character(0)  # Replace NA with an empty character vector
    } else {
        as.character(na.omit(x))  # Remove NA values and convert to character
    }
})
g <- graph_from_adj_list(idxssna, mode="all")
ajm<-as_adjacency_matrix(g, sparse = TRUE) ## this needs to be checked to see if there is any missing cell
dim(ajm)
colPair(spe) <- ajm

dim(spe)

## loading lig-recep and creating genes lig-rec graph
lr <- read.csv("/Users/inzirio/Downloads/human_lr_pair.txt",sep="\t")
# Filtra il dataframe lr
filtered_lr <- lr[lr$ligand_gene_symbol %in% rownames(spe) & lr$receptor_gene_symbol %in% rownames(spe), ]

# Visualizza il risultato
head(filtered_lr)
dim(filtered_lr)

gg <- graph_from_data_frame(filtered_lr[,c(2,3)])
ajmg <- as_adjacency_matrix(gg, sparse=FALSE)
dim(ajmg)
dim(spe)
m <- matrix(nrow=dim(spe)[1], ncol=dim(spe)[1], data=0)
rownames(m) <- colnames(m) <- rownames(spe)


# Controlla l'intersezione dei geni per allineare le matrici
common_genes <- intersect(rownames(ajmg), rownames(m))

# Assicurati che ajmg sia ordinata in base ai geni comuni
ajmg <- ajmg[common_genes, common_genes]

# Inserisci i valori di `1` dalla matrice ajmg nella matrice `m`
m[common_genes, common_genes] <- ajmg

m1 <- apply(m, 2, as.integer)
rowPair(spe) <- m1

###################
# trovare quelle coppie di lig-rec espressi dal grafo dei geni
# vedere il liv di espressione dalla matrice dei conteggi
# trovare le coppie di cellule in colpair che esprimono lig-rec
# farne un plot sull'immagine di riferimento, per singolo ligando
#

############################################### nnon funziona un cazzo
# Estrai la matrice di espressione
expression_matrix <- counts(spe)

# Ottieni le coppie di geni e cellule
row_pairs <- rowPair(spe)
col_pairs <- colPair(spe)

# Prepara un data frame per salvare i risultati
results <- data.frame(
    gene_from = rownames(expression_matrix)[row_pairs@from],
    gene_to = rownames(expression_matrix)[row_pairs@to],
    cell_pairs = I(vector("list", length(row_pairs)))  # Colonna per le coppie di indici delle cellule
)

# Itera su ogni coppia di geni
for (i in seq_along(row_pairs)) {
    # Geni della coppia
    gene_from <- row_pairs@from[i]
    gene_to <- row_pairs@to[i]

    # Trova le coppie di cellule in colPair
    cell_from_indices <- col_pairs@from
    cell_to_indices <- col_pairs@to

    # Filtra le coppie dove i geni sono espressi
    expressed_from <- expression_matrix[gene_from, cell_from_indices] > 0
    expressed_to <- expression_matrix[gene_to, cell_to_indices] > 0

    # Trova le coppie di indici dove entrambe le condizioni sono vere
    valid_pairs <- which(expressed_from & expressed_to)

    # Salva le coppie di indici nella colonna cell_pairs
    results$cell_pairs[[i]] <- cbind(cell_from = cell_from_indices[valid_pairs],
                                     cell_to = cell_to_indices[valid_pairs])
}

# Visualizza i risultati
head(results)

plot_lr_network <- function(spe, ligand, results) {
    # Ottieni la matrice di espressione
    expression_matrix <- counts(spe)

    # Identifica i cell_id che esprimono il ligand
    cells_expressing_ligand <- which(expression_matrix[ligand, ] > 0)

    # Identifica i gene_to associati al ligand in results
    related_rows <- results[results$gene_from == ligand, ]

    # Identifica i cell_id che esprimono i gene_to
    cells_expressing_receptor <- unique(unlist(related_rows$cell_pairs))

    # Crea un vettore per colorare i cell_id
    cell_colors <- rep("red", ncol(expression_matrix))  # Default rosso
    cell_colors[cells_expressing_ligand] <- "yellow"   # Giallo per ligand
    cell_colors[cells_expressing_receptor] <- "green"   # Verde per gene_to

    # Ottieni le coordinate delle cellule
    cell_coords <- data.frame(
        x = colData(spe)$CenterX_local_px,
        y = -colData(spe)$CenterY_local_px,
        color = cell_colors
    )

    # Plot
    library(ggplot2)
    ggplot(cell_coords, aes(x = x, y = y)) +
        geom_point(aes(color = color), size = 3) +
        scale_color_identity() +
        theme_void() +
        ggtitle(paste("Expression Network of", ligand))
}


ligand <- "ADM2"
library(ggplot2)
plot_lr_network(spe, ligand, results)


plot_lr_network_with_image <- function(spe, ligand, results, sample_id = unique(spe$sample_id), image_id) {
    # Assicurati che spe sia un oggetto SpatialExperiment
    stopifnot(is(spe, "SpatialExperiment"))

    # Ottieni la matrice di espressione
    expression_matrix <- counts(spe)

    # Identifica i cell_id che esprimono il ligand
    cells_expressing_ligand <- which(expression_matrix[ligand, ] > 0)

    # Identifica i gene_to associati al ligand in results
    related_rows <- results[results$gene_from == ligand, ]

    # Identifica i cell_id che esprimono i gene_to
    cells_expressing_receptor <- unique(unlist(related_rows$cell_pairs))

    # Crea un vettore per colorare i cell_id
    cell_colors <- rep("red", ncol(expression_matrix))  # Default rosso
    cell_colors[cells_expressing_ligand] <- "yellow"    # Giallo per ligand
    cell_colors[cells_expressing_receptor] <- "green"   # Verde per gene_to

    # Ottieni le coordinate delle cellule
    data <- as.data.frame(colData(spe)) %>%
        mutate(
            CenterY_local_px = -CenterY_local_px,
            color = cell_colors
        )

    # Carica l'immagine di sfondo
    image_raster <- imgRaster(spe, sample_id, image_id)

    # Plot
    library(ggplot2)
    ggplot() +
        annotation_raster(image_raster, xmin = 0, xmax = 4256, ymin = -4256, ymax = 0) +
        geom_text(
            data = data,
            aes(
                x = CenterX_local_px,
                y = CenterY_local_px,
                label = cell_id,
                color = color
            ),
            size = 2
        ) +
        scale_color_identity() +
        theme_void() +
        coord_fixed() +
        ggtitle(paste("Expression Network of", ligand))
}

plot_lr_network_with_image(spe, ligand, results, image_id="fov146")
















