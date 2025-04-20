## creating data for cosmx starting from DBKero

f=11
c=c(1:50)
gc<- 1:12

countmat <- data.table::fread(countmat_file, showProgress=FALSE) # cell count matrix
c1 <- countmat[countmat$fov==f,]
c1 <- c1[c1$cell_ID%in%c,1:12]
write.csv(x=c1, file="/Users/inzirio/Library/CloudStorage/GoogleDrive-drighelli@gmail.com/My\ Drive/works/coding/SpaceTrooper/inst/extdata/CosMx_DBKero_Tiny/DBKero_exprMat_file.csv", row.names=FALSE)

metadata <- data.table::fread(metadata_file, showProgress=FALSE) # cell metadata
m1 <- metadata[metadata$fov==f,]
m1 <- m1[m1$cell_ID%in%c,]
write.csv(x=m1, file="/Users/inzirio/Library/CloudStorage/GoogleDrive-drighelli@gmail.com/My\ Drive/works/coding/SpaceTrooper/inst/extdata/CosMx_DBKero_Tiny/DBKero_metadata_file.csv", row.names=FALSE)

fov_positions <- as.data.frame(data.table::fread(fovpos_file, header=TRUE))
f1 <- fov_positions[fov_positions$fov==f,]
write.csv(x=f1, file="/Users/inzirio/Library/CloudStorage/GoogleDrive-drighelli@gmail.com/My\ Drive/works/coding/SpaceTrooper/inst/extdata/CosMx_DBKero_Tiny/DBKero_fov_positions_file.csv", row.names=FALSE)


spat_obj <- switch(type, csv=fread(polygonsFile),
                   parquet=read_parquet(polygonsFile))
s1 <- spat_obj[spat_obj$fov==f, ]
s1 <- s1[s1$cellID%in%c,]
write.csv(x=s1, file="/Users/inzirio/Library/CloudStorage/GoogleDrive-drighelli@gmail.com/My\ Drive/works/coding/SpaceTrooper/inst/extdata/CosMx_DBKero_Tiny/DBKero-polygons.csv", row.names=FALSE)

## data creation for merscope

library(devtools)
load_all()
# folder <- "/Users/inzirio/Downloads/MERSCOPE_data/Human_brain"
folder <- "/Users/inzirio/Downloaspeds/Merfish_data/human_uterine_cancer_patient2"
spe <- readMerfishSPE(folder, compute_missing_metrics=FALSE, keep_polygons=FALSE)
pols <- readPolygonsMerfish(metadata(spe)$polygons, type="parquet")
spe <- addPolygonsToSPE(spe, pols) # fov not present in polygons!!!!


