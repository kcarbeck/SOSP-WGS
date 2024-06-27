# create species range mask file
# based on https://www.mountainmanmaier.com/software/pop_genom/ & https://nbviewer.org/github/brandonlind/offset_validation/blob/main/07_fit_gradient_forest_models_to_climate_data.ipynb#fit


# create a single .tif file based on the species range shapefile
# then save rasters as RDS, where all values within shapefile boundaries are set to 0 and those outside are NA
R

library(raster)
library(rgdal)
library(reshape2)
library(sf)

setwd("/workdir/kcarbeck/data")

# load shapefile of range
sh <- readOGR("/workdir/kcarbeck/data/eBird/merged_range_shapefile.shp")
    # OGR data source with driver: ESRI Shapefile
    # Source: "/local/workdir/kcarbeck/data/eBird", layer: "merged_range_shapefile"
    # with 1 features
    # It has 1 fields
    # Integer64 fields read as strings:  FID

# load raster 
r <- raster("/workdir/kcarbeck/data/2001_2011/standardized_MAT.asc")
    # class      : RasterLayer
    # dimensions : 1086, 3026, 3286236  (nrow, ncol, ncell)
    # resolution : 0.04166667, 0.04166667  (x, y)
    # extent     : -178.6875, -52.60417, 18.84866, 64.09866  (xmin, xmax, ymin, ymax)
    # crs        : +proj=longlat +datum=WGS84 +no_defs
    # source     : standardized_MAT.asc
    # names      : standardized_MAT


# project shapefile to match raster
projection(r) <- projection(sh)
masked_raster <- mask(r, sh)

# rasterize shapefile
mask <- rasterize(sh, masked_raster, field = 0)
cellStats(mask, min)

# save mask
writeRaster(mask, "/workdir/kcarbeck/data/maskdir/MAT_WGS84_mask_range.tif", format = "GTiff", overwrite = TRUE)
# save mask
saveRDS(mask, "/workdir/kcarbeck/data/maskdir/rangeWGS84_mask.RDS")






