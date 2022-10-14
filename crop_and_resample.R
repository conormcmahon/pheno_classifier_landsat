
library(tidyverse)
library(raster)
library(here)

# to do - fix messiness with choosing whether or not to reproject

directory <- "D:/landsat_vafb/temp/"
files <- list.files(path=directory, pattern="*_B4.TIF")

# Load an image, crop it coarsely, reproject to a new CRS, and then crop more finely
source(here::here("aggregate_custom.R"))
load_crop_reproject <- function(filename, extent1, interp_method="bilinear")
{
  temp <- raster(filename)
  # Only reproject if new raster is in a different CRS
  if(!compareCRS(temp,target_projection_raster))
  {
    return(aggregate_custom(temp,target_projection_raster,reproj_method=interp_method))
  }
  else
  {
    return(crop(temp,extent1))
  }
}

# Set up reprojections and cropping of Landsat scenes
# Huachuca / San Pedro Extent:
#subset_extent <- extent(547429,616400,3436600,3551845)
# Vandenber Extents: 
#subset_extent <- extent(710679, 739881, 3828117, 3881215)
#subset_extent_second <- extent(158406, 195962, 3827024, 3886773)# If multiple Landsat scenes with more than one projection are used, choose ONE scene with the desired projection:
# Santa Ynez Front Range
subset_extent <- extent(186696,281906,3798792,3840797)
files <- files[(grep("042",files))] # remove western tiles which don't wholly cover the target area 

# Huachuca Target Scene ID
#target_scene_ID <- "035038"
# Vandenberg Target Scene ID
#target_scene_ID <- "043036"
# Santa Ynez Front Range Target Scene ID
target_scene_ID <- "042036"
files_in_target_scene <- grep(target_scene_ID, files)
first_file_in_scene <- files[[files_in_target_scene[[1]]]]
target_projection_raster <- raster(paste(directory,"/",first_file_in_scene,sep=""))
target_projection_raster <- load_crop_reproject(paste(directory,"/",first_file_in_scene,sep=""), subset_extent)

generate_indices <- function(filename)
{
  filename_base <- substring(filename, 1, nchar(filename)-6)
  mission <- substr(filename,1,4)
  print(paste("Beginning work on scene ",filename_base," from mission ",mission,sep=""))
  if(mission %in% c("LT05","LT04"))
  {
    print("  Loading and cropping images!")
    blue <- load_crop_reproject(paste(directory,"/",filename_base, "B1.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected blue image.")
    green <- load_crop_reproject(paste(directory,"/",filename_base, "B2.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected green image.")
    red <- load_crop_reproject(paste(directory,"/",filename_base, "B3.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected red image.")
    NIR <- load_crop_reproject(paste(directory,"/",filename_base, "B4.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected NIR image.")
    SWIR <- load_crop_reproject(paste(directory,"/",filename_base, "B5.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected SWIR image.")
    SWIR2 <- load_crop_reproject(paste(directory,"/",filename_base, "B7.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected SWIR 2 image.")
    
    output_raster <- stack(blue,green,red,NIR,SWIR,SWIR2)
    output_raster <- output_raster
    writeRaster(output_raster, paste(directory,"/stacked_scenes_south/",filename_base,"stacked.TIF",sep=""),overwrite=TRUE)
  }
  if(mission=="LC08")
  {
    print("  Loading and cropping images!")
    ultrablue <- load_crop_reproject(paste(directory,"/",filename_base, "B1.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected ultrablue image.")
    blue <- load_crop_reproject(paste(directory,"/",filename_base, "B2.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected blue image.")
    green <- load_crop_reproject(paste(directory,"/",filename_base, "B3.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected green image.")
    red <- load_crop_reproject(paste(directory,"/",filename_base, "B4.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected red image.")
    NIR <- load_crop_reproject(paste(directory,"/",filename_base, "B5.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected NIR image.")
    SWIR <- load_crop_reproject(paste(directory,"/",filename_base, "B6.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected SWIR image.")
    SWIR2 <- load_crop_reproject(paste(directory,"/",filename_base, "B7.TIF",sep=""), subset_extent)
    print("  Cropped and reprojected SWIR 2 image.")
    
    output_raster <- stack(ultrablue,blue,green,red,NIR,SWIR,SWIR2)
    # Scale down reflectance values for MESMA
    output_raster <- output_raster
    writeRaster(output_raster, paste(directory,"/stacked_scenes_south/",filename_base,"stacked.TIF",sep=""),overwrite=TRUE)
  }
}

all_scenes <- lapply(files, generate_indices)