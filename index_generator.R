
library(tidyverse)
library(raster)
library(here)

# to do - fix messiness with choosing whether or not to reproject

directory <- "D:/landsat_vafb/data/"
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
subset_extent_second <- extent(710679, 739881, 3828117, 3881215)
subset_extent <- extent(158406, 195962, 3827024, 3886773)
# Santa Ynez Front Range
#subset_extent <- extent(186696,281906,3798792,3840797) # for eastern tile, 042036 in UTM 11N
#files <- files[(grep("042",files))] # remove western tiles which don't wholly cover the target area 
#subset_extent <- extent(738349,833736,3801988,3836759) # for eastern tile, 043036 in UTM 10N
# If multiple Landsat scenes with more than one projection are used, choose ONE scene with the desired projection:

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
    # print("  Loading and cropping images!")
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
    if(file.exists(paste(directory,"/",filename_base, "BQA.TIF",sep="")))
    {
     quality <- load_crop_reproject(paste(directory,"/",filename_base, "BQA.TIF",sep=""), subset_extent, interp_method="ngb")
    }
    else if(file.exists(paste(directory,"/",filename_base, "CLOUD_QA.TIF",sep="")))
    {
     quality <- load_crop_reproject(paste(directory,"/",filename_base, "CLOUD_QA.TIF",sep=""), subset_extent, interp_method="ngb")
    }
    else
    {
     print("ERROR: Could not find a quality image for file ", filename_base, " with expected ending QBA.TIF or CLOUD_QA.TIF", sep="")
    }
    quality <- load_crop_reproject(paste(directory,"/",substr(filename_base,1,nchar(filename_base)-3), "QA_Pixel.tif",sep=""), subset_extent, interp_method="ngb")
    print("  Cropped and reprojected Quality image.")
    # 
    print("  Making Tasseled Cap values")
    tc_bright <- 0.2043*blue + 0.4158*green + 0.5524*red + 0.5741*NIR + 0.3124*SWIR + 0.2303*SWIR2
    tc_green <- -0.1603*blue + -0.2819*green + -0.4934*red + 0.7940*NIR + -0.0002*SWIR + -0.1446*SWIR2
    tc_wet <- 0.0315*blue + 0.2021*green + 0.3102*red + 0.1594*NIR + -0.6806*SWIR + -0.6109*SWIR2
    tc_4 <- -0.2117*blue + -0.0284*green + 0.1302*red + -0.1007*NIR + 0.6529*SWIR + -0.7078*SWIR2
    tc_5 <- -0.8669*blue + -0.1835*green + 0.3856*red + 0.0408*NIR + -0.1132*SWIR + 0.2272*SWIR2
    tc_6 <- 0.3677*blue + -0.8200*green + 0.4354*red + 0.0518*NIR + -0.0066*SWIR + -0.0104*SWIR2
    # 
    output_raster <- stack(blue,green,red,NIR,SWIR,SWIR2)
    writeRaster(output_raster, paste(directory,"/stacked_scenes/",filename_base,"stacked.TIF",sep=""),overwrite=TRUE)
  }
  if(mission=="LC08")
  {
    print(filename_base)
    print("  Loading and cropping images!")
    ultrablue <- load_crop_reproject(paste(directory,"/",filename_base, "B1.TIF",sep=""), subset_extent) * 0.0000275 - 0.2
    print("  Loading and cropping images!")
    blue <- load_crop_reproject(paste(directory,"/",filename_base, "B2.TIF",sep=""), subset_extent) * 0.0000275 - 0.2
    print("  Cropped and reprojected blue image.")
    green <- load_crop_reproject(paste(directory,"/",filename_base, "B3.TIF",sep=""), subset_extent) * 0.0000275 - 0.2
    print("  Cropped and reprojected green image.")
    red <- load_crop_reproject(paste(directory,"/",filename_base, "B4.TIF",sep=""), subset_extent) * 0.0000275 - 0.2
    print("  Cropped and reprojected red image.")
    NIR <- load_crop_reproject(paste(directory,"/",filename_base, "B5.TIF",sep=""), subset_extent) * 0.0000275 - 0.2
    print("  Cropped and reprojected NIR image.")
    SWIR <- load_crop_reproject(paste(directory,"/",filename_base, "B6.TIF",sep=""), subset_extent) * 0.0000275 - 0.2
    print("  Cropped and reprojected SWIR image.")
    SWIR2 <- load_crop_reproject(paste(directory,"/",filename_base, "B7.TIF",sep=""), subset_extent) * 0.0000275 - 0.2
    print("  Cropped and reprojected SWIR 2 image.")
    quality <- load_crop_reproject(paste(directory,"/",substr(filename_base,1,nchar(filename_base)-3), "QA_Pixel.tif",sep=""), subset_extent, interp_method="ngb")
#    quality <- load_crop_reproject(paste(directory,"/",filename_base, "BQA.TIF",sep=""), subset_extent, interp_method="ngb")
    print("  Cropped and reprojected Quality image.")
    
    tc_bright <- 0.3029*blue + 0.2786*green + 0.4733*red + 0.5599*NIR + 0.508*SWIR + 0.1872*SWIR2
    tc_green <- -0.2941*blue + -0.243*green + -0.5424*red + 0.7276*NIR + 0.0713*SWIR + -0.1608*SWIR2
    tc_wet <- 0.1511*blue + 0.1973*green + 0.3283*red + 0.3407*NIR + -0.7117*SWIR + -0.4559*SWIR2
    tc_4 <- -0.8239*blue + 0.0849*green + 0.4396*red + -0.058*NIR + 0.2013*SWIR + -0.2773*SWIR2
    tc_5 <- -0.3294*blue + 0.0557*green + 0.1056*red + 0.1855*NIR + -0.4349*SWIR + 0.8085*SWIR2
    tc_6 <- 0.1079*blue + -0.9023*green + 0.4119*red + 0.0575*NIR + -0.0259*SWIR + 0.0252*SWIR2
    
    output_raster <- stack(ultrablue,blue,green,red,NIR,SWIR,SWIR2)
    writeRaster(output_raster, paste(directory,"/stacked_scenes/",filename_base,"stacked.TIF",sep=""),overwrite=TRUE)
  }
  print("  Making NDVI and NDWI images")
  NDVI <- (NIR-red)/(NIR+red)
  NDWI <- (NIR-SWIR)/NIR+SWIR

  print("  Writing rasters to disk")
  writeRaster(NDVI, paste(directory,"/stacked_scenes/",filename_base,"NDVI.TIF",sep=""),overwrite=TRUE)
  print("  Saved one file")
  writeRaster(NDWI, paste(directory,"/stacked_scenes/",filename_base,"NDWI.TIF",sep=""),overwrite=TRUE)
  writeRaster(tc_bright, paste(directory,"/stacked_scenes/",filename_base,"tc_bright.TIF",sep=""),overwrite=TRUE)
  writeRaster(tc_green, paste(directory,"/stacked_scenes/",filename_base,"tc_green.TIF",sep=""),overwrite=TRUE)
  writeRaster(tc_wet, paste(directory,"/stacked_scenes/",filename_base,"tc_wet.TIF",sep=""),overwrite=TRUE)
  writeRaster(tc_4, paste(directory,"/stacked_scenes/",filename_base,"tc_4.TIF",sep=""),overwrite=TRUE)
  writeRaster(tc_5, paste(directory,"/stacked_scenes/",filename_base,"tc_5.TIF",sep=""),overwrite=TRUE)
  writeRaster(tc_6, paste(directory,"/stacked_scenes/",filename_base,"tc_6.TIF",sep=""),overwrite=TRUE)
  writeRaster(quality, paste(directory,"/stacked_scenes/",filename_base,"quality.TIF",sep=""),overwrite=TRUE)
}

all_scenes <- lapply(files, generate_indices)