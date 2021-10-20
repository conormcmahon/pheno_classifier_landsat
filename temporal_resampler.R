
library(tidyverse)
library(raster)

# TODO 
#   Something is wrong with December interpolation for old files (Landsat 4/5)

directory <- "D:/landsat_vafb/data/stacked_scenes"

files <- list.files(path=directory, pattern="*NDVI.tif")
files <- substr(files, 1, nchar(files)-8)

target_years <- seq(1984,2010,2)

dates <- data.frame(year = as.numeric(as.character(substr(files,18,21))),
                    month = as.numeric(as.character(substr(files,22,23))),
                    day = as.numeric(as.character(substr(files,24,25))),
                    index = 1:length(files))
dates <- dates %>%
  mutate(doy = lubridate::yday(paste(year,month,day,sep="-")))

load_file <- function(index, type)
{
  # Get Dataset
  collection <- substr(files[[index]],6,9)
  dataset <- substr(files[[index]],1,4)
  # Load scene 
  scene <- raster(paste(directory,"/",files[[index]],type,".tif",sep=""))
  quality <- raster(paste(directory,"/",files[[index]],"quality.tif",sep=""))
#  if(collection=="L2SP")
#  {
#    good <- (quality==0)
#  }
  if(dataset %in% c("LT05","LT04"))
  {
    good <- (quality %in% c(0, 672,676,680,684, 5440))
  }
  else if(dataset=="LC08")
  {
    good <- (quality %in% c(21824))
  }
  else
  {
    print(paste("ERROR: dataset ",dataset," not in expected values."),sep="")
  }
  
  # Return raster as a data frame 
  data <- as.data.frame(scene, xy=TRUE)
  names(data) <- c("easting","northing","value")
  data$doy <- rep(dates[index,]$doy, nrow(data))
  
  goodness <- as.data.frame(good, xy=TRUE)
  names(goodness) <- c("easting","northing","good")
  bad_points <- sum(goodness$good==0)
  print(paste("  For scene ",paste(directory,"/",files[[index]],type,".tif",sep=""),sep=""))
  print(paste("  And quality raster ", paste(directory,"/",files[[index]],"quality.tif",sep=""),sep=""))
  print(paste("  Dropping ",
              bad_points," 
              points out of ",
              nrow(goodness),
              " total points, or ",
              bad_points/nrow(goodness)*100,"%.",
              sep=""))
  print("")
    
  temp <- cbind(data,goodness)
  temp <- temp[temp$good > 0,]
  data <- temp[,c(1,2,3,4,7)] %>%
    drop_na()
  
  return(data)
}

aggregate_scenes <- function(indices, type)
{
  scenes <- lapply(indices, load_file, type=type)
  year_data <- bind_rows(scenes)
}

# Get raster size, list of all Northing/Easting values, and map between pixel indices and coordinates
scene_for_crs <- raster(paste(directory,"/",files[[1]],"NDVI.tif",sep=""))
scene_for_crs_df <- as.data.frame(scene_for_crs, xy=TRUE)
northings <- sort(unique(scene_for_crs_df$y))
num_rows <- length(northings)
eastings <- sort(unique(scene_for_crs_df$x))
num_cols <- length(eastings)
output <- array(rep(0,num_rows*num_cols*12), c(num_rows,num_cols,12))
north_map <- 1:length(northings)
names(north_map) <- northings
east_map <- 1:length(eastings)
names(east_map) <- eastings

generate_pheno_raster <- function(type)
{
  for(current_year in target_years)
  {
    # Get dates for all scenes within this or the next year
    year_dates <- dates %>% filter(year>=current_year & year<=current_year+1)
    
    year_data <- aggregate_scenes(year_dates$index, type)
    
    output_filename <- paste(directory,"/pheno_",type,"_",current_year,".tif",sep="")
    
    print(paste("Starting work on ", output_filename, ", which is ",num_cols,"x",num_rows,"in size.", sep=""))
    
    year_data$northing <- as.character(year_data$northing)
    year_data$easting <- as.character(year_data$easting)
    ind_i <- 0
    for(northing_i in unique(northings))
    {
      ind_j <- 0
      data_i <- filter(year_data, northing==northing_i)
      data_i_list <- data_i %>% group_split(easting) #split(data_i, data_i[1])
      for(pixel in data_i_list)
      {
        if(nrow(pixel) < 2)
        {
          ind_j <- ind_j+1
          next
        }
        new_series <- approx(pixel$doy,pixel$value,xout=((1:12)*30))
        output[north_map[[pixel[1,]$northing]],east_map[[pixel[1,]$easting]],1:12] <- new_series$y
        ind_j <- ind_j+1
      }
      print(paste(ind_i/num_rows*100, "% done with file ", output_filename, sep=""))
      ind_i <- ind_i+1
    }
    
    pheno_raster <- flip(brick(output),direction=2)
    crs(pheno_raster) <- crs(scene_for_crs)
    origin(pheno_raster) <- origin(scene_for_crs)
    extent(pheno_raster) <- extent(scene_for_crs)
    writeRaster(pheno_raster, filename=output_filename,overwrite=TRUE)
  }
}

generate_pheno_raster("NDVI")
generate_pheno_raster("tc_bright")
generate_pheno_raster("tc_wet")
#generate_pheno_raster("tc_4")
#generate_pheno_raster("tc_5")
#generate_pheno_raster("tc_6")